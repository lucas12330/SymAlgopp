/**
 * @file EquationDifferentielle.cpp
 * @author Lucas Bezanilla
 * @date 2026-04-15
 * @brief Fichier d'en-tête pour la sous-classe EquationDifferentielle.
 *
 * Projet SymAlgo++ :
 * Cette classe hérite de la classe de base Equation. Elle est dédiée
 * à la représentation et la résolution (numérique ou symbolique) des
 * équations différentielles.
 */

#include "EquationDifferentielle.hpp"
#include "EquationClassique.hpp"
#include "ASTNode.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

// Racine du polynôme caractéristique avec sa multiplicité
struct RacineCaracteristique {
    std::complex<double> valeur;
    int multiplicite;
};

/*
 * Nom : estRacineDeMultiplicite
 * Description : Vérifie que mu annule le polynôme de coefficients c (c[i] devant X^i)
 *               ainsi que ses m-1 premières dérivées, à une tolérance relative près
 *               (rapportée à la somme des modules des termes, comme une erreur inverse).
 */
bool estRacineDeMultiplicite(std::vector<double> c, std::complex<double> mu, int m) {
    for (int k = 0; k < m; ++k) {
        std::complex<double> valeur = 0.0, puissance = 1.0;
        double echelle = 0.0;
        for (double ci : c) {
            valeur += ci * puissance;
            echelle += std::abs(ci) * std::abs(puissance);
            puissance *= mu;
        }
        if (std::abs(valeur) > 1e-6 * echelle) return false;
        // Dérivée du polynôme
        for (size_t i = 1; i < c.size(); ++i) c[i - 1] = c[i] * static_cast<double>(i);
        c.pop_back();
    }
    return true;
}

/*
 * Nom : regrouperRacines
 * Description : Regroupe les racines numériques (Eigen) en racines multiples. Une racine
 *               de multiplicité m est calculée avec une erreur en eps^(1/m) (jusqu'à 1e-2
 *               pour m = 6) mais la moyenne du groupe reste exacte : on regroupe avec une
 *               tolérance large, puis on valide chaque groupe par les dérivées du polynôme.
 *               Un groupe non validé est traité comme des racines simples distinctes.
 */
std::vector<RacineCaracteristique> regrouperRacines(const std::vector<double>& coeffs,
                                                   const Eigen::VectorXcd& racines) {
    std::vector<RacineCaracteristique> resultat;
    std::vector<bool> utilisee(racines.size(), false);
    for (Eigen::Index i = 0; i < racines.size(); ++i) {
        if (utilisee[i]) continue;
        std::vector<Eigen::Index> groupe = {i};
        const double tolerance = 0.05 * std::max(1.0, std::abs(racines[i]));
        for (Eigen::Index j = i + 1; j < racines.size(); ++j) {
            if (!utilisee[j] && std::abs(racines[j] - racines[i]) <= tolerance) groupe.push_back(j);
        }
        std::complex<double> moyenne = 0.0;
        for (Eigen::Index j : groupe) moyenne += racines[j];
        moyenne /= static_cast<double>(groupe.size());

        const int m = static_cast<int>(groupe.size());
        if (m > 1 && !estRacineDeMultiplicite(coeffs, moyenne, m)) {
            utilisee[i] = true;
            resultat.push_back({racines[i], 1});
            continue;
        }
        for (Eigen::Index j : groupe) utilisee[j] = true;
        // Nettoie les résidus d'arrondi (ex. partie imaginaire 1e-17 d'une racine réelle)
        const double seuil = 1e-9 * std::max(1.0, std::abs(moyenne));
        if (std::abs(moyenne.imag()) < seuil) moyenne.imag(0.0);
        if (std::abs(moyenne.real()) < seuil) moyenne.real(0.0);
        resultat.push_back({moyenne, m});
    }
    return resultat;
}

// x^k (ou nullptr si k = 0)
ExprPtr puissanceDeX(const ExprPtr& X, int k) {
    if (k == 0) return nullptr;
    if (k == 1) return X;
    return ast_pow(X, static_cast<double>(k));
}

// Produit des facteurs non nuls
ExprPtr produit(std::initializer_list<ExprPtr> facteurs) {
    ExprPtr r = nullptr;
    for (const ExprPtr& f : facteurs) {
        if (f) r = r ? r * f : f;
    }
    return r ? r : cst(1.0);
}

} // namespace

EquationDifferentielle::EquationDifferentielle() {}

void EquationDifferentielle::ajouterTerme(unsigned int rang, double coeff) {
    double& total = m_terme[rang];
    total += coeff;
    // Un terme qui s'annule disparaît : le terme dominant n'est jamais nul, ce qui
    // évite une division par zéro dans la matrice compagnon
    if (total == 0.0) m_terme.erase(rang);
}

void EquationDifferentielle::setConditionsInitiales(const std::vector<double>& ci) {
    m_conditions_initiales = ci;
}

void EquationDifferentielle::afficher() const {
    bool first = true;
    for (auto it = m_terme.rbegin(); it != m_terme.rend(); ++it) {
        if (it->second == 0.0) continue;
        
        if (!first) {
            if (it->second > 0) std::cout << " + ";
            else std::cout << " - ";
        } else {
            if (it->second < 0) std::cout << "-";
        }
        
        double abs_coeff = std::abs(it->second);
        if (abs_coeff != 1.0) std::cout << abs_coeff << "*";

        if (it->first == 0) {
            std::cout << "y";
        } else if (it->first == 1) {
            std::cout << "y'";
        } else if (it->first == 2) {
            std::cout << "y''";
        } else {
            std::cout << "y^(" << it->first << ")";
        }
        first = false;
    }
    if (first) std::cout << "0";
    std::cout << " = 0" << std::endl;
}

Eigen::MatrixXd EquationDifferentielle::getMatriceCompagnon() const {
    if (m_terme.empty()) return Eigen::MatrixXd::Zero(1, 1);
    unsigned int n = m_terme.rbegin()->first;
    if (n == 0) return Eigen::MatrixXd::Zero(1, 1);
    
    double an = m_terme.rbegin()->second;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    for (unsigned int i = 0; i < n - 1; ++i) {
        A(i, i + 1) = 1.0;
    }
    for (unsigned int j = 0; j < n; ++j) {
        double aj = 0.0;
        if (m_terme.count(j)) aj = m_terme.at(j);
        A(n - 1, j) = -aj / an;
    }
    return A;
}

std::vector<ExprPtr> EquationDifferentielle::baseDeSolutions() const {
    std::vector<ExprPtr> base;
    if (m_terme.empty() || m_terme.rbegin()->first == 0) return base;

    const unsigned int n = m_terme.rbegin()->first;
    std::vector<double> coeffs(n + 1, 0.0); // polynôme caractéristique sum a_i r^i
    for (const auto& [rang, coeff] : m_terme) coeffs[rang] = coeff;

    const Eigen::EigenSolver<Eigen::MatrixXd> solveur(getMatriceCompagnon(), false);
    const auto X = var("x");
    for (const RacineCaracteristique& r : regrouperRacines(coeffs, solveur.eigenvalues())) {
        const double alpha = r.valeur.real();
        const double beta = r.valeur.imag();
        if (beta < 0.0) continue; // traitée avec sa conjuguée
        const ExprPtr exponentielle = alpha == 0.0 ? nullptr : ast_exp(cst(alpha) * X);
        for (int k = 0; k < r.multiplicite; ++k) {
            // Racine de multiplicité m : x^k e^(alpha x) [cos(beta x), sin(beta x)], k < m
            if (beta == 0.0) {
                base.push_back(produit({puissanceDeX(X, k), exponentielle}));
            } else {
                base.push_back(produit({puissanceDeX(X, k), exponentielle, ast_cos(cst(beta) * X)}));
                base.push_back(produit({puissanceDeX(X, k), exponentielle, ast_sin(cst(beta) * X)}));
            }
        }
    }
    return base;
}

EquationClassique* EquationDifferentielle::resoudreLitteral() const {
    ExprPtr solution = cst(0.0);
    int indice = 1;
    for (const ExprPtr& phi : baseDeSolutions()) {
        solution = solution + param("C" + std::to_string(indice++)) * phi;
    }
    return new EquationClassique(solution->simplifier());
}

EquationClassique* EquationDifferentielle::resoudreProblemeCauchy() const {
    const std::vector<ExprPtr> base = baseDeSolutions();
    const Eigen::Index n = static_cast<Eigen::Index>(base.size());
    if (n == 0) return new EquationClassique(cst(0.0));

    // M(i, j) = phi_j^(i)(0) ; les constantes c vérifient M c = (y(0), y'(0), ...)
    Eigen::MatrixXd M(n, n);
    for (Eigen::Index j = 0; j < n; ++j) {
        ExprPtr derivee = base[j];
        for (Eigen::Index i = 0; i < n; ++i) {
            M(i, j) = derivee->eval(0.0);
            derivee = derivee->derivee()->simplifier();
        }
    }
    Eigen::VectorXd conditions = Eigen::VectorXd::Zero(n);
    for (Eigen::Index i = 0; i < n && i < static_cast<Eigen::Index>(m_conditions_initiales.size()); ++i) {
        conditions(i) = m_conditions_initiales[i];
    }
    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(M);
    if (qr.rank() < n) {
        throw std::runtime_error("resoudreProblemeCauchy : base de solutions degeneree");
    }
    const Eigen::VectorXd c = qr.solve(conditions);

    const double echelle = std::max(1.0, c.cwiseAbs().maxCoeff());
    ExprPtr solution = cst(0.0);
    for (Eigen::Index j = 0; j < n; ++j) {
        if (std::abs(c(j)) > 1e-12 * echelle) solution = solution + cst(c(j)) * base[j];
    }
    return new EquationClassique(solution->simplifier());
}

double EquationDifferentielle::eval(double x) const {
    unsigned int n = m_terme.empty() ? 0 : m_terme.rbegin()->first;
    if (n == 0) return 0.0;
    
    Eigen::VectorXd Y = Eigen::VectorXd::Zero(n);
    for (size_t i = 0; i < std::min((size_t)n, m_conditions_initiales.size()); ++i) {
        Y(i) = m_conditions_initiales[i];
    }
    
    if (std::abs(x) < 1e-9) return Y(0);
    
    Eigen::MatrixXd A = getMatriceCompagnon();
    int steps = std::max(100, (int)(std::abs(x) / 0.01));
    double h = x / steps;
    
    for (int i = 0; i < steps; ++i) {
        Eigen::VectorXd k1 = A * Y;
        Eigen::VectorXd k2 = A * (Y + 0.5 * h * k1);
        Eigen::VectorXd k3 = A * (Y + 0.5 * h * k2);
        Eigen::VectorXd k4 = A * (Y + h * k3);
        Y += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    
    return Y(0);
}