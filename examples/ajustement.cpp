/**
 * @file ajustement.cpp
 * @brief Cas d'usage : ajustement d'un modèle non linéaire par la méthode de Newton, avec un
 *        gradient et une hessienne exacts.
 *
 * Décharge d'un condensateur mesurée avec du bruit : V(t) = a exp(-t/b) + c (amplitude a,
 * constante de temps b, tension résiduelle c). On minimise la somme des carrés des écarts
 * S(a, b, c). Le gradient et la hessienne de S sont dérivés symboliquement, compilés dans un
 * seul programme, puis Newton amorti converge en quelques itérations ; la hessienne donne
 * aussi l'incertitude sur chaque paramètre. Le programme renvoie un code non nul si
 * l'ajustement ne retrouve pas les paramètres ayant servi à générer les mesures.
 */

#include <symalgopp>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace symalgo;

int main() {
    // Mesures simulées : paramètres vrais (4.2, 1.8, 0.35), bruit de 0.05 V
    const double aVrai = 4.2, bVrai = 1.8, cVrai = 0.35, bruit = 0.05;
    const int nombreMesures = 40;
    std::mt19937_64 alea(7);
    std::normal_distribution<double> loi(0.0, bruit);
    std::vector<double> temps, tension;
    for (int i = 0; i < nombreMesures; ++i) {
        const double t = 0.25 * i;
        temps.push_back(t);
        tension.push_back(aVrai * std::exp(-t / bVrai) + cVrai + loi(alea));
    }

    // Modèle lu depuis du texte, une équation par mesure, sommées en un critère S
    const ExprPtr t = var("t");
    const ExprPtr modele = lire("a*exp(-t/b) + c", {"t", {"a", "b", "c"}});
    ExprPtr S = nombre(Nombre(0));
    for (int i = 0; i < nombreMesures; ++i) {
        const ExprPtr ecart = substituer(modele, t, cst(temps[i])) - cst(tension[i]);
        S = S + ecart * ecart;
    }

    const std::vector<ExprPtr> parametres = variables(S); // a, b, c
    const std::size_t n = parametres.size();

    const auto debut = std::chrono::steady_clock::now();
    const std::vector<ExprPtr> g = gradient(S, parametres);
    const auto H = hessienne(S, parametres);
    std::vector<ExprPtr> sorties = {S};
    sorties.insert(sorties.end(), g.begin(), g.end());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) sorties.push_back(H[i][j]);
    }
    const ProgrammeEvaluation programme(sorties, parametres);
    const double msDerivation =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - debut).count();

    // Évalue S, son gradient et sa hessienne en un point
    struct Etat {
        double S;
        Eigen::VectorXd grad;
        Eigen::MatrixXd hess;
    };
    const auto evaluerEn = [&](const Eigen::VectorXd& p) {
        const std::vector<double> v = programme.evaluerEn(std::vector<double>(p.data(), p.data() + n));
        Etat e{v[0], Eigen::VectorXd(n), Eigen::MatrixXd(n, n)};
        std::size_t k = 1;
        for (std::size_t i = 0; i < n; ++i) e.grad[i] = v[k++];
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i; j < n; ++j) e.hess(i, j) = e.hess(j, i) = v[k++];
        }
        return e;
    };

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Modèle : V(t) = " << modele << "   (" << nombreMesures << " mesures)\n";
    std::cout << "Dérivation et compilation de S, gradient et hessienne : " << std::setprecision(1) << msDerivation
              << " ms (" << programme.nombreInstructions() << " instructions)\n\n"
              << std::setprecision(5);

    // Newton amorti : pas de Newton (hessienne régularisée si nécessaire), divisé par 2 tant
    // que le critère ne baisse pas
    Eigen::VectorXd p(n);
    p << 3.0, 1.0, 0.0; // a, b, c : point de départ grossier
    Etat etat = evaluerEn(p);
    int iterations = 0;
    std::cout << "iter        S          |gradient|      a        b        c\n";
    for (; iterations < 50; ++iterations) {
        std::cout << std::setw(3) << iterations << std::setw(13) << etat.S << std::setw(14) << std::scientific
                  << std::setprecision(2) << etat.grad.norm() << std::fixed << std::setprecision(5) << std::setw(10)
                  << p[0] << std::setw(9) << p[1] << std::setw(9) << p[2] << "\n";
        if (etat.grad.norm() < 1e-9) break;
        double lambda = 0.0;
        Eigen::VectorXd pas;
        for (;; lambda = lambda == 0.0 ? 1e-3 : 10.0 * lambda) {
            const Eigen::LDLT<Eigen::MatrixXd> facto(etat.hess + lambda * Eigen::MatrixXd::Identity(n, n));
            pas = facto.solve(-etat.grad);
            if (facto.info() == Eigen::Success && facto.isPositive() && pas.dot(etat.grad) < 0.0) break;
        }
        double echelle = 1.0;
        Etat suivant = evaluerEn(p + pas);
        while (!(suivant.S < etat.S) && echelle > 1e-10) {
            echelle *= 0.5;
            suivant = evaluerEn(p + echelle * pas);
        }
        if (!(suivant.S < etat.S)) break;
        p += echelle * pas;
        etat = suivant;
    }

    // Incertitudes : covariance = 2 s² H^-1 avec s² = S / (mesures - paramètres)
    const double s2 = etat.S / (nombreMesures - static_cast<int>(n));
    const Eigen::MatrixXd covariance = 2.0 * s2 * etat.hess.inverse();
    const double vrais[] = {aVrai, bVrai, cVrai};
    const char* noms[] = {"a", "b", "c"};
    std::cout << "\nRésultat (" << iterations << " itérations), bruit estimé = " << std::sqrt(s2) << " V :\n";
    bool ok = iterations < 50;
    for (std::size_t i = 0; i < n; ++i) {
        const double sigma = std::sqrt(covariance(i, i));
        std::cout << "  " << noms[i] << " = " << p[i] << " ± " << sigma << "   (valeur vraie " << vrais[i] << ")\n";
        ok = ok && std::abs(p[i] - vrais[i]) < 4.0 * sigma;
    }
    std::cout << (ok ? "\nLes paramètres vrais sont retrouvés à 4 sigma près.\n"
                     : "\nÉCART : l'ajustement ne retrouve pas les paramètres vrais.\n");
    return ok ? 0 : 1;
}
