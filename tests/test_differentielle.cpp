/**
 * @file test_differentielle.cpp
 * @brief Tests unitaires de EquationDifferentielle (affichage, matrice
 *        compagnon, résolution numérique RK4 et résolution littérale).
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

// Capture ce qu'une fonction écrit sur std::cout
template <typename F>
std::string capturerSortie(F&& fonction) {
    std::ostringstream tampon;
    std::streambuf* ancien = std::cout.rdbuf(tampon.rdbuf());
    fonction();
    std::cout.rdbuf(ancien);
    return tampon.str();
}

EquationDifferentielle oscillateurHarmonique() {
    EquationDifferentielle eq; // y'' + 4y = 0
    eq.ajouterTerme(2, 1.0);
    eq.ajouterTerme(0, 4.0);
    return eq;
}

} // namespace

// ============================================================================
// Affichage
// ============================================================================

TEST_CASE(affichage_oscillateur) {
    const auto eq = oscillateurHarmonique();
    CHECK_EQ(capturerSortie([&] { eq.afficher(); }), std::string("y'' + 4*y = 0\n"));
}

TEST_CASE(affichage_equation_complexe) {
    EquationDifferentielle eq; // 3y''' - 2y' + 5y = 0
    eq.ajouterTerme(3, 3.0);
    eq.ajouterTerme(1, -2.0);
    eq.ajouterTerme(0, 5.0);
    CHECK_EQ(capturerSortie([&] { eq.afficher(); }), std::string("3*y^(3) - 2*y' + 5*y = 0\n"));
}

TEST_CASE(affichage_cumul_des_termes) {
    EquationDifferentielle eq;
    eq.ajouterTerme(1, 1.0);
    eq.ajouterTerme(1, 2.0); // cumulé : 3y'
    eq.ajouterTerme(0, -1.0);
    CHECK_EQ(capturerSortie([&] { eq.afficher(); }), std::string("3*y' - y = 0\n"));
}

// ============================================================================
// Matrice compagnon
// ============================================================================

TEST_CASE(matrice_compagnon_ordre_3) {
    EquationDifferentielle eq; // 3y''' - 2y' + 5y = 0  =>  y''' = -5/3 y + 2/3 y'
    eq.ajouterTerme(3, 3.0);
    eq.ajouterTerme(1, -2.0);
    eq.ajouterTerme(0, 5.0);
    const Eigen::MatrixXd A = eq.getMatriceCompagnon();
    CHECK_EQ(A.rows(), 3);
    CHECK_EQ(A.cols(), 3);
    CHECK_NEAR(A(0, 1), 1.0, 1e-15);
    CHECK_NEAR(A(1, 2), 1.0, 1e-15);
    CHECK_NEAR(A(2, 0), -5.0 / 3.0, 1e-15);
    CHECK_NEAR(A(2, 1), 2.0 / 3.0, 1e-15);
    CHECK_NEAR(A(2, 2), 0.0, 1e-15);
    CHECK_NEAR(A(0, 0), 0.0, 1e-15);
}

// ============================================================================
// Résolution numérique (RK4)
// ============================================================================

TEST_CASE(rk4_oscillateur_contre_solution_exacte) {
    auto eq = oscillateurHarmonique();
    eq.setConditionsInitiales({1.0, 0.0}); // solution exacte : cos(2x)
    for (double x : {0.0, 0.25, 3.14159265358979 / 4.0, 1.0, 2.5, -1.0}) {
        CHECK_NEAR(eq.eval(x), std::cos(2.0 * x), 1e-7);
    }
}

TEST_CASE(rk4_decroissance_exponentielle) {
    EquationDifferentielle eq; // y' + 2y = 0, y(0) = 3  =>  y = 3 e^{-2x}
    eq.ajouterTerme(1, 1.0);
    eq.ajouterTerme(0, 2.0);
    eq.setConditionsInitiales({3.0});
    for (double x : {0.5, 1.0, 2.0}) {
        CHECK_NEAR(eq.eval(x), 3.0 * std::exp(-2.0 * x), 1e-8);
    }
}

TEST_CASE(rk4_sans_conditions_initiales) {
    auto eq = oscillateurHarmonique(); // conditions initiales nulles => y = 0
    CHECK_NEAR(eq.eval(1.0), 0.0, 1e-15);
}

// ============================================================================
// Résolution littérale
// ============================================================================

TEST_CASE(solution_litterale_oscillateur) {
    const auto eq = oscillateurHarmonique();
    const EquationClassique sol = eq.resoudreLitteral();
    CHECK_EQ(capturerSortie([&] { sol.afficher(); }),
             std::string("C1*cos(2*x) + C2*sin(2*x) = 0\n"));
    // Les constantes sont des paramètres symboliques : pas d'évaluation silencieuse
    CHECK_THROWS(sol.eval(1.0), std::logic_error);
}

namespace {

// EDO dont le polynôme caractéristique a pour coefficients c[i] devant r^i
EquationDifferentielle depuisCoefficients(const std::vector<double>& c) {
    EquationDifferentielle eq;
    for (size_t i = 0; i < c.size(); ++i) {
        if (c[i] != 0.0) eq.ajouterTerme(static_cast<unsigned int>(i), c[i]);
    }
    return eq;
}

// Vérifie que y satisfait sum c_i y^(i) = 0 en quelques points
void verifierSatisfaitEquation(const std::vector<double>& c, const EquationClassique& y) {
    for (double x : {-0.7, 0.4, 1.3}) {
        double somme = 0.0, echelle = 0.0;
        EquationClassique d = y;
        for (size_t i = 0; i < c.size(); ++i) {
            const double v = d.eval(x);
            somme += c[i] * v;
            echelle += std::abs(c[i] * v);
            d = d.derivee();
        }
        CHECK_NEAR(somme, 0.0, 1e-9 * std::max(1.0, echelle));
    }
}

// Compare la solution exacte du problème de Cauchy à RK4 et vérifie l'équation
void verifierCauchy(const std::vector<double>& c, const std::vector<double>& ci) {
    auto eq = depuisCoefficients(c);
    eq.setConditionsInitiales(ci);
    const EquationClassique y = eq.resoudreProblemeCauchy();
    for (double x : {0.0, 0.5, 1.0, -0.8}) {
        const double exacte = y.eval(x);
        CHECK_NEAR(exacte, eq.eval(x), 1e-6 * std::max(1.0, std::abs(exacte)));
    }
    verifierSatisfaitEquation(c, y);
}

std::string solutionGenerale(const std::vector<double>& c) {
    const auto eq = depuisCoefficients(c);
    const EquationClassique sol = eq.resoudreLitteral();
    return capturerSortie([&] { sol.afficher(); });
}

} // namespace

TEST_CASE(racine_double_solution_generale) {
    // y'' + 2y' + y = 0 : r = -1 double => (C1 + C2 x) e^(-x). Bug corrigé : C1 e^-x + C2 e^-x
    CHECK_EQ(solutionGenerale({1.0, 2.0, 1.0}),
             std::string("C1*exp(-x) + C2*x*exp(-x) = 0\n"));
    // y'' = 0 : r = 0 double => C1 + C2 x
    CHECK_EQ(solutionGenerale({0.0, 0.0, 1.0}), std::string("C2*x + C1 = 0\n"));
}

TEST_CASE(cauchy_racines_simples) {
    verifierCauchy({4.0, 0.0, 1.0}, {1.0, 0.0});          // oscillateur : cos(2x)
    verifierCauchy({4.0, 0.2, 1.0}, {1.0, 0.0});          // ressort amorti
    verifierCauchy({5.0, -2.0, 0.0, 3.0}, {1.0, -1.0, 0.5}); // 3y''' - 2y' + 5y
    verifierCauchy({-6.0, 1.0, 1.0}, {2.0, 1.0});         // racines 2 et -3
    verifierCauchy({2.0, 1.0}, {3.0});                    // ordre 1
}

TEST_CASE(cauchy_racines_multiples) {
    verifierCauchy({1.0, 2.0, 1.0}, {1.0, 1.0});                  // (D+1)^2
    verifierCauchy({-8.0, 12.0, -6.0, 1.0}, {1.0, 0.0, -1.0});   // (D-2)^3
    verifierCauchy({1.0, 4.0, 6.0, 4.0, 1.0}, {1.0, -1.0, 2.0, 0.5}); // (D+1)^4
    verifierCauchy({1.0, 0.0, 2.0, 0.0, 1.0}, {0.0, 1.0, 0.0, 0.0}); // (D^2+1)^2 : complexes doubles
    verifierCauchy({0.0, 0.0, 1.0}, {2.0, -3.0});                 // y'' = 0 : 2 - 3x
}

TEST_CASE(racines_proches_mais_distinctes) {
    // r = 1 et r = 1.02 : ne doivent pas être fusionnées en racine double
    verifierCauchy({1.02, -2.02, 1.0}, {1.0, 0.0});
}

TEST_CASE(cauchy_solution_exacte_connue) {
    auto eq = oscillateurHarmonique();
    eq.setConditionsInitiales({1.0, 0.0});
    const EquationClassique y = eq.resoudreProblemeCauchy();
    CHECK_EQ(capturerSortie([&] { y.afficher(); }), std::string("cos(2*x) = 0\n"));
}

TEST_CASE(annulation_du_terme_dominant) {
    EquationDifferentielle eq; // 1*y'' - 1*y'' + y' + 2y = 0 : équation d'ordre 1
    eq.ajouterTerme(2, 1.0);
    eq.ajouterTerme(2, -1.0);
    eq.ajouterTerme(1, 1.0);
    eq.ajouterTerme(0, 2.0);
    const Eigen::MatrixXd A = eq.getMatriceCompagnon();
    CHECK_EQ(A.rows(), 1);
    CHECK_NEAR(A(0, 0), -2.0, 1e-15);
    CHECK_EQ(capturerSortie([&] { eq.afficher(); }), std::string("y' + 2*y = 0\n"));
}

TEST_CASE(derivee_de_l_equation_differentielle) {
    // La dérivée de la solution de Cauchy, calculée par RK4 sur l'EDO dérivée,
    // doit coïncider avec la dérivée symbolique de la solution exacte
    for (const auto& c : {std::vector<double>{4.0, 0.2, 1.0}, std::vector<double>{5.0, -2.0, 0.0, 3.0}}) {
        auto eq = depuisCoefficients(c);
        eq.setConditionsInitiales({1.0, -0.5, 0.25});
        const EquationClassique yPrime = eq.resoudreProblemeCauchy().derivee();
        const EquationDifferentielle dEq = eq.derivee();
        for (double x : {0.0, 0.6, 1.5}) {
            CHECK_NEAR(dEq.eval(x), yPrime.eval(x), 1e-6);
        }
    }
}

TEST_CASE(derivee_generique_polymorphe) {
    auto edo = oscillateurHarmonique();
    edo.setConditionsInitiales({1.0, 0.0}); // y = cos(2x), y' = -2 sin(2x)
    const EquationClassique classique(ast_sin(var("x")));
    const std::vector<const Equation*> equations = {&edo, &classique};
    const std::unique_ptr<Equation> d0 = equations[0]->deriveeGenerique();
    const std::unique_ptr<Equation> d1 = equations[1]->deriveeGenerique();
    CHECK_NEAR(d0->eval(0.3), -2.0 * std::sin(0.6), 1e-7);
    CHECK_NEAR(d1->eval(0.3), std::cos(0.3), 1e-15);
}

int main() { return test::executerTous(); }
