/**
 * @file test_differentielle.cpp
 * @brief Tests unitaires de EquationDifferentielle (affichage, matrice
 *        compagnon, résolution numérique RK4 et résolution littérale).
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"
#include "test_framework.hpp"

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
    std::unique_ptr<EquationClassique> sol(eq.resoudreLitteral());
    CHECK_EQ(capturerSortie([&] { sol->afficher(); }),
             std::string("(C1 * cos(2 * x) + C2 * sin(2 * x)) = 0\n"));
}

int main() { return test::executerTous(); }
