/**
 * @file test_solveur.cpp
 * @brief Tests du solveur d'équations : polynômes, isolement, familles trigonométriques,
 *        intervalles, recherche numérique, complétude.
 */

#include <cmath>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "Polynome.hpp"
#include "Solveur.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");
constexpr double PI = 3.14159265358979323846;

// "valeur1 ; valeur2 ; ..." (familles suivies de [k])
std::string texte(const Solutions& s) {
    std::string t;
    for (const Solution& sol : s.liste) {
        if (!t.empty()) t += " ; ";
        t += sol.exacte ? sol.valeur->texte() : "~" + std::to_string(sol.approximation);
        for (const std::string& k : sol.entiers) t += " [" + k + "]";
        if (sol.multiplicite > 1) t += " (x" + std::to_string(sol.multiplicite) + ")";
    }
    return t;
}

// Chaque solution (hors famille) vérifie l'équation
void verifierSolutions(const Solutions& s, const ExprPtr& f) {
    for (const Solution& sol : s.liste) {
        if (sol.estFamille()) continue;
        CHECK_NEAR(f->eval(sol.approximation), 0.0, 1e-9 * (1.0 + std::abs(sol.approximation)));
    }
}

} // namespace

TEST_CASE(equations_polynomiales) {
    CHECK_EQ(texte(resoudre(ast_pow(X, 2.0) - 3.0 * X + 2.0)), std::string("1 ; 2"));
    CHECK_EQ(texte(resoudre(cst(2.0) * X + 1.0, cst(7.0))), std::string("3"));
    CHECK_EQ(texte(resoudre(ast_pow(X, 2.0), cst(2.0))), std::string("-2^(1/2) ; 2^(1/2)"));
    CHECK_EQ(texte(resoudre(ast_pow(X - 1.0, 3.0) * (X + 2.0))), std::string("-2 ; 1 (x3)"));
    const Solutions s = resoudre(ast_pow(X, 5.0) - X - 1.0); // degré 5 : numérique certifié
    CHECK(s.complet);
    CHECK_EQ(s.liste.size(), std::size_t(1));
    CHECK(!s.liste[0].exacte);
    verifierSolutions(s, ast_pow(X, 5.0) - X - 1.0);
    CHECK(resoudre(ast_pow(X, 2.0) + 1.0).liste.empty());
    CHECK(resoudre(ast_pow(X, 2.0) + 1.0).complet);
}

TEST_CASE(isolement_fonctions) {
    CHECK_EQ(texte(resoudre(ast_exp(cst(2.0) * X - 1.0), cst(3.0))), std::string("ln(3)/2 + 1/2"));
    CHECK_EQ(texte(resoudre(ast_ln(X), cst(2.0))), std::string("exp(2)"));
    CHECK_EQ(texte(resoudre(ast_pow(cst(2.0), X), cst(8.0))), std::string("3"));
    CHECK_EQ(texte(resoudre(ast_pow(X, frac(1, 2)), cst(3.0))), std::string("9"));
    CHECK_EQ(texte(resoudre(ast_pow(X, 3.0), cst(-8.0))), std::string("-2"));
    CHECK_EQ(texte(resoudre(ast_asin(X), pi() / 6.0)), std::string("1/2"));
    CHECK_EQ(texte(resoudre(ast_atan(cst(2.0) * X), pi() / 4.0)), std::string("1/2"));
    // Pas de solution réelle, et c'est certain
    const Solutions aucune = resoudre(ast_exp(X), cst(-1.0));
    CHECK(aucune.liste.empty() && aucune.complet);
    CHECK(resoudre(ast_sin(X), cst(2.0)).liste.empty());
}

TEST_CASE(produit_nul_et_domaine) {
    CHECK_EQ(texte(resoudre(X * (X - 2.0) * ast_exp(X))), std::string("0 ; 2"));
    // x = 0 annule x mais pas x*ln(x), non définie en 0
    CHECK_EQ(texte(resoudre(X * ast_ln(X))), std::string("1"));
}

TEST_CASE(substitution) {
    // sin(x)^2 - sin(x) = 0 : sin(x) = 0 ou sin(x) = 1
    const Solutions s = resoudre(ast_pow(ast_sin(X), 2.0) - ast_sin(X));
    CHECK_EQ(s.liste.size(), std::size_t(3));
    CHECK(s.complet);
    // exp(2x) - 3 exp(x) + 2 = 0 : exp(x) = 1 ou exp(x) = 2
    CHECK_EQ(texte(resoudre(ast_exp(cst(2.0) * X) - 3.0 * ast_exp(X) + 2.0)), std::string("0 ; ln(2)"));
}

TEST_CASE(familles_trigonometriques) {
    CHECK_EQ(texte(resoudre(ast_sin(X), frac(1, 2))), std::string("pi/6 + 2*pi*k [k] ; 5*pi/6 + 2*pi*k [k]"));
    CHECK_EQ(texte(resoudre(ast_cos(cst(2.0) * X), cst(0.0))), std::string("pi/4 + pi*k [k] ; -pi/4 + pi*k [k]"));
    CHECK_EQ(texte(resoudre(ast_tan(X), cst(1.0))), std::string("pi/4 + pi*k [k]"));
    CHECK_EQ(texte(resoudre(ast_sin(X), cst(1.0))), std::string("pi/2 + 2*pi*k [k]"));
    // Chaque membre de la famille est solution
    const Solutions s = resoudre(ast_sin(cst(3.0) * X + 1.0), frac(1, 3));
    for (const Solution& sol : s.liste) {
        CHECK(sol.estFamille());
        for (int n : {-2, 0, 3}) {
            const double v = substituer(sol.valeur, param(sol.entiers[0]), cst(n))->eval(0.0);
            CHECK_NEAR(std::sin(3.0 * v + 1.0), 1.0 / 3.0, 1e-12);
        }
    }
}

TEST_CASE(resolution_sur_intervalle) {
    // Familles dépliées : valeurs exactes dans [0, 10]
    CHECK_EQ(texte(resoudreSurIntervalle(ast_sin(X) - frac(1, 2), 0.0, 10.0)),
             std::string("pi/6 ; 5*pi/6 ; 13*pi/6 ; 17*pi/6"));
    const Solutions s = resoudreSurIntervalle(ast_cos(X), -PI, PI);
    CHECK_EQ(texte(s), std::string("-pi/2 ; pi/2"));
    // Forme non résolue exactement : complétée numériquement, résultat non garanti complet
    const Solutions w = resoudreSurIntervalle(ast_exp(X) + X, -5.0, 5.0);
    CHECK(!w.complet);
    CHECK_EQ(w.liste.size(), std::size_t(1));
    CHECK_NEAR(w.liste[0].approximation, -0.5671432904097838, 1e-15); // constante oméga
}

TEST_CASE(resolution_numerique) {
    // Racine double sans changement de signe, racine simple, et un pôle écarté
    const ExprPtr f = ast_pow(X - 1.0, 2.0) * (X + 2.0) / (X - 3.0);
    const Solutions s = resoudreNumerique(f, -5.0, 5.0);
    CHECK(!s.complet);
    CHECK_EQ(s.liste.size(), std::size_t(2));
    CHECK_NEAR(s.liste[0].approximation, -2.0, 1e-12);
    CHECK_NEAR(s.liste[1].approximation, 1.0, 1e-7);
    CHECK(resoudreNumerique(cst(1.0) / X, -1.0, 1.0).liste.empty()); // pôle, pas racine
    CHECK(resoudreNumerique(param("a") * X, -1.0, 1.0).liste.empty()); // non évaluable
}

TEST_CASE(cas_particuliers) {
    CHECK(resoudre(ast_pow(X + 1.0, 2.0), developper(ast_pow(X + 1.0, 2.0))).toutReel);
    const Solutions incomplet = resoudre(ast_exp(X) + X);
    CHECK(!incomplet.complet);
    CHECK(resoudre(cst(3.0)).liste.empty());
    CHECK(resoudre(cst(0.0)).toutReel);
}

TEST_CASE(equation_classique_resoudre) {
    const EquationClassique eq(ast_pow(X, 3.0) - X);
    CHECK_EQ(texte(eq.resoudre()), std::string("-1 ; 0 ; 1"));
    CHECK_EQ(eq.factoriser().getExpression()->texte(), std::string("x*(x - 1)*(x + 1)"));
    CHECK_EQ(EquationClassique(ast_pow(X + 2.0, 2.0)).developper().getExpression()->texte(), std::string("x^2 + 4*x + 4"));
    CHECK_EQ(texte(EquationClassique(ast_sin(X)).resoudre(-1.0, 4.0)), std::string("0 ; pi"));
}

int main() { return test::executerTous(); }
