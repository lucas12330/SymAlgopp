/**
 * @file test_evaluateur.cpp
 * @brief Tests de l'évaluation compilée : identité avec l'évaluation directe de l'arbre,
 *        évaluation par blocs, partage des sous-expressions, erreurs.
 */

#include <cmath>
#include <stdexcept>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "Evaluateur.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");

std::vector<ExprPtr> expressionsDeTest() {
    return {
        cst(3.5),
        X,
        ast_pow(X, 2.0) + X * 5.0 + 6.0,
        ast_sin(X) * ast_cos(X) + ast_pow(X, 3.0) - cst(2.0) * X,
        ast_exp(ast_sin(X)) * ast_pow(X, 2.0),
        ast_ln(ast_pow(X, 2.0) + 1.0) / (X + 3.0),
        ast_pow(X, frac(1, 2)) + ast_pow(X, -3.0),       // exposants non entier et négatif
        ast_pow(cst(2.0), X) + ast_pow(X, X),           // exposants variables
        ast_tan(X / 3.0) - frac(2, 7) * ast_exp(-X),
        (ast_sin(X) + 1.0) * (ast_sin(X) + 2.0) / ast_pow(ast_sin(X) + 3.0, 2.0), // partage
        ast_exp(ast_sin(X))->derivee()->derivee()->derivee(),                      // gros graphe
        ast_asin(X / 10.0) + ast_acos(X / 8.0) * ast_atan(X) + pi() * X,           // réciproques, pi
    };
}

// Égalité à l'arrondi près (l'ordre des opérations peut différer légèrement)
void verifierProche(double obtenu, double attendu) {
    if (std::isnan(attendu)) {
        CHECK(std::isnan(obtenu));
    } else {
        CHECK_NEAR(obtenu, attendu, 1e-12 * std::max(1.0, std::abs(attendu)));
    }
}

} // namespace

TEST_CASE(compile_identique_a_l_arbre) {
    for (const ExprPtr& e : expressionsDeTest()) {
        const ProgrammeEvaluation p(e);
        CHECK(p.estValide());
        for (double x : {0.3, 1.0, 2.5, 7.0}) verifierProche(p.evaluer(x), e->eval(x));
    }
}

TEST_CASE(evaluation_par_blocs) {
    // Plus d'un bloc (256 points) et un bloc incomplet
    std::vector<double> xs;
    for (int i = 0; i < 1000; ++i) xs.push_back(0.01 + i * 0.013);
    for (const ExprPtr& e : expressionsDeTest()) {
        const ProgrammeEvaluation p(e);
        const std::vector<double> ys = p.evaluer(xs);
        CHECK_EQ(ys.size(), xs.size());
        for (std::size_t i = 0; i < xs.size(); i += 97) verifierProche(ys[i], e->eval(xs[i]));
        verifierProche(ys.back(), e->eval(xs.back()));
    }
    CHECK(ProgrammeEvaluation(X).evaluer(std::vector<double>{}).empty());
}

TEST_CASE(sous_expressions_partagees_calculees_une_fois) {
    // sin(x) apparaît trois fois mais n'est compilé qu'une fois
    const ExprPtr s = ast_sin(X);
    const ExprPtr e = (s + 1.0) * (s + 2.0) * (s + 4.0);
    const ProgrammeEvaluation p(e);
    // x, sin(x), trois sommes (Echelle/Axpy déjà fusionnés), deux produits
    CHECK(p.nombreInstructions() <= 12);
    // Réutilisation des registres : bien moins de registres que d'instructions
    const ProgrammeEvaluation gros(ast_exp(ast_sin(X))->derivee()->derivee()->derivee()->derivee());
    CHECK(gros.nombreRegistres() < gros.nombreInstructions() / 2);
}

TEST_CASE(valeurs_speciales) {
    const ProgrammeEvaluation inverse(cst(1.0) / X);
    CHECK(std::isinf(inverse.evaluer(0.0)));
    const ProgrammeEvaluation logarithme(ast_ln(X));
    CHECK(std::isnan(logarithme.evaluer(-1.0)));
}

TEST_CASE(expressions_non_evaluables) {
    const ProgrammeEvaluation p(param("C1") * X);
    CHECK(!p.estValide());
    CHECK_THROWS(p.evaluer(1.0), std::logic_error);
    CHECK_THROWS(p.evaluer(std::vector<double>{1.0, 2.0}), std::logic_error);
    const ProgrammeEvaluation q((X * ast_sin(X))->integrer());
    CHECK_THROWS(q.evaluer(1.0), std::logic_error);
}

TEST_CASE(equation_compilation_automatique) {
    const ExprPtr f = ast_sin(X) * ast_exp(X);
    EquationClassique eq(f);
    // Les premières évaluations passent par l'arbre, les suivantes par le programme compilé
    for (int i = 0; i < 20; ++i) verifierProche(eq.eval(0.1 * i), f->eval(0.1 * i));
    const std::vector<double> ys = eq.eval(std::vector<double>{0.5, 1.5});
    verifierProche(ys[1], f->eval(1.5));
    // simplifier() invalide le programme compilé
    EquationClassique g(cst(2.0) * X);
    for (int i = 0; i < 10; ++i) g.eval(1.0);
    g.simplifier();
    verifierProche(g.eval(3.0), 6.0);
    // Paramètre : même erreur, compilé ou non
    EquationClassique h(param("C") + X);
    for (int i = 0; i < 10; ++i) CHECK_THROWS(h.eval(1.0), std::logic_error);
}

int main() { return test::executerTous(); }
