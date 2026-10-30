/**
 * @file test_evaluateur.cpp
 * @brief Tests de l'évaluation compilée : identité avec l'évaluation directe de l'arbre,
 *        évaluation par blocs, partage des sous-expressions, erreurs.
 */

#include <algorithm>
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

TEST_CASE(abs_sign_hyperboliques_compiles) {
    const std::vector<ExprPtr> expressions = {
        ast_abs(X - 1.0),
        ast_sign(X - 2.0) * X,
        ast_sinh(X) + ast_cosh(X / 2.0),
        ast_tanh(ast_sin(X)) * ast_abs(ast_cos(X)),
        ast_sinh(X)->derivee() - ast_cosh(X) + ast_tanh(X)->derivee() * ast_abs(X),
        ast_abs(ast_sign(X - 1.5) - X),
    };
    std::vector<double> xs;
    for (int i = 0; i < 700; ++i) xs.push_back(-4.0 + i * 0.0123); // plus d'un bloc, passe par 0
    xs.push_back(0.0);
    xs.push_back(2.0);
    for (const ExprPtr& e : expressions) {
        const ProgrammeEvaluation p(e);
        CHECK(p.estValide());
        for (double x : {-3.2, -0.5, 0.0, 1.5, 2.0, 7.0}) verifierProche(p.evaluer(x), e->eval(x)); // scalaire
        const std::vector<double> ys = p.evaluer(xs);                                              // par blocs
        for (std::size_t i = 0; i < xs.size(); ++i) verifierProche(ys[i], e->eval(xs[i]));
    }
    // Signe à zéro et valeurs spéciales
    const ProgrammeEvaluation s(ast_sign(X));
    CHECK_NEAR(s.evaluer(0.0), 0.0, 0.0);
    CHECK_NEAR(s.evaluer(-0.0), 0.0, 0.0);
    CHECK_NEAR(s.evaluer(-5.0), -1.0, 0.0);
    CHECK_NEAR(ProgrammeEvaluation(ast_tanh(X)).evaluer(1e3), 1.0, 0.0);
    CHECK(std::isinf(ProgrammeEvaluation(ast_cosh(X)).evaluer(1e3)));
}

namespace {

// Écart en ulp (de la valeur de référence) entre une valeur calculée et une référence std::
double ecartUlp(double obtenu, double ref) {
    if (std::isnan(ref)) return std::isnan(obtenu) ? 0.0 : 1e300;
    if (std::isinf(ref) || ref == 0.0) return obtenu == ref ? 0.0 : 1e300;
    const double ulp = std::nextafter(std::fabs(ref), HUGE_VAL) - std::fabs(ref);
    return std::fabs(obtenu - ref) / ulp;
}

} // namespace

/*
 * Les fonctions hyperboliques maison (une seule exp/expm1) sont comparées à std:: de la glibc,
 * correctement arrondie (0,5 ulp). Écart mesuré hors tests contre long double sur 6 millions de
 * points : au plus 2,1 ulp ; vis-à-vis de std:: il ne dépasse donc pas 2,1 + 0,5 < 3 ulp.
 * Les chemins scalaire, arbre et par blocs doivent donner exactement les mêmes bits.
 */
TEST_CASE(hyperboliques_precision_et_coherence) {
    struct Fonction { ExprPtr e; double (*ref)(double); };
    const Fonction fonctions[] = {{ast_sinh(X), [](double v) { return std::sinh(v); }},
                                  {ast_cosh(X), [](double v) { return std::cosh(v); }},
                                  {ast_tanh(X), [](double v) { return std::tanh(v); }}};
    std::vector<double> xs;
    for (double x = -720.0; x <= 720.0; x += 0.0731) xs.push_back(x);        // grille large
    for (double x = -3.0; x <= 3.0; x += 0.00123) xs.push_back(x);           // autour de 0
    for (double m = 1e-300; m < 1e3; m *= 1.37) { xs.push_back(m); xs.push_back(-m); }
    for (double x : {0.0, -0.0, 5e-324, 1e-300, 3.7e-9, 0.55, 1.0, 21.99, 22.0, 709.0, 709.5, 710.4758600739439,
                     710.48, -710.4, 1e10, HUGE_VAL, -HUGE_VAL, std::nan("")})
        xs.push_back(x);
    for (const Fonction& f : fonctions) {
        const ProgrammeEvaluation p(f.e);
        const std::vector<double> ys = p.evaluer(xs);
        double pire = 0.0;
        for (std::size_t i = 0; i < xs.size(); ++i) {
            pire = std::max(pire, ecartUlp(ys[i], f.ref(xs[i])));
            const double s = p.evaluer(xs[i]);
            CHECK(s == ys[i] || (std::isnan(s) && std::isnan(ys[i])));      // scalaire == blocs
            const double d = f.e->eval(xs[i]);
            CHECK(d == ys[i] || (std::isnan(d) && std::isnan(ys[i])));      // arbre == blocs
        }
        CHECK(pire < 3.0);
    }
    // Signe de zéro et conservation de la parité
    CHECK(std::signbit(ProgrammeEvaluation(ast_sinh(X)).evaluer(-0.0)));
    CHECK(std::signbit(ProgrammeEvaluation(ast_tanh(X)).evaluer(-1e-20)));
    CHECK(std::isinf(ProgrammeEvaluation(ast_sinh(X)).evaluer(-1e3)) && ProgrammeEvaluation(ast_sinh(X)).evaluer(-1e3) < 0);
}

int main() { return test::executerTous(); }
