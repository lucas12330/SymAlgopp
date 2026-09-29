/**
 * @file test_ast.cpp
 * @brief Tests unitaires de l'AST et de EquationClassique.
 */

#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "test_framework.hpp"

namespace {

// Représentation textuelle d'une expression (via afficher)
std::string texte(const ExprPtr& e) {
    std::ostringstream os;
    e->afficher(os);
    return os.str();
}

// Dérivée numérique par différence centrée, pour valider les dérivées formelles
double deriveeNumerique(const ExprPtr& f, double x, double h = 1e-5) {
    return (f->eval(x + h) - f->eval(x - h)) / (2.0 * h);
}

// Vérifie que d/dx[F] == f en plusieurs points (F primitive de f)
void verifierPrimitive(const ExprPtr& f, const std::vector<double>& points) {
    const ExprPtr F = f->integrer()->simplifier();
    for (double x : points) {
        CHECK_NEAR(deriveeNumerique(F, x), f->eval(x), 1e-6);
    }
}

const ExprPtr X = var("x");

} // namespace

// ============================================================================
// Évaluation
// ============================================================================

TEST_CASE(eval_polynome) {
    EquationClassique eq(ast_pow(X, 2) + (X * 5) + 6);
    CHECK_NEAR(eq.eval(2.0), 20.0, 1e-12);
    CHECK_NEAR(eq.eval(-3.0), 0.0, 1e-12);
}

TEST_CASE(eval_fonctions_usuelles) {
    CHECK_NEAR(ast_sin(X)->eval(0.5), std::sin(0.5), 1e-15);
    CHECK_NEAR(ast_cos(X)->eval(0.5), std::cos(0.5), 1e-15);
    CHECK_NEAR(ast_tan(X)->eval(0.5), std::tan(0.5), 1e-15);
    CHECK_NEAR(ast_exp(X)->eval(0.5), std::exp(0.5), 1e-15);
    CHECK_NEAR(ast_ln(X)->eval(0.5), std::log(0.5), 1e-15);
    CHECK_NEAR((cst(1.0) / X)->eval(4.0), 0.25, 1e-15);
    CHECK_NEAR((X - 3.0)->eval(1.0), -2.0, 1e-15);
}

// ============================================================================
// Fractions exactes
// ============================================================================

TEST_CASE(fraction_reduction_et_signe) {
    auto f = std::dynamic_pointer_cast<Fraction>(frac(2, 6));
    CHECK(f != nullptr);
    CHECK_EQ(f->getNum(), 1);
    CHECK_EQ(f->getDen(), 3);
    CHECK_EQ(texte(frac(2, 6)), std::string("(1/3)"));

    auto g = std::dynamic_pointer_cast<Fraction>(frac(1, -2));
    CHECK_EQ(g->getNum(), -1);
    CHECK_EQ(g->getDen(), 2);
    CHECK_NEAR(g->eval(0.0), -0.5, 1e-15);
}

TEST_CASE(fraction_denominateur_nul) {
    CHECK_THROWS(frac(1, 0), std::invalid_argument);
}

TEST_CASE(fraction_dans_expression) {
    EquationClassique eq(frac(2, 6) * X);
    CHECK_NEAR(eq.eval(3.0), 1.0, 1e-15);
    std::unique_ptr<EquationClassique> d(eq.derivee());
    CHECK_NEAR(d->eval(42.0), 1.0 / 3.0, 1e-15);
}

// ============================================================================
// Dérivation
// ============================================================================

TEST_CASE(derivee_polynome) {
    EquationClassique eq(ast_pow(X, 2) + (X * 5) + 6);
    std::unique_ptr<EquationClassique> d(eq.derivee());
    CHECK_NEAR(d->eval(2.0), 9.0, 1e-12);
}

TEST_CASE(derivee_comparee_a_la_derivee_numerique) {
    const std::vector<ExprPtr> expressions = {
        ast_sin(X) * ast_cos(X),
        ast_pow(X, 3) - cst(2.0) * X,
        ast_exp(cst(2.0) * X) / (X + 1.0),
        ast_ln(ast_pow(X, 2) + 1.0),
        ast_tan(X),
        ast_pow(X, X),                  // u(x)^v(x)
        ast_pow(cst(2.0), X),           // a^x
        ast_sin(ast_cos(X)),            // composition
    };
    for (const auto& e : expressions) {
        const ExprPtr d = e->derivee()->simplifier();
        for (double x : {0.3, 0.9, 1.7}) {
            CHECK_NEAR(d->eval(x), deriveeNumerique(e, x), 1e-6);
        }
    }
}

TEST_CASE(derivee_tangente) {
    EquationClassique eq(ast_tan(X));
    std::unique_ptr<EquationClassique> d(eq.derivee());
    CHECK_NEAR(eq.eval(0.0), 0.0, 1e-15);
    CHECK_NEAR(d->eval(0.0), 1.0, 1e-15);
}

// ============================================================================
// Paramètres symboliques (constantes C1, C2...)
// ============================================================================

TEST_CASE(parametre_derivee_nulle) {
    const ExprPtr C = param("C1");
    CHECK_EQ(texte(C->derivee()->simplifier()), std::string("0"));
    // d/dx (C1 * x) = C1 et non x + C1
    CHECK_EQ(texte((C * X)->derivee()->simplifier()), std::string("C1"));
    // d/dx (C1 * sin x) = C1 * cos x
    CHECK_EQ(texte((C * ast_sin(X))->derivee()->simplifier()), std::string("C1 * cos(x)"));
}

TEST_CASE(parametre_evaluation_impossible) {
    CHECK_THROWS(param("C1")->eval(5.0), std::logic_error);
    CHECK_THROWS((param("C1") * X)->eval(5.0), std::logic_error);
}

TEST_CASE(parametre_egalite) {
    CHECK(param("C1")->estEgal(*param("C1")));
    CHECK(!param("C1")->estEgal(*param("C2")));
    CHECK(!param("x")->estEgal(*var("x")));
    CHECK(!var("x")->estEgal(*param("x")));
}

TEST_CASE(variable_de_nom_quelconque) {
    // La bibliothèque est à une variable : var("v") est la variable d'évaluation
    const ExprPtr V = var("v");
    CHECK_NEAR((cst(5.0) * ast_pow(V, 2))->eval(10.0), 500.0, 1e-12);
    CHECK_NEAR((cst(5.0) * ast_pow(V, 2))->derivee()->simplifier()->eval(10.0), 100.0, 1e-12);
}

// ============================================================================
// Simplification
// ============================================================================

TEST_CASE(simplification_factorisation) {
    EquationClassique eq((cst(2) * ast_sin(X)) + (cst(3) * ast_sin(X)));
    eq.simplifier();
    CHECK_NEAR(eq.eval(1.0), 5.0 * std::sin(1.0), 1e-15);
    CHECK_EQ(texte((cst(2) * ast_sin(X) + cst(3) * ast_sin(X))->simplifier()),
             std::string("5 * sin(x)"));
}

TEST_CASE(simplification_elements_neutres) {
    CHECK_EQ(texte((X + 0.0)->simplifier()), std::string("x"));
    CHECK_EQ(texte((cst(0.0) + X)->simplifier()), std::string("x"));
    CHECK_EQ(texte((X * 1.0)->simplifier()), std::string("x"));
    CHECK_EQ(texte((X * 0.0)->simplifier()), std::string("0"));
    CHECK_EQ(texte((X - X)->simplifier()), std::string("0"));
    CHECK_EQ(texte((X / X)->simplifier()), std::string("1"));
    CHECK_EQ(texte(ast_pow(X, 1.0)->simplifier()), std::string("x"));
    CHECK_EQ(texte(ast_pow(X, 0.0)->simplifier()), std::string("1"));
    CHECK_EQ(texte(ast_exp(ast_ln(X))->simplifier()), std::string("x"));
}

// ============================================================================
// Intégration
// ============================================================================

TEST_CASE(integrale_polynome) {
    verifierPrimitive(ast_pow(X, 2) + (X * 5) + 6, {-1.0, 0.5, 2.0});
}

TEST_CASE(integrales_fonctions_usuelles) {
    verifierPrimitive(ast_sin(X), {0.2, 1.1});
    verifierPrimitive(ast_cos(X), {0.2, 1.1});
    verifierPrimitive(ast_tan(X), {0.2, 1.1});
    verifierPrimitive(ast_exp(X), {0.2, 1.1});
    verifierPrimitive(ast_exp(cst(3.0) * X), {0.2, 1.1});
    verifierPrimitive(ast_ln(X), {0.2, 1.1});
    verifierPrimitive(ast_pow(X, -1.0), {0.2, 1.1});
    verifierPrimitive(ast_sin(X) / 4.0, {0.2, 1.1});
}

TEST_CASE(integrales_substitution_lineaire) {
    verifierPrimitive(ast_sin(cst(2.0) * X), {0.2, 1.1});
    verifierPrimitive(ast_cos(X / 3.0 + 1.0), {0.2, 1.1});
    verifierPrimitive(ast_exp(cst(3.0) * X + 1.0), {0.2, 1.1});
    verifierPrimitive(ast_exp(cst(-0.5) * X), {0.2, 1.1});
    verifierPrimitive(ast_tan(cst(0.5) * X), {0.2, 1.1});
    verifierPrimitive(ast_ln(cst(2.0) * X + 1.0), {0.2, 1.1});
    verifierPrimitive(ast_pow(cst(2.0) * X + 1.0, 3.0), {0.2, 1.1});
    verifierPrimitive(ast_pow(cst(2.0), X), {0.2, 1.1});
    verifierPrimitive(ast_pow(cst(2.0), cst(3.0) * X), {0.2, 1.1});
}

TEST_CASE(integrales_quotients_et_produits) {
    verifierPrimitive(cst(1.0) / X, {0.2, 1.1});
    verifierPrimitive(cst(3.0) / (cst(2.0) * X + 1.0), {0.2, 1.1});
    verifierPrimitive(cst(1.0) / ast_pow(X, 2.0), {0.2, 1.1});
    verifierPrimitive(X * X, {0.2, 1.1});
    verifierPrimitive(ast_cos(X) * 4.0, {0.2, 1.1});
    // Facteur constant composé (ne dépend pas de x)
    verifierPrimitive(ast_sin(cst(1.0)) * X, {0.2, 1.1});
}

TEST_CASE(integrale_avec_parametre) {
    // ∫ C1 * cos(x) dx = C1 * sin(x)
    CHECK_EQ(texte((param("C1") * ast_cos(X))->integrer()->simplifier()),
             std::string("C1 * sin(x)"));
    CHECK_EQ(texte(param("C1")->integrer()->simplifier()), std::string("C1 * x"));
}

TEST_CASE(integrale_non_evaluee) {
    // x * sin(x) demanderait une intégration par parties : pas de résultat faux
    const ExprPtr f = X * ast_sin(X);
    const ExprPtr F = f->integrer();
    CHECK(std::dynamic_pointer_cast<IntegraleNonEvaluee>(F) != nullptr);
    CHECK_EQ(texte(F), std::string("integrale(x * sin(x))"));
    CHECK_THROWS(F->eval(1.0), std::logic_error);
    // (∫f)' = f
    CHECK(F->derivee()->estEgal(*f));
    // Une somme dont un terme n'est pas intégrable garde le reste calculé
    const ExprPtr G = (ast_cos(X) + ast_sin(ast_pow(X, 2.0)))->integrer()->simplifier();
    CHECK_EQ(texte(G), std::string("(sin(x) + integrale(sin((x)^(2))))"));
    // Plus aucun 0 silencieux
    CHECK(std::dynamic_pointer_cast<IntegraleNonEvaluee>(ast_exp(ast_pow(X, 2.0))->integrer()) != nullptr);
    CHECK(std::dynamic_pointer_cast<IntegraleNonEvaluee>((ast_ln(X) / X)->integrer()) != nullptr);
}

// ============================================================================
// Limites
// ============================================================================

TEST_CASE(limite_continue) {
    CHECK_NEAR((ast_pow(X, 2) + 1.0)->limite(3.0)->simplifier()->eval(0.0), 10.0, 1e-12);
}

TEST_CASE(limite_hopital) {
    EquationClassique eq(ast_sin(X) / X);
    std::unique_ptr<EquationClassique> l(eq.limite(0.0));
    CHECK_NEAR(l->eval(0.0), 1.0, 1e-12);

    // Exemple du README : (x^2 - 1) / (x - 1) en 1 -> 2
    CHECK_NEAR(((ast_pow(X, 2) - 1.0) / (X - 1.0))->limite(1.0)->simplifier()->eval(0.0),
               2.0, 1e-12);
    // (1 - cos x) / x^2 en 0 -> 1/2 (L'Hôpital appliqué deux fois)
    CHECK_NEAR(((cst(1.0) - ast_cos(X)) / ast_pow(X, 2))->limite(0.0)->simplifier()->eval(0.0),
               0.5, 1e-12);
}

namespace {

// Valeur numérique d'une limite déterminée
double valeurLimite(const ExprPtr& f, double a) { return f->limite(a)->simplifier()->eval(0.0); }

bool limiteNonDeterminee(const ExprPtr& f, double a) {
    return std::dynamic_pointer_cast<LimiteNonEvaluee>(f->limite(a)->simplifier()) != nullptr;
}

const double INF = std::numeric_limits<double>::infinity();

} // namespace

TEST_CASE(limite_zero_fois_infini) {
    // Bug corrigé : donnait NaN
    CHECK_NEAR(valeurLimite(X * (cst(1.0) / X), 0.0), 1.0, 1e-12);
    CHECK_NEAR(valeurLimite(ast_sin(X) * (cst(1.0) / X), 0.0), 1.0, 1e-12);
}

TEST_CASE(limite_signe_de_l_infini) {
    CHECK_EQ(valeurLimite(cst(1.0) / ast_pow(X, 2.0), 0.0), INF);
    CHECK_EQ(valeurLimite(cst(-1.0) / ast_pow(X, 2.0), 0.0), -INF);
    CHECK_EQ(valeurLimite(cst(1.0) / ast_pow(X - 2.0, 4.0), 2.0), INF);
    CHECK_EQ(valeurLimite(cst(3.0) / (cst(1.0) - ast_cos(X)), 0.0), INF);
    CHECK_EQ(valeurLimite(ast_pow(X, -2.0), 0.0), INF);
}

TEST_CASE(limite_inexistante) {
    // Bug corrigé : -1/x en 0 donnait +inf ; à gauche -inf, à droite +inf
    CHECK(limiteNonDeterminee(cst(-1.0) / X, 0.0));
    CHECK(limiteNonDeterminee(cst(1.0) / X, 0.0));
    CHECK(limiteNonDeterminee(ast_pow(X, -1.0), 0.0));
    CHECK(limiteNonDeterminee(ast_tan(X), 3.14159265358979323846 / 2.0));
    CHECK(limiteNonDeterminee(ast_sin(cst(1.0) / X), 0.0));
    CHECK_THROWS(((cst(1.0) / X)->limite(0.0)->eval(0.0)), std::logic_error);
}

TEST_CASE(limite_logarithme) {
    CHECK_EQ(valeurLimite(ast_ln(X), 0.0), -INF);
    CHECK_NEAR(valeurLimite(ast_ln(X), 1.0), 0.0, 1e-15);
    // Bug corrigé : ln d'un argument négatif donnait -inf
    CHECK(limiteNonDeterminee(ast_ln(X - 2.0), 1.0));
}

TEST_CASE(limite_hopital_avances) {
    // Résidu d'arrondi : sin(pi) ne vaut pas exactement 0 en double
    const double pi = 3.14159265358979323846;
    CHECK_NEAR(valeurLimite(ast_sin(X) / (X - pi), pi), -1.0, 1e-12);
    CHECK_NEAR(valeurLimite((ast_exp(X) - 1.0) / X, 0.0), 1.0, 1e-12);
    CHECK_NEAR(valeurLimite((X - ast_sin(X)) / ast_pow(X, 3.0), 0.0), 1.0 / 6.0, 1e-12);
    CHECK_NEAR(valeurLimite(ast_ln(X) / (X - 1.0), 1.0), 1.0, 1e-12);
}

TEST_CASE(limite_formes_exponentielles) {
    // (1 + x)^(1/x) -> e (forme 1^inf)
    CHECK_NEAR(valeurLimite(ast_pow(X + 1.0, cst(1.0) / X), 0.0), std::exp(1.0), 1e-12);
    // (1 + 2x)^(3/x) -> e^6
    CHECK_NEAR(valeurLimite(ast_pow(cst(2.0) * X + 1.0, cst(3.0) / X), 0.0), std::exp(6.0), 1e-9);
}

TEST_CASE(limite_avec_parametre) {
    // Les paramètres symboliques traversent le calcul
    CHECK_EQ(texte((param("C1") * X + 1.0)->limite(2.0)->simplifier()),
             std::string("(2 * C1 + 1)"));
}

TEST_CASE(limite_sans_recursion_infinie) {
    // L'Hôpital tournerait en rond (exp(1/x) sur x) : la profondeur est bornée
    const ExprPtr f = ast_exp(cst(-1.0) / ast_pow(X, 2.0)) / ast_pow(X, 2.0);
    const ExprPtr l = f->limite(0.0);
    CHECK(l != nullptr);
}

// ============================================================================
// Développements limités
// ============================================================================

TEST_CASE(dl_polynome) {
    EquationClassique eq(ast_pow(X, 3) + ast_pow(X, 2) + X + 1);
    std::unique_ptr<EquationClassique> dl(eq.DL(0.0, 2));
    for (double x : {-0.5, 0.1, 2.0}) {
        CHECK_NEAR(dl->eval(x), 1.0 + x + x * x, 1e-12);
    }
}

TEST_CASE(dl_exponentielle_en_1) {
    // DL de exp en a = 1 à l'ordre 4 : e * sum (x-1)^k / k!
    const ExprPtr dl = ast_exp(X)->DL(1.0, 4);
    const double h = 0.1;
    const double attendu = std::exp(1.0) * (1 + h + h * h / 2 + h * h * h / 6 + h * h * h * h / 24);
    CHECK_NEAR(dl->eval(1.0 + h), attendu, 1e-12);
}

// ============================================================================
// Tracé adaptatif
// ============================================================================

TEST_CASE(trace_points_ordonnes_et_exacts) {
    EquationClassique eq(ast_pow(X, 3) + ast_pow(X, 2) + X + 1);
    const auto points = eq.genererPointsTrace(0.0, 1.0, 0.1);
    CHECK(points.size() >= 11);
    CHECK_NEAR(points.front().first, 0.0, 1e-15);
    CHECK_NEAR(points.back().first, 1.0, 1e-12);
    for (size_t i = 0; i < points.size(); ++i) {
        CHECK_NEAR(points[i].second, eq.eval(points[i].first), 1e-15);
        if (i > 0) CHECK(points[i].first > points[i - 1].first);
    }
}

TEST_CASE(trace_raffine_les_zones_courbes) {
    EquationClassique droite(cst(2.0) * X + 1.0);
    EquationClassique courbe(ast_sin(cst(10.0) * X));
    const auto pDroite = droite.genererPointsTrace(0.0, 1.0, 1e-3);
    const auto pCourbe = courbe.genererPointsTrace(0.0, 1.0, 1e-3);
    CHECK_EQ(pDroite.size(), size_t(11));   // aucun raffinement sur une droite
    CHECK(pCourbe.size() > pDroite.size());
}

TEST_CASE(trace_intervalle_vide) {
    EquationClassique eq(X);
    CHECK(eq.genererPointsTrace(1.0, 1.0).empty());
    CHECK(eq.genererPointsTrace(2.0, 1.0).empty());
}

int main() { return test::executerTous(); }
