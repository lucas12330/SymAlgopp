/**
 * @file test_ast.cpp
 * @brief Tests unitaires de l'AST et de EquationClassique.
 */

#include <chrono>
#include <limits>
#include <memory>
#include <type_traits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "test_framework.hpp"

using namespace symalgo;

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
    const ExprPtr ef = frac(2, 6);
    const Fraction* f = comme<Fraction>(ef);
    CHECK(f != nullptr);
    CHECK_EQ(f->getNum(), 1);
    CHECK_EQ(f->getDen(), 3);
    CHECK_EQ(texte(frac(2, 6)), std::string("(1/3)"));

    const ExprPtr eg = frac(1, -2);
    const Fraction* g = comme<Fraction>(eg);
    CHECK_EQ(g->getNum(), -1);
    CHECK_EQ(g->getDen(), 2);
    CHECK_NEAR(g->eval(0.0), -0.5, 1e-15);
}

TEST_CASE(fraction_denominateur_nul) {
    CHECK_THROWS(frac(1, 0), std::invalid_argument);
}

TEST_CASE(noeuds_crees_uniquement_par_les_helpers) {
    // Un noeud sur la pile serait libéré par erreur par son premier ExprPtr :
    // la construction directe est interdite à la compilation (clé CleFabrique)
    static_assert(!std::is_constructible_v<Constante, double>);
    static_assert(!std::is_constructible_v<Sinus, ExprPtr>);
    // Via les helpers, le partage fonctionne sans copie
    const ExprPtr e = cst(5.0);
    CHECK(e->clone().get() == e.get());
    CHECK_EQ(e->nombreReferences(), 1u);
    {
        const ExprPtr copie = e;
        CHECK_EQ(e->nombreReferences(), 2u);
    }
    CHECK_EQ(e->nombreReferences(), 1u);
}

TEST_CASE(equation_expression_nulle) {
    CHECK_THROWS(EquationClassique(nullptr), std::invalid_argument);
}

TEST_CASE(types_des_noeuds) {
    const ExprPtr x = var("x");
    CHECK(comme<Constante>(cst(1.0)) != nullptr);
    CHECK(comme<Fraction>(frac(1, 3)) != nullptr);
    CHECK(comme<Variable>(x) != nullptr);
    CHECK(comme<Parametre>(param("C")) != nullptr);
    CHECK(comme<Addition>(x + x) != nullptr);
    CHECK(comme<Soustraction>(x - x) != nullptr);
    CHECK(comme<Multiplication>(x * x) != nullptr);
    CHECK(comme<Division>(x / x) != nullptr);
    CHECK(comme<Puissance>(ast_pow(x, 2.0)) != nullptr);
    CHECK(comme<Sinus>(ast_sin(x)) != nullptr);
    CHECK(comme<Cosinus>(ast_cos(x)) != nullptr);
    CHECK(comme<Tangente>(ast_tan(x)) != nullptr);
    CHECK(comme<Exponentielle>(ast_exp(x)) != nullptr);
    CHECK(comme<Logarithme>(ast_ln(x)) != nullptr);
    // Groupes
    CHECK(comme<OperateurBinaire>(x * x) != nullptr);
    CHECK(comme<OperateurBinaire>(ast_sin(x)) == nullptr);
    CHECK(comme<FonctionUnaire>(ast_ln(x)) != nullptr);
    CHECK(comme<FonctionUnaire>(x) == nullptr);
    // Mauvais type
    CHECK(comme<Sinus>(ast_cos(x)) == nullptr);
    CHECK(comme<Constante>(frac(1, 3)) == nullptr);
    CHECK(comme<Sinus>(ExprPtr()) == nullptr);
}

TEST_CASE(fraction_dans_expression) {
    EquationClassique eq(frac(2, 6) * X);
    CHECK_NEAR(eq.eval(3.0), 1.0, 1e-15);
    const EquationClassique d = eq.derivee();
    CHECK_NEAR(d.eval(42.0), 1.0 / 3.0, 1e-15);
}

// ============================================================================
// Dérivation
// ============================================================================

TEST_CASE(derivee_polynome) {
    EquationClassique eq(ast_pow(X, 2) + (X * 5) + 6);
    const EquationClassique d = eq.derivee();
    CHECK_NEAR(d.eval(2.0), 9.0, 1e-12);
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
    const EquationClassique d = eq.derivee();
    CHECK_NEAR(eq.eval(0.0), 0.0, 1e-15);
    CHECK_NEAR(d.eval(0.0), 1.0, 1e-15);
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

TEST_CASE(simplification_petites_constantes) {
    // Bug corrigé : toute constante < 1e-9 était traitée comme nulle
    const double G = 6.674e-11;
    const ExprPtr e = (cst(G) * X)->simplifier();
    CHECK_NEAR(e->eval(2.0), 2.0 * G, 1e-25);
    CHECK_NEAR((X * cst(1e-12))->simplifier()->eval(1.0), 1e-12, 1e-27);
}

TEST_CASE(simplification_repli_des_constantes) {
    CHECK_EQ(texte((cst(2.0) * (cst(3.0) * X))->simplifier()), std::string("6 * x"));
    CHECK_EQ(texte(((cst(2.0) * X) * cst(3.0))->simplifier()), std::string("6 * x"));
    CHECK_EQ(texte((X * (cst(3.0) * ast_sin(X)))->simplifier()), std::string("3 * x * sin(x)"));
    CHECK_EQ(texte(((cst(6.0) * X) / 3.0)->simplifier()), std::string("2 * x"));
    CHECK_EQ(texte((cst(0.0) - X)->simplifier()), std::string("-1 * x"));
    // Primitive de sin(2x) : plus de « 0.5 * -1 * cos(2 * x) »
    CHECK_EQ(texte(ast_sin(cst(2.0) * X)->integrer()->simplifier()), std::string("-0.5 * cos(2 * x)"));
}

TEST_CASE(simplification_puissances) {
    CHECK_EQ(texte((X * X)->simplifier()), std::string("(x)^(2)"));
    CHECK_EQ(texte((ast_pow(X, 2.0) * X)->simplifier()), std::string("(x)^(3)"));
    CHECK_EQ(texte((ast_pow(X, 2.0) * ast_pow(X, -2.0))->simplifier()), std::string("1"));
    CHECK_EQ(texte((ast_sin(X) * ast_sin(X))->simplifier()), std::string("(sin(x))^(2)"));
}

TEST_CASE(fractions_arithmetique_exacte) {
    // Bug corrigé : 1/3 + 1/3 donnait 0.666667 (exactitude perdue)
    CHECK_EQ(texte((frac(1, 3) + frac(1, 3))->simplifier()), std::string("(2/3)"));
    CHECK_EQ(texte((frac(1, 3) + frac(1, 6))->simplifier()), std::string("(1/2)"));
    CHECK_EQ(texte((frac(1, 3) - frac(1, 3))->simplifier()), std::string("0"));
    CHECK_EQ(texte((frac(2, 3) * frac(3, 4))->simplifier()), std::string("(1/2)"));
    CHECK_EQ(texte((frac(2, 3) / frac(4, 9))->simplifier()), std::string("(3/2)"));
    CHECK_EQ(texte((frac(1, 3) * cst(2.0))->simplifier()), std::string("(2/3)"));
    CHECK_EQ(texte((frac(1, 3) * X + frac(1, 6) * X)->simplifier()), std::string("(1/2) * x"));
    // Entre deux Constante, le calcul reste en double
    CHECK_NEAR((cst(1.0) / cst(3.0))->simplifier()->eval(0.0), 1.0 / 3.0, 1e-16);
    // Constante non entière avec une fraction : double
    CHECK_NEAR((frac(1, 3) + cst(0.5))->simplifier()->eval(0.0), 1.0 / 3.0 + 0.5, 1e-15);
}

TEST_CASE(fractions_debordement) {
    // Le produit déborderait int64 : repli en double plutôt qu'un résultat faux
    const int64_t grand = int64_t(1) << 40;
    const ExprPtr r = (frac(1, grand) * frac(1, grand))->simplifier();
    CHECK_NEAR(r->eval(0.0), 1.0 / (double(grand) * double(grand)), 1e-40);
}

TEST_CASE(simplification_division_par_zero_conservee) {
    // 0/0 ne doit pas devenir 0
    const ExprPtr e = (cst(0.0) / cst(0.0))->simplifier();
    CHECK(std::isnan(e->eval(0.0)));
    CHECK(std::isinf((X / cst(0.0))->simplifier()->eval(1.0)));
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
    CHECK(comme<IntegraleNonEvaluee>(F) != nullptr);
    CHECK_EQ(texte(F), std::string("integrale(x * sin(x))"));
    CHECK_THROWS(F->eval(1.0), std::logic_error);
    // (∫f)' = f
    CHECK(F->derivee()->estEgal(*f));
    // Une somme dont un terme n'est pas intégrable garde le reste calculé
    const ExprPtr G = (ast_cos(X) + ast_sin(ast_pow(X, 2.0)))->integrer()->simplifier();
    CHECK_EQ(texte(G), std::string("(sin(x) + integrale(sin((x)^(2))))"));
    // Plus aucun 0 silencieux
    CHECK(comme<IntegraleNonEvaluee>(ast_exp(ast_pow(X, 2.0))->integrer()) != nullptr);
    CHECK(comme<IntegraleNonEvaluee>((ast_ln(X) / X)->integrer()) != nullptr);
}

// ============================================================================
// Limites
// ============================================================================

TEST_CASE(limite_continue) {
    CHECK_NEAR((ast_pow(X, 2) + 1.0)->limite(3.0)->simplifier()->eval(0.0), 10.0, 1e-12);
}

TEST_CASE(limite_hopital) {
    EquationClassique eq(ast_sin(X) / X);
    const EquationClassique l = eq.limite(0.0);
    CHECK_NEAR(l.eval(0.0), 1.0, 1e-12);

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
    return comme<LimiteNonEvaluee>(f->limite(a)->simplifier()) != nullptr;
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
    const EquationClassique dl = eq.DL(0.0, 2);
    for (double x : {-0.5, 0.1, 2.0}) {
        CHECK_NEAR(dl.eval(x), 1.0 + x + x * x, 1e-12);
    }
}

TEST_CASE(dl_exponentielle_en_1) {
    // DL de exp en a = 1 à l'ordre 4 : e * sum (x-1)^k / k!
    const ExprPtr dl = ast_exp(X)->DL(1.0, 4);
    const double h = 0.1;
    const double attendu = std::exp(1.0) * (1 + h + h * h / 2 + h * h * h / 6 + h * h * h * h / 24);
    CHECK_NEAR(dl->eval(1.0 + h), attendu, 1e-12);
}

TEST_CASE(dl_precision_ordre_eleve) {
    // |DL_n(a + h) - f(a + h)| = O(h^(n+1)) : avec h = 0.01 et n = 6, erreur ~1e-14
    struct Cas { ExprPtr f; double a; };
    const std::vector<Cas> cas = {
        {ast_exp(ast_sin(X)), 0.0},
        {ast_exp(ast_sin(X)), 0.7},
        {ast_ln(X + 1.0) / (X + 2.0), 0.0},
        {ast_tan(X) * ast_cos(cst(2.0) * X), 0.3},
        {ast_pow(X + 1.0, 0.5), 0.0},          // binôme généralisé
        {ast_pow(X, 3.0) - cst(2.0) * X, 0.0}, // puissance entière en un zéro de la base
        {ast_pow(X, X), 1.5},                  // u^v par exp(v ln u)
        {cst(1.0) / (cst(1.0) - X), 0.0},
        {ast_pow(ast_sin(X), 2.0), 0.0},
    };
    const double h = 0.01;
    for (const auto& c : cas) {
        const ExprPtr dl = c.f->DL(c.a, 6);
        CHECK_NEAR(dl->eval(c.a + h), c.f->eval(c.a + h), 1e-11);
        CHECK_NEAR(dl->eval(c.a - h), c.f->eval(c.a - h), 1e-11);
    }
}

TEST_CASE(dl_coefficients_exacts) {
    // exp : tous les coefficients 1/k! jusqu'à l'ordre 15 (1/13! = 1.6e-10 était perdu)
    const ExprPtr dl = ast_exp(X)->DL(0.0, 15);
    const double x = 1.0;
    double somme = 0.0, terme = 1.0;
    for (int k = 0; k <= 15; ++k) {
        somme += terme;
        terme /= (k + 1);
    }
    CHECK_NEAR(dl->eval(x), somme, 1e-15);
    // Série géométrique : 1 + x + ... + x^5
    CHECK_EQ(texte((cst(1.0) / (cst(1.0) - X))->DL(0.0, 3)),
             std::string("(((1 + x) + (x)^(2)) + (x)^(3))"));
    // sin : les coefficients pairs sont exactement nuls
    CHECK_EQ(texte(ast_sin(X)->DL(0.0, 4)), std::string("(x + -0.166667 * (x)^(3))"));
}

TEST_CASE(dl_ordre_eleve_rapide) {
    // Bug corrigé : l'ordre 10 prenait ~3 s (arbre de dérivées exponentiel)
    const auto debut = std::chrono::steady_clock::now();
    const ExprPtr dl = ast_exp(ast_sin(X))->DL(0.0, 20);
    const double duree = std::chrono::duration<double>(std::chrono::steady_clock::now() - debut).count();
    CHECK(duree < 0.1);
    CHECK_NEAR(dl->eval(0.1), std::exp(std::sin(0.1)), 1e-15);
}

TEST_CASE(dl_cas_limites) {
    CHECK_THROWS(ast_ln(X)->DL(0.0, 2), std::domain_error);          // ln non défini en 0
    CHECK_THROWS(ast_pow(X, 0.5)->DL(0.0, 2), std::domain_error);    // sqrt non dérivable en 0
    CHECK_THROWS(ast_exp(X)->DL(0.0, -1), std::invalid_argument);
    CHECK_EQ(texte(ast_cos(X)->DL(0.0, 0)), std::string("1"));
    CHECK_EQ(texte(ast_sin(X)->DL(0.0, 0)), std::string("0"));
    // Paramètre symbolique : pas de valeur numérique, erreur explicite
    CHECK_THROWS((param("C1") * ast_sin(X))->DL(0.0, 3), std::logic_error);
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
