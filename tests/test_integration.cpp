/**
 * @file test_integration.cpp
 * @brief Tests de l'intégrale définie : valeur exacte, quadrature adaptative, bornes
 *        infinies, singularités (jamais de valeur exacte fausse), paramètres, orientation.
 */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "Integration.hpp"
#include "Multivariable.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");
constexpr double PI = 3.14159265358979323846;
constexpr double INF = std::numeric_limits<double>::infinity();

OptionsIntegration numerique() {
    OptionsIntegration o;
    o.exacteSiPossible = false;
    return o;
}

// Calcule par les deux chemins et vérifie qu'ils s'accordent
void accord(const ExprPtr& f, const ExprPtr& a, const ExprPtr& b, const std::string& attendu, double tolerance = 1e-9) {
    const ResultatIntegrale e = integrerDefinie(f, a, b);
    const ResultatIntegrale n = integrerDefinie(f, a, b, numerique());
    CHECK(e.estExacte());
    CHECK(!n.estExacte());
    if (e.estExacte()) CHECK_EQ(e.exacte->texte(), attendu);
    CHECK_NEAR(e.valeur, n.valeur, tolerance);
    CHECK(n.erreurEstimee <= 1e-9 * (1.0 + std::abs(n.valeur)));
    CHECK_EQ(e.erreurEstimee, 0.0);
}

} // namespace

TEST_CASE(exact_polynomes_et_trigonometrie) {
    accord(ast_pow(X, 2.0), cst(0.0), cst(1.0), "1/3");
    accord(ast_sin(X), cst(0.0), pi(), "2");
    accord(ast_cos(X), cst(0.0), pi() / 2.0, "1");
    accord(ast_exp(X), cst(0.0), cst(1.0), "exp(1) - 1");
    accord(cst(3.0) * ast_pow(X, 2.0) + cst(2.0) * X + 1.0, cst(-1.0), cst(2.0), "15");
    accord(ast_pow(X, 3.0), frac(1, 2), frac(3, 2), "5/4");
    // Bornes du type ExprPtr : de 1 à e
    accord(1.0 / X, cst(1.0), ast_exp(cst(1.0)), "1");
    accord(ast_exp(X), cst(0.0), ast_ln(cst(2.0)), "1");
    // sin entre pi/6 et pi/2 : cos(pi/6) = 3^(1/2)/2
    accord(ast_sin(X), pi() / 6.0, pi() / 2.0, "3^(1/2)/2");
}

TEST_CASE(surcharge_double_bornes_entieres) {
    const ResultatIntegrale r = integrerDefinie(ast_pow(X, 2.0), 0.0, 1.0);
    CHECK(r.estExacte());
    CHECK_EQ(r.exacte->texte(), std::string("1/3"));
    CHECK_NEAR(r.valeur, 1.0 / 3.0, 1e-15);
    // Borne non entière : numérique seulement (cst(0.5) serait un réel)
    const ResultatIntegrale s = integrerDefinie(ast_pow(X, 2.0), 0.0, 0.5);
    CHECK(!s.estExacte());
    CHECK_NEAR(s.valeur, 1.0 / 24.0, 1e-13);
    CHECK(s.erreurEstimee < 1e-12);
}

TEST_CASE(numerique_primitive_inconnue) {
    // La primitive de 1/(1+x^2) n'est pas connue du moteur symbolique : quadrature
    const ResultatIntegrale r = integrerDefinie(1.0 / (X * X + 1.0), 0.0, 1.0);
    CHECK(!r.estExacte());
    CHECK_NEAR(r.valeur, PI / 4.0, 1e-12);
    CHECK(r.erreurEstimee <= 1e-10);
    CHECK(std::abs(r.valeur - PI / 4.0) <= r.erreurEstimee + 1e-15);

    // exp(-x^2) sur [0, 1] : erf(1) * sqrt(pi) / 2
    CHECK_NEAR(integrerDefinie(ast_exp(-X * X), 0.0, 1.0).valeur, std::sqrt(PI) / 2.0 * std::erf(1.0), 1e-12);
    // Oscillations : sin(x^2) sur [0, 10]
    CHECK_NEAR(integrerDefinie(ast_sin(X * X), 0.0, 10.0).valeur, 0.5836708999296392, 1e-8);
    // Fonction discontinue de dérivée : |x| écrit avec une racine, sur [-1, 1]
    CHECK_NEAR(integrerDefinie(ast_pow(X * X, frac(1, 2)), -1.0, 1.0).valeur, 1.0, 1e-9);
}

TEST_CASE(bornes_infinies) {
    const ResultatIntegrale a = integrerDefinie(ast_exp(-X), 0.0, INF);
    CHECK(!a.estExacte());
    CHECK_NEAR(a.valeur, 1.0, 1e-9);
    CHECK(a.erreurEstimee < 1e-8);

    CHECK_NEAR(integrerDefinie(ast_exp(-X * X), -INF, INF).valeur, std::sqrt(PI), 1e-9);
    CHECK_NEAR(integrerDefinie(ast_exp(-X * X), 0.0, INF).valeur, std::sqrt(PI) / 2.0, 1e-9);
    CHECK_NEAR(integrerDefinie(ast_exp(X), -INF, 0.0).valeur, 1.0, 1e-9);
    CHECK_NEAR(integrerDefinie(ast_exp(X), -INF, 3.0).valeur, std::exp(3.0), 1e-8);
    CHECK_NEAR(integrerDefinie(1.0 / (X * X + 1.0), -INF, INF).valeur, PI, 1e-8);
    CHECK_NEAR(integrerDefinie(1.0 / (X * X), 1.0, INF).valeur, 1.0, 1e-8);
    // Orientation : de +inf à 0
    CHECK_NEAR(integrerDefinie(ast_exp(-X), INF, 0.0).valeur, -1.0, 1e-9);
    // Bornes infinies égales
    CHECK_EQ(integrerDefinie(ast_exp(-X), INF, INF).valeur, 0.0);
    // Divergentes : 1/x sur [1, inf[, constante sur R
    const ResultatIntegrale d = integrerDefinie(1.0 / X, 1.0, INF);
    CHECK(!std::isfinite(d.valeur));
    CHECK(std::isinf(d.erreurEstimee));
    CHECK(!integrerDefinie(cst(1.0), -INF, INF).estFinie());
}

TEST_CASE(singularites_jamais_de_valeur_exacte_fausse) {
    // Piège classique : la primitive -1/x donnerait -2 sur [-1, 1]
    for (const OptionsIntegration& o : {OptionsIntegration{}, numerique()}) {
        const ResultatIntegrale r = integrerDefinie(1.0 / (X * X), cst(-1.0), cst(1.0), o);
        CHECK(!r.estExacte());
        CHECK(!std::isfinite(r.valeur));
        CHECK(std::isinf(r.erreurEstimee));
    }
    CHECK(!integrerDefinie(1.0 / (X * X), -1.0, 1.0).estFinie());
    // Singularité en une borne, non intégrable
    CHECK(!integrerDefinie(1.0 / (X * X), 0.0, 1.0).estFinie());
    CHECK(!integrerDefinie(1.0 / X, 0.0, 1.0).estFinie());
    // 1/x sur [-1, 1] : valeur principale non définie ici
    CHECK(!integrerDefinie(1.0 / X, -1.0, 1.0).estFinie());
    // Pôle intérieur de 1/(x-1)
    CHECK(!integrerDefinie(1.0 / (X - 1.0), 0.0, 2.0).estFinie());
    CHECK(!integrerDefinie(1.0 / (X - 1.0), 0.0, 2.0).estExacte());
    // Mais l'exact reste exact hors de la singularité
    accord(1.0 / (X * X), cst(1.0), cst(2.0), "1/2");
    accord(1.0 / (X - 1.0), cst(2.0), cst(3.0), "ln(2)");
    // tan : pôle en pi/2
    CHECK(!integrerDefinie(ast_tan(X), 0.0, 2.0).estFinie());
    CHECK(!integrerDefinie(ast_tan(X), 0.0, 2.0).estExacte());
    CHECK(!integrerDefinie(ast_tan(X), cst(0.0), pi() / 2.0).estFinie());
    accord(ast_tan(X), cst(0.0), cst(1.0), "-ln(cos(1))");
    // ln d'un argument négatif ou nul : jamais d'exact (la primitive ln(x) n'a pas de valeur en -2)
    const ResultatIntegrale n = integrerDefinie(1.0 / X, -3.0, -2.0);
    CHECK(!n.estExacte());
    CHECK_NEAR(n.valeur, -std::log(1.5), 1e-10);
    // Puissance non entière d'une base négative : non définie sur les réels
    CHECK(!integrerDefinie(ast_pow(X, frac(1, 2)), -1.0, 1.0).estFinie());
}

TEST_CASE(singularites_integrables) {
    // 1/sqrt(x) sur [0, 1] : la primitive 2*sqrt(x) est continue en 0, résultat exact
    accord(ast_pow(X, frac(-1, 2)), cst(0.0), cst(1.0), "2", 1e-8);
    const ResultatIntegrale e = integrerDefinie(ast_pow(X, frac(-1, 2)), 0.0, 1.0);
    CHECK(e.estExacte());
    CHECK_EQ(e.valeur, 2.0);
    // sqrt(x) en 0
    accord(ast_pow(X, frac(1, 2)), cst(0.0), cst(4.0), "16/3");
    // ln(x) en 0 : intégrable, la primitive x*ln(x) - x n'est pas évaluable en 0 -> numérique
    const ResultatIntegrale l = integrerDefinie(ast_ln(X), 0.0, 1.0);
    CHECK(!l.estExacte());
    CHECK_NEAR(l.valeur, -1.0, 1e-8);
    // 1/sqrt(1-x^2) : singularités aux deux bornes, résultat pi
    const ResultatIntegrale a = integrerDefinie(ast_pow(1.0 - X * X, frac(-1, 2)), -1.0, 1.0);
    CHECK_NEAR(a.valeur, PI, 1e-6); // pas d'extrapolation : la queue en 1e-16 d'une borne non nulle reste inconnue
    CHECK(std::abs(a.valeur - PI) <= a.erreurEstimee);
    // Singularité intérieure amovible : sin(x)/x, coupure en 0 (2*Si(1))
    const ResultatIntegrale s = integrerDefinie(ast_sin(X) / X, -1.0, 1.0);
    CHECK(!s.estExacte());
    CHECK_NEAR(s.valeur, 1.8921661407343662, 1e-10);
    // 1/sqrt(|x|) : singularité intégrable à l'intérieur, coupure en 0, résultat 4
    CHECK_NEAR(integrerDefinie(ast_pow(X * X, frac(-1, 4)), -1.0, 1.0).valeur, 4.0, 1e-6);
}

TEST_CASE(parametres) {
    const ExprPtr p = param("p");
    // Résultat exact symbolique : p/3, bornes paramétrées
    const ResultatIntegrale r = integrerDefinie(p * X * X, cst(0.0), cst(1.0));
    CHECK(r.estExacte());
    CHECK_EQ(r.exacte, p / 3.0);
    CHECK(std::isnan(r.valeur)); // pas de valeur numérique
    const ResultatIntegrale s = integrerDefinie(X, cst(0.0), p);
    CHECK(s.estExacte());
    CHECK_EQ(s.exacte, p * p / 2.0);
    const ResultatIntegrale t = integrerDefinie(ast_sin(X), cst(0.0), p);
    CHECK(t.estExacte());
    CHECK_EQ(t.exacte, 1.0 - ast_cos(p));
    // Primitive inconnue et paramètre : le numérique est impossible
    CHECK_THROWS(integrerDefinie(ast_pow(X, p), cst(0.0), cst(1.0)), std::domain_error);
    CHECK_THROWS(integrerDefinie(p / X, cst(0.0), cst(1.0)), std::domain_error);
    CHECK_THROWS(integrerDefinie(1.0 / (X * X + p), 0.0, 1.0), std::domain_error);
    // Singularité possible avec des bornes paramétrées : on ne certifie rien
    CHECK_THROWS(integrerDefinie(1.0 / X, cst(1.0), p), std::domain_error);
    // Le numérique forcé lève aussi
    CHECK_THROWS(integrerDefinie(p * X, cst(0.0), cst(1.0), numerique()), std::domain_error);
}

TEST_CASE(refus_et_orientation) {
    const ExprPtr y = var("y");
    CHECK_THROWS(integrerDefinie(X * y, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS(integrerDefinie(X * y, cst(0.0), cst(1.0)), std::invalid_argument);
    // Borne contenant une variable
    CHECK_THROWS(integrerDefinie(X, cst(0.0), X), std::invalid_argument);
    CHECK_THROWS(integrerDefinie(ExprPtr(), 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS(integrerDefinie(X, 0.0, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);

    // a > b : opposé ; a == b : 0 exact
    const ResultatIntegrale direct = integrerDefinie(ast_pow(X, 2.0), cst(0.0), cst(3.0));
    const ResultatIntegrale inverse = integrerDefinie(ast_pow(X, 2.0), cst(3.0), cst(0.0));
    CHECK(direct.estExacte() && inverse.estExacte());
    CHECK_EQ(direct.exacte->texte(), std::string("9"));
    CHECK_EQ(inverse.exacte->texte(), std::string("-9"));
    CHECK_EQ(inverse.valeur, -9.0);
    const ResultatIntegrale nul = integrerDefinie(1.0 / X, cst(2.0), cst(2.0));
    CHECK(nul.estExacte());
    CHECK_EQ(nul.exacte->texte(), std::string("0"));
    CHECK_EQ(nul.valeur, 0.0);
    CHECK_EQ(nul.erreurEstimee, 0.0);
    const ResultatIntegrale nul2 = integrerDefinie(1.0 / (X * X), 1.5, 1.5);
    CHECK_EQ(nul2.valeur, 0.0);
    // Numérique inversé
    const ResultatIntegrale n = integrerDefinie(1.0 / (X * X + 1.0), 1.0, 0.0);
    CHECK_NEAR(n.valeur, -PI / 4.0, 1e-12);
}

TEST_CASE(tolerances_et_budget) {
    OptionsIntegration serre;
    serre.toleranceRelative = 1e-13;
    serre.toleranceAbsolue = 1e-14;
    CHECK_NEAR(integrerDefinie(1.0 / (X * X + 1.0), 0.0, 1.0, serre).valeur, PI / 4.0, 1e-13);

    // Budget épuisé : pas d'exception, erreur honnête
    OptionsIntegration court = numerique();
    court.subdivisionsMax = 3;
    court.toleranceRelative = 1e-13;
    const ResultatIntegrale r = integrerDefinie(ast_ln(X), 0.0, 1.0, court);
    CHECK(std::isfinite(r.valeur));
    CHECK(r.erreurEstimee > 1e-13);
    CHECK(std::abs(r.valeur + 1.0) <= r.erreurEstimee);

    // Tolérance relative lâche : moins d'intervalles, erreur sous la tolérance
    OptionsIntegration lache;
    lache.toleranceRelative = 1e-4;
    lache.toleranceAbsolue = 1e-6;
    const ResultatIntegrale l = integrerDefinie(ast_sin(X * X), 0.0, 10.0, lache);
    CHECK(l.erreurEstimee <= 1e-4 * std::abs(l.valeur) + 1e-6);
    CHECK_NEAR(l.valeur, 0.5836708999296392, 1e-3);
}

TEST_CASE(equation_classique_integrer_borne) {
    EquationClassique eq("x^2");
    const ResultatIntegrale r = eq.integrer(0.0, 3.0);
    CHECK(r.estExacte());
    CHECK_EQ(r.exacte->texte(), std::string("9"));
    const ResultatIntegrale s = eq.integrer(cst(0.0), frac(1, 2));
    CHECK(s.estExacte());
    CHECK_EQ(s.exacte->texte(), std::string("1/24"));

    EquationClassique arctan("1/(1+x^2)");
    const ResultatIntegrale a = arctan.integrer(0.0, 1.0);
    CHECK(!a.estExacte());
    CHECK_NEAR(a.valeur, PI / 4.0, 1e-12);
    CHECK_NEAR(arctan.integrer(-INF, INF).valeur, PI, 1e-8);
    OptionsIntegration o = numerique();
    CHECK_NEAR(eq.integrer(cst(0.0), cst(3.0), o).valeur, 9.0, 1e-12);

    // La primitive de l'équation intégrée en x donne la même valeur que integrer(a, b)
    const EquationClassique F = EquationClassique("cos(x)").integrer();
    CHECK_NEAR(F.eval(1.0) - F.eval(0.0), EquationClassique("cos(x)").integrer(0.0, 1.0).valeur, 1e-14);
}

TEST_CASE(coherence_exact_numerique) {
    // Une série de fonctions classiques : les deux chemins s'accordent
    const std::vector<ExprPtr> fs = {
        ast_exp(X) * X, ast_sin(X) * ast_cos(X), ast_pow(X, 5.0) - 3.0 * X, ast_ln(X) / X, ast_exp(cst(2.0) * X) + 1.0,
        ast_pow(X, frac(3, 2)), ast_atan(X)};
    for (const ExprPtr& f : fs) {
        const ResultatIntegrale e = integrerDefinie(f, 1.0, 2.0);
        const ResultatIntegrale n = integrerDefinie(f, 1.0, 2.0, numerique());
        CHECK(std::isfinite(e.valeur));
        CHECK_NEAR(e.valeur, n.valeur, 1e-10 * (1.0 + std::abs(n.valeur)));
    }
}

int main() { return test::executerTous(); }
