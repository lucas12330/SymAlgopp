/**
 * @file test_multivariable.cpp
 * @brief Tests des fonctions de plusieurs variables : variables d'une expression, dérivées
 *        partielles, gradient, jacobienne, hessienne, laplacien, divergence, évaluation
 *        numérique et refus des opérations à une variable.
 */

#include <cmath>
#include <string>
#include <vector>

#include <symalgopp>

#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");
const ExprPtr Y = var("y");
const ExprPtr Z = var("z");

ExprPtr l(const std::string& texte) { return lire(texte, {"x", {"y", "z"}}); }

// Dérivée partielle par différences finies centrées (référence indépendante des règles)
double diffFinies(const ExprPtr& f, const std::string& nom, Valeurs v) {
    const double h = 1e-6;
    v[nom] += h;
    const double haut = evaluer(f, v);
    v[nom] -= 2 * h;
    return (haut - evaluer(f, v)) / (2 * h);
}

} // namespace

TEST_CASE(variables_de_l_expression) {
    const ExprPtr f = l("y*x + sin(z) + a*x");
    const std::vector<ExprPtr> v = variables(f);
    CHECK_EQ(v.size(), std::size_t(3));
    CHECK(v[0] == X);
    CHECK(v[1] == Y);
    CHECK(v[2] == Z);
    CHECK(variables(l("3 + a")).empty());
    CHECK(variables(l("x^2 + x*sin(x)")).size() == 1);
    CHECK(f->plusieursVariables());
    CHECK(!l("x^2")->plusieursVariables());
    CHECK(!l("a")->contientVariable());
}

TEST_CASE(derivee_partielle_exacte) {
    const ExprPtr f = l("x^2*y + 3*y^2");
    CHECK(f->derivee("x")->simplifier() == l("2*x*y"));
    CHECK(f->derivee(Y)->simplifier() == l("x^2 + 6*y"));
    CHECK(f->derivee("z")->simplifier() == nombre(Nombre(0))); // ne dépend pas de z
    CHECK(l("x*y")->derivee("y") == X);
    CHECK(l("x^y")->derivee("y")->simplifier() == l("x^y*ln(x)"));
    CHECK(l("x^y")->derivee("x")->simplifier() == l("y*x^(y - 1)"));
    CHECK(l("sin(x*y)")->derivee("x")->simplifier() == l("y*cos(x*y)"));
    CHECK(l("a*x*y")->derivee("y") == l("a*x")); // un paramètre reste une constante
}

TEST_CASE(derivee_sans_ambiguite_a_une_variable) {
    CHECK(l("x^3")->derivee() == l("3*x^2"));
    // le nom de la variable n'a pas d'importance : t^3 est une fonction de t
    const ExprPtr t = lire("t^3", {"t"});
    CHECK(t->derivee() == lire("3*t^2", {"t"}));
    CHECK(t->derivee("t") == t->derivee());
    CHECK_THROWS(l("x*y")->derivee(), std::invalid_argument);
    CHECK_THROWS(l("x*y")->derivee(cst(2.0)), std::invalid_argument);
    // dériver par rapport à une variable absente donne 0
    CHECK(t->derivee("x") == nombre(Nombre(0)));
}

TEST_CASE(derivees_partielles_comparees_aux_differences_finies) {
    const std::vector<std::string> textes = {
        "x^2*y + sin(x*y)",       "exp(x*y)/(1 + z^2)",       "atan(x/y) + ln(x^2 + y^2 + z^2)",
        "x^y + y^z + z^x",        "(x + y*z)^3*cos(x - z)",   "sin(x)*cos(y)*tan(z/4)",
        "x/(y^2 + 1) - y/(x^2 + 2)",
    };
    const Valeurs point = {{"x", 0.7}, {"y", 1.3}, {"z", 0.4}};
    for (const std::string& texte : textes) {
        const ExprPtr f = l(texte);
        for (const char* nom : {"x", "y", "z"}) {
            const double exacte = evaluer(f->derivee(nom), point);
            CHECK_NEAR(exacte, diffFinies(f, nom, point), 1e-6 * (1.0 + std::abs(exacte)));
        }
    }
}

TEST_CASE(sous_expressions_partagees_derivees_une_fois) {
    // u apparaît sous deux formes : le cache de dérivation ne doit pas mélanger les variables
    const ExprPtr u = l("sin(x*y)");
    const ExprPtr f = u * u + ast_exp(u) * X;
    const Valeurs point = {{"x", 0.9}, {"y", 0.5}};
    for (const char* nom : {"x", "y"}) {
        CHECK_NEAR(evaluer(f->derivee(nom), point), diffFinies(f, nom, point), 1e-6);
    }
    // deux dérivées successives avec des cibles différentes (chaque appel a son cache)
    CHECK_NEAR(evaluer(f->derivee("x")->derivee("y"), point), evaluer(f->derivee("y")->derivee("x"), point), 1e-12);
}

TEST_CASE(derivees_croisees_symetriques) {
    const ExprPtr f = l("exp(x*y) + x^3*y^2 + sin(y*z)*x");
    CHECK(deriveeMixte(f, {X, Y}) == deriveeMixte(f, {Y, X}));
    CHECK(deriveeMixte(f, {X, Y, Z}) == deriveeMixte(f, {Z, Y, X}));
    CHECK(deriveeMixte(l("x^3*y^2"), {X, X, Y}) == l("12*x*y"));
    CHECK(deriveeMixte(f, {}) == f->simplifier());
}

TEST_CASE(gradient_et_hessienne) {
    const ExprPtr f = l("x^2*y + y^3");
    const std::vector<ExprPtr> g = gradient(f);
    CHECK_EQ(g.size(), std::size_t(2));
    CHECK(g[0] == l("2*x*y"));
    CHECK(g[1] == l("x^2 + 3*y^2"));

    const auto H = hessienne(f);
    CHECK(H[0][0] == l("2*y"));
    CHECK(H[0][1] == l("2*x"));
    CHECK(H[1][0] == H[0][1]); // même noeud : symétrique
    CHECK(H[1][1] == l("6*y"));

    // ordre imposé par l'appelant, variable absente de f
    const auto ordre = gradient(f, {Y, X, Z});
    CHECK(ordre[0] == g[1]);
    CHECK(ordre[1] == g[0]);
    CHECK(ordre[2] == nombre(Nombre(0)));
}

TEST_CASE(jacobienne_laplacien_divergence) {
    const auto J = jacobienne({l("x*y"), l("x + sin(y)")}, {X, Y});
    CHECK(J[0][0] == Y);
    CHECK(J[0][1] == X);
    CHECK(J[1][0] == cst(1.0));
    CHECK(J[1][1] == l("cos(y)"));

    CHECK(laplacien(l("x^2 + y^2 + z^2")) == cst(6.0));
    CHECK(laplacien(l("x^3 - 3*x*y^2")) == nombre(Nombre(0))); // fonction harmonique
    CHECK(divergence({l("x*y"), l("y*z"), l("z*x")}, {X, Y, Z}) == l("x + y + z"));
    CHECK_THROWS(divergence({X}, {X, Y}), std::invalid_argument);
}

TEST_CASE(evaluation_numerique) {
    const ExprPtr f = l("x^2*y + z");
    CHECK_NEAR(evaluer(f, {{"x", 3.0}, {"y", 2.0}, {"z", 1.0}}), 19.0, 1e-12);
    CHECK_NEAR(evaluer(f, {{"x", 3.0}, {"y", 2.0}, {"z", 1.0}, {"inutile", 9.0}}), 19.0, 1e-12);
    CHECK_THROWS(evaluer(f, {{"x", 3.0}, {"y", 2.0}}), std::invalid_argument);
    CHECK_THROWS(evaluer(l("a*x"), {{"x", 1.0}}), std::logic_error); // paramètre sans valeur
    const std::vector<double> g = evaluer(gradient(f), {{"x", 3.0}, {"y", 2.0}, {"z", 1.0}});
    CHECK_NEAR(g[0], 12.0, 1e-12);
    CHECK_NEAR(g[1], 9.0, 1e-12);
    CHECK_NEAR(g[2], 1.0, 1e-12);
}

TEST_CASE(programme_a_plusieurs_sorties_partage_les_calculs) {
    const ExprPtr f = l("exp(x*y)*sin(x + y)");
    const std::vector<ExprPtr> sorties = {f, f->derivee("x"), f->derivee("y"), deriveeMixte(f, {X, Y})};
    const ProgrammeEvaluation ensemble(sorties, {X, Y});
    CHECK(ensemble.estValide());
    CHECK_EQ(ensemble.nombreEntrees(), std::size_t(2));
    CHECK_EQ(ensemble.nombreSorties(), std::size_t(4));

    std::size_t separees = 0;
    for (const ExprPtr& s : sorties) separees += ProgrammeEvaluation({s}, {X, Y}).nombreInstructions();
    CHECK(ensemble.nombreInstructions() < separees);

    const std::vector<double> v = ensemble.evaluerEn({0.6, 1.1});
    for (std::size_t i = 0; i < sorties.size(); ++i) {
        CHECK_NEAR(v[i], evaluer(sorties[i], {{"x", 0.6}, {"y", 1.1}}), 1e-12);
    }
    CHECK_THROWS(ensemble.evaluerEn({1.0}), std::invalid_argument);
    CHECK_THROWS(ensemble.evaluer(1.0), std::logic_error);

    // l'ordre des entrées est celui de la liste
    CHECK_NEAR(ProgrammeEvaluation({l("x - y")}, {Y, X}).evaluerEn({1.0, 5.0})[0], 4.0, 1e-15);
    // variable sans entrée, entrée qui n'est pas une variable, entrée en double
    CHECK(!ProgrammeEvaluation({l("x*y")}, {X}).estValide());
    CHECK_THROWS(ProgrammeEvaluation({X}, {cst(1.0)}), std::invalid_argument);
    CHECK_THROWS(ProgrammeEvaluation({X}, {X, X}), std::invalid_argument);
}

TEST_CASE(lecture_a_plusieurs_variables) {
    const ExprPtr f = lire("x*y + a", {"x", {"y"}});
    CHECK(f == X * Y + param("a"));
    // sans la déclaration, y est un paramètre
    CHECK(lire("x*y") == X * param("y"));
    // aller-retour de l'affichage
    for (const char* texte : {"x^2*y - z/x", "sin(x*y)^z", "1/(x + y*z)", "-x*y + 2*z^3"}) {
        const ExprPtr e = l(texte);
        CHECK(l(e->texte()) == e);
    }
}

TEST_CASE(operations_a_une_variable_refusees) {
    const ExprPtr f = l("x*y");
    CHECK(f->integrer()->type() == TypeNoeud::IntegraleNonEvaluee);
    CHECK(f->limite(0.0)->type() == TypeNoeud::LimiteNonEvaluee);
    CHECK_THROWS(f->DL(0.0, 3), std::domain_error);
    CHECK_THROWS(resoudre(f, cst(1.0)), std::invalid_argument);
    CHECK_THROWS(resoudreNumerique(f, 0.0, 1.0), std::invalid_argument);
    CHECK_THROWS(resoudreSurIntervalle(f, 0.0, 1.0), std::invalid_argument);
    CHECK(!ProgrammeEvaluation(f).estValide());
    CHECK_THROWS(ProgrammeEvaluation(f).evaluer(1.0), std::logic_error);
    // avec une constante à la place de y, la résolution fonctionne comme avant
    CHECK(resoudre(lire("2*x"), lire("1")).liste.size() == 1);
}

TEST_CASE(equation_classique_plusieurs_variables) {
    const EquationClassique eq("x^2*y = 1", {"x", {"y"}});
    CHECK_EQ(eq.variables().size(), std::size_t(2));
    CHECK_NEAR(eq.eval({{"x", 2.0}, {"y", 0.5}}), 1.0, 1e-12);
    CHECK_NEAR(eq.eval({{"x", 2.0}, {"y", 0.5}}), 1.0, 1e-12); // programme réutilisé
    CHECK_THROWS(eq.eval({{"x", 2.0}}), std::invalid_argument);
    CHECK_THROWS(eq.eval(2.0), std::logic_error);
    CHECK_THROWS(eq.derivee(), std::invalid_argument);
    CHECK(eq.derivee("y").getExpression() == l("x^2"));
    const std::vector<EquationClassique> g = eq.gradient();
    CHECK_EQ(g.size(), std::size_t(2));
    CHECK(g[0].getExpression() == l("2*x*y"));
    CHECK(g[1].getExpression() == l("x^2"));
}

TEST_CASE(une_variable_inchangee) {
    // la dérivation classique, l'intégration et les limites ne changent pas
    const ExprPtr f = lire("x^2*sin(x)");
    CHECK(f->derivee() == f->derivee("x"));
    CHECK(lire("sin(x)/x")->limite(0.0) == cst(1.0));
    CHECK(lire("2*x")->integrer()->simplifier() == lire("x^2"));
    const EquationClassique eq("x^2 = 2");
    CHECK_NEAR(eq.eval(3.0), 7.0, 1e-12);
    CHECK_EQ(eq.resoudre().liste.size(), std::size_t(2));
}

int main() { return test::executerTous(); }
