/**
 * @file test_simplification.cpp
 * @brief Tests de simplifier() : options (OptionsSimplification) et réduction des fractions
 *        rationnelles par PGCD.
 */

#include <cmath>
#include <string>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "Lecture.hpp"
#include "Simplification.hpp"
#include "Solveur.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

// Forme simplifiée (texte), après vérification des invariants : relecture du texte en le
// même noeud, simplification idempotente, mêmes valeurs que l'expression d'origine
std::string simplifie(const std::string& texte) {
    const ExprPtr e = lire(texte);
    const ExprPtr s = e->simplifier();
    CHECK(lire(s->texte()).get() == s.get());
    CHECK(s->simplifier().get() == s.get());
    if (!e->plusieursVariables()) {
        for (double x : {-2.7, -0.4, 0.3, 2.2, 3.9}) {
            double attendu = 0.0;
            try {
                attendu = e->eval(x);
            } catch (const std::logic_error&) {
                continue; // paramètre : pas de valeur numérique
            }
            CHECK_NEAR(s->eval(x), attendu, 1e-9 * (1.0 + std::abs(attendu)));
        }
    }
    return s->texte();
}

// L'expression est rendue telle quelle (même noeud)
bool inchangee(const std::string& texte) {
    const ExprPtr e = lire(texte, OptionsLecture{"x", {"y"}});
    return e->simplifier().get() == e.get();
}

} // namespace

TEST_CASE(rationnel_quotients_de_polynomes) {
    CHECK_EQ(simplifie("(x^2 - 1)/(x - 1)"), std::string("x + 1"));
    CHECK_EQ(simplifie("(x^3 - 1)/(x^2 - 1)"), std::string("(x^2 + x + 1)/(x + 1)"));
    CHECK_EQ(simplifie("(2*x + 2)/(x + 1)"), std::string("2"));
    CHECK_EQ(simplifie("x^3/(x^2 + x)"), std::string("x^2/(x + 1)"));
    CHECK_EQ(simplifie("(x^2 - 1)/(1 - x)"), std::string("-x - 1"));
    CHECK_EQ(simplifie("(x^2 - 1)/(2*x - 2)"), std::string("x/2 + 1/2"));
    CHECK_EQ(simplifie("(x^2 - 1/4)/(x - 1/2)"), std::string("x + 1/2"));
}

TEST_CASE(rationnel_puissances_et_multiplicites) {
    // Les facteurs ne sont pas développés : P^a / Q^b par PGCD des bases
    CHECK_EQ(simplifie("(x^2 - 1)^2/(x - 1)"), std::string("(x - 1)*(x + 1)^2"));
    CHECK_EQ(simplifie("(x - 1)/(x^2 - 2*x + 1)"), std::string("1/(x - 1)"));
    CHECK_EQ(simplifie("(x^2 - 1)^3/((x - 1)^5*(x + 1))"), std::string("(x + 1)^2/(x - 1)^2"));
}

TEST_CASE(rationnel_facteurs_non_polynomiaux_gardes) {
    CHECK_EQ(simplifie("sin(x)*(x^2 - 1)/(x - 1)"), std::string("(x + 1)*sin(x)"));
    CHECK_EQ(simplifie("3*exp(x)*(x^2 - 4)/(x + 2)"), std::string("3*(x - 2)*exp(x)"));
    // Dans les sous-expressions : sommes, arguments de fonctions
    CHECK_EQ(simplifie("1 + (x^2 - 1)/(x - 1)"), std::string("x + 2"));
    CHECK_EQ(simplifie("exp((x^2 - 1)/(x - 1))"), std::string("exp(x + 1)"));
}

TEST_CASE(rationnel_cas_laisses_intacts) {
    CHECK(inchangee("(x^2 + 1)/(x + 1)"));         // premiers entre eux
    CHECK(inchangee("x/(x + 1)"));                  // rejet sans allocation
    CHECK(inchangee("sin(x)/x"));
    CHECK(inchangee("(x^2 - a^2)/(x - a)"));        // paramètre
    CHECK(inchangee("(x^2 - y^2)/(x - y)"));        // plusieurs variables
    CHECK(inchangee("y*(x^2 - 1)/(x - 1)"));
    CHECK(inchangee("(x^2 - 1)/(x - 1)^(1/2)"));    // exposant non entier
    CHECK(inchangee("x^2 + 3*x + 2"));
    // Coefficients réels (non exacts ; lire("0.5") serait le rationnel exact 1/2) : sans exception
    const ExprPtr x = var("x");
    const ExprPtr reel = (cst(0.5) * ast_pow(x, 2.0) - cst(0.5)) / (x - 1.0);
    CHECK(reel->simplifier().get() == reel.get());
    // Coefficient réel d'un produit de polynômes exacts : la réduction s'applique
    const ExprPtr coefficientReel = cst(0.5) * ast_sin(x) * (ast_pow(x, 2.0) - 1.0) / (x - 1.0);
    CHECK_EQ(coefficientReel->simplifier()->texte(), std::string("0.5*(x + 1)*sin(x)"));
}

TEST_CASE(options_desactivent_la_reduction) {
    const ExprPtr e = lire("(x^2 - 1)/(x - 1)");
    CHECK(e->simplifier(OptionsSimplification::aucune()).get() == e.get());
    OptionsSimplification sansRationnel;
    sansRationnel.rationnel = false;
    CHECK(e->simplifier(sansRationnel).get() == e.get());
    CHECK_EQ(e->simplifier(OptionsSimplification())->texte(), std::string("x + 1"));
}

TEST_CASE(options_ne_polluent_pas_le_memo) {
    // Autres options d'abord : la simplification par défaut est ensuite calculée
    const ExprPtr e = lire("(x^3 - x)/(x^2 + x)");
    CHECK(e->simplifier(OptionsSimplification::aucune()).get() == e.get());
    CHECK_EQ(e->simplifier()->texte(), std::string("x - 1"));
    // Défaut d'abord (mémorisé) : aucune() rend toujours l'expression d'origine
    CHECK(e->simplifier(OptionsSimplification::aucune()).get() == e.get());
    CHECK_EQ(e->simplifier()->texte(), std::string("x - 1"));
}

TEST_CASE(equation_classique_simplifier_avec_options) {
    EquationClassique eq("(x^2 - 1)/(x - 1)");
    eq.simplifier(OptionsSimplification::aucune());
    CHECK_EQ(eq.getExpression()->texte(), std::string("(x^2 - 1)/(x - 1)"));
    eq.simplifier();
    CHECK_EQ(eq.getExpression()->texte(), std::string("x + 1"));
    CHECK_NEAR(eq.eval(1.0), 2.0, 1e-12); // singularité levée (domaine élargi, comme cancel)
}

TEST_CASE(resolution_garde_le_domaine) {
    // Le solveur ne simplifie pas l'équation : x = 1 (hors domaine) n'est pas une solution
    const Solutions s = resoudre(lire("(x^2 - 1)/(x - 1)"));
    CHECK_EQ(s.liste.size(), std::size_t(1));
    if (s.liste.size() == 1) CHECK_EQ(s.liste[0].valeur->texte(), std::string("-1"));
    const Solutions t = EquationClassique("(x^2 - 1)/(x - 1) = 0").resoudre();
    CHECK_EQ(t.liste.size(), std::size_t(1));
}

TEST_CASE(derivee_et_primitive_simplifiees) {
    // derivee() simplifie : chaque quotient de la règle du quotient est réduit
    const EquationClassique d = EquationClassique("(x^2 - 1)/(x - 1)").derivee();
    for (double x : {-1.5, 0.5, 2.5}) CHECK_NEAR(d.eval(x), 1.0, 1e-12);
    const ExprPtr f = lire("x^2/(x^3 + x^2)");
    CHECK_EQ(f->simplifier()->texte(), std::string("1/(x + 1)"));
    EquationClassique eq(f); // l'intégration d'un quotient demande sa forme réduite
    eq.simplifier();
    CHECK_EQ(eq.integrer().getExpression()->texte(), std::string("ln(x + 1)"));
}

int main() { return test::executerTous(); }
