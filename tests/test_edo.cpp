/**
 * @file test_edo.cpp
 * @brief Tests de la résolution symbolique des équations différentielles non linéaires :
 *        classes du premier ordre, second ordre réductible, solutions singulières,
 *        équations hors de portée.
 *
 * Chaque solution est contrôlée deux fois : par verifierEDO (réinjection symbolique ou
 * numérique), et, pour les formes explicites, par le résidu de l'équation calculé ici en des
 * points choisis du domaine, avec des valeurs numériques des constantes.
 */

#include <cmath>
#include <string>
#include <vector>

#include "EDO.hpp"
#include "Journal.hpp"
#include "Lecture.hpp"
#include "Multivariable.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const InconnuesEDO V; // x, y, dy, d2y

// « y'' » et « y' » remplacés par les variables d2y et dy, puis gauche - droite
ExprPtr equationDe(std::string texte) {
    for (const auto& [motif, nom] : {std::pair<std::string, std::string>{"y''", "d2y"}, {"y'", "dy"}}) {
        for (std::size_t i = texte.find(motif); i != std::string::npos; i = texte.find(motif, i)) texte.replace(i, motif.size(), nom);
    }
    const EgaliteLue e = lireEquation(texte, OptionsLecture{"x", {"y", "dy", "d2y"}});
    return e.gauche - e.droite;
}

// « methode : solution | solution », les singulières entre crochets
std::string resume(const SolutionsEDO& s) {
    std::string t = s.resolue() ? s.methode + " : " : "non résolue";
    for (std::size_t i = 0; i < s.liste.size(); ++i) {
        const std::string sol = s.liste[i].texte();
        t += (i ? " | " : "") + (s.liste[i].singuliere ? "[" + sol + "]" : sol);
    }
    return t;
}

// Résidu de l'équation pour y = Y(x), en x0, les constantes C1, C2 valant c1, c2
double residuExplicite(const ExprPtr& equation, const ExprPtr& Y, double x0, double c1, double c2) {
    ExprPtr y = substituer(substituer(Y, param("C1"), cst(c1)), param("C2"), cst(c2));
    const ExprPtr dy = y->derivee(V.x), d2y = dy->derivee(V.x);
    ExprPtr r = substituer(equation, V.d2y, d2y);
    r = substituer(substituer(r, V.dy, dy), V.y, y);
    return evaluer(r, {{"x", x0}});
}

struct Cas {
    const char* equation;
    const char* attendu; // résumé attendu (vide : non comparé)
    double x0, c1, c2;   // point de contrôle des solutions explicites
};

void controler(const Cas& c) {
    const SolutionsEDO s = resoudreEDO(c.equation);
    CHECK(s.resolue());
    if (!s.resolue()) {
        std::cerr << "      non résolue : " << c.equation << "\n";
        return;
    }
    if (c.attendu[0] != '\0') CHECK_EQ(resume(s), std::string(c.attendu));
    const ExprPtr equation = equationDe(c.equation);
    for (const SolutionEDO& sol : s.liste) {
        CHECK(verifierEDO(equation, V, sol));
        if (!sol.y) continue;
        const double r = residuExplicite(equation, sol.y, c.x0, c.c1, c.c2);
        if (!std::isfinite(r)) continue; // branche non définie en ce point
        CHECK_NEAR(r, 0.0, 1e-8);
        if (std::fabs(r) > 1e-8) std::cerr << "      " << c.equation << " : " << sol.texte() << ", résidu " << r << "\n";
    }
}

} // namespace

// ============================================================================
// Premier ordre
// ============================================================================

TEST_CASE(quadrature_et_lineaire) {
    controler({"y' = cos(x)", "quadrature : y = C1 + sin(x)", 0.7, 1.3, 0});
    controler({"y' = x*y", "lineaire : y = C1*exp(x^2/2)", 0.7, 1.3, 0});
    controler({"y' + y/x = x^2", "lineaire : y = x^3/4 + C1/x", 0.7, 1.3, 0});
    controler({"y' - 2*x*y = x", "lineaire : y = C1*exp(x^2) - 1/2", 0.7, 1.3, 0});
    controler({"y' + y*tan(x) = cos(x)", "lineaire : y = C1*cos(x) + x*cos(x)", 0.7, 1.3, 0});
    controler({"y' = y + x", "lineaire : y = -x + C1*exp(x) - 1", 0.7, 1.3, 0});
    // Somme canonique qui est (x + 1)(y + 1) : linéaire, la constante absorbe y = -1
    controler({"y' = x*y + x + y + 1", "lineaire : y = C1*exp(x^2/2 + x) - 1", 0.7, 1.3, 0});
    // Coefficient de y' variable
    controler({"x*y' + 2*y = sin(x)", "", 0.7, 1.3, 0});
}

TEST_CASE(bernoulli) {
    controler({"y' = x*y^2", "bernoulli : y = 1/(-x^2/2 + C1) | [y = 0]", 0.7, 1.3, 0});
    controler({"y' + y/x = x*y^2", "bernoulli : y = 1/(-x^2 + C1*x) | [y = 0]", 0.7, 1.3, 0});
    controler({"x*y' + y = y^2*ln(x)", "bernoulli : y = 1/(C1*x + ln(x) + 1) | [y = 0]", 0.7, 1.3, 0});
    controler({"y' = y*(1 - y)", "bernoulli : y = 1/(C1*exp(-x) + 1) | [y = 0]", 0.7, 1.3, 0});
    // Deux branches pour un exposant impair
    controler({"y' = y - y^3", "bernoulli : y = 1/(C1*exp(-2*x) + 1)^(1/2) | y = -1/(C1*exp(-2*x) + 1)^(1/2) | [y = 0]", 0.7, 1.3, 0});
    // Exposant fractionnaire : solution valable là où C exp(x/2) >= 1 ; non-unicité en y = 0
    controler({"y' = y + y^(1/2)", "", 0.7, 2.6, 0});
    controler({"y' = 3*y^(2/3)", "bernoulli : y = (x + C1)^3 | [y = 0]", 0.7, 1.3, 0});
    // Exposant négatif : y' = y/x + x/y, homogène aussi
    controler({"y' = y/x + x/y", "", 0.7, 1.3, 0});
    controler({"(x^2 + y^2) - 2*x*y*y' = 0", "bernoulli : y = (x^2 + C1*x)^(1/2) | y = -(x^2 + C1*x)^(1/2)", 0.7, 1.3, 0});
}

TEST_CASE(variables_separables) {
    controler({"y' = (1 + y^2)*exp(x)", "separable : y = tan(C1 + exp(x))", 0.3, 0.2, 0});
    controler({"y' = y^2 + 1", "separable : y = tan(x + C1)", 0.3, 0.2, 0});
    // Séparable sans produit apparent : exp(x + y) = exp(x) exp(y)
    controler({"y' = exp(x + y)", "separable : y = -ln(-C1 - exp(x))", 0.3, -5.0, 0});
    controler({"y' = exp(y)*sin(x)", "separable : y = -ln(-C1 + cos(x))", 0.3, -5.0, 0});
    // Constante absorbée par le logarithme, et équilibres perdus par la division
    controler({"x*y' = y*ln(y)", "separable : y = exp(C1*x) | [y = 1]", 0.7, 1.3, 0});
    controler({"y' = (y - 1)*(y - 2)", "separable : -x - C1 + ln(abs((y - 2)/(y - 1))) = 0 | [y = 1] | [y = 2]", 0.7, 1.3, 0});
}

TEST_CASE(fonction_d_une_combinaison_affine) {
    controler({"y' = (x + y)^2", "affine : y = -x + tan(x + C1)", 0.3, 0.2, 0});
    controler({"y' = sin(x + y)", "", 0.3, 0.2, 0});
    controler({"y' = (y - x)^2 + 1", "affine : y = x + 1/(-x - C1) | [y = x]", 0.3, 0.2, 0});
    controler({"y' = 1/(x + y)", "affine : y - C1 - ln(abs(x + y + 1)) = 0 | [y = -x - 1]", 0.3, 0.2, 0});
    controler({"y' = (x + y + 1)/(x + y - 1)", "", 0.3, 0.2, 0});
}

TEST_CASE(homogenes) {
    controler({"y' = (y^2 + x*y)/x^2", "", 0.7, 1.3, 0});
    controler({"y' = y/x + tan(y/x)", "homogene : y = x*asin(C1*x)", 0.3, 0.7, 0});
    controler({"y' = (x + y)/(x - y)", "homogene : -C1 - ln(y^2/x^2 + 1)/2 - ln(abs(x)) + atan(y/x) = 0", 0.7, 1.3, 0});
    // Réductible à une homogène par translation de l'origine en (1, 1)
    controler({"y' = (x + y - 2)/(x - y)",
               "homogene : -C1 - ln((y - 1)^2/(x - 1)^2 + 1)/2 - ln(abs(x - 1)) + atan((y - 1)/(x - 1)) = 0", 0.7, 1.3, 0});
}

TEST_CASE(exactes_et_facteurs_integrants) {
    controler({"(2*x*y + 3) + (x^2 + 4*y)*y' = 0", "exacte : x^2*y + 2*y^2 + 3*x - C1 = 0", 0.7, 1.3, 0});
    controler({"(3*x*y + y^2) + (x^2 + x*y)*y' = 0", "facteur_integrant : x^2*y^2/2 + x^3*y - C1 = 0", 0.7, 1.3, 0}); // mu = x
    controler({"y + (2*x - y*exp(y))*y' = 0", "facteur_integrant : x*y^2 - C1 - (y^2 - 2*y + 2)*exp(y) = 0", 0.7, 1.3, 0}); // mu = y
    controler({"y*(1 + x*y) + x*(1 - x*y)*y' = 0", "", 0.7, 1.3, 0}); // mu(x y)
    controler({"(x^2 + y^2 + x) + x*y*y' = 0", "", 0.7, 1.3, 0});
    // Exacte dont y s'isole
    controler({"(y + 2*x) + x*y' = 0", "", 0.7, 1.3, 0});
}

TEST_CASE(riccati_et_clairaut) {
    controler({"y' = y^2 - 2/x^2", "riccati : y = 1/x + 1/(-x/3 + C1/x^2) | [y = 1/x]", 0.7, 1.3, 0});
    controler({"y' = 1 + x^2 - 2*x*y + y^2", "", 0.3, 0.2, 0});
    // Clairaut : droites et enveloppe
    controler({"y = x*y' + y'^2", "clairaut : y = C1*x + C1^2 | [y = -x^2/4]", 0.7, 1.3, 0});
    controler({"y = x*y' - ln(y')", "clairaut : y = C1*x - ln(C1) | [y = -ln(1/x) + 1]", 0.7, 1.3, 0});
    controler({"y = x*y' + (1 + y'^2)^(1/2)", "clairaut : y = C1*x + (C1^2 + 1)^(1/2)", 0.7, 1.3, 0});
}

// ============================================================================
// Second ordre réductible
// ============================================================================

TEST_CASE(second_ordre_sans_y) {
    controler({"y'' = x", "ordre2_sans_y : y = x^3/6 + C1*x + C2", 0.7, 1.3, 0.4});
    controler({"y'' = y'/x", "ordre2_sans_y : y = C1*x^2/2 + C2", 0.7, 1.3, 0.4});
    controler({"x*y'' + y' = x", "ordre2_sans_y : y = x^2/4 + C2 + C1*ln(abs(x))", 0.7, 1.3, 0.4});
    controler({"y'' + y'^2 = 0", "ordre2_sans_y : y = C2 + ln(abs(x + C1)) | [y = C2]", 0.7, 1.3, 0.4});
    controler({"y'' = 1 + y'^2", "", 0.3, 0.2, 0.4});
}

TEST_CASE(second_ordre_sans_x) {
    controler({"y'' = y'^2/y", "ordre2_sans_x : y = exp(C1*(x + C2)) | y = -exp(C1*(x + C2))", 0.7, 1.3, 0.4});
    // p dp/dy = -y : p = +-sqrt(C1 - y^2), primitive à paramètre (asin)
    controler({"y'' = -y", "ordre2_sans_x : y = (4*C1)^(1/2)*sin(x + C2)/2 | y = (4*C1)^(1/2)*sin(-x - C2)/2", 0.7, 1.3, 0.4});
    // Solution implicite : int dy/sqrt(C1 y^2 - 1) = x + C2
    const SolutionsEDO s = resoudreEDO("y''*y = y'^2 + 1");
    CHECK(s.resolue());
    CHECK_EQ(s.methode, std::string("ordre2_sans_x"));
    for (const SolutionEDO& sol : s.liste) {
        CHECK(!sol.estExplicite());
        CHECK(verifierEDO(equationDe("y''*y = y'^2 + 1"), V, sol));
    }
}

// ============================================================================
// Hors de portée, interface
// ============================================================================

TEST_CASE(equations_non_resolues) {
    // Pas de solution élémentaire (Riccati y' = x^2 + y^2, erfi), forme paramétrée (Lagrange),
    // primitive dont la forme dépend du signe d'une constante : résultat vide, jamais faux
    for (const char* texte : {"y' = x^2 + y^2", "y' = y^2 - x*y + 1", "y = 2*x*y' + y'^2", "y'' = 2*y*y'", "y'' = exp(y)",
                              "y' = sin(x*y)", "y'' = x*y", "y'^3 + y^3 = x"}) {
        const SolutionsEDO s = resoudreEDO(texte);
        CHECK(!s.resolue());
        CHECK(s.methode.empty());
        if (s.resolue()) std::cerr << "      " << texte << " : " << resume(s) << "\n";
    }
    // Pas une équation différentielle
    CHECK(!resoudreEDO("y = x").resolue());
    CHECK_THROWS(resoudreEDO("y' = "), ErreurLecture);
}

TEST_CASE(interface) {
    // Expressions : les dérivées sont des variables
    const SolutionsEDO s = resoudreEDO(V.dy - V.x * V.y, V);
    CHECK_EQ(resume(s), std::string("lineaire : y = C1*exp(x^2/2)"));
    CHECK(s.liste[0].estExplicite());
    // Autres noms de variable et de fonction
    const SolutionsEDO t = resoudreEDO("u' = t*u^2", "t", "u");
    CHECK_EQ(resume(t), std::string("bernoulli : y = 1/(-t^2/2 + C1) | [y = 0]"));
    // verifierEDO refuse une fonction qui n'est pas solution
    const ExprPtr equation = equationDe("y' = x*y");
    CHECK(verifierEDO(equation, V, {lire("3*exp(x^2/2)"), ExprPtr(), false}));
    CHECK(!verifierEDO(equation, V, {lire("exp(x^2)"), ExprPtr(), false}));
    CHECK(!verifierEDO(equation, V, {ExprPtr(), ExprPtr(), false}));
    // Forme implicite : x^2 + y^2 = C pour y' = -x/y
    CHECK(verifierEDO(equationDe("y' = -x/y"), V, {ExprPtr(), lire("x^2 + y^2 - C1", OptionsLecture{"x", {"y"}}), false}));
    CHECK(!verifierEDO(equationDe("y' = -x/y"), V, {ExprPtr(), lire("x^2 - y^2 - C1", OptionsLecture{"x", {"y"}}), false}));
}

// ============================================================================
// Conditions initiales
// ============================================================================

namespace {

std::string cauchy(const char* texte, const char* x0, const char* y0, const char* dy0 = nullptr) {
    return resume(resoudreCauchy(texte, lire(x0), lire(y0), dy0 ? lire(dy0) : ExprPtr()));
}

} // namespace

TEST_CASE(probleme_de_cauchy_premier_ordre) {
    CHECK_EQ(cauchy("y' = x*y", "0", "2"), std::string("lineaire : y = 2*exp(x^2/2)"));
    CHECK_EQ(cauchy("y' + y/x = x^2", "1", "1"), std::string("lineaire : y = x^3/4 + 3/(4*x)"));
    CHECK_EQ(cauchy("y' = x*y^2", "0", "1"), std::string("bernoulli : y = 1/(-x^2/2 + 1)"));
    CHECK_EQ(cauchy("y' = (1 + y^2)*exp(x)", "0", "1"), std::string("separable : y = tan(pi/4 + exp(x) - 1)"));
    // La famille générale ne passe pas par (0, 0) (il faudrait C infini) : seule la solution nulle
    CHECK_EQ(cauchy("y' = x*y^2", "0", "0"), std::string("bernoulli : [y = 0]"));
    // La condition choisit la branche
    CHECK_EQ(cauchy("y' = y - y^3", "0", "-1/2"), std::string("bernoulli : y = -1/(3*exp(-2*x) + 1)^(1/2)"));
    CHECK_EQ(cauchy("y' = y/x + x/y", "1", "-2"), std::string("bernoulli : y = -(4*x^2 + 2*x^2*ln(abs(x)))^(1/2)"));
    // Perte d'unicité : deux solutions par (0, 0) pour y' = 3 y^(2/3), deux droites de Clairaut
    CHECK_EQ(cauchy("y' = 3*y^(2/3)", "0", "0"), std::string("bernoulli : y = x^3 | [y = 0]"));
    CHECK_EQ(cauchy("y = x*y' + y'^2", "0", "1"), std::string("clairaut : y = x + 1 | y = -x + 1"));
    // Forme implicite : la constante vaut Phi(x0, y0)
    CHECK_EQ(cauchy("(2*x*y + 3) + (x^2 + 4*y)*y' = 0", "1", "1"), std::string("exacte : x^2*y + 2*y^2 + 3*x - 6 = 0"));
    // Riccati : la solution particulière passe par (1, 1), pas la famille
    CHECK_EQ(cauchy("y' = y^2 - 2/x^2", "1", "1"), std::string("riccati : [y = 1/x]"));
}

TEST_CASE(probleme_de_cauchy_second_ordre) {
    CHECK_EQ(cauchy("y'' = x", "0", "1", "2"), std::string("ordre2_sans_y : y = x^3/6 + 2*x + 1"));
    CHECK_EQ(cauchy("y'' = y'/x", "1", "0", "3"), std::string("ordre2_sans_y : y = 3*x^2/2 - 3/2"));
    CHECK_EQ(cauchy("y'' + y'^2 = 0", "0", "0", "1"), std::string("ordre2_sans_y : y = ln(abs(x + 1))"));
    CHECK_EQ(cauchy("y'' = -y", "0", "0", "1"), std::string("ordre2_sans_x : y = sin(x)"));
    CHECK_EQ(cauchy("y'' = -y", "0", "1", "0"), std::string("ordre2_sans_x : y = sin(x + pi/2)")); // cos(x), écrit une seule fois
    CHECK_EQ(cauchy("y'' = y'^2/y", "0", "1", "2"), std::string("ordre2_sans_x : y = exp(2*x)"));
    // Second ordre sans y'(x0) : pas de solution unique, rien n'est renvoyé
    CHECK(!resoudreCauchy("y'' = x", cst(0.0), cst(1.0)).resolue());
    // Chaque solution particulière vérifie l'équation
    const SolutionsEDO s = resoudreCauchy("y'' = -y", cst(0.0), cst(0.0), cst(1.0));
    CHECK(s.resolue() && verifierEDO(equationDe("y'' = -y"), V, s.liste[0]));
}

TEST_CASE(journal_des_etapes) {
    Journal j;
    const SolutionsEDO s = resoudreEDO("y' + y/x = x*y^2", j);
    CHECK(s.resolue());
    CHECK(journal::actif() == nullptr);
    std::vector<std::string> regles;
    for (const Etape& e : j.etapes()) regles.push_back(e.regle);
    CHECK(!regles.empty() && regles.front() == "edo.forme");
    bool bernoulli = false;
    for (const std::string& r : regles) bernoulli = bernoulli || r == "edo.bernoulli";
    CHECK(bernoulli);
    // Sans journal : rien n'est noté
    resoudreEDO("y' + y/x = x*y^2");
    CHECK_EQ(j.etapes().size(), regles.size());
}

int main() { return test::executerTous(); }
