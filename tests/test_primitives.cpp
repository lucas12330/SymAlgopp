/**
 * @file test_primitives.cpp
 * @brief Tests de l'intégration symbolique : table étendue, ln|u|, fractions rationnelles,
 *        changement de variable « la dérivée divise », classe exponentielle-polynôme,
 *        règles de Bioche, équation de Risch, intégration par parties, radicaux,
 *        intégrales hors de portée et journal des étapes.
 *
 * Chaque primitive F de f est vérifiée par dérivation : F' - f est nul exactement quand la
 * forme canonique le permet, sinon F' et f sont comparés en plusieurs points (écart relatif
 * < 1e-10). Les textes attendus des cas phares sont aussi vérifiés.
 */

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "Journal.hpp"
#include "Lecture.hpp"
#include "Primitives.hpp"
#include "Solveur.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");

bool contientIntegrale(const ExprPtr& e) { return e->texte().find("integrale") != std::string::npos; }

/*
 * Cas : intégrande, texte attendu de la primitive (vide : non vérifié), F' - f nul
 * exactement exigé, points de vérification numérique.
 */
struct Cas {
    const char* f;
    const char* attendu;
    bool exact;
    std::vector<double> points;
};

void verifier(const Cas& c, const OptionsLecture& options = OptionsLecture{}) {
    const ExprPtr f = lire(c.f, options);
    const ExprPtr F = f->integrer();
    if (contientIntegrale(F)) {
        CHECK(!contientIntegrale(F));
        std::cerr << "      integrale(" << c.f << ") = " << F->texte() << "\n";
        return;
    }
    if (c.attendu[0] != '\0') {
        CHECK_EQ(F->texte(), std::string(c.attendu));
        // Le texte se relit en la même expression (invariant de l'affichage)
        CHECK(lire(F->texte(), options).get() == F.get());
    }
    const ExprPtr d = F->derivee();
    const ExprPtr ecart = (d - f)->simplifier();
    const bool nul = ecart->estConstante() && ecart->getValeurConstante() == 0.0;
    if (c.exact) {
        CHECK(nul);
        if (!nul) std::cerr << "      " << c.f << " : F' - f = " << ecart->texte() << "\n";
    }
    for (double x : c.points) {
        const double attendu = f->eval(x);
        CHECK_NEAR(d->eval(x), attendu, 1e-10 * (1.0 + std::fabs(attendu)));
    }
}

const std::vector<double> POINTS = {0.3, 0.7, 1.3};

} // namespace

// ============================================================================
// Table étendue
// ============================================================================

TEST_CASE(table_atan) {
    verifier({"1/(x^2 + 1)", "atan(x)", true, POINTS});
    verifier({"1/(x^2 + 4)", "atan(x/2)/2", false, {-1.2, 0.3, 2.5}});
    verifier({"1/(x^2 + 2)", "2^(1/2)*atan(2^(1/2)*x/2)/2", false, {-1.2, 0.3, 2.5}});
    verifier({"1/(1 + (2*x + 1)^2)", "atan(2*x + 1)/2", true, {-1.2, 0.3, 2.5}});
    verifier({"3/(2*x^2 + 3)", "", false, {-1.2, 0.3, 2.5}});
    verifier({"1/(x^2 + 2*x + 5)", "atan(x/2 + 1/2)/2", false, {-3.0, 0.3, 2.5}});
    verifier({"1/(x^2 + 2*x + 1)", "-1/(x + 1)", false, {-3.0, 0.3, 2.5}});
}

TEST_CASE(table_racines_carrees) {
    verifier({"1/sqrt(1 - x^2)", "asin(x)", true, {-0.6, 0.2, 0.9}});
    verifier({"1/sqrt(4 - x^2)", "asin(x/2)", false, {-1.6, 0.2, 1.9}});
    verifier({"1/sqrt(3 - 2*x - x^2)", "", false, {-2.5, -0.5, 0.6}});
    verifier({"1/sqrt(x^2 + 1)", "ln(x + (x^2 + 1)^(1/2))", false, {-2.0, 0.3, 4.0}});
    verifier({"1/sqrt(x^2 - 1)", "ln(abs(x + (x^2 - 1)^(1/2)))", false, {-2.5, 1.5, 4.0}});
}

TEST_CASE(table_trigonometrique_hyperbolique) {
    verifier({"1/cos(x)^2", "tan(x)", false, {-1.0, 0.3, 2.0}});
    verifier({"1/sin(x)^2", "-cos(x)/sin(x)", false, {-1.0, 0.3, 2.0}});
    verifier({"1/cos(3*x + 1)^2", "tan(3*x + 1)/3", false, {-1.0, 0.3, 2.0}});
    verifier({"tan(x)^2", "-x + tan(x)", true, {-1.0, 0.3, 2.0}});
    verifier({"tanh(x)^2", "x - tanh(x)", true, {-1.0, 0.3, 2.0}});
    verifier({"1/cosh(x)^2", "tanh(x)", false, {-1.0, 0.3, 2.0}});
    verifier({"1/sinh(2*x)^2", "-cosh(2*x)/(2*sinh(2*x))", false, {-1.0, 0.3, 2.0}});
    // tan -> -ln|cos| : juste aussi là où cos < 0
    verifier({"tan(x)", "-ln(abs(cos(x)))", false, {-1.0, 0.3, 2.0, 4.0}});
    verifier({"tan(2*x + 1)", "", false, {-1.0, 0.3, 2.0}});
    // Déjà présentes avant : argument affine
    verifier({"sinh(x)", "cosh(x)", true, POINTS});
    verifier({"cosh(2*x)", "sinh(2*x)/2", true, POINTS});
}

TEST_CASE(table_puissances_et_exponentielles) {
    verifier({"2^x", "2^x/ln(2)", true, POINTS});
    verifier({"pi^(3*x)", "pi^(3*x)/(3*ln(pi))", true, POINTS});
    verifier({"x^5", "x^6/6", true, POINTS});
    verifier({"(2*x + 1)^(-3)", "", true, POINTS});
    verifier({"sqrt(x)", "2*x^(3/2)/3", true, POINTS});
}

// ============================================================================
// 1/u = ln|u|
// ============================================================================

TEST_CASE(logarithme_valeur_absolue) {
    verifier({"1/x", "ln(abs(x))", true, {-2.0, -0.5, 0.5, 3.0}});
    verifier({"1/(3*x - 2)", "ln(abs(3*x - 2))/3", true, {-1.0, 0.2, 2.0}});
    verifier({"x/(x^2 - 4)", "ln(abs(x^2 - 4))/2", true, {-3.0, 0.5, 3.0}});
    // Argument manifestement positif : ln(u) sans valeur absolue
    verifier({"2*x/(x^2 + 1)", "ln(x^2 + 1)", true, {-2.0, 0.5, 3.0}});
    verifier({"exp(x)/(exp(x) + 1)", "ln(exp(x) + 1)", true, POINTS});
    // ln|u|' = u'/u directement : derivee(integrer(1/x)) redonne exactement 1/x
    const ExprPtr f = 1.0 / X;
    CHECK(f->integrer()->derivee().get() == f.get());
    CHECK(ast_ln(ast_abs(ast_pow(X, 2.0) - 1.0))->derivee().get() == (2.0 * X / (ast_pow(X, 2.0) - 1.0)).get());
    // Le solveur inverse ln(abs(u)) = R (exp puis abs)
    const Solutions s = resoudre(ast_ln(ast_abs(X - 1.0)), cst(0.0));
    CHECK_EQ(s.liste.size(), std::size_t(2));
    if (s.liste.size() == 2) {
        CHECK_EQ(s.liste[0].valeur->texte(), std::string("0"));
        CHECK_EQ(s.liste[1].valeur->texte(), std::string("2"));
    }
    const Solutions t = resoudre(lire("1/x")->integrer(), cst(0.0));
    CHECK_EQ(t.liste.size(), std::size_t(2)); // ln|x| = 0 : x = -1 ou x = 1
}

// ============================================================================
// Changement de variable « la dérivée divise »
// ============================================================================

TEST_CASE(substitution_cas_phares) {
    verifier({"x*cos(x^2)", "sin(x^2)/2", true, POINTS});
    verifier({"x*exp(x^2)", "exp(x^2)/2", true, POINTS});
    verifier({"cos(x)*sin(x)^3", "sin(x)^4/4", true, POINTS});
    verifier({"ln(x)/x", "ln(x)^2/2", true, POINTS});
    verifier({"exp(x)/(1 + exp(2*x))", "atan(exp(x))", false, {-1.0, 0.3, 2.0}});
    verifier({"tan(x)", "-ln(abs(cos(x)))", false, {0.3, 2.0}});
    verifier({"(2*x + 1)/(x^2 + x + 1)", "ln(x^2 + x + 1)", true, {-2.0, 0.3, 2.0}});
    // Coefficient non trivial : (4x + 2)/(2x + 1) ne se réduit que par simplifier()
    verifier({"(4*x + 2)/(x^2 + x + 1)", "2*ln(x^2 + x + 1)", false, {-2.0, 0.3, 2.0}});
    verifier({"(x + 1)/(x^2 + 2*x + 3)", "ln(x^2 + 2*x + 3)/2", false, {-2.0, 0.3, 2.0}});
    verifier({"x/sqrt(1 - x^2)", "-(-x^2 + 1)^(1/2)", true, {-0.6, 0.2, 0.9}});
}

TEST_CASE(substitution_autres_cas) {
    verifier({"x^2*(x^3 + 1)^5", "(x^3 + 1)^6/18", true, POINTS});
    verifier({"exp(sqrt(x))/sqrt(x)", "2*exp(x^(1/2))", true, POINTS});
    verifier({"1/(x*ln(x))", "ln(abs(ln(x)))", true, {0.5, 2.0, 3.0}});
    verifier({"cos(x)/(1 + sin(x)^2)", "atan(sin(x))", true, POINTS});
    verifier({"x*2^(x^2)", "2^(x^2)/(2*ln(2))", true, POINTS});
    verifier({"sin(x)/cos(x)^3", "", false, {-1.0, 0.3, 1.2}});
    verifier({"sin(2*x)*exp(cos(2*x))", "-exp(cos(2*x))/2", true, POINTS});
    verifier({"3*x^2/(x^3 + 1)", "ln(abs(x^3 + 1))", true, {-3.0, 0.3, 2.0}});
    verifier({"exp(x)*exp(exp(x))", "exp(exp(x))", true, POINTS});
    // Substitution dans une somme (linéarité) et facteur constant
    verifier({"3*x^2 + 4*x*cos(x^2) + 5", "x^3 + 5*x + 2*sin(x^2)", true, POINTS});
}

TEST_CASE(substitution_constantes_symboliques) {
    // Après u = g(x), un terme constant (pi, racine, paramètre) s'intègre en c u, pas en c x
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    verifier({"cos(x)*(sin(x)^2 + pi)", "sin(x)^3/3 + pi*sin(x)", false, points});
    verifier({"(ln(x)^2 + pi)/x", "ln(x)^3/3 + pi*ln(x)", false, points});
    verifier({"x*((x^2 + 1)^(1/2) + pi)", "", false, points});
    verifier({"sin(x)*(cos(x)^2 + 3^(1/2))", "", false, points});
    verifier({"x*sin(x^2)*(cos(x^2) + pi)", "", false, points});
    verifier({"exp(x)*sin(exp(x))*(cos(exp(x)) + pi)", "", false, points});
    CHECK_EQ(lire("cos(x)*(sin(x)^3 + a)")->integrer()->texte(), std::string("sin(x)^4/4 + a*sin(x)"));
    // Une somme en une autre variable que x
    CHECK_EQ(lire("u^2 + pi", OptionsLecture{"u"})->integrer()->texte(), std::string("u^3/3 + pi*u"));
    CHECK_EQ(lire("y + 1", OptionsLecture{"y"})->integrer()->texte(), std::string("y^2/2 + y"));
}

TEST_CASE(substitution_logarithme_et_exponentielle) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    // u = ln(x) : x = exp(u), puis retour exp(r ln(x)) = x^r
    verifier({"sin(ln(x))", "x*(sin(ln(x))/2 - cos(ln(x))/2)", false, points});
    verifier({"x*cos(ln(x))", "x^2*(sin(ln(x))/5 + 2*cos(ln(x))/5)", false, points});
    // u = exp(x) : exp(x + exp(x)) = u exp(u)
    verifier({"exp(x + exp(x))", "exp(exp(x))", false, {-1.3, 0.31, 0.77}});
    verifier({"exp(2*x + exp(x))", "(exp(x) - 1)*exp(exp(x))", false, {-1.3, 0.31, 0.77}});
    // Base symbolique : a^x = exp(x ln(a))
    CHECK_EQ(lire("a^x")->integrer()->texte(), std::string("a^x/ln(a)"));
    CHECK_EQ(lire("x*a^x")->integrer()->texte(), std::string("a^x*(-1/ln(a)^2 + x/ln(a))"));
}

TEST_CASE(substitution_autre_variable) {
    // Variable nommée y, puis u (le nom de la variable interne du changement de variable)
    verifier({"2*y/(y^2 + 1)", "ln(y^2 + 1)", true, POINTS}, OptionsLecture{"y"});
    verifier({"y + 1", "y^2/2 + y", true, POINTS}, OptionsLecture{"y"});
    verifier({"u*cos(u^2)", "sin(u^2)/2", true, POINTS}, OptionsLecture{"u"});
    verifier({"t*exp(t^2)", "exp(t^2)/2", true, POINTS}, OptionsLecture{"t"});
}

// ============================================================================
// Fractions rationnelles : division, Hermite, éléments simples
// ============================================================================

TEST_CASE(rationnelle_logarithmes) {
    const std::vector<double> points = {-3.7, -0.45, 0.31, 4.3};
    verifier({"1/(x^2 - 1)", "ln(abs((x - 1)/(x + 1)))/2", false, points});
    verifier({"1/(1 - x^2)", "ln(abs((x + 1)/(x - 1)))/2", false, points});
    verifier({"1/(4*x^2 - 1)", "ln(abs((2*x - 1)/(2*x + 1)))/4", false, points});
    // Logarithmes de même coefficient regroupés
    verifier({"x/(x^2 - 4)", "ln(abs(x^2 - 4))/2", true, points});
    verifier({"1/((x - 1)*(x - 2)*(x - 3))", "ln(abs(x^2 - 4*x + 3))/2 - ln(abs(x - 2))", false, {-3.7, 0.31, 2.5, 4.3}});
    verifier({"(x^5 + 2)/(x^3 - x)", "x^3/3 + x - 2*ln(abs(x)) + 3*ln(abs(x - 1))/2 + ln(abs(x + 1))/2", false, points});
    // Dérivée logarithmique d'un dénominateur irréductible de degré 3 : pas de factorisation
    verifier({"(3*x^2 + 1)/(x^3 + x + 1)", "ln(abs(x^3 + x + 1))", true, points});
    verifier({"x^3/(x^4 + 1)", "ln(x^4 + 1)/4", true, points});
    verifier({"1/(x^5 + x)", "-ln(x^4 + 1)/4 + ln(abs(x))", false, points});
    // Racines irrationnelles de degré 2
    verifier({"1/(x^2 - 2)", "2^(1/2)*ln(abs((x - 2^(1/2))/(x + 2^(1/2))))/4", false, points});
    verifier({"1/(x^2 - 3*x + 1)", "", false, points});
}

TEST_CASE(rationnelle_arctangentes) {
    const std::vector<double> points = {-3.7, -0.45, 0.31, 4.3};
    verifier({"1/(x^2 + x + 1)", "2*3^(1/2)*atan(2*3^(1/2)*(x + 1/2)/3)/3", false, points});
    verifier({"(x^3 + 1)/(x^2 + 1)", "x^2/2 - ln(x^2 + 1)/2 + atan(x)", false, points});
    verifier({"1/(x^3 + 1)", "3^(1/2)*atan(2*3^(1/2)*(x - 1/2)/3)/3 - ln(x^2 - x + 1)/6 + ln(abs(x + 1))/3", false, points});
    verifier({"1/(x^4 - 1)", "ln(abs((x - 1)/(x + 1)))/4 - atan(x)/2", false, points});
    // x^4 + x^2 + 1 = (x^2 - x + 1)(x^2 + x + 1) : factorisation complète sur Q nécessaire
    verifier({"1/(x^4 + x^2 + 1)", "", false, points});
    verifier({"1/(x^6 - 1)", "", false, points});
    verifier({"1/(2*x^2 + 3)", "6^(1/2)*atan(6^(1/2)*x/3)/6", false, points});
}

TEST_CASE(rationnelle_hermite) {
    const std::vector<double> points = {-3.7, -0.45, 0.31, 4.3};
    verifier({"1/(x^2 + 1)^2", "x/2/(x^2 + 1) + atan(x)/2", false, points});
    verifier({"x^2/(x^2 + 1)^2", "-x/2/(x^2 + 1) + atan(x)/2", false, points});
    verifier({"1/(x^2 + 1)^3", "(3*x^3 + 5*x)/(8*(x^2 + 1)^2) + 3*atan(x)/8", false, points});
    verifier({"1/(x*(x + 1)^2)", "1/(x + 1) + ln(abs(x/(x + 1)))", false, points});
    // Purement rationnelle : aucun logarithme
    verifier({"(2*x + 3)/(x^2 + 3*x + 5)^2", "-1/(x^2 + 3*x + 5)", true, points});
    verifier({"x/(1 - x)^3", "(2*x - 1)/(2*(x - 1)^2)", false, points});
    verifier({"x^4/(x - 1)^3", "x^2/2 + 3*x - (8*x - 7)/(2*(x - 1)^2) + 6*ln(abs(x - 1))", false, {-3.7, -0.45, 0.31, 4.3}});
    verifier({"1/(x^2*(x^2 + 2)^2)", "", false, points});
}

TEST_CASE(rationnelle_grandes_multiplicites) {
    // Vérification exacte : F' - f est la fraction nulle. La réduction d'Hermite accumule sa
    // partie rationnelle sur un dénominateur commun (une seule réduction, à la fin).
    for (const char* texte : {"1/(x^2 + 1)^6", "1/(x^2 + 1)^32", "1/((x - 1)^9*(x + 2)^7*(x^2 + x + 1)^3)", "x^3/(2*x^2 - 3)^12"}) {
        const ExprPtr f = lire(texte);
        const ExprPtr F = f->integrer();
        CHECK(!contientIntegrale(F));
        Fraction ecart;
        const bool rationnelle = versFraction(F->derivee() - f, X, ecart);
        CHECK(rationnelle);
        CHECK(rationnelle && ecart.estNulle());
    }
    verifier({"1/(x^2 + 1)^6", "(315*x^9 + 1470*x^7 + 2688*x^5 + 2370*x^3 + 965*x)/(1280*(x^2 + 1)^5) + 63*atan(x)/256", false, {}});
    verifier({"(x^2 + 2*x - 1)/(x^3 - x)", "ln(abs((x^2 - x)/(x + 1)))", false, {-3.7, -0.45, 0.31, 4.3}});
}

TEST_CASE(rationnelle_facteurs_de_degre_4) {
    // Facteur irréductible sur Q scindé en deux trinômes réels conjugués sur Q(sqrt(m))
    const std::vector<double> points = {-3.7, -0.45, 0.31, 1.9, 4.3};
    verifier({"1/(x^4 + 1)",
              "2^(1/2)*ln((x^2 + 2^(1/2)*x + 1)/(x^2 - 2^(1/2)*x + 1))/8 + 2^(1/2)*atan(2^(1/2)*(x - 2^(1/2)/2))/4 + "
              "2^(1/2)*atan(2^(1/2)*(x + 2^(1/2)/2))/4",
              false, points});
    verifier({"(x^2 + 1)/(x^4 + 1)", "2^(1/2)*atan(2^(1/2)*(x - 2^(1/2)/2))/2 + 2^(1/2)*atan(2^(1/2)*(x + 2^(1/2)/2))/2", false, points});
    verifier({"x^2/(x^4 + 1)", "", false, points});
    verifier({"(x^3 + 2*x^2 - x + 5)/(x^4 + 1)", "", false, points});
    verifier({"1/(x^4 + 1)^2", "", false, points});
    verifier({"1/((x^4 + 1)*(x - 1))", "", false, points});
    // x^6 + 1 = (x^2 + 1)(x^4 - x^2 + 1), second facteur scindé sur Q(sqrt(3))
    verifier({"1/(x^6 + 1)", "3^(1/2)*ln((x^2 + 3^(1/2)*x + 1)/(x^2 - 3^(1/2)*x + 1))/12 + atan(x)/3 + atan(2*x - 3^(1/2))/6 + atan(2*x + 3^(1/2))/6",
              false, points});
    // Terme en x^3 (décalage), résolvante cubique à racine rationnelle (polynômes cyclotomiques)
    verifier({"1/(x^4 + 4*x^3 + 6*x^2 + 4*x + 2)", "", false, points});
    verifier({"1/(x^4 + x^3 + x^2 + x + 1)", "", false, points});
    verifier({"(x^3 + 1)/(x^4 - x^3 + x^2 - x + 1)", "", false, points});
    // Trinômes à racines réelles : logarithmes, radicaux emboîtés
    verifier({"1/(x^4 - 2)", "", false, points});
    verifier({"1/(x^4 - 10*x^2 + 1)", "", false, {-0.45, 0.2, 0.31, 1.9, 4.3}});
    verifier({"1/(x^4 - 5*x^2 + 5)", "", false, {-0.45, 0.2, 0.31, 2.6, 4.3}});
    // Par ricochet : par parties, puis t = exp(x)
    verifier({"atan(x^2)", "", false, points});
    verifier({"tanh(2*x)*sinh(x)", "", false, {-1.3, 0.31, 0.77, 1.9}});
    // Factorisation réelle hors des corps quadratiques : laissé non évalué
    for (const char* texte : {"1/(x^4 + x + 1)", "1/(x^4 + 3)", "1/(x^4 + x^2 + 2)", "1/(x^8 + 1)"}) {
        CHECK(contientIntegrale(lire(texte)->integrer()));
    }
}

TEST_CASE(rationnelle_substitution_monomiale) {
    // u = x^k quand x f(x) est une fraction en x^k
    const std::vector<double> points = {-1.7, -0.45, 0.31, 2.3};
    verifier({"x/(x^4 + 1)", "atan(x^2)/2", true, points});
    verifier({"x^2/(x^6 + 1)", "atan(x^3)/3", true, points});
    verifier({"x^3/(x^8 + 1)", "atan(x^4)/4", true, points});
    verifier({"x/((x^4 + 1)*(x^2 + 1))", "ln(x^2 + 1)/4 - ln(x^4 + 1)/8 + atan(x^2)/4", false, points});
    // Au journal : la tentative en x, abandonnée, laisse la place au changement de variable
    Journal j;
    lire("x/((x^4 + 1)*(x^2 + 1))")->integrer(j);
    std::vector<std::string> regles;
    for (const Etape& e : j.etapes()) regles.push_back(e.regle);
    CHECK(regles == std::vector<std::string>({"integration.substitution", "integration.elements_simples", "integration.rationnelle"}));
    if (!j.vide()) CHECK_EQ(j.etapes()[0].description, std::string("Changement de variable u = x^2"));
}

TEST_CASE(rationnelle_outils) {
    // Fraction : toujours réduite, dénominateur unitaire
    const Fraction x = Fraction::variable();
    const Fraction un(Nombre(1));
    const Fraction f = (x * x - un) / (x - un); // x + 1
    CHECK(f.estPolynome());
    CHECK(f == x + un);
    CHECK_EQ(((x + un) / (x * x * Nombre(2) + x * Nombre(2))).versExpression(X)->texte(), std::string("1/(2*x)"));
    CHECK((un / x).puissance(-3) == x * x * x);
    CHECK((x / (x + un)).derivee() == (x + un).puissance(-2));
    // Composition : f(g) avec f = x/(x + 1), g = 1/x donne 1/(x + 1)
    CHECK((x / (x + un)).composer(un / x) == un / (x + un));
    CHECK_THROWS(un / Fraction(), std::domain_error);
    // versFraction : refuse ce qui n'est pas rationnel exact
    Fraction g;
    CHECK(versFraction(lire("(x^2 + 1)/(x - 3)^2 + 1/x"), X, g));
    CHECK(!versFraction(lire("sin(x)/x"), X, g));
    CHECK(!versFraction(lire("x^(1/2)"), X, g));
    CHECK(!versFraction(lire("a*x", {"x", {}}), X, g)); // paramètre
    CHECK(!versFraction(X * cst(0.5), X, g));              // coefficient réel
    // Résultat partiel signalé
    bool complete = true;
    CHECK(versFraction(lire("1/(x^3 + x + 1) + 1/x"), X, g));
    const ExprPtr F = primitiveFraction(g, X, &complete);
    CHECK(!complete);
    CHECK(contientIntegrale(F));
    // Produit de polynômes : changement de variable d'abord, développement sinon
    CHECK_EQ(lire("x^2*(x^3 + 1)^5")->integrer()->texte(), std::string("(x^3 + 1)^6/18"));
    CHECK_EQ(lire("(x + 1)^2*(x - 3)")->integrer()->texte(), std::string("x^4/4 - x^3/3 - 5*x^2/2 - 3*x"));
}

// ============================================================================
// Classe exponentielle-polynôme : polynômes, exp, sin, cos, sinh, cosh d'arguments affines
// ============================================================================

TEST_CASE(exponentielle_polynome) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"x*exp(x)", "(x - 1)*exp(x)", false, points});
    verifier({"x^2*exp(x)", "(x^2 - 2*x + 2)*exp(x)", false, points});
    verifier({"(x + 1)^2*exp(3*x)", "(x^2/3 + 4*x/9 + 5/27)*exp(3*x)", false, points});
    verifier({"x/exp(x)^2", "(-x/2 - 1/4)*exp(-2*x)", false, points});
    verifier({"x*sin(x)", "-x*cos(x) + sin(x)", false, points});
    verifier({"x^2*sin(2*x)", "x*sin(2*x)/2 - x^2*cos(2*x)/2 + cos(2*x)/4", false, points});
    verifier({"exp(x)*sin(x)", "(sin(x)/2 - cos(x)/2)*exp(x)", false, points});
    verifier({"exp(2*x)*sin(3*x)", "(2*sin(3*x)/13 - 3*cos(3*x)/13)*exp(2*x)", false, points});
    verifier({"exp(-x)*cos(2*x + 1)", "", false, points});
    verifier({"x*exp(x)*sin(x)", "", false, points});
    verifier({"x^2*exp(x)*cos(x)", "", false, points});
    // Pentes symboliques et bases constantes
    verifier({"x*sin(pi*x)", "sin(pi*x)/pi^2 - x*cos(pi*x)/pi", false, points});
    verifier({"x*2^x", "2^x*(-1/ln(2)^2 + x/ln(2))", false, points});
    verifier({"2^x*sin(x)", "", false, points});
    // Un facteur constant symbolique sort de l'intégrale
    CHECK_EQ(lire("a*x*exp(x)")->integrer()->texte(), std::string("a*(x - 1)*exp(x)"));
}

TEST_CASE(linearisation_trigonometrique) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"sin(x)^2", "x/2 - sin(2*x)/4", false, points});
    verifier({"cos(x)^2", "x/2 + sin(2*x)/4", false, points});
    verifier({"sin(x)^4", "3*x/8 - sin(2*x)/4 + sin(4*x)/32", false, points});
    verifier({"sin(x)^2*cos(x)^2", "x/8 - sin(4*x)/32", false, points});
    verifier({"sin(2*x)*cos(3*x)", "cos(x)/2 - cos(5*x)/10", false, points});
    verifier({"cos(3*x)*cos(5*x)", "sin(2*x)/4 + sin(8*x)/16", false, points});
    verifier({"sin(x)*sin(2*x)*sin(3*x)", "", false, points});
    verifier({"x*sin(x)^2", "", false, points});
    // Phases : un argument constant passe dans le coefficient
    verifier({"sin(x + 1)*cos(x - 1)", "x*sin(2)/2 - cos(2*x)/4", false, points});
}

TEST_CASE(hyperboliques_produits) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"sinh(x)^2", "-x/2 + sinh(2*x)/4", false, points});
    verifier({"x*cosh(x)", "x*sinh(x) - cosh(x)", false, points});
    verifier({"x^2*cosh(x)", "-2*x*cosh(x) + x^2*sinh(x) + 2*sinh(x)", false, points});
    verifier({"sinh(2*x)*cosh(3*x)", "-cosh(x)/2 + cosh(5*x)/10", false, points});
    verifier({"cos(x)*cosh(x)", "sin(x)*cosh(x)/2 + cos(x)*sinh(x)/2", false, points});
    verifier({"x*cos(x)*cosh(x)", "", false, points});
    verifier({"exp(x)*cosh(x)", "x/2 + exp(2*x)/4", false, points});
    verifier({"x*exp(x)*cosh(2*x)", "", false, points});
}

// ============================================================================
// Fractions rationnelles de sinus et cosinus : règles de Bioche, tan(u/2)
// ============================================================================

TEST_CASE(trigonometrie_puissances_impaires) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"sin(x)^3", "cos(x)^3/3 - cos(x)", false, points});
    verifier({"cos(x)^3", "-sin(x)^3/3 + sin(x)", false, points});
    verifier({"sin(x)^3*cos(x)^2", "-cos(x)^3/3 + cos(x)^5/5", false, points});
    verifier({"sin(x)^2*cos(x)^3", "sin(x)^3/3 - sin(x)^5/5", false, points});
    // Les deux règles s'appliquent : la plus courte
    verifier({"sin(x)^3*cos(x)", "sin(x)^4/4", false, points});
    verifier({"sin(2*x)^3", "cos(2*x)^3/6 - cos(2*x)/2", false, points});
}

TEST_CASE(trigonometrie_bioche) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    verifier({"1/sin(x)", "ln(abs((cos(x) - 1)/(cos(x) + 1)))/2", false, points});
    verifier({"1/cos(x)", "ln(abs((sin(x) + 1)/(sin(x) - 1)))/2", false, points});
    verifier({"sin(x)/cos(x)^2", "1/cos(x)", false, points});
    verifier({"sin(x)/(1 + cos(x)^2)", "-atan(cos(x))", false, points});
    verifier({"cos(x)/(1 + sin(x)^2)", "atan(sin(x))", false, points});
    verifier({"1/cos(x)^3", "", false, points});
    verifier({"cos(2*x)/cos(x)", "", false, points});
    // t = tan(u) ; atan(tan(u)) est récrit u
    verifier({"tan(x)^4", "x + tan(x)^3/3 - tan(x)", false, points});
    verifier({"1/cos(x)^4", "tan(x)^3/3 + tan(x)", false, points});
    verifier({"1/(1 + sin(x)^2)", "2^(1/2)*atan(2^(1/2)*tan(x))/2", false, points});
    verifier({"1/(4*sin(x)^2 + 9*cos(x)^2)", "atan(2*tan(x)/3)/6", false, points});
    verifier({"sin(x)/(sin(x) + cos(x))", "", false, points});
    verifier({"tan(3*x)^2", "-x + tan(3*x)/3", false, {0.1, 0.3, 0.7}});
}

TEST_CASE(trigonometrie_tangente_de_l_angle_moitie) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    verifier({"1/(1 + cos(x))", "tan(x/2)", false, points});
    verifier({"1/(1 + sin(x))", "-2/(tan(x/2) + 1)", false, points});
    verifier({"1/(2 + cos(x))", "2*3^(1/2)*atan(3^(1/2)*tan(x/2)/3)/3", false, points});
    verifier({"1/(3 + 5*cos(x))", "ln(abs((tan(x/2) + 2)/(tan(x/2) - 2)))/4", false, points});
    verifier({"1/(5 + 4*sin(x))", "2*atan(5*tan(x/2)/3 + 4/3)/3", false, points});
    verifier({"1/(1 + sin(x) + cos(x))", "ln(abs(tan(x/2) + 1))", false, points});
    verifier({"(1 - cos(x))/(1 + cos(x))", "-x + 2*tan(x/2)", false, points});
    verifier({"1/(sin(x) + cos(x))", "", false, points});
}

// ============================================================================
// Fractions rationnelles d'exponentielles et de fonctions hyperboliques : t = exp(u)
// ============================================================================

TEST_CASE(exponentielle_rationnelle) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"1/(exp(x) + 1)", "x - ln(exp(x) + 1)", false, points});
    verifier({"1/(1 + exp(2*x))", "x - ln(exp(2*x) + 1)/2", false, points});
    verifier({"(exp(x) - 1)/(exp(x) + 1)", "-x + 2*ln(exp(x) + 1)", false, points});
    verifier({"1/cosh(x)", "2*atan(exp(x))", false, points});
    verifier({"1/sinh(x)", "ln(abs((exp(x) - 1)/(exp(x) + 1)))", false, points});
    verifier({"1/(sinh(x) + 2*cosh(x))", "2*3^(1/2)*atan(3^(1/2)*exp(x))/3", false, points});
    verifier({"1/(exp(x/2) + exp(x))", "", false, points});
    verifier({"exp(3*x)/(exp(x) + 2)", "", false, points});
    verifier({"1/(exp(x) - 1)^2", "", false, points});
}

// ============================================================================
// Équation de Risch : R(x) exp(W(x))
// ============================================================================

TEST_CASE(equation_de_risch) {
    const std::vector<double> points = {-1.3, 0.31, 0.77, 1.9};
    verifier({"exp(x)*(x - 1)/x^2", "exp(x)/x", false, points});
    verifier({"x*exp(x)/(x + 1)^2", "exp(x)/(x + 1)", false, points});
    verifier({"(x^2 - 2*x + 2)*exp(x)/x^3", "(x - 1)*exp(x)/x^2", false, points});
    verifier({"(2*x^2 + 1)*exp(x^2)", "x*exp(x^2)", false, points});
    verifier({"(1 + 2*x)*exp(x + x^2)", "exp(x^2 + x)", true, points});
    verifier({"x^3*exp(x^2)", "(x^2/2 - 1/2)*exp(x^2)", false, points});
    // Somme dont aucun terme ne s'intègre seul
    verifier({"exp(x)/x - exp(x)/x^2", "exp(x)/x", false, points});
    verifier({"exp(x)/x - exp(x)/x^2 + cos(x)", "exp(x)/x + sin(x)", false, points});
    // Pas de solution rationnelle : primitive non élémentaire, laissée non évaluée
    for (const char* texte : {"exp(x^2)", "x^2*exp(x^2)", "exp(x)/x", "exp(x^2)*(2*x^3 + 3*x^2 - x + 1)", "exp(x)/x + exp(x)/x^2"}) {
        CHECK(contientIntegrale(lire(texte)->integrer()));
    }
}

// ============================================================================
// Intégration par parties
// ============================================================================

TEST_CASE(integration_par_parties) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    verifier({"x*ln(x)", "-x^2/4 + x^2*ln(x)/2", false, points});
    verifier({"x^2*ln(x)", "-x^3/9 + x^3*ln(x)/3", false, points});
    verifier({"ln(x)^2", "2*x - 2*x*ln(x) + x*ln(x)^2", false, points});
    verifier({"ln(x)^3", "-6*x + 6*x*ln(x) - 3*x*ln(x)^2 + x*ln(x)^3", false, points});
    verifier({"ln(x)/x^2", "-1/x - ln(x)/x", false, points});
    verifier({"ln(x^2 + 1)", "-2*x + x*ln(x^2 + 1) + 2*atan(x)", false, points});
    verifier({"x*atan(x)", "-x/2 + x^2*atan(x)/2 + atan(x)/2", false, points});
    verifier({"atan(x)/x^2", "-atan(x)/x - ln(x^2 + 1)/2 + ln(abs(x))", false, points});
    verifier({"x*ln(x)^2", "", false, points});
    verifier({"x*ln(x^2 + x + 1)", "", false, points});
    verifier({"ln(x + (x^2 + 1)^(1/2))", "-(x^2 + 1)^(1/2) + x*ln(x + (x^2 + 1)^(1/2))", false, points});
    verifier({"atan(x^(1/2))", "", false, points});
    // Arcs : le reste est une intégrale à radical
    const std::vector<double> arcs = {-0.6, 0.2, 0.5, 0.8};
    verifier({"x*asin(x)", "x*(-x^2 + 1)^(1/2)/4 + x^2*asin(x)/2 - asin(x)/4", false, arcs});
    verifier({"x*acos(x)", "", false, arcs});
    verifier({"asin(x)^2", "-2*x + x*asin(x)^2 + 2*(-x^2 + 1)^(1/2)*asin(x)", false, arcs});
    verifier({"asin(x)/x^2", "", false, arcs});
}

// ============================================================================
// Radicaux : (a x + b)^(p/q) et sqrt(a x^2 + b x + c)
// ============================================================================

TEST_CASE(radicaux_d_une_fonction_affine) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    verifier({"x*(x + 1)^(1/2)", "-2*(x + 1)^(3/2)/3 + 2*(x + 1)^(5/2)/5", false, points});
    verifier({"x^2*(2*x + 1)^(1/2)", "", false, points});
    verifier({"1/(1 + x^(1/2))", "2*x^(1/2) - 2*ln(x^(1/2) + 1)", false, points});
    verifier({"x^(1/2)/(1 + x)", "2*x^(1/2) - 2*atan(x^(1/2))", false, points});
    verifier({"1/(x^(1/2) + x^(1/3))", "2*x^(1/2) - 3*x^(1/3) + 6*x^(1/6) - 6*ln(x^(1/6) + 1)", false, points});
    verifier({"1/(x*(x + 1)^(1/2))", "ln(abs(((x + 1)^(1/2) - 1)/((x + 1)^(1/2) + 1)))", false, points});
    verifier({"(x + 1)^(1/3)*x", "-3*(x + 1)^(4/3)/4 + 3*(x + 1)^(7/3)/7", false, points});
    verifier({"x^(1/2)/(x - 1)", "", false, {0.31, 0.77, 1.3, 2.6}});
}

TEST_CASE(radicaux_d_un_trinome) {
    const std::vector<double> points = {0.31, 0.77, 1.3, 2.6};
    const std::vector<double> dansLeDisque = {-0.6, 0.2, 0.5, 0.8};
    // Coefficients indéterminés : S(x) y + lambda int dx/y
    verifier({"(x^2 + 1)^(1/2)", "x*(x^2 + 1)^(1/2)/2 + ln(x + (x^2 + 1)^(1/2))/2", false, points});
    verifier({"(1 - x^2)^(1/2)", "x*(-x^2 + 1)^(1/2)/2 + asin(x)/2", false, dansLeDisque});
    verifier({"(x^2 - 1)^(1/2)", "x*(x^2 - 1)^(1/2)/2 - ln(abs(x + (x^2 - 1)^(1/2)))/2", false, {1.3, 1.9, 2.6}});
    verifier({"x^2/(1 - x^2)^(1/2)", "-x*(-x^2 + 1)^(1/2)/2 + asin(x)/2", false, dansLeDisque});
    verifier({"x^2*(x^2 + 1)^(1/2)", "(x^3/4 + x/8)*(x^2 + 1)^(1/2) - ln(x + (x^2 + 1)^(1/2))/8", false, points});
    verifier({"x^3/(x^2 + 1)^(1/2)", "(x^2/3 - 2/3)*(x^2 + 1)^(1/2)", false, points});
    verifier({"(x^2 + 2*x + 5)^(1/2)", "", false, points});
    verifier({"(3 - 2*x - x^2)^(1/2)", "(-x^2 - 2*x + 3)^(1/2)*(x/2 + 1/2) + 2*asin(x/2 + 1/2)", false, dansLeDisque});
    verifier({"(x^2 + 1)^(3/2)", "", false, points});
    verifier({"(2*x^2 + 3)^(1/2)", "", false, points});
    // Substitutions d'Euler ; la partie rationnelle revient sous la forme A(x) + B(x) y
    verifier({"1/(x*(x^2 + 1)^(1/2))", "ln(abs((x + (x^2 + 1)^(1/2) - 1)/(x + (x^2 + 1)^(1/2) + 1)))", false, points});
    verifier({"1/(x^2*(x^2 + 1)^(1/2))", "-(x^2 + 1)^(1/2)/x", false, points});
    verifier({"1/(x^2 + 1)^(3/2)", "x/(x^2 + 1)^(1/2)", false, points});
    verifier({"1/(x + (x^2 + 1)^(1/2))", "-x^2/2 + x*(x^2 + 1)^(1/2)/2 + ln(x + (x^2 + 1)^(1/2))/2", false, points});
    verifier({"1/((x + 1)*(x^2 + 1)^(1/2))", "", false, points});
    verifier({"1/(x*(2*x^2 + 3)^(1/2))", "", false, points});
    // Racines rationnelles, coefficient dominant négatif : t = y/(x - alpha)
    verifier({"1/(x*(1 - x^2)^(1/2))", "", false, {0.2, 0.5, 0.8}});
    verifier({"1/(x*(4 - x^2)^(1/2))", "", false, {0.2, 0.5, 0.8, 1.5}});
    verifier({"1/((x - 1)*(x^2 - 3*x + 2)^(1/2))", "2*(x^2 - 3*x + 2)^(1/2)/(x - 1)", false, {2.3, 2.6, 3.5}});
    // Trinôme carré parfait : y = sqrt(a) |x + h|, primitive écrite avec sign(x + h)
    const std::vector<double> deuxCotes = {-2.4, -0.45, 0.31, 0.77, 1.9, 3.7};
    verifier({"(x^2)^(1/2)", "x^2*sign(x)/2", false, deuxCotes});
    verifier({"1/(x*(x^2 - 2*x + 1)^(1/2))", "ln(abs((x - 1)/x))*sign(x - 1)", false, deuxCotes});
    verifier({"(x^2 - 2*x + 1)^(1/2)/x", "(x - ln(abs(x)))*sign(x - 1)", false, deuxCotes});
    verifier({"1/((x^2 - 2*x + 1)^(1/2) + x)", "", false, deuxCotes});
    verifier({"1/((x^2 + 1)*(x^2 - 2*x + 1)^(1/2))", "", false, deuxCotes});
    verifier({"x/((x + 3)*(4*x^2 + 4*x + 1)^(1/2))", "", false, deuxCotes});
    verifier({"1/(x^2*(x^2)^(1/2))", "-sign(x)/(2*x^2)", false, deuxCotes});
    verifier({"1/(x*(x^2 - 2*x + 1)^(3/2))", "", false, deuxCotes});
    verifier({"1/(x^2 + 2*x + 1)^(1/2)", "ln(abs(x + 1))*sign(x + 1)", false, deuxCotes});
    // Jamais réel : laissé non évalué
    CHECK(contientIntegrale(lire("(-x^2 - 1)^(1/2)")->integrer()));
    // Hors de portée : racine d'un polynôme de degré 3 (intégrale elliptique)
    CHECK(contientIntegrale(lire("(x^3 + 1)^(1/2)")->integrer()));
    CHECK(contientIntegrale(lire("1/(x^4 + 1)^(1/2)")->integrer()));
}

TEST_CASE(outils_des_changements_de_variable) {
    const Fraction x = Fraction::variable(), un(Nombre(1));
    // Parité et fraction en x^k
    CHECK_EQ((x * x / (x * x + un)).parite(), 1);
    CHECK_EQ((x / (x * x + un)).parite(), -1);
    CHECK_EQ((x / (x + un)).parite(), 0);
    CHECK_EQ(Fraction().parite(), 1);
    CHECK((x * x / (x * x * x * x + un)).enPuissance(2) == x / (x * x + un));
    // Extension quadratique y^2 = x^2 + 1 : (x + y)(y - x) = 1, 1/(x + y) = y - x
    const Fraction carre = x * x + un;
    const ElementQuadratique t{x, un, &carre};
    const ElementQuadratique inverse = t.puissance(-1);
    CHECK(inverse.a == -x);
    CHECK(inverse.b == un);
    const ElementQuadratique produit = t * inverse;
    CHECK(produit.a == un);
    CHECK(produit.b.estNulle());
    // Évaluation d'une fraction : (t^2 - 1)/(2 t) = x pour t = x + y
    const ElementQuadratique valeur = t.evaluer((x * x - un) / (x * Nombre(2)));
    CHECK(valeur.a == x);
    CHECK(valeur.b.estNulle());
    CHECK_THROWS((ElementQuadratique{Fraction(), Fraction(), &carre}).puissance(-1), std::domain_error);
    // Racine d'ordre pair : positive là où elle est définie
    CHECK(estPositifOuNul(lire("x^(1/2) + 1")));
    CHECK(!estPositifOuNul(lire("x^(1/3) + 1")));
}

// ============================================================================
// Hors de portée : intégrale non évaluée, jamais de résultat faux
// ============================================================================

TEST_CASE(integrales_hors_de_portee) {
    // Primitives non élémentaires (exp(x^2), sin(x)/x, 1/ln(x), x^x...) ou demandant des
    // nombres algébriques de degré > 2
    for (const char* texte : {"exp(x^2)", "sin(x^2)", "exp(x)/x", "sin(x)/x", "1/ln(x)", "x^x", "sin(sin(x))", "ln(ln(x))",
                              "x^2*exp(x^2)", "exp(x)*ln(x)", "1/(x^3 + x + 1)", "1/(x^4 + x + 1)", "(x^3 + 1)^(1/2)",
                              "exp(x^2)/x", "tan(x)/x"}) {
        const ExprPtr f = lire(texte);
        const ExprPtr F = f->integrer();
        CHECK(comme<IntegraleNonEvaluee>(F) != nullptr);
        if (!comme<IntegraleNonEvaluee>(F)) std::cerr << "      integrale(" << texte << ") = " << F->texte() << "\n";
        CHECK(F->derivee().get() == f.get());
    }
    // Plusieurs variables : refusé
    const ExprPtr g = lire("x*cos(x*y)", {"x", {"y"}})->integrer();
    CHECK(comme<IntegraleNonEvaluee>(g) != nullptr);
}

// ============================================================================
// Journal des étapes
// ============================================================================

TEST_CASE(journal_substitution) {
    const ExprPtr f = lire("2*x*cos(x^2)");
    Journal j;
    const ExprPtr F = f->integrer(j);
    CHECK(F.get() == lire("sin(x^2)").get());
    CHECK(journal::actif() == nullptr); // désactivé à la sortie
    const std::vector<Etape>& e = j.etapes();
    CHECK_EQ(e.size(), std::size_t(3));
    if (e.size() == 3) {
        CHECK_EQ(e[0].regle, std::string("integration.substitution"));
        CHECK_EQ(e[0].description, std::string("Changement de variable u = x^2, du = 2*x dx"));
        CHECK(e[0].avant.get() == f.get());
        CHECK_EQ(e[0].apres->texte(), std::string("cos(u)"));
        CHECK_EQ(e[0].niveau, 0);
        CHECK_EQ(e[1].regle, std::string("integration.table"));
        CHECK_EQ(e[1].apres->texte(), std::string("sin(u)"));
        CHECK_EQ(e[1].niveau, 1);
        CHECK_EQ(e[2].regle, std::string("integration.retour"));
        CHECK(e[2].apres.get() == F.get());
        CHECK_EQ(e[2].niveau, 0);
    }
    // Sans journal : même résultat, rien n'est noté
    CHECK(f->integrer().get() == F.get());
    CHECK_EQ(j.etapes().size(), std::size_t(3));
}

TEST_CASE(journal_fraction_rationnelle) {
    // Dérivée logarithmique : une seule étape
    Journal j;
    const ExprPtr F = lire("2*x/(x^2 + 1)")->integrer(j);
    CHECK(F.get() == lire("ln(x^2 + 1)").get());
    CHECK_EQ(j.etapes().size(), std::size_t(1));
    if (!j.etapes().empty()) CHECK_EQ(j.etapes()[0].regle, std::string("integration.rationnelle"));
    // Pôles multiples puis éléments simples
    Journal k;
    const ExprPtr f = lire("1/(x^2*(x - 1)*(x + 2))");
    const ExprPtr G = f->integrer(k);
    const std::vector<Etape>& e = k.etapes();
    CHECK_EQ(e.size(), std::size_t(3));
    if (e.size() == 3) {
        CHECK_EQ(e[0].regle, std::string("integration.hermite"));
        CHECK_EQ(e[1].regle, std::string("integration.elements_simples"));
        CHECK_EQ(e[2].regle, std::string("integration.rationnelle"));
        CHECK(e[2].avant.get() == f.get());
        CHECK(e[2].apres.get() == G.get());
    }
    // Échec : les étapes de la tentative sont retirées
    Journal l;
    lire("1/(x^3 + x + 1)")->integrer(l);
    CHECK_EQ(l.etapes().size(), std::size_t(1));
    if (!l.etapes().empty()) CHECK_EQ(l.etapes()[0].regle, std::string("integration.non_evaluee"));
}

TEST_CASE(journal_linearite_et_echec) {
    Journal j;
    lire("3*x^2 + pi*exp(x)/(1 + exp(2*x))")->integrer(j);
    std::vector<std::string> regles;
    std::vector<int> niveaux;
    for (const Etape& e : j.etapes()) {
        regles.push_back(e.regle);
        niveaux.push_back(e.niveau);
    }
    // Somme terme à terme ; x^2 (le 3 est le coefficient du terme) ; pi sort de l'intégrale ;
    // t = exp(x) ; fraction rationnelle 1/(t^2 + 1)
    const std::vector<std::string> attendues = {"integration.linearite", "integration.table", "integration.linearite",
                                                "integration.substitution", "integration.rationnelle"};
    const std::vector<int> niveauxAttendus = {0, 1, 1, 2, 2};
    CHECK(regles == attendues);
    CHECK(niveaux == niveauxAttendus);
    if (regles != attendues || niveaux != niveauxAttendus) {
        for (const Etape& e : j.etapes()) std::cerr << "      " << e.niveau << " " << e.regle << " : " << e.description << "\n";
    }
    // Un changement de variable abandonné (u = x^2 : exp(u)/(2u) n'a pas de primitive élémentaire)
    // ne laisse pas d'étape
    Journal k;
    lire("exp(x^2)/x")->integrer(k);
    CHECK_EQ(k.etapes().size(), std::size_t(1));
    if (!k.vide()) CHECK_EQ(k.etapes()[0].regle, std::string("integration.non_evaluee"));
}

TEST_CASE(journal_nouvelles_methodes) {
    auto regles = [](const char* texte) {
        Journal j;
        lire(texte)->integrer(j);
        std::string r;
        for (const Etape& e : j.etapes()) r += (r.empty() ? "" : " ") + e.regle.substr(std::string("integration.").size());
        return r;
    };
    CHECK_EQ(regles("x*exp(x)*sin(x)"), std::string("exponentielle_polynome"));
    CHECK_EQ(regles("sin(x)^2"), std::string("linearisation exponentielle_polynome"));
    CHECK_EQ(regles("sin(x)^3*cos(x)^2"), std::string("bioche rationnelle"));
    CHECK_EQ(regles("1/(2 + cos(x))"), std::string("bioche rationnelle"));
    CHECK_EQ(regles("1/cosh(x)"), std::string("substitution rationnelle"));
    // Par parties : la primitive du cofacteur ne laisse pas d'étape
    CHECK_EQ(regles("x^2*ln(x)"), std::string("parties linearite table"));
    CHECK_EQ(regles("x^2/(1 - x^2)^(1/2)"), std::string("radicaux"));
    CHECK_EQ(regles("x*(x + 1)^(1/2)"), std::string("substitution radicaux"));
    CHECK_EQ(regles("1/(x*(x^2 + 1)^(1/2))"), std::string("euler elements_simples radicaux"));
    CHECK_EQ(regles("exp(x)/x - exp(x)/x^2"), std::string("linearite non_evaluee non_evaluee somme"));
    CHECK_EQ(regles("exp(x)*(x - 1)/x^2"), std::string("risch"));
    // Une tentative abandonnée ne laisse rien : par parties (exp(x)/x n'a pas de primitive élémentaire)
    CHECK_EQ(regles("exp(x)*ln(x)"), std::string("non_evaluee"));
}

int main() { return test::executerTous(); }
