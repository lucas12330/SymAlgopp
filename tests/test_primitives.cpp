/**
 * @file test_primitives.cpp
 * @brief Tests de l'intégration symbolique : table étendue, ln|u|, fractions rationnelles,
 *        changement de variable « la dérivée divise », intégrales hors de portée et journal
 *        des étapes.
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

TEST_CASE(rationnelle_substitution_monomiale) {
    // u = x^k quand x f(x) est une fraction en x^k
    const std::vector<double> points = {-1.7, -0.45, 0.31, 2.3};
    verifier({"x/(x^4 + 1)", "atan(x^2)/2", true, points});
    verifier({"x^2/(x^6 + 1)", "atan(x^3)/3", true, points});
    verifier({"x^3/(x^8 + 1)", "atan(x^4)/4", true, points});
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
// Hors de portée : intégrale non évaluée, jamais de résultat faux
// ============================================================================

TEST_CASE(integrales_hors_de_portee) {
    for (const char* texte : {"x*exp(x)", "1/(x^3 + x + 1)", "sin(x)^2", "exp(x)*sin(x)", "exp(x^2)", "sin(x^2)",
                              "x*sin(x)", "x^3*exp(x^2)", "1/(x^4 + 1)", "sin(sin(x))", "ln(ln(x))"}) {
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
    // Somme terme à terme ; x^2 (le 3 est le coefficient du terme) ; u = exp(x) ; pi sort de
    // l'intégrale de pi/(u^2 + 1) ; 1/(u^2 + 1) ; retour à x
    const std::vector<std::string> attendues = {"integration.linearite",    "integration.table",
                                                "integration.substitution", "integration.linearite",
                                                "integration.table",        "integration.retour"};
    CHECK(regles == attendues);
    CHECK(niveaux == std::vector<int>({0, 1, 1, 2, 3, 1}));
    if (regles != attendues || niveaux != std::vector<int>({0, 1, 1, 2, 3, 1})) {
        for (const Etape& e : j.etapes()) std::cerr << "      " << e.niveau << " " << e.regle << " : " << e.description << "\n";
    }
    // Un changement de variable abandonné (u = x^2 : u exp(u) / 2 n'est pas intégrable ici)
    // ne laisse pas d'étape
    Journal k;
    lire("x^3*exp(x^2)")->integrer(k);
    CHECK_EQ(k.etapes().size(), std::size_t(1));
    if (!k.vide()) CHECK_EQ(k.etapes()[0].regle, std::string("integration.non_evaluee"));
}

int main() { return test::executerTous(); }
