/**
 * @file test_latex.cpp
 * @brief Tests du rendu LaTeX (expressions, solutions, journal des étapes) et du journal de
 *        résolution : formes attendues, accolades équilibrées sur des expressions
 *        aléatoires, expression à fort partage, et équivalence résolution avec / sans journal.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <symalgopp> // l'en-tête unique suffit (et vérifie qu'il compile seul)

#include "test_framework.hpp"

using namespace symalgo;

namespace {

// Créées dans cet ordre : x avant y (l'ordre canonique des facteurs en dépend)
const ExprPtr X = var("x");
const ExprPtr Y = var("y");

std::string L(const std::string& texte, const OptionsLatex& o = {}) { return latex(lire(texte, {"x", {"y"}}), o); }

// Accolades et \left / \right équilibrés
bool equilibre(const std::string& s) {
    int accolades = 0, gauches = 0, droits = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\\' && i + 1 < s.size()) {
            if (s.compare(i, 5, "\\left") == 0) ++gauches;
            if (s.compare(i, 6, "\\right") == 0) ++droits;
            if (s[i + 1] == '{' || s[i + 1] == '}') ++i; // accolade échappée
            continue;
        }
        if (s[i] == '{') ++accolades;
        if (s[i] == '}' && --accolades < 0) return false;
    }
    return accolades == 0 && gauches == droits;
}

bool contient(const std::string& s, const std::string& morceau) { return s.find(morceau) != std::string::npos; }

} // namespace

TEST_CASE(latex_polynomes_et_signes) {
    CHECK_EQ(L("x^2 + 5*x + 6"), "x^{2} + 5 x + 6");
    CHECK_EQ(L("-x"), "-x");
    CHECK_EQ(L("x - 1"), "x - 1");
    CHECK_EQ(L("2*x*y"), "2 x y");
    CHECK_EQ(L("6 + 5*x + x^2"), "x^{2} + 5 x + 6");
    CHECK_EQ(L("-x^2 - x - 1"), "-x^{2} - x - 1");
    CHECK_EQ(L("3"), "3");
    CHECK_EQ(L("-3"), "-3");
    CHECK_EQ(L("(x + 1)^2"), "\\left(x + 1\\right)^{2}");
    CHECK_EQ(L("(-2)^x"), "\\left(-2\\right)^{x}");
    CHECK_EQ(L("x^x"), "x^{x}");
    CHECK_EQ(L("x^(-x)"), "x^{-x}");
    CHECK_EQ(L("x^2*y^3"), "x^{2} y^{3}");
}

TEST_CASE(latex_fractions) {
    CHECK_EQ(L("(x+1)/(x-1)"), "\\frac{x + 1}{x - 1}");
    CHECK_EQ(L("5/7"), "\\frac{5}{7}");
    CHECK_EQ(L("-5/7"), "-\\frac{5}{7}");
    CHECK_EQ(L("-3/4*x"), "-\\frac{3 x}{4}");
    CHECK_EQ(L("1/(x+1)"), "\\frac{1}{x + 1}");
    CHECK_EQ(L("sin(x)/x"), "\\frac{\\sin\\left(x\\right)}{x}");
    CHECK_EQ(L("3/(2*x)"), "\\frac{3}{2 x}");
    CHECK_EQ(L("x/2 + 1/3"), "\\frac{x}{2} + \\frac{1}{3}");
    CHECK_EQ(L("2*(x+1)/(x-1)"), "\\frac{2 \\left(x + 1\\right)}{x - 1}");
    OptionsLatex sans;
    sans.fractions = false;
    CHECK_EQ(L("5/7", sans), "5/7");
    CHECK_EQ(L("(x+1)/(x-1)", sans), "\\left(x + 1\\right)/\\left(x - 1\\right)");
    CHECK_EQ(L("3/(2*x)", sans), "3/\\left(2 x\\right)");
}

TEST_CASE(latex_racines_et_puissances) {
    CHECK_EQ(L("x^(1/2)"), "\\sqrt{x}");
    CHECK_EQ(L("x^(1/3)"), "\\sqrt[3]{x}");
    CHECK_EQ(L("2^(1/2)/2"), "\\frac{\\sqrt{2}}{2}");
    CHECK_EQ(L("x^(-1/2)"), "\\frac{1}{\\sqrt{x}}");
    CHECK_EQ(L("(x^2+1)^(-1/2)"), "\\frac{1}{\\sqrt{x^{2} + 1}}");
    CHECK_EQ(L("x^(2/3)"), "x^{\\frac{2}{3}}");
    OptionsLatex sans;
    sans.racines = false;
    CHECK_EQ(L("x^(1/2)", sans), "x^{\\frac{1}{2}}");
    CHECK_EQ(L("x^(-1/2)", sans), "\\frac{1}{x^{\\frac{1}{2}}}");
}

TEST_CASE(latex_fonctions) {
    CHECK_EQ(L("sin(x)^2"), "\\sin^{2}\\left(x\\right)");
    CHECK_EQ(L("sin(x)^(1/2)"), "\\sqrt{\\sin\\left(x\\right)}");
    CHECK_EQ(L("sin(x)^x"), "\\sin\\left(x\\right)^{x}");
    CHECK_EQ(L("1/cos(x)^2"), "\\frac{1}{\\cos^{2}\\left(x\\right)}");
    CHECK_EQ(L("exp(x)"), "e^{x}");
    CHECK_EQ(L("exp(x)^2"), "\\left(e^{x}\\right)^{2}");
    CHECK_EQ(L("2*exp(-x)"), "2 e^{-x}");
    CHECK_EQ(L("ln(x)"), "\\ln\\left(x\\right)");
    CHECK_EQ(L("ln(x)^2"), "\\ln^{2}\\left(x\\right)");
    CHECK_EQ(L("asin(x)"), "\\arcsin\\left(x\\right)");
    CHECK_EQ(L("acos(x)"), "\\arccos\\left(x\\right)");
    CHECK_EQ(L("atan(2*x)"), "\\arctan\\left(2 x\\right)");
    CHECK_EQ(L("tan(x)"), "\\tan\\left(x\\right)");
    CHECK_EQ(L("abs(x - 1)"), "\\left|x - 1\\right|");
    CHECK_EQ(L("abs(x)^2 * 3"), "3 \\left|x\\right|^{2}");
    CHECK_EQ(L("sign(x)"), "\\operatorname{sgn}\\left(x\\right)");
    CHECK_EQ(L("sinh(x)"), "\\sinh\\left(x\\right)");
    CHECK_EQ(L("cosh(2*x)"), "\\cosh\\left(2 x\\right)");
    CHECK_EQ(L("tanh(x)^2"), "\\tanh^{2}\\left(x\\right)");
}

TEST_CASE(latex_noms_et_constantes) {
    CHECK_EQ(latex(pi()), "\\pi");
    CHECK_EQ(L("2*pi*x"), "2 x \\pi");
    CHECK_EQ(latex(param("alpha")), "\\alpha");
    CHECK_EQ(latex(param("alpha") * param("beta")), "\\alpha \\beta");
    CHECK_EQ(latex(param("theta1")), "\\theta_{1}");
    CHECK_EQ(latex(param("C1") * X), "C_{1} x");
    CHECK_EQ(latex(param("k_2")), "k_{2}");
    CHECK_EQ(latex(param("vitesse")), "\\mathrm{vitesse}");
    CHECK_EQ(latex(param("Delta")), "\\Delta");
    CHECK_EQ(latex(param("a")), "a");
    // Un nombre devant un nombre : \cdot obligatoire
    CHECK_EQ(L("2*3^x"), "2 \\cdot 3^{x}");
    OptionsLatex explicite;
    explicite.multiplicationExplicite = true;
    CHECK_EQ(L("2*x*y", explicite), "2 \\cdot x \\cdot y");
    CHECK_EQ(L("(x+1)/(x-1)", explicite), "\\frac{x + 1}{x - 1}");
    CHECK_EQ(L("3*x/(2*y)", explicite), "\\frac{3 \\cdot x}{2 \\cdot y}");
}

TEST_CASE(latex_nombres_reels) {
    auto reel = [](double v) { return latex(nombre(Nombre::reel(v))); };
    CHECK_EQ(reel(1.5), "1.5");
    CHECK_EQ(reel(-1.5), "-1.5");
    CHECK_EQ(reel(0.1), "0.1");
    CHECK_EQ(reel(1e-10), "10^{-10}");
    CHECK_EQ(reel(2.5e20), "2.5 \\times 10^{20}");
    CHECK_EQ(reel(-3e-7), "-3 \\times 10^{-7}");
    CHECK_EQ(reel(std::numeric_limits<double>::infinity()), "\\infty");
    CHECK_EQ(reel(-std::numeric_limits<double>::infinity()), "-\\infty");
    CHECK_EQ(reel(std::numeric_limits<double>::quiet_NaN()), "\\mathrm{NaN}");
    // La notation se relit au même double
    for (double v : {0.1, 1.0 / 3.0, 6.02214076e23, 1.7976931348623157e308, 5e-324, 123456.789}) {
        const std::string t = reel(v);
        std::string c;
        // 10^{e} / m \times 10^{e} -> m e e
        const std::size_t p = t.find("10^{");
        if (p == std::string::npos) {
            c = t;
        } else {
            std::string m = p == 0 ? "1" : t.substr(0, t.find(" \\times"));
            c = m + "e" + t.substr(p + 4, t.size() - p - 5);
        }
        CHECK(std::strtod(c.c_str(), nullptr) == v);
    }
    CHECK_EQ(latex(nombre(Nombre::reel(0.5)) * X), "0.5 x");
    CHECK_EQ(latex(nombre(Nombre::reel(1e-10)) * X), "10^{-10} x");
}

TEST_CASE(latex_noeuds_non_evalues) {
    const ExprPtr integrale = ast_pow(X, 2.0)->integrer();
    CHECK_EQ(latex(integrale), "\\frac{x^{3}}{3}");
    // Primitive qu'aucune règle ne connaît : intégrale non évaluée
    const ExprPtr inconnue = ast_sin(ast_sin(X))->integrer();
    if (inconnue->type() == TypeNoeud::IntegraleNonEvaluee) {
        CHECK_EQ(latex(inconnue), "\\int \\sin\\left(\\sin\\left(x\\right)\\right)\\,dx");
    }
    const ExprPtr limite = (ast_sin(X) / X)->limite(0.0);
    CHECK(limite.get() == cst(1.0).get());
    const ExprPtr nonDefinie = (cst(1.0) / X)->limite(0.0);
    if (nonDefinie->type() == TypeNoeud::LimiteNonEvaluee) {
        CHECK_EQ(latex(nonDefinie), "\\lim_{x \\to 0} \\frac{1}{x}");
    }
}

TEST_CASE(latex_solutions) {
    CHECK_EQ(latex(resoudre(ast_pow(X, 2.0) - 5.0 * X + 6.0)), "x = 2 \\text{ ou } x = 3");
    CHECK_EQ(latex(resoudre(ast_pow(X - 1.0, 2.0))), "x = 1 \\; (\\text{double})");
    CHECK_EQ(latex(resoudre(ast_pow(X - 1.0, 3.0))), "x = 1 \\; (\\text{triple})");
    CHECK_EQ(latex(resoudre(ast_pow(X - 1.0, 5.0))), "x = 1 \\; (\\text{multiplicit\xC3\xA9 } 5)");
    CHECK_EQ(latex(resoudre(ast_sin(X), frac(1, 2))),
             "x = \\frac{\\pi}{6} + 2 k \\pi,\\ k \\in \\mathbb{Z} \\text{ ou } x = \\frac{5 \\pi}{6} + 2 k \\pi,\\ k \\in \\mathbb{Z}");
    CHECK_EQ(latex(resoudre(ast_tan(X), cst(1.0))), "x = \\frac{\\pi}{4} + k \\pi,\\ k \\in \\mathbb{Z}");
    CHECK_EQ(latex(resoudre(ast_pow(X, 2.0) + 1.0)), "x \\in \\emptyset");
    CHECK_EQ(latex(resoudre(X - X)), "x \\in \\mathbb{R}");
    CHECK_EQ(latex(resoudre(ast_pow(X, 2.0) - 2.0)), "x = -\\sqrt{2} \\text{ ou } x = \\sqrt{2}");
    CHECK_EQ(latex(resoudre(ast_exp(X), cst(2.0)), "t"), "t = \\ln\\left(2\\right)");
    // Solution approchée (racine irrationnelle d'un polynôme non résoluble exactement)
    const Solutions cubique = resoudre(ast_pow(X, 3.0) - X - 1.0);
    CHECK_EQ(latex(cubique), "x \\approx 1.3247");
    // Résultat incomplet signalé
    Solutions incomplet;
    incomplet.complet = false;
    CHECK(contient(latex(incomplet), "incomplet"));
    incomplet.liste = resoudre(ast_pow(X, 2.0) - 5.0 * X + 6.0).liste;
    CHECK(contient(latex(incomplet), "x = 2 \\text{ ou } x = 3 \\text{ ou d'autres solutions"));
    // Plusieurs entiers libres
    CHECK(equilibre(latex(resoudre(ast_sin(ast_cos(X)), frac(1, 2)))));
}

TEST_CASE(latex_document) {
    const std::string d = documentLatex("\\[ x^{2} \\]");
    CHECK(contient(d, "\\documentclass{article}"));
    CHECK(contient(d, "\\usepackage[utf8]{inputenc}"));
    CHECK(contient(d, "\\usepackage{amsmath, amssymb}"));
    CHECK(contient(d, "\\begin{document}\n\\[ x^{2} \\]\n\\end{document}"));
}

namespace {

// Expression aléatoire de profondeur bornée, à une variable x et des paramètres
ExprPtr aleatoire(std::mt19937& g, int profondeur) {
    auto entre = [&](int a, int b) { return std::uniform_int_distribution<int>(a, b)(g); };
    if (profondeur <= 0 || entre(0, 9) == 0) {
        switch (entre(0, 7)) {
            case 0: return X;
            case 1: return X;
            case 2: return cst(entre(-9, 9));
            case 3: return frac(entre(-9, 9), entre(1, 9));
            case 4: return param(entre(0, 1) ? "a" : "alpha");
            case 5: return pi();
            case 6: return nombre(Nombre::reel(entre(-50, 50) * 0.37));
            default: return param("C" + std::to_string(entre(1, 3)));
        }
    }
    const ExprPtr a = aleatoire(g, profondeur - 1);
    switch (entre(0, 19)) {
        case 0: case 1: return a + aleatoire(g, profondeur - 1);
        case 2: return a - aleatoire(g, profondeur - 1);
        case 3: case 4: return a * aleatoire(g, profondeur - 1);
        case 5: return a / aleatoire(g, profondeur - 1);
        case 6: return ast_pow(a, entre(-4, 5));
        case 7: return ast_pow(a, frac(entre(-5, 5), entre(2, 4)));
        case 8: return ast_pow(a, aleatoire(g, profondeur - 2));
        case 9: return ast_sin(a);
        case 10: return ast_cos(a);
        case 11: return ast_tan(a);
        case 12: return ast_exp(a);
        case 13: return ast_ln(a);
        case 14: return ast_asin(a);
        case 15: return ast_atan(a);
        case 16: return ast_abs(a);
        case 17: return ast_sign(a);
        case 18: return entre(0, 1) ? ast_sinh(a) : ast_cosh(a);
        default: return -a;
    }
}

} // namespace

TEST_CASE(latex_accolades_equilibrees_aleatoire) {
    std::mt19937 g(20260502);
    for (int i = 0; i < 2000; ++i) {
        const ExprPtr e = aleatoire(g, 5);
        for (int variante = 0; variante < 3; ++variante) {
            OptionsLatex o;
            o.fractions = variante != 1;
            o.racines = variante != 2;
            o.multiplicationExplicite = variante == 1;
            const std::string t = latex(e, o);
            if (!equilibre(t) || t.empty() || contient(t, "nan") || contient(t, "e+") ||
                contient(t, "e-") || contient(t, "\\left()")) {
                CHECK(false);
                std::cerr << "      " << e->texte() << "  ->  " << t << "\n";
                return;
            }
        }
    }
    // Les dérivées successives restent bien formées
    ExprPtr d = ast_exp(ast_sin(X)) * ast_pow(X, 2.0);
    for (int k = 0; k < 5; ++k) {
        d = d->derivee();
        CHECK(equilibre(latex(d)));
    }
}

TEST_CASE(latex_fort_partage) {
    // Dérivée 7e : un graphe très partagé, dont le texte est grand ; le rendu reste linéaire
    ExprPtr d = ast_tan(ast_exp(ast_sin(X))) * ast_ln(X + 1.0);
    for (int k = 0; k < 7; ++k) d = d->derivee();
    std::string t;
    double meilleur = 1e9; // meilleur temps sur 5 rendus (le premier paie les défauts de page)
    for (int i = 0; i < 5; ++i) {
        const auto debut = std::chrono::steady_clock::now();
        t = latex(d);
        meilleur = std::min(meilleur, std::chrono::duration<double>(std::chrono::steady_clock::now() - debut).count());
    }
    CHECK(equilibre(t));
    CHECK(t.size() > 20000);
    CHECK(meilleur < 1.0);
    // Débit : au moins 1 Mo/s (un rendu quadratique en la profondeur ou le partage en serait très loin)
    CHECK(static_cast<double>(t.size()) / std::max(meilleur, 1e-9) > 1e6);
    std::cout << "    7e derivee : " << t.size() << " caracteres en " << meilleur * 1000 << " ms\n";
}

// ============================================================================
// Journal de résolution
// ============================================================================

namespace {

bool aRegle(const Journal& j, const std::string& regle) {
    for (const Etape& e : j.etapes()) {
        if (e.regle == regle) return true;
    }
    return false;
}

std::string visible(const std::vector<Solution>& l) {
    std::string r;
    for (const Solution& s : l) r += s.valeur->texte() + (s.exacte ? "" : "~") + "x" + std::to_string(s.multiplicite) + ";";
    return r;
}

} // namespace

TEST_CASE(journal_inactif_ne_change_rien) {
    CHECK(journal::actif() == nullptr);
    const ExprPtr a = ast_pow(X, 2.0) - 5.0 * X + 6.0;
    const ExprPtr b = ast_exp(2.0 * X) - 3.0 * ast_exp(X) + 2.0;
    const ExprPtr c = ast_sin(X) - frac(1, 2);
    const ExprPtr d = X * ast_ln(X);
    const ExprPtr e = ast_pow(X, 3.0) - X - 1.0;
    const ExprPtr f = ast_abs(X - 1.0) - 2.0;
    for (const ExprPtr& eq : {a, b, c, d, e, f}) {
        const Solutions sans = resoudre(eq);
        Journal j;
        const Solutions avec = resoudre(eq, j);
        CHECK(!j.vide());
        CHECK_EQ(visible(avec.liste), visible(sans.liste));
        CHECK(avec.complet == sans.complet);
        CHECK(avec.toutReel == sans.toutReel);
        CHECK(journal::actif() == nullptr); // l'activation est restaurée
    }
    // Un journal vide tant que rien n'est demandé
    Journal vide;
    resoudre(a);
    CHECK(vide.vide());
    // Intervalle et EquationClassique
    Journal j;
    const Solutions s = resoudreSurIntervalle(ast_sin(X), 0.0, 7.0, j);
    CHECK_EQ(visible(s.liste), visible(resoudreSurIntervalle(ast_sin(X), 0.0, 7.0).liste));
    CHECK(aRegle(j, "resolution.intervalle"));
    Journal k;
    CHECK_EQ(visible(EquationClassique(a).resoudre(k).liste), visible(resoudre(a).liste));
    CHECK(!k.vide());
}

TEST_CASE(journal_polynome_du_second_degre) {
    Journal j;
    resoudre(ast_pow(X, 2.0) - 5.0 * X + 6.0, cst(0.0), j);
    CHECK(aRegle(j, "resolution.equation"));
    CHECK(aRegle(j, "resolution.polynome"));
    CHECK(aRegle(j, "resolution.discriminant"));
    CHECK(aRegle(j, "resolution.racine"));
    CHECK(aRegle(j, "resolution.solution"));
    const std::string t = latex(j);
    CHECK(equilibre(t));
    CHECK(contient(t, "\\begin{itemize}"));
    CHECK(contient(t, "\\item \xC3\x89quation \xC3\xA0 r\xC3\xA9soudre : $\\displaystyle x^{2} - 5 x + 6 = 0$"));
    CHECK(contient(t, "Discriminant $\\Delta$ = b$^{2}$ $-$ 4ac, avec a = 1, b = -5, c = 6"));
    CHECK(aRegle(j, "resolution.discriminant_signe"));
    CHECK(contient(t, "\\Delta = 1$"));
    CHECK(contient(t, "$\\displaystyle x = 2$"));
    CHECK(contient(t, "$\\displaystyle x = 3$"));
}

TEST_CASE(journal_changement_de_variable) {
    Journal j;
    const Solutions s = resoudre(ast_exp(2.0 * X) - 3.0 * ast_exp(X) + 2.0, cst(0.0), j);
    CHECK_EQ(s.liste.size(), static_cast<std::size_t>(2));
    CHECK(aRegle(j, "resolution.regroupement_exponentielles"));
    CHECK(aRegle(j, "resolution.changement_variable"));
    CHECK(aRegle(j, "resolution.polynome_en_t"));
    CHECK(aRegle(j, "resolution.retour_variable"));
    CHECK(aRegle(j, "resolution.logarithme"));
    bool retrait = false;
    for (const Etape& e : j.etapes()) retrait = retrait || e.niveau > 0;
    CHECK(retrait);
    const std::string t = latex(j);
    CHECK(equilibre(t));
    CHECK(contient(t, "t = e^{x}"));
    CHECK(contient(t, "t^{2} - 3 t + 2 = 0"));
    CHECK(contient(t, "e^{x} = 2"));
    CHECK(contient(t, "x = \\ln\\left(2\\right)"));
    CHECK(contient(t, "\n  \\begin{itemize}\n    \\item")); // sous-niveau imbriqué
}

TEST_CASE(journal_famille_trigonometrique) {
    Journal j;
    resoudre(ast_sin(X), frac(1, 2), j);
    CHECK(aRegle(j, "resolution.famille_trigonometrique"));
    const std::string t = latex(j);
    CHECK(equilibre(t));
    CHECK(contient(t, "\\sin\\left(x\\right) = \\frac{1}{2}"));
    CHECK(contient(t, "x = \\frac{\\pi}{6} + 2 k \\pi"));
    CHECK(contient(t, "x = \\frac{5 \\pi}{6} + 2 k \\pi"));
}

TEST_CASE(journal_solution_rejetee_hors_domaine) {
    Journal j;
    const Solutions s = resoudre(X * ast_ln(X), cst(0.0), j);
    CHECK_EQ(s.liste.size(), static_cast<std::size_t>(1));
    CHECK(aRegle(j, "resolution.produit_nul"));
    CHECK(aRegle(j, "resolution.solution_rejetee"));
    const std::string t = latex(j);
    CHECK(equilibre(t));
    CHECK(contient(t, "hors du domaine"));
    CHECK(contient(t, "$\\displaystyle x = 0$"));
    // La solution restante est x = 1
    CHECK(contient(t, "Solution : $\\displaystyle x = 1$"));
    int rejetees = 0;
    for (const Etape& e : j.etapes()) rejetees += e.regle == "resolution.solution_rejetee";
    CHECK_EQ(rejetees, 1);
}

TEST_CASE(journal_domaine_et_incomplet) {
    Journal j;
    resoudre(ast_exp(X), cst(-1.0), j);
    CHECK(aRegle(j, "resolution.hors_domaine"));
    CHECK(aRegle(j, "resolution.bilan"));
    Journal k;
    const Solutions s = resoudre(X + ast_sin(X) * ast_exp(X), cst(1.0), k);
    CHECK(!s.complet);
    CHECK(aRegle(k, "resolution.incomplet"));
    CHECK(equilibre(latex(k)));
    // Profondeur d'imbrication limitée à 4 listes
    Journal profond;
    for (int i = 0; i < 9; ++i) {
        profond.ajouter("test", "niveau " + std::to_string(i));
        profond.entrer();
    }
    const std::string t = latex(profond);
    CHECK_EQ(std::count(t.begin(), t.end(), '\n') > 0, true);
    std::size_t ouvertes = 0, maximum = 0;
    for (std::size_t p = 0; (p = t.find("\\begin{itemize}", p)) != std::string::npos; p += 5) maximum = std::max(maximum, ++ouvertes);
    CHECK_EQ(maximum, static_cast<std::size_t>(4));
    CHECK(equilibre(t));
}

TEST_CASE(journal_echappement_du_texte) {
    Journal j;
    j.ajouter("test", "50 % de x_1 & y # {z} < 3", X, X + 1.0);
    const std::string t = latex(j);
    CHECK(contient(t, "50 \\% de x\\_1 \\& y \\# \\{z\\} $<$ 3"));
    CHECK(contient(t, "$\\displaystyle x \\;\\longrightarrow\\; x + 1$"));
    CHECK(equilibre(t));
}

TEST_CASE(latex_expression_nulle) {
    CHECK_THROWS(latex(ExprPtr()), std::invalid_argument);
}

TEST_CASE(journal_racine_approchee_notee_environ) {
    Journal j;
    resoudre(ast_pow(X, 3.0) - X - 1.0, cst(0.0), j);
    const std::string t = latex(j);
    CHECK(equilibre(t));
    CHECK(contient(t, "$\\displaystyle x \\approx 1.32"));
    // Une racine exacte reste notée avec « = »
    Journal k;
    resoudre(ast_pow(X, 2.0) - 5.0 * X + 6.0, cst(0.0), k);
    CHECK(!contient(latex(k), "\\approx"));
}

TEST_CASE(latex_noms_sans_lettre_initiale) {
    // « _1 » (lu comme paramètre) n'est pas le nombre 1 : x + _1 ne doit pas s'écrire x + 1
    CHECK(latex(lire("_1 + x")) != latex(lire("1 + x")));
    CHECK_EQ(latex(param("_1")), "\\mathrm{\\_1}");
    CHECK_EQ(latex(param("2a") * X), "\\mathrm{2a} x");
    // Un paramètre « _a » n'est pas le paramètre « a »
    CHECK(latex(lire("_a*x")) != latex(lire("a*x")));
    // Variable interne du solveur : toujours écrite t
    CHECK_EQ(latex(var("_t")), "t");
}

int main() { return test::executerTous(); }
