/**
 * @file test_lecture.cpp
 * @brief Tests de la lecture d'expressions : nombres exacts, précédences, multiplication
 *        implicite, fonctions, paramètres, Unicode, aller-retour avec l'affichage, équations
 *        et messages d'erreur.
 */

#include <string>
#include <vector>

#include <symalgopp> // l'en-tête unique suffit (et vérifie qu'il compile seul)

#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");

// Même noeud (hash-consing) : égalité exacte de la forme canonique
bool pareil(const std::string& texte, const ExprPtr& attendu) { return lire(texte).get() == attendu.get(); }

// Position signalée pour un texte invalide (0 si aucune erreur)
std::size_t positionErreur(const std::string& texte) {
    try {
        lireEquation(texte);
    } catch (const ErreurLecture& e) {
        return e.position();
    }
    return 0;
}

} // namespace

TEST_CASE(nombres_exacts) {
    CHECK(pareil("42", frac(42)));
    CHECK(pareil("0.1", frac(1, 10)));
    CHECK(pareil(".5", frac(1, 2)));
    CHECK(pareil("1.", frac(1)));
    CHECK(pareil("2.5e-3", frac(1, 400)));
    CHECK(pareil("1E3", frac(1000)));
    CHECK(pareil("12.50e+1", frac(125)));
    CHECK_EQ(lire("123456789012345678901234567890")->texte(), std::string("123456789012345678901234567890"));
    CHECK_EQ(lire("0.1 + 0.2")->texte(), std::string("3/10")); // exact, contrairement au double
    CHECK(lire("0.1")->getValeurConstante() == 0.1);
}

TEST_CASE(precedences) {
    CHECK(pareil("1 + 2*3", frac(7)));
    CHECK(pareil("(1 + 2)*3", frac(9)));
    CHECK(pareil("2^3^2", frac(512)));   // associative à droite
    CHECK(pareil("-x^2", -ast_pow(X, 2.0)));
    CHECK(pareil("2^-1", frac(1, 2)));
    CHECK(pareil("8/4/2", frac(1)));     // associative à gauche
    CHECK(pareil("10 - 4 - 3", frac(3)));
    CHECK(pareil("--x", X));
    CHECK(pareil("+x", X));
    CHECK(pareil("2**3", frac(8)));
    CHECK(pareil("x ** 2", ast_pow(X, 2.0)));
    CHECK(pareil("x - y - z", X - param("y") - param("z")));
    CHECK(pareil("2*x^2 + 3*x - 1", 2.0 * ast_pow(X, 2.0) + 3.0 * X - 1.0));
}

TEST_CASE(multiplication_implicite) {
    CHECK(pareil("2x", 2.0 * X));
    CHECK(pareil("2x^2", 2.0 * ast_pow(X, 2.0)));
    CHECK(pareil("2(x + 1)", 2.0 * X + 2.0));
    CHECK(pareil("(x + 1)(x - 1)", (X + 1.0) * (X - 1.0)));
    CHECK(pareil("x sin(x)", X * ast_sin(X)));
    CHECK(pareil("3 sin(x) cos(x)", 3.0 * ast_sin(X) * ast_cos(X)));
    CHECK(pareil("1/2x", X / 2.0));      // même précédence que *
    CHECK(pareil("2pi", 2.0 * pi()));
    CHECK(pareil("2e", 2.0 * ast_exp(frac(1))));
    CHECK(pareil("a (x + 1)", param("a") * (X + 1.0)));
    CHECK(pareil("x(x + 1)", X * (X + 1.0)));
}

TEST_CASE(fonctions_et_constantes) {
    CHECK(pareil("sin(x) + cos(x) + tan(x)", ast_sin(X) + ast_cos(X) + ast_tan(X)));
    CHECK(pareil("exp(x)*ln(x)", ast_exp(X) * ast_ln(X)));
    CHECK(pareil("log(x)", ast_ln(X)));
    CHECK(pareil("asin(x) + acos(x) + atan(x)", ast_asin(X) + ast_acos(X) + ast_atan(X)));
    CHECK(pareil("arcsin(x) + arccos(x) + arctan(x)", ast_asin(X) + ast_acos(X) + ast_atan(X)));
    CHECK(pareil("sqrt(x)", ast_pow(X, frac(1, 2))));
    CHECK(pareil("sqrt(8)", 2.0 * ast_pow(cst(2.0), frac(1, 2))));
    CHECK(pareil("sin(pi/6)", frac(1, 2)));
    CHECK(pareil("acos(1/2)", pi() / 3.0));
    CHECK(pareil("e", ast_exp(frac(1))));
    CHECK(pareil("exp(sin(x)^2)", ast_exp(ast_pow(ast_sin(X), 2.0))));
    CHECK(pareil("ln(8)", 3.0 * ast_ln(cst(2.0))));
}

TEST_CASE(variable_et_parametres) {
    const ExprPtr a = param("a"), b = param("b");
    CHECK(pareil("a*x + b", a * X + b));
    CHECK(lire("a*x^2")->derivee().get() == (2.0 * a * X).get());
    CHECK(lire("C1*cos(2*x)")->derivee().get() == (-2.0 * param("C1") * ast_sin(2.0 * X)).get());

    const ExprPtr t = lire("t^2 - x", {"t"});
    CHECK(t.get() == (ast_pow(var("t"), 2.0) - param("x")).get());
    CHECK_NEAR(t->derivee()->eval(3.0), 6.0, 1e-15);

    // Les noms de paramètres peuvent contenir chiffres et soulignés
    CHECK(pareil("k_1 + k2", param("k_1") + param("k2")));
}

TEST_CASE(symboles_unicode) {
    CHECK(pareil("2×x", 2.0 * X));
    CHECK(pareil("2·x", 2.0 * X));
    CHECK(pareil("6÷3", frac(2)));
    CHECK(pareil("−x", -X));
    CHECK(pareil("x² + x³", ast_pow(X, 2.0) + ast_pow(X, 3.0)));
    CHECK(pareil("π", pi()));
    CHECK(pareil("2π", 2.0 * pi()));
    CHECK(pareil("sin(x)²", ast_pow(ast_sin(X), 2.0)));
}

TEST_CASE(aller_retour_avec_affichage) {
    const ExprPtr k = param("k");
    const std::vector<ExprPtr> expressions = {
        ast_pow(X, 3.0) / 3.0 + 5.0 * ast_pow(X, 2.0) / 2.0 + 6.0 * X,
        ast_sin(X) / X,
        1.0 / (X + 1.0),
        -X - 1.0,
        X - ast_pow(cst(2.0), frac(1, 2)),
        ast_pow(cst(2.0), frac(1, 2)) / 2.0,
        ast_pow(cst(-8.0), frac(1, 3)),
        ast_pow(X, -X),
        ast_pow(X, frac(-1, 2)),
        pi() / 6.0 + 2.0 * pi() * k,
        ast_exp(2.0 * X) * ast_sin(3.0 * X),
        ast_pow(ast_ln(X), 2.0),
        ast_atan(X / 2.0),
        param("C1") * ast_cos(2.0 * X) + param("C2") * ast_sin(2.0 * X),
        (param("C1") + param("C2") * X) * ast_exp(-X),
        nombre(Nombre(2).puissanceEntiere(100) / Nombre(3)) * X,
        ast_pow(X, 2.0) * ast_pow(X + 1.0, -2.0) - ast_ln(ast_cos(X)),
        factoriser(ast_pow(X, 3.0) - ast_pow(X, 2.0) - 2.0 * X + 2.0),
        developper(ast_pow(X + 2.0, 5.0)),
        (ast_exp(ast_sin(X)) * ast_pow(X, 2.0))->derivee()->derivee()->derivee()->derivee(),
        ast_asin(X)->derivee(),
        (ast_pow(X, 2.0) + ast_cos(3.0 * X) + ast_exp(X / 2.0))->integrer(),
    };
    for (const ExprPtr& e : expressions) {
        const std::string texte = e->texte();
        const ExprPtr relue = lire(texte, {"x"});
        CHECK(relue.get() == e.get());
        if (relue.get() != e.get()) std::cerr << "      " << texte << "  ->  " << relue->texte() << "\n";
    }
    // Les solutions d'une équation se relisent aussi
    for (const Solution& s : resoudre(ast_sin(X), frac(1, 2)).liste) {
        CHECK(lire(s.valeur->texte()).get() == s.valeur.get());
    }
}

TEST_CASE(equations) {
    const EgaliteLue eq = lireEquation("x^2 = 2");
    CHECK(eq.gauche.get() == ast_pow(X, 2.0).get());
    CHECK(eq.droite.get() == frac(2).get());
    CHECK(lireEquation("x + 1").droite.get() == frac(0).get());

    const Solutions s = resoudre(eq.gauche, eq.droite);
    CHECK_EQ(s.liste.size(), 2u);
    CHECK_EQ(s.liste[1].valeur->texte(), std::string("2^(1/2)"));

    const Solutions familles = EquationClassique("sin(x) = 1/2").resoudre();
    CHECK_EQ(familles.liste.size(), 2u);
    CHECK_EQ(familles.liste[0].valeur->texte(), std::string("pi/6 + 2*pi*k"));

    const Solutions surIntervalle = EquationClassique("exp(2x) - 3exp(x) + 2 = 0").resoudre(-10.0, 10.0);
    CHECK_EQ(surIntervalle.liste.size(), 2u);

    CHECK_NEAR(EquationClassique("x^2 = 2x").eval(3.0), 3.0, 1e-15);

    using namespace symalgo::litteraux;
    CHECK("x^2 + 1"_expr.get() == (ast_pow(X, 2.0) + 1.0).get());
}

TEST_CASE(erreurs) {
    CHECK_THROWS(lire(""), ErreurLecture);
    CHECK_THROWS(lire("   "), ErreurLecture);
    CHECK_THROWS(lire("2 +"), ErreurLecture);
    CHECK_THROWS(lire("(x + 1"), ErreurLecture);
    CHECK_THROWS(lire("x + 1)"), ErreurLecture);
    CHECK_THROWS(lire("()"), ErreurLecture);
    CHECK_THROWS(lire("sin x"), ErreurLecture);
    CHECK_THROWS(lire("sin()"), ErreurLecture);
    CHECK_THROWS(lire("sinh(x)"), ErreurLecture);
    CHECK_THROWS(lire("f(x)"), ErreurLecture);
    CHECK_THROWS(lire("2 3"), ErreurLecture);
    CHECK_THROWS(lire("x = 2"), ErreurLecture);
    CHECK_THROWS(lireEquation("x = 1 = 2"), ErreurLecture);
    CHECK_THROWS(lireEquation("x ="), ErreurLecture);
    CHECK_THROWS(lireEquation("= 2"), ErreurLecture);
    CHECK_THROWS(lire("x $ 2"), ErreurLecture);
    CHECK_THROWS(lire("1e99999"), ErreurLecture);
    CHECK_THROWS(lire("0." + std::string(5000, '1')), ErreurLecture);
    CHECK_THROWS(lire("2^"), ErreurLecture);
    CHECK_THROWS(EquationClassique("x +"), std::invalid_argument); // ErreurLecture en dérive

    // Une imbrication démesurée est refusée proprement, sans épuiser la pile
    CHECK_THROWS(lire(std::string(100000, '(') + "x" + std::string(100000, ')')), ErreurLecture);
    CHECK_THROWS(lire("x^" + std::string(100000, '(') + "x" + std::string(100000, ')')), ErreurLecture);
    CHECK(pareil(std::string(100000, '-') + "x", X));  // les signes se lisent sans récursion
    CHECK(pareil("2*" + std::string(100001, '-') + "x", -2.0 * X));
    CHECK(pareil(std::string(200, '(') + "x" + std::string(200, ')'), X));

    // Positions (en caractères, à partir de 1)
    CHECK_EQ(positionErreur("x + $"), 5u);
    CHECK_EQ(positionErreur("π + $"), 5u);        // π compte pour un caractère
    CHECK_EQ(positionErreur("(x + 1"), 7u);       // fin du texte
    CHECK_EQ(positionErreur("x + 1)"), 6u);
    CHECK_EQ(positionErreur("2 + sinh(x)"), 5u);
    CHECK_EQ(positionErreur("x = 1 = 2"), 7u);

    // Le message montre le texte et un repère sous la position
    try {
        lire("x + * 2");
        CHECK(false);
    } catch (const ErreurLecture& e) {
        CHECK_EQ(e.position(), 5u);
        CHECK_EQ(std::string(e.what()), std::string("lecture, position 5 : expression attendue\n  x + * 2\n      ^"));
    }
}

int main() { return test::executerTous(); }
