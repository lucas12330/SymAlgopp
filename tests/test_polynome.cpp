/**
 * @file test_polynome.cpp
 * @brief Tests du développement des expressions et des polynômes exacts.
 */

#include <cmath>
#include <string>
#include <vector>

#include "ASTNode.hpp"
#include "Polynome.hpp"
#include "test_framework.hpp"

using namespace symalgo;

namespace {

const ExprPtr X = var("x");

// La forme développée doit prendre les mêmes valeurs que l'expression d'origine
void verifierMemesValeurs(const ExprPtr& e) {
    const ExprPtr d = developper(e);
    for (double x : {-1.7, -0.3, 0.4, 1.9, 3.2}) {
        const double attendu = e->eval(x);
        CHECK_NEAR(d->eval(x), attendu, 1e-9 * std::max(1.0, std::abs(attendu)));
    }
}

} // namespace

TEST_CASE(developper_identites_remarquables) {
    CHECK_EQ(developper(ast_pow(X + 1.0, 2.0))->texte(), std::string("x^2 + 2*x + 1"));
    CHECK_EQ(developper((X + 1.0) * (X - 1.0))->texte(), std::string("x^2 - 1"));
    CHECK_EQ(developper(ast_pow(X - 2.0, 3.0))->texte(), std::string("x^3 - 6*x^2 + 12*x - 8"));
    CHECK_EQ(developper((cst(2.0) * X + 3.0) * (X - 4.0))->texte(), std::string("2*x^2 - 5*x - 12"));
}

TEST_CASE(developper_coefficients_binomiaux) {
    // (x + 1)^10 : coefficient de x^5 = C(10, 5) = 252
    const ExprPtr d = developper(ast_pow(X + 1.0, 10.0));
    const Somme* s = comme<Somme>(d);
    CHECK(s != nullptr);
    CHECK_EQ(s->getTermes().size(), std::size_t(10)); // x^1..x^10, constante à part
    bool trouve = false;
    for (const Terme& t : s->getTermes()) {
        if (t.expression.get() == ast_pow(X, 5.0).get()) {
            CHECK_EQ(t.coefficient.texte(), std::string("252"));
            trouve = true;
        }
    }
    CHECK(trouve);
    // Coefficients exacts : (x/3 + 1/2)^2 = x^2/9 + x/3 + 1/4
    CHECK_EQ(developper(ast_pow(X / 3.0 + frac(1, 2), 2.0))->texte(), std::string("x^2/9 + x/3 + 1/4"));
}

TEST_CASE(developper_en_profondeur) {
    CHECK_EQ(developper(ast_sin(ast_pow(X + 1.0, 2.0)))->texte(), std::string("sin(x^2 + 2*x + 1)"));
    CHECK_EQ(developper(cst(1.0) / ast_pow(X + 1.0, 2.0))->texte(), std::string("1/(x^2 + 2*x + 1)"));
    // Les numérateurs sont distribués, le dénominateur reste un facteur
    CHECK_EQ(developper((X + 1.0) * ast_exp(X))->texte(), std::string("x*exp(x) + exp(x)"));
}

TEST_CASE(developper_conserve_les_valeurs) {
    verifierMemesValeurs(ast_pow(X + 1.0, 2.0) * (X - 3.0));
    verifierMemesValeurs(ast_pow(cst(2.0) * X - 1.0, 5.0) - ast_pow(X, 5.0));
    verifierMemesValeurs((X + ast_sin(X)) * (X - ast_cos(X)));
    verifierMemesValeurs(ast_pow(X + 1.0, -2.0) * (X + 2.0));
    verifierMemesValeurs(ast_exp(X) * ast_pow(X + ast_exp(X), 3.0));
}

TEST_CASE(developper_expression_deja_developpee) {
    const ExprPtr e = ast_pow(X, 2.0) + 3.0 * X + 1.0;
    CHECK(developper(e).get() == e.get()); // hash-consing : même noeud
    CHECK(developper(cst(5.0)).get() == cst(5.0).get());
}

namespace {

Polynome polynome(const ExprPtr& e) {
    Polynome p;
    const bool ok = Polynome::depuisExpression(developper(e), p);
    CHECK(ok);
    return p;
}

std::string texteRacines(const std::vector<RacineReelle>& racines) {
    std::string t;
    for (const RacineReelle& r : racines) {
        if (!t.empty()) t += " ; ";
        t += r.exacte ? r.exacte->texte() : std::string("~");
        if (r.multiplicite > 1) t += " (x" + std::to_string(r.multiplicite) + ")";
    }
    return t;
}

} // namespace

TEST_CASE(polynome_reconnaissance) {
    Polynome p;
    ExprPtr x;
    CHECK(Polynome::depuisExpression(developper(ast_pow(X + 1.0, 3.0)), p, &x));
    CHECK_EQ(p.degre(), 3);
    CHECK_EQ(p.coefficient(1).texte(), std::string("3"));
    CHECK(x.get() == X.get());
    CHECK(Polynome::depuisExpression(cst(7.0), p) && p.degre() == 0);
    CHECK(!Polynome::depuisExpression(ast_sin(X), p));
    CHECK(!Polynome::depuisExpression(cst(1.0) / X, p));
    CHECK(!Polynome::depuisExpression(X + var("y"), p)); // deux variables
    CHECK(!Polynome::depuisExpression(param("a") * X, p)); // coefficient symbolique
    CHECK(Polynome::depuisExpression(cst(0.5) * X, p) && !p.estExact());
}

TEST_CASE(polynome_arithmetique_exacte) {
    const Polynome a = polynome(ast_pow(X, 3.0) - 1.0), b = polynome(X - 1.0);
    Polynome q, r;
    a.diviser(b, q, r);
    CHECK_EQ(q.versExpression(X)->texte(), std::string("x^2 + x + 1"));
    CHECK(r.estNul());
    CHECK_EQ((a * b).versExpression(X)->texte(), std::string("x^4 - x^3 - x + 1"));
    CHECK_EQ(a.derivee().versExpression(X)->texte(), std::string("3*x^2"));
    // pgcd((x-1)^2 (x+2), (x-1)(x+3)) = x - 1
    const Polynome g = Polynome::pgcd(polynome(ast_pow(X - 1.0, 2.0) * (X + 2.0)), polynome((X - 1.0) * (X + 3.0)));
    CHECK_EQ(g.versExpression(X)->texte(), std::string("x - 1"));
    CHECK_EQ(polynome(ast_pow(X, 2.0) / 3.0 + frac(1, 7)).evaluer(Nombre(2)).texte(), std::string("31/21"));
}

TEST_CASE(polynome_sans_carre) {
    // (x - 1)^3 (x + 2)^2 (x - 5) : facteurs de multiplicités 1, 2 et 3
    const auto facteurs = polynome(ast_pow(X - 1.0, 3.0) * ast_pow(X + 2.0, 2.0) * (X - 5.0)).sansCarre();
    CHECK_EQ(facteurs.size(), std::size_t(3));
    CHECK_EQ(facteurs[0].first.versExpression(X)->texte(), std::string("x - 5"));
    CHECK_EQ(facteurs[1].first.versExpression(X)->texte(), std::string("x + 2"));
    CHECK_EQ(facteurs[1].second, 2);
    CHECK_EQ(facteurs[2].first.versExpression(X)->texte(), std::string("x - 1"));
    CHECK_EQ(facteurs[2].second, 3);
}

TEST_CASE(polynome_comptage_de_sturm) {
    const Polynome p = polynome(ast_pow(X, 5.0) - X); // racines -1, 0, 1
    CHECK_EQ(p.nombreRacinesReelles(Nombre(-10), Nombre(10)), 3);
    CHECK_EQ(p.nombreRacinesReelles(Nombre(0), Nombre(10)), 1);   // ]0, 10] : seulement 1
    CHECK_EQ(p.nombreRacinesReelles(Nombre(-1), Nombre(0)), 1);   // ]-1, 0] : seulement 0
    CHECK_EQ(polynome(ast_pow(X, 2.0) + 1.0).nombreRacinesReelles(Nombre(-100), Nombre(100)), 0);
    CHECK_EQ(polynome(ast_pow(X - 1.0, 4.0)).nombreRacinesReelles(Nombre(0), Nombre(2)), 1); // distinctes
}

TEST_CASE(racines_exactes) {
    CHECK_EQ(texteRacines(polynome(ast_pow(X, 2.0) - 3.0 * X + 2.0).racinesReelles()), std::string("1 ; 2"));
    CHECK_EQ(texteRacines(polynome(cst(6.0) * ast_pow(X, 3.0) - cst(11.0) * ast_pow(X, 2.0) + cst(6.0) * X - 1.0).racinesReelles()),
             std::string("1/3 ; 1/2 ; 1"));
    CHECK_EQ(texteRacines(polynome(ast_pow(X - 1.0, 3.0) * ast_pow(X + 2.0, 2.0)).racinesReelles()),
             std::string("-2 (x2) ; 1 (x3)"));
    // Nombre d'or : racines de degré 2 par radicaux
    CHECK_EQ(texteRacines(polynome(ast_pow(X, 2.0) - X - 1.0).racinesReelles()),
             std::string("-5^(1/2)/2 + 1/2 ; 5^(1/2)/2 + 1/2"));
    CHECK(polynome(ast_pow(X, 2.0) + 1.0).racinesReelles().empty());
}

TEST_CASE(racines_numeriques_certifiees) {
    // x^3 - 2 : racine irrationnelle de degré 3, au double le plus proche
    const auto r = polynome(ast_pow(X, 3.0) - 2.0).racinesReelles();
    CHECK_EQ(r.size(), std::size_t(1));
    CHECK(!r[0].exacte);
    CHECK_EQ(r[0].valeur, std::cbrt(2.0)); // au double près
    CHECK(r[0].gauche < r[0].droite);
    // Polynôme de Wilkinson (degré 10, très mal conditionné) : les 10 racines exactes
    ExprPtr w = cst(1.0);
    for (int k = 1; k <= 10; ++k) w = w * (X - cst(k));
    CHECK_EQ(texteRacines(polynome(w).racinesReelles()), std::string("1 ; 2 ; 3 ; 4 ; 5 ; 6 ; 7 ; 8 ; 9 ; 10"));
    // Coefficient réel : racines certifiées, sans forme exacte
    const auto s = polynome(ast_pow(X, 2.0) - cst(0.1)).racinesReelles();
    CHECK_EQ(s.size(), std::size_t(2));
    CHECK(!s[0].exacte && !s[1].exacte);
    CHECK_EQ(s[1].valeur, std::sqrt(0.1)); // racine de la valeur exacte du double 0.1
    CHECK_EQ(s[0].valeur, -std::sqrt(0.1));
}

TEST_CASE(factorisation) {
    CHECK_EQ(factoriser(ast_pow(X, 3.0) - X)->texte(), std::string("x*(x - 1)*(x + 1)"));
    CHECK_EQ(factoriser(ast_pow(X, 5.0) - X)->texte(), std::string("x*(x - 1)*(x + 1)*(x^2 + 1)"));
    CHECK_EQ(factoriser(cst(6.0) * ast_pow(X, 3.0) - cst(11.0) * ast_pow(X, 2.0) + cst(6.0) * X - 1.0)->texte(),
             std::string("(x - 1)*(2*x - 1)*(3*x - 1)"));
    CHECK_EQ(factoriser(developper(ast_pow(X - 1.0, 3.0) * ast_pow(X + 2.0, 2.0)))->texte(),
             std::string("(x - 1)^3*(x + 2)^2"));
    CHECK_EQ(factoriser(cst(2.0) * ast_pow(X, 2.0) - 2.0)->texte(), std::string("2*(x - 1)*(x + 1)"));
    CHECK_EQ(factoriser(ast_pow(X, 2.0) - X - 1.0)->texte(), std::string("x^2 - x - 1")); // irréductible sur Q
    // Non polynomial ou approché : inchangé
    const ExprPtr e = ast_sin(X) + 1.0;
    CHECK(factoriser(e).get() == e.get());
    const ExprPtr f = ast_pow(X, 2.0) - cst(0.25);
    CHECK(factoriser(f).get() == f.get());
    // La factorisation redonne bien le polynôme
    const ExprPtr g = cst(6.0) * ast_pow(X, 4.0) - cst(5.0) * ast_pow(X, 3.0) - cst(2.0) * X;
    CHECK(developper(factoriser(g)).get() == developper(g).get());
}

int main() { return test::executerTous(); }
