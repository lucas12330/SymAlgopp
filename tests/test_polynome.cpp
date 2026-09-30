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

int main() { return test::executerTous(); }
