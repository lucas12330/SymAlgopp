/**
 * @file Integration.hpp
 * @brief Intégrale définie d'une expression à une variable : valeur exacte quand la
 *        primitive le permet, quadrature adaptative sinon.
 *
 * Chemin exact : F = f->integrer(), puis F(b) - F(a) par substitution (bornes rationnelles,
 * pi, e...). Il n'est pris que si F et f sont certifiées sans singularité sur [a, b]
 * (dénominateurs, ln, puissances non entières, tan, asin/acos, vérifiés avec le solveur) :
 * dans le doute, on bascule en numérique plutôt que de renvoyer une valeur exacte fausse
 * (l'intégrale de 1/x^2 sur [-1, 1] n'est pas -2).
 *
 * Chemin numérique : Gauss-Kronrod G7-K15 adaptatif (bissection des intervalles de plus
 * grande erreur, évaluation compilée par blocs), bornes infinies par changement de variable,
 * coupure aux singularités intérieures repérées. Une intégrale divergente est rendue avec
 * une valeur NaN et une erreur infinie.
 */

#pragma once

#include <cmath>
#include <limits>

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : OptionsIntegration
 * Description : Réglages de l'intégration. L'arrêt a lieu quand l'erreur estimée est
 *               inférieure à max(toleranceAbsolue, toleranceRelative * |I|), ou après
 *               subdivisionsMax bissections (le résultat porte alors son erreur estimée).
 *               exacteSiPossible = faux force le chemin numérique.
 */
struct OptionsIntegration {
    double toleranceRelative = 1e-10;
    double toleranceAbsolue = 1e-12;
    int subdivisionsMax = 200;
    bool exacteSiPossible = true;
};

/*
 * Nom : ResultatIntegrale
 * Description : exacte est la valeur symbolique (1/3, pi/4, ln(2)...), nulle si elle n'a pas
 *               été obtenue ; valeur est son évaluation, ou l'approximation numérique (NaN
 *               si l'intégrale diverge, ou si une valeur exacte contient des paramètres) ;
 *               erreurEstimee vaut 0 pour une valeur exacte, +infini si l'intégrale diverge.
 */
struct ResultatIntegrale {
    ExprPtr exacte;
    double valeur = std::numeric_limits<double>::quiet_NaN();
    double erreurEstimee = 0.0;

    bool estExacte() const { return static_cast<bool>(exacte); }
    bool estFinie() const { return std::isfinite(valeur); }
};

/*
 * Nom : integrerDefinie
 * Description : Intégrale de f sur [a, b] (a > b donne l'opposé de l'intégrale sur [b, a],
 *               a == b donne 0 exact). Les bornes ExprPtr doivent être constantes (pi(),
 *               frac(1, 2), param("p")...) ; un paramètre n'autorise que le résultat exact
 *               symbolique. Avec des bornes double, ±infini est permis (le résultat est alors
 *               numérique) et le chemin exact n'est tenté que pour des bornes entières.
 *               Lève std::invalid_argument pour une expression à plusieurs variables, une
 *               borne contenant une variable ou NaN ; std::domain_error si l'intégrale ne
 *               peut être calculée que numériquement et que f contient un paramètre ou un
 *               noeud non évaluable.
 * Utilisation : ResultatIntegrale r = integrerDefinie(ast_sin(var("x")), cst(0.0), pi());
 *               r.exacte  // 2
 *               ResultatIntegrale g = integrerDefinie(ast_exp(-var("x") * var("x")), -INFINITY, INFINITY);
 *               g.valeur  // sqrt(pi) = 1.7724538509...
 */
ResultatIntegrale integrerDefinie(const ExprPtr& f, const ExprPtr& a, const ExprPtr& b,
                                  const OptionsIntegration& options = {});
ResultatIntegrale integrerDefinie(const ExprPtr& f, double a, double b, const OptionsIntegration& options = {});

} // namespace symalgo
