/**
 * @file Primitives.hpp
 * @brief Outils internes de l'intégration symbolique (src/Primitives.cpp), partagés avec
 *        les règles par noeud de src/Regles.cpp. En-tête interne : il ne fait pas partie de
 *        l'API publique (symalgopp.hpp).
 */

#pragma once

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : coefficientLineaire
 * Description : Si u est une fonction affine de sa variable (u = a*x + b avec a non nul,
 *               a pouvant être symbolique comme pi ou C1), renvoie vrai et écrit a : la
 *               dérivée de u ne dépend alors pas de x. Sert à la règle F(ax+b)/a.
 */
bool coefficientLineaire(const ExprPtr& u, ExprPtr& a);

/*
 * Nom : variableDe
 * Description : Le noeud Variable de l'expression (la première rencontrée en descendant),
 *               ou un pointeur nul si elle n'en contient pas. Coût : la profondeur de
 *               l'expression (une seule branche est suivie).
 */
ExprPtr variableDe(const ASTNode& e);

/*
 * Nom : estPositifOuNul
 * Description : Vrai si u >= 0 partout où u est défini, reconnu sur sa forme seulement
 *               (constante positive, pi, exp, cosh, abs, puissance paire, somme de termes
 *               positifs, trinôme du second degré sans racine réelle...). Faux ne signifie
 *               pas que u prend des valeurs négatives.
 */
bool estPositifOuNul(const ExprPtr& u);

/*
 * Nom : lnAbs
 * Description : ln|u|, écrit ln(u) quand u est manifestement positif (estPositifOuNul) :
 *               ln(x^2 + 1), ln(exp(x) + 1), mais ln(abs(x)).
 */
ExprPtr lnAbs(const ExprPtr& u);

} // namespace symalgo
