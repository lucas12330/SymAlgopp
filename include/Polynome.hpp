/**
 * @file Polynome.hpp
 * @brief Développement des expressions et polynômes à coefficients exacts.
 */

#pragma once

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : developper
 * Description : Développe l'expression : distribue les produits sur les sommes et
 *               développe les puissances entières positives de sommes, à tous les niveaux
 *               (y compris dans les arguments des fonctions). Le résultat est canonique :
 *               les termes semblables sont regroupés. Les dénominateurs sont développés
 *               mais restent au dénominateur.
 * Utilisation : developper(ast_pow(x + 1.0, 2.0))  // x^2 + 2*x + 1
 */
ExprPtr developper(const ExprPtr& e);

} // namespace symalgo
