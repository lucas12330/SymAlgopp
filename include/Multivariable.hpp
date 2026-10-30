/**
 * @file Multivariable.hpp
 * @brief Fonctions de plusieurs variables : variables d'une expression, dérivées partielles,
 *        gradient, jacobienne, hessienne, laplacien, divergence et évaluation numérique.
 *
 * Toute `var("nom")` est une variable distincte ; les paramètres (`param("a")`) restent des
 * constantes symboliques. Dériver par rapport à une variable traite les autres comme des
 * constantes :
 *
 *   ExprPtr f = lire("x^2*y + sin(x*y)", {"x", {"y"}});
 *   f->derivee("y")            // x^2 + x*cos(x*y)
 *   gradient(f)                // { dérivée en x, dérivée en y } (variables triées par nom)
 *
 * Les résultats de ce fichier sont simplifiés. Les opérations à une variable (intégrale,
 * limite, développement limité, résolution, évaluation en un réel) ne s'appliquent pas aux
 * expressions à plusieurs variables.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : Valeurs
 * Description : Valeur numérique de chaque variable, par nom.
 */
using Valeurs = std::map<std::string, double>;

/*
 * Nom : variables
 * Description : Variables de l'expression (noeuds Variable), triées par nom, sans doublon.
 * Utilisation : for (const ExprPtr& v : variables(f)) std::cout << v << '\n';
 */
std::vector<ExprPtr> variables(const ExprPtr& e);
std::vector<ExprPtr> variables(const std::vector<ExprPtr>& expressions);

/*
 * Nom : gradient
 * Description : Vecteur des dérivées partielles de f, dans l'ordre des variables données
 *               (par défaut, celles de f triées par nom).
 * Utilisation : std::vector<ExprPtr> g = gradient(f, {var("x"), var("y")});
 */
std::vector<ExprPtr> gradient(const ExprPtr& f, const std::vector<ExprPtr>& vars = {});

/*
 * Nom : jacobienne
 * Description : Matrice des dérivées partielles J[i][j] = d fs[i] / d vars[j].
 */
std::vector<std::vector<ExprPtr>> jacobienne(const std::vector<ExprPtr>& fs, const std::vector<ExprPtr>& vars = {});

/*
 * Nom : hessienne
 * Description : Matrice des dérivées secondes H[i][j] = d²f / d vars[i] d vars[j] (symétrique
 *               pour une fonction de classe C², seul le triangle supérieur est calculé).
 */
std::vector<std::vector<ExprPtr>> hessienne(const ExprPtr& f, const std::vector<ExprPtr>& vars = {});

/*
 * Nom : laplacien
 * Description : Somme des dérivées secondes pures, d²f/dx² + d²f/dy² + ...
 */
ExprPtr laplacien(const ExprPtr& f, const std::vector<ExprPtr>& vars = {});

/*
 * Nom : divergence
 * Description : Divergence d'un champ de vecteurs : somme des d champ[i] / d vars[i] (autant
 *               de composantes que de variables, sinon std::invalid_argument).
 */
ExprPtr divergence(const std::vector<ExprPtr>& champ, const std::vector<ExprPtr>& vars);

/*
 * Nom : deriveeMixte
 * Description : Dérivées partielles successives, dans l'ordre donné (la première variable
 *               est dérivée en premier).
 * Utilisation : ExprPtr fxy = deriveeMixte(f, {var("x"), var("y")});   // d²f/dy dx
 */
ExprPtr deriveeMixte(const ExprPtr& f, const std::vector<ExprPtr>& ordre);

/*
 * Nom : evaluer
 * Description : Valeur numérique de l'expression (ou de chaque expression, en partageant les
 *               sous-expressions communes) avec les valeurs données. Lève std::invalid_argument
 *               si une variable n'a pas de valeur, et std::logic_error si l'expression
 *               contient un paramètre ou un noeud non évalué.
 * Utilisation : double v = evaluer(f, {{"x", 1.0}, {"y", 2.0}});
 *               std::vector<double> g = evaluer(gradient(f), {{"x", 1.0}, {"y", 2.0}});
 */
double evaluer(const ExprPtr& e, const Valeurs& valeurs);
std::vector<double> evaluer(const std::vector<ExprPtr>& expressions, const Valeurs& valeurs);

} // namespace symalgo
