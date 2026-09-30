/**
 * @file Solveur.hpp
 * @brief Résolution d'équations réelles à une inconnue.
 *
 * La résolution exacte isole l'inconnue en inversant les opérations une à une (sommes,
 * produits, puissances, exp, ln, fonctions trigonométriques et réciproques) et traite les
 * polynômes par leurs racines certifiées (suites de Sturm). Les équations trigonométriques
 * donnent des familles de solutions paramétrées par un entier (pi/6 + 2*k*pi).
 * La résolution sur un intervalle déplie ces familles et complète, si besoin, par une
 * recherche numérique.
 */

#pragma once

#include <string>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : Solution
 * Description : Une solution réelle : valeur exacte quand elle est connue, approximation,
 *               multiplicité (racines de polynômes), et entiers libres pour une famille
 *               (par exemple k dans pi/6 + 2*k*pi, k parcourant les entiers relatifs).
 */
struct Solution {
    ExprPtr valeur;                     // exacte ; approchée (constante réelle) sinon
    double approximation = 0.0;         // NaN pour une famille (dépend de k)
    int multiplicite = 1;
    bool exacte = true;
    std::vector<std::string> entiers;   // paramètres entiers de la famille (vide sinon)

    bool estFamille() const { return !entiers.empty(); }
};

/*
 * Nom : Solutions
 * Description : Ensemble des solutions trouvées. complet est vrai si toutes les solutions
 *               réelles sont listées (sinon, certaines formes n'ont pas pu être résolues) ;
 *               toutReel si l'équation est vraie pour tout x.
 */
struct Solutions {
    std::vector<Solution> liste;
    bool complet = true;
    bool toutReel = false;
};

/*
 * Nom : resoudre
 * Description : Résout gauche = droite (ou expression = 0) sur les réels.
 * Utilisation : Solutions s = resoudre(ast_pow(x, 2.0) - 3.0 * x + 2.0);  // x = 1, x = 2
 *               Solutions t = resoudre(ast_sin(x), frac(1, 2));          // familles avec k
 */
Solutions resoudre(const ExprPtr& gauche, const ExprPtr& droite);
Solutions resoudre(const ExprPtr& expression);

/*
 * Nom : resoudreSurIntervalle
 * Description : Solutions de expression = 0 dans [a, b] : les familles sont dépliées pour
 *               chaque entier k donnant une solution de l'intervalle, et une recherche
 *               numérique complète les formes non résolues exactement.
 */
Solutions resoudreSurIntervalle(const ExprPtr& expression, double a, double b);

/*
 * Nom : resoudreNumerique
 * Description : Racines de f dans [a, b] par recherche numérique seule : échantillonnage,
 *               méthode de Brent sur les changements de signe, et racines doubles repérées
 *               aux points critiques. Le résultat n'est jamais garanti complet.
 */
Solutions resoudreNumerique(const ExprPtr& f, double a, double b);

} // namespace symalgo
