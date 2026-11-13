/**
 * @file Latex.hpp
 * @brief Rendu LaTeX des expressions, des solutions d'une équation et du journal d'un calcul.
 *
 *   std::string l = latex(e);                       // "\frac{x + 1}{x - 1}"
 *   std::string s = latex(resoudre(ast_sin(x), frac(1, 2)));
 *   Journal j; resoudre(g, d, j);
 *   std::string doc = documentLatex(latex(j));      // document compilable tel quel
 *
 * Le rendu suit la précédence de l'affichage texte (somme < produit < puissance < atome) :
 * sommes écrites par degré décroissant, quotients regroupés en une seule fraction, racines
 * pour les exposants 1/n, parenthèses \left( \right) seulement quand elles sont nécessaires.
 * Le temps est linéaire en la taille du texte produit, écrit dans une seule chaîne.
 *
 * Choix de rendu :
 *   - un coefficient rationnel et les facteurs forment une seule fraction : -3/4*x -> -\frac{3 x}{4} ;
 *   - x^(2/3) -> x^{\frac{2}{3}} (la racine n'est utilisée que pour les exposants 1/n) ;
 *   - sin(x)^2 -> \sin^{2}\left(x\right) (idem cos, tan, sinh, cosh, tanh, ln) ;
 *   - exp(u) -> e^{u} ; abs(u) -> \left|u\right| ; sign(u) -> \operatorname{sgn}\left(u\right).
 */

#pragma once

#include <string>

#include "ASTNode.hpp"
#include "Journal.hpp"
#include "Solveur.hpp"

namespace symalgo {

/*
 * Nom : OptionsLatex
 * Description : fractions : \frac{a}{b} (sinon a/b) ; racines : u^(1/2) -> \sqrt{u},
 *               u^(1/3) -> \sqrt[3]{u} (sinon u^{\frac{1}{2}}) ; multiplicationExplicite :
 *               2 \cdot x au lieu de 2 x (le \cdot reste obligatoire entre deux nombres).
 */
struct OptionsLatex {
    bool fractions = true;
    bool racines = true;
    bool multiplicationExplicite = false;
};

/*
 * Nom : latex
 * Description : Texte LaTeX (mode mathématique, sans $) de l'expression. Les nombres réels
 *               sont écrits dans la notation la plus courte qui se relit au même double
 *               (1.5, 10^{-10}, 2.5 \times 10^{20}) ; les noms de lettres grecques deviennent
 *               \alpha..., les noms longs \mathrm{nom}, k2 devient k_{2}.
 *               Lève std::invalid_argument si l'expression est nulle.
 */
std::string latex(const ExprPtr& e, const OptionsLatex& options = {});

/*
 * Nom : latex (solutions)
 * Description : Ensemble de solutions de l'équation en `inconnue` : « x = 2 \text{ ou } x = 3 »,
 *               familles « x = \frac{\pi}{6} + 2 k \pi,\ k \in \mathbb{Z} », multiplicité
 *               « (\text{double}) », valeur approchée « x \approx 1.3247 », « x \in \mathbb{R} »
 *               (identité), « x \in \emptyset » (aucune solution), résultat incomplet signalé.
 */
std::string latex(const Solutions& solutions, const std::string& inconnue = "x", const OptionsLatex& options = {});

/*
 * Nom : latex (journal)
 * Description : Les étapes d'un calcul en liste itemize imbriquée (une ligne par étape, retrait
 *               selon le niveau, au plus 4 niveaux de liste) : description, puis expressions
 *               en mode mathématique ($avant \longrightarrow apres$ ou $gauche = droite$).
 *               À inclure dans un document (documentLatex) ; nécessite amsmath.
 */
std::string latex(const Journal& journal, const OptionsLatex& options = {});

/*
 * Nom : documentLatex
 * Description : Document autonome (article, utf8, T1, amsmath, amssymb) dont le corps est donné
 *               tel quel (par exemple latex(journal), ou des \[ ... \] autour de latex(e)).
 */
std::string documentLatex(const std::string& corps);

} // namespace symalgo
