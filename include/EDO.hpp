/**
 * @file EDO.hpp
 * @brief Résolution symbolique d'équations différentielles ordinaires non linéaires.
 *
 * Premier ordre F(x, y, y') = 0 : quadrature, linéaire, Bernoulli, variables séparables,
 * y' = F(a x + b y), homogène, exacte et facteurs intégrants mu(x) et mu(y), Riccati avec
 * solution particulière, Clairaut. Second ordre réductible : y absente (y'' = f(x, y')) ou
 * x absente (y'' = f(y, y')).
 *
 * Les équations linéaires à coefficients constants d'ordre quelconque restent le domaine de
 * EquationDifferentielle.
 *
 * Toute solution renvoyée a été réinjectée dans l'équation : une famille qui ne la vérifie
 * pas est écartée, et une équation hors de portée donne un résultat vide, jamais une
 * solution fausse.
 */

#pragma once

#include "ASTNode.hpp"
#include "Journal.hpp"

#include <string>
#include <vector>

namespace symalgo {

/*
 * Nom : InconnuesEDO
 * Description : Les variables d'une équation différentielle : x, la fonction inconnue y et
 *               ses dérivées, représentées par des variables distinctes.
 * Utilisation : InconnuesEDO v;   // x, y, dy, d2y
 *               resoudreEDO(v.dy - v.x * v.y, v);   // y' = x y
 */
struct InconnuesEDO {
    ExprPtr x = var("x");
    ExprPtr y = var("y");
    ExprPtr dy = var("dy");   // y'
    ExprPtr d2y = var("d2y"); // y''
};

/*
 * Nom : SolutionEDO
 * Description : Une famille de solutions, les constantes d'intégration étant les paramètres
 *               C1, C2... Elle est explicite (y = y(x, C1...)) quand y a pu être isolée,
 *               sinon implicite : relation(x, y, C1...) = 0. Une solution singulière n'est
 *               obtenue pour aucune valeur des constantes de la famille générale (équilibre
 *               perdu par une division, enveloppe d'une famille de Clairaut).
 */
struct SolutionEDO {
    ExprPtr y;        // forme explicite, ou pointeur nul
    ExprPtr relation; // forme implicite (nulle sur les solutions), ou pointeur nul
    bool singuliere = false;

    bool estExplicite() const { return static_cast<bool>(y); }
    std::string texte() const; // « y = ... » ou « ... = 0 »
};

/*
 * Nom : SolutionsEDO
 * Description : Résultat de resoudreEDO. liste est vide si l'équation n'a pas été résolue.
 *               methode nomme la classe reconnue : quadrature, lineaire, bernoulli,
 *               separable, affine, homogene, exacte, facteur_integrant, riccati, clairaut,
 *               ordre2_sans_y, ordre2_sans_x.
 */
struct SolutionsEDO {
    std::vector<SolutionEDO> liste;
    std::string methode;

    bool resolue() const { return !liste.empty(); }
};

/*
 * Nom : resoudreEDO
 * Description : Résout equation = 0, où equation est une expression de x, y, dy (et d2y).
 *               Avec un Journal, les étapes y sont notées (identifiants « edo.* »).
 * Utilisation : InconnuesEDO v;
 *               SolutionsEDO s = resoudreEDO(v.dy + v.y / v.x - ast_pow(v.x, 2.0), v);
 *               s.liste[0].y;   // x^3/4 + C1/x
 */
SolutionsEDO resoudreEDO(const ExprPtr& equation, const InconnuesEDO& inconnues = InconnuesEDO());
SolutionsEDO resoudreEDO(const ExprPtr& equation, const InconnuesEDO& inconnues, Journal& journal);

/*
 * Nom : resoudreEDO (texte)
 * Description : Lit l'équation : y' et y'' désignent les dérivées de la fonction inconnue.
 *               Lève ErreurLecture si le texte est mal formé.
 * Utilisation : resoudreEDO("y' = x*y^2")             // y = -2/(x^2 + C1), et y = 0
 *               resoudreEDO("y'' = y'/x")             // y = C1*x^2/2 + C2
 *               resoudreEDO("u' = t*u", "t", "u")
 */
SolutionsEDO resoudreEDO(const std::string& texte, const std::string& variable = "x", const std::string& fonction = "y");
SolutionsEDO resoudreEDO(const std::string& texte, Journal& journal, const std::string& variable = "x",
                         const std::string& fonction = "y");

/*
 * Nom : verifierEDO
 * Description : Réinjecte une solution dans equation = 0 : vrai si le résidu est nul, par
 *               simplification ou, à défaut, numériquement en plusieurs points tirés dans le
 *               domaine (valeurs arbitraires de x, de y pour une forme implicite, et des
 *               constantes). Faux si le résidu est non nul ou si trop peu de points sont
 *               évaluables.
 */
bool verifierEDO(const ExprPtr& equation, const InconnuesEDO& inconnues, const SolutionEDO& solution);

} // namespace symalgo
