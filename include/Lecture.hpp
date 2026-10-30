/**
 * @file Lecture.hpp
 * @brief Lecture d'expressions et d'équations depuis du texte.
 *
 * Grammaire (précédences usuelles, puissance associative à droite) :
 *
 *   equation   := expression [ '=' expression ]
 *   expression := terme { ('+' | '-') terme }
 *   terme      := unaire { ('*' | '/') unaire | puissance }   // 2x, 3(x + 1), x sin(x)
 *   unaire     := ('-' | '+') unaire | puissance
 *   puissance  := primaire [ ('^' | '**') unaire | '²' | '³' ]
 *   primaire   := nombre | identifiant | fonction '(' expression ')' | '(' expression ')'
 *
 * - Les nombres sont exacts, décimaux compris : « 0.1 » est le rationnel 1/10, « 2.5e-3 »
 *   vaut 1/400, et les entiers sont de taille arbitraire.
 * - La multiplication implicite ne s'applique que devant un identifiant ou une parenthèse
 *   (« 2x », « 2(x + 1) », « (x + 1)(x - 1) ») ; elle a la même précédence que « * », donc
 *   « 1/2x » vaut x/2. « -x^2 » vaut -(x^2) et « 2^-1 » vaut 1/2.
 * - Fonctions : sin, cos, tan, exp, ln (et log), sqrt, asin, acos, atan (et arcsin, arccos,
 *   arctan). Constantes : pi (et π), e (= exp(1)).
 * - L'identifiant choisi comme variable (« x » par défaut) devient la variable, de même que
 *   ceux d'OptionsLecture::autresVariables (fonctions de plusieurs variables) ; tout autre
 *   identifiant devient un paramètre symbolique (« a*x + b »). Un paramètre collé à une
 *   parenthèse (« f(x) », « sinh(x) ») est refusé comme fonction inconnue ; « a (x + 1) »
 *   ou « a*(x + 1) » sont des produits.
 * - Symboles Unicode acceptés : π, ×, ·, ÷, − (signe moins) et les exposants ², ³.
 *
 * Le texte produit par l'affichage d'une expression exacte se relit à l'identique :
 * lire(e->texte()) est le même noeud que e.
 */

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : OptionsLecture
 * Description : Réglages de la lecture : nom de la variable et noms des variables
 *               supplémentaires d'une fonction de plusieurs variables (les autres
 *               identifiants sont des paramètres).
 * Utilisation : lire("x^2*y + a", {"x", {"y"}});   // x et y variables, a paramètre
 */
struct OptionsLecture {
    OptionsLecture() = default;
    OptionsLecture(std::string variable, std::vector<std::string> autresVariables = {})
        : variable(std::move(variable)), autresVariables(std::move(autresVariables)) {}

    std::string variable = "x";
    std::vector<std::string> autresVariables;
};

/*
 * Nom : ErreurLecture
 * Description : Texte mal formé. what() contient le message, le texte et un repère sous la
 *               position fautive ; position() est l'indice (à partir de 1) du caractère
 *               fautif, raison() le message seul.
 */
class ErreurLecture : public std::invalid_argument {
public:
    ErreurLecture(const std::string& texte, std::size_t position, const std::string& raison);

    std::size_t position() const { return m_position; }
    const std::string& raison() const { return m_raison; }

private:
    std::size_t m_position;
    std::string m_raison;
};

/*
 * Nom : lire
 * Description : Expression écrite dans le texte, sous forme canonique. Lève ErreurLecture
 *               si le texte est mal formé ou contient un « = ».
 * Utilisation : ExprPtr e = lire("x^2 + 3*sin(x)");
 *               ExprPtr f = lire("t^2 - 1", {"t"});
 */
ExprPtr lire(const std::string& texte, const OptionsLecture& options = {});

/*
 * Nom : EgaliteLue
 * Description : Les deux membres d'une équation lue (droite vaut 0 sans « = »).
 */
struct EgaliteLue {
    ExprPtr gauche;
    ExprPtr droite;
};

/*
 * Nom : lireEquation
 * Description : Équation « gauche = droite », ou expression seule (= 0).
 * Utilisation : EgaliteLue eq = lireEquation("sin(x) = 1/2");
 *               Solutions s = resoudre(eq.gauche, eq.droite);
 */
EgaliteLue lireEquation(const std::string& texte, const OptionsLecture& options = {});

namespace litteraux {

/*
 * Nom : _expr
 * Description : Littéral d'expression (variable x).
 * Utilisation : using namespace symalgo::litteraux;
 *               ExprPtr e = "x^2 + 1"_expr;
 */
inline ExprPtr operator""_expr(const char* texte, std::size_t taille) { return lire(std::string(texte, taille)); }

} // namespace litteraux

} // namespace symalgo
