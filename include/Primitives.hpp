/**
 * @file Primitives.hpp
 * @brief Outils internes de l'intégration symbolique (src/Primitives*.cpp), partagés avec
 *        les règles par noeud de src/Regles.cpp. En-tête interne : il ne fait pas partie de
 *        l'API publique (symalgopp.hpp).
 */

#pragma once

#include "ASTNode.hpp"
#include "Polynome.hpp"

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

/*
 * Nom : contientNonEvaluee
 * Description : Vrai si e contient une intégrale non évaluée (primitive seulement partielle).
 */
bool contientNonEvaluee(const ASTNode& e);

// ============================================================================
// Fractions rationnelles (src/PrimitivesRationnelles.cpp)
// ============================================================================

/*
 * Nom : Fraction
 * Description : Fraction rationnelle exacte n/d d'une variable, toujours réduite : n et d
 *               premiers entre eux, d unitaire (0 s'écrit 0/1). C'est la forme de calcul
 *               des méthodes d'intégration : après un changement de variable (règles de
 *               Bioche, t = exp(x), racine d'une fonction affine...), l'intégrande y est
 *               construite directement, sans passer par un arbre d'expression.
 *               Les opérations lèvent std::domain_error pour une division par 0.
 */
class Fraction {
public:
    Fraction();                                           // 0
    explicit Fraction(const Nombre& constante);
    explicit Fraction(Polynome numerateur);
    Fraction(Polynome numerateur, Polynome denominateur);

    static Fraction variable();                           // x

    const Polynome& numerateur() const { return m_n; }
    const Polynome& denominateur() const { return m_d; }
    bool estNulle() const { return m_n.estNul(); }
    bool estPolynome() const { return m_d.degre() == 0; }

    Fraction operator+(const Fraction& b) const;
    Fraction operator-(const Fraction& b) const;
    Fraction operator*(const Fraction& b) const;
    Fraction operator/(const Fraction& b) const;
    Fraction operator*(const Nombre& k) const;
    Fraction operator-() const;
    bool operator==(const Fraction& b) const;
    bool operator!=(const Fraction& b) const { return !(*this == b); }

    Fraction puissance(long long k) const;                // k de signe quelconque
    Fraction derivee() const;
    Fraction composer(const Fraction& g) const;           // f(g(x))

    /*
     * Nom : versExpression
     * Description : Expression de la fraction : numérateur primitif sur les facteurs sans
     *               carré primitifs du dénominateur, la constante devant (x/(2*(x^2 + 1)),
     *               1/(x + 1)^2).
     */
    ExprPtr versExpression(const ExprPtr& x) const;

private:
    struct DejaReduite {};
    Fraction(Polynome n, Polynome d, DejaReduite) : m_n(std::move(n)), m_d(std::move(d)) {}
    void reduire();

    Polynome m_n, m_d;
};

/*
 * Nom : versFraction
 * Description : Écrit e comme fraction rationnelle exacte de la variable x (sommes,
 *               produits, puissances entières, constantes exactes). Renvoie faux si e n'en
 *               est pas une (fonction, paramètre, coefficient réel, autre variable, exposant
 *               non entier) ou si elle est trop grande (degré > 256).
 */
bool versFraction(const ExprPtr& e, const ExprPtr& x, Fraction& f);

/*
 * Nom : primitiveFraction
 * Description : Primitive de f en la variable x : partie polynomiale, partie rationnelle par
 *               la réduction d'Hermite, puis logarithmes et arctangentes des facteurs
 *               irréductibles de degré 1 et 2 du dénominateur sans carré (les logarithmes de
 *               même coefficient sont regroupés : ln|x^2 - 4|/2, ln|(x - 1)/(x + 1)|/2).
 *               Un facteur irréductible de degré >= 3 ne s'intègre que si son numérateur est
 *               proportionnel à sa dérivée ; sinon sa contribution reste une intégrale non
 *               évaluée et complete reçoit faux.
 */
ExprPtr primitiveFraction(const Fraction& f, const ExprPtr& x, bool* complete = nullptr);

// ============================================================================
// Méthodes d'intégration : chacune renvoie un pointeur nul si elle ne s'applique pas
// ============================================================================

/*
 * Nom : parFractionRationnelle
 * Description : Fraction rationnelle à coefficients rationnels qui n'est pas un polynôme :
 *               primitiveFraction. Pointeur nul si rien n'a pu être intégré.
 */
ExprPtr parFractionRationnelle(const ExprPtr& f);

/*
 * Nom : parDeveloppement
 * Description : Dernier recours pour un produit : s'il se développe en une somme, celle-ci
 *               est intégrée terme à terme ((x + 1)^2 (x - 3)). Pointeur nul si un terme
 *               reste non évalué.
 */
ExprPtr parDeveloppement(const ExprPtr& f);

} // namespace symalgo
