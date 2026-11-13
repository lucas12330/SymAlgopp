/**
 * @file Primitives.hpp
 * @brief Outils internes de l'intégration symbolique (src/Primitives*.cpp), partagés avec
 *        les règles par noeud de src/Regles.cpp. En-tête interne : il ne fait pas partie de
 *        l'API publique (symalgopp.hpp).
 */

#pragma once

#include "ASTNode.hpp"
#include "Canonique.hpp"
#include "Polynome.hpp"

#include <utility>
#include <vector>

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
     * Nom : parite
     * Description : 1 si f(-x) = f(x), -1 si f(-x) = -f(x), 0 sinon (la fraction nulle est paire).
     */
    int parite() const;

    /*
     * Nom : enPuissance
     * Description : La fraction g telle que f(x) = g(x^k), pour f ne contenant que des
     *               puissances de x multiples de k (une fraction paire pour k = 2).
     */
    Fraction enPuissance(int k) const;

    // Vrai si le degré du numérateur ou du dénominateur dépasse degreMax
    bool depasse(int degreMax) const { return m_n.degre() > degreMax || m_d.degre() > degreMax; }

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

/*
 * Nom : primitiveFractionSeparee
 * Description : Comme primitiveFraction, mais la partie rationnelle de la primitive (partie
 *               polynomiale et réduction d'Hermite) est rendue à part, sous forme de fraction ;
 *               seuls les logarithmes et arctangentes sont renvoyés en expression de x. Sert
 *               aux changements de variable dont la partie rationnelle se récrit ensuite dans
 *               la variable de départ. Renvoie faux si la primitive n'est pas complète.
 */
bool primitiveFractionSeparee(const Fraction& f, const ExprPtr& x, Fraction& rationnelle, ExprPtr& transcendante);

// ============================================================================
// Extension quadratique : A + B y avec y^2 = carre (trigonométrie : s^2 = 1 - c^2 ;
// radicaux : y^2 = a x^2 + b x + c)
// ============================================================================

/*
 * Nom : ElementQuadratique
 * Description : Élément A + B y de l'extension K(y), K le corps des fractions rationnelles et
 *               y^2 = *carre (qui n'est pas un carré dans K). Tous les éléments d'un calcul
 *               partagent le même carré. Forme de calcul des fractions rationnelles en
 *               (sin, cos) et en (x, sqrt(a x^2 + b x + c)).
 */
struct ElementQuadratique {
    Fraction a, b;
    const Fraction* carre = nullptr;

    ElementQuadratique operator+(const ElementQuadratique& e) const { return {a + e.a, b + e.b, carre}; }
    ElementQuadratique operator*(const ElementQuadratique& e) const {
        return {a * e.a + b * e.b * *carre, a * e.b + b * e.a, carre};
    }
    ElementQuadratique operator*(const Nombre& k) const { return {a * k, b * k, carre}; }
    ElementQuadratique puissance(long long k) const; // std::domain_error pour 0^(-k)
    bool estNul() const { return a.estNulle() && b.estNulle(); }
    bool depasse(int degreMax) const { return a.depasse(degreMax) || b.depasse(degreMax); }

    /*
     * Nom : evaluer
     * Description : f(*this) pour une fraction rationnelle f (schéma de Horner dans l'extension).
     */
    ElementQuadratique evaluer(const Fraction& f) const;
};

// ============================================================================
// Conversion d'une expression en fraction rationnelle de « générateurs »
// ============================================================================

enum class Feuille { Convertie, Refusee, Descendre };

/*
 * Nom : convertirRationnel
 * Description : Évalue e comme fraction rationnelle à coefficients rationnels exacts de
 *               générateurs, dans un corps dont Element est le type (operator+, operator*,
 *               operator*(Nombre), puissance(long long), depasse(int)). Pour chaque noeud,
 *               regles.feuille(noeud, r) répond Convertie (r est sa valeur), Refusee (échec)
 *               ou Descendre (constantes, sommes, produits et puissances entières sont alors
 *               parcourus) ; regles.puissance(base, exposant, r) convertit une puissance
 *               d'exposant non entier, ou renvoie faux. unite est l'élément 1.
 *               Renvoie faux si e n'est pas de cette forme, ou si un degré dépasse degreMax.
 */
template <class Element, class Regles>
bool convertirRationnel(const ASTNode& e, Element& r, Regles& regles, const Element& unite, int degreMax);

namespace detail {

template <class Element, class Regles>
bool convertirPuissance(const ASTNode& base, const ASTNode& exposant, Element& r, Regles& regles, const Element& unite,
                        int degreMax) {
    long long k;
    const Constante* c = comme<Constante>(&exposant);
    if (!c || !c->getNombre().versEntier(k)) return regles.puissance(base, exposant, r);
    if (k > degreMax || k < -degreMax) return false;
    Element b = unite;
    if (!convertirRationnel(base, b, regles, unite, degreMax)) return false;
    if (k < 0 && b.estNul()) return false;
    r = b.puissance(k);
    return !r.depasse(degreMax);
}

} // namespace detail

template <class Element, class Regles>
bool convertirRationnel(const ASTNode& e, Element& r, Regles& regles, const Element& unite, int degreMax) {
    switch (regles.feuille(e, r)) {
        case Feuille::Convertie: return true;
        case Feuille::Refusee: return false;
        case Feuille::Descendre: break;
    }
    switch (e.type()) {
        case TypeNoeud::Constante: {
            const Nombre& c = static_cast<const Constante&>(e).getNombre();
            if (!c.estExact()) return false;
            r = unite * c;
            return true;
        }
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(e);
            if (!s.getConstante().estExact()) return false;
            Element total = unite * s.getConstante();
            for (const Terme& t : s.getTermes()) {
                Element terme = unite;
                if (!t.coefficient.estExact() || !convertirRationnel(*t.expression, terme, regles, unite, degreMax)) return false;
                total = total + terme * t.coefficient;
                if (total.depasse(degreMax)) return false;
            }
            r = std::move(total);
            return true;
        }
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(e);
            if (!p.getCoefficient().estExact()) return false;
            Element total = unite * p.getCoefficient();
            for (const Facteur& facteur : p.getFacteurs()) {
                Element terme = unite;
                if (!detail::convertirPuissance(*facteur.base, *facteur.exposant, terme, regles, unite, degreMax)) return false;
                total = total * terme;
                if (total.depasse(degreMax)) return false;
            }
            r = std::move(total);
            return true;
        }
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return detail::convertirPuissance(*p.getBase(), *p.getExposant(), r, regles, unite, degreMax);
        }
        default: return false;
    }
}

// ============================================================================
// Fonctions d'arguments commensurables : sin(2x), cos(3x)... multiples d'un même angle
// ============================================================================

/*
 * Nom : Multiples
 * Description : Arguments affines des fonctions d'une expression, écrits comme multiples
 *               entiers d'un même argument de base u (2x et 3x : u = x ; x/2 et x : u = x/2).
 */
struct Multiples {
    std::vector<const ASTNode*> arguments;
    std::vector<long long> multiples; // arguments[i] = multiples[i] * u
    ExprPtr u, pente;                 // u affine, pente = u'

    bool ajouter(const ExprPtr& argument);   // faux si l'argument n'est pas affine
    bool resoudre(long long multipleMax);    // calcule u ; faux si non commensurables
    long long multipleDe(const ASTNode* argument) const;

private:
    ExprPtr m_premier;
};

/*
 * Nom : collecter
 * Description : Vérifie que e est une fraction rationnelle à coefficients exacts des
 *               fonctions acceptées (accepte(type)) d'arguments affines, sans x isolé, et
 *               collecte ces arguments dans m.
 */
template <class Accepte>
bool collecter(const ASTNode& e, Multiples& m, const Accepte& accepte) {
    if (!e.contientVariable()) {
        const Constante* c = comme<Constante>(&e);
        return c && c->getNombre().estExact();
    }
    switch (e.type()) {
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) {
                if (!collecter(*t.expression, m, accepte)) return false;
            }
            return true;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                if (f.exposant->contientVariable() || !collecter(*f.base, m, accepte)) return false;
            }
            return true;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return !p.getExposant()->contientVariable() && collecter(*p.getBase(), m, accepte);
        }
        default:
            return accepte(e.type()) && m.ajouter(static_cast<const FonctionUnaire&>(e).m_argument);
    }
}

// Élément du corps Q(t) pour convertirRationnel
struct ElementFraction {
    Fraction f;
    ElementFraction operator+(const ElementFraction& e) const { return {f + e.f}; }
    ElementFraction operator*(const ElementFraction& e) const { return {f * e.f}; }
    ElementFraction operator*(const Nombre& k) const { return {f * k}; }
    ElementFraction puissance(long long k) const { return {f.puissance(k)}; }
    bool estNul() const { return f.estNulle(); }
    bool depasse(int degreMax) const { return f.depasse(degreMax); }
};

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
 * Nom : tableEtendue
 * Description : Primitives usuelles d'une puissance : 1/(a x^2 + b x + c) -> atan,
 *               1/sqrt(a x^2 + b x + c) -> asin ou ln, 1/cos^2 -> tan, tan^2, c^(a x + b)...
 */
ExprPtr tableEtendue(const ExprPtr& f);

/*
 * Nom : parSubstitution
 * Description : Changement de variable « la dérivée divise » : u = g(x) pris parmi les
 *               sous-expressions de f ; si f/g' ne s'écrit qu'en fonction de u, on intègre
 *               en u puis on revient à x.
 */
ExprPtr parSubstitution(const ExprPtr& f);

/*
 * Nom : parExponentiellePolynome
 * Description : Produit de polynômes en x, d'exponentielles, de sinus, cosinus, sinh et cosh
 *               d'arguments affines (puissances entières positives) : linéarisation en
 *               somme de x^n exp(W) cos(V)..., chaque terme ayant une primitive explicite.
 *               x exp(x), exp(x) sin(x), sin(x)^2, x^2 cos(3x), sin(2x) cos(5x), x cosh(x).
 */
ExprPtr parExponentiellePolynome(const ExprPtr& f);

/*
 * Nom : parExponentielleRationnelle
 * Description : Fraction rationnelle de exp(k u), sinh(k u), cosh(k u), tanh(k u) (u affine, k
 *               entiers) sans x isolé : t = exp(u) ramène à une fraction rationnelle de t.
 */
ExprPtr parExponentielleRationnelle(const ExprPtr& f);

/*
 * Nom : parEquationDeRisch
 * Description : R(x) exp(W(x)), R fraction rationnelle et W polynôme : la primitive est
 *               S(x) exp(W(x)) où S est solution rationnelle de S' + W' S = R, si elle existe
 *               (sinon la primitive n'est pas de cette forme : pointeur nul). Accepte aussi une
 *               somme de tels termes : exp(x)/x - exp(x)/x^2 = (exp(x)/x)'.
 */
ExprPtr parEquationDeRisch(const ExprPtr& f);

/*
 * Nom : parTrigonometrieImpaire, parTrigonometrie
 * Description : Fraction rationnelle de sin(k u), cos(k u), tan(k u) (u affine, k entiers),
 *               ramenée à la forme A(c) + s B(c) puis à une fraction rationnelle par les
 *               règles de Bioche : t = cos u, t = sin u, t = tan u, sinon t = tan(u/2).
 *               parTrigonometrieImpaire ne traite que les polynômes en sin u et cos u impairs
 *               en l'un des deux (sin^3 cos^2 : primitive polynomiale en cos u) ; les autres
 *               polynômes sont linéarisés par parExponentiellePolynome.
 */
ExprPtr parTrigonometrieImpaire(const ExprPtr& f);
ExprPtr parTrigonometrie(const ExprPtr& f);

/*
 * Nom : parParties
 * Description : f = g L avec L = ln(v)^n, atan(v), asin(v) ou acos(v) : G L - int G L', où G
 *               est une primitive de g (x ln(x), ln(x)^2, x atan(x), ln(x^2 + 1), asin(x)/x^2).
 */
ExprPtr parParties(const ExprPtr& f);

/*
 * Nom : parRadicaux
 * Description : Fractions rationnelles de x et d'un radical :
 *                 (a x + b)^(p/q) : t = (a x + b)^(1/q) ;
 *                 sqrt(a x^2 + b x + c) : forme A(x) + B(x) y, puis coefficients indéterminés
 *                 (int P/y = S y + lambda int dx/y) ou substitution d'Euler.
 */
ExprPtr parRadicaux(const ExprPtr& f);

/*
 * Nom : secondeChanceSomme
 * Description : Pour une somme dont aucun terme ne s'intègre seul : méthodes qui traitent la
 *               somme entière (équation de Risch). Pointeur nul en cas d'échec.
 */
ExprPtr secondeChanceSomme(const ExprPtr& somme);

/*
 * Nom : parDeveloppement
 * Description : Dernier recours pour un produit : s'il se développe en une somme, celle-ci
 *               est intégrée terme à terme ((x + 1)^2 (x - 3)). Pointeur nul si un terme
 *               reste non évalué.
 */
ExprPtr parDeveloppement(const ExprPtr& f);

} // namespace symalgo
