/**
 * @file PrimitivesRationnelles.cpp
 * @brief Fractions rationnelles exactes et leur intégration.
 *
 * int n/d se calcule en trois temps, tous exacts sur Q :
 *   1. division euclidienne : la partie polynomiale s'intègre terme à terme ;
 *   2. réduction d'Hermite (version linéaire de Mack) : int a/d = g + int r/d*, où g est une
 *      fraction rationnelle et d* le produit des facteurs sans carré de d. Elle ne demande
 *      que des PGCD : aucune factorisation ;
 *   3. partie logarithmique de r/d* : si r = c d*', c'est c ln|d*| ; sinon d* est factorisé
 *      sur Q (Polynome::facteursIrreductibles) et r/d* décomposée en éléments simples.
 *      Un facteur de degré 1 donne un logarithme, un facteur de degré 2 un logarithme et une
 *      arctangente (discriminant négatif) ou un second logarithme (discriminant positif).
 *      Un facteur irréductible de degré 4 est scindé en deux trinômes réels quand leurs
 *      coefficients sont dans un corps quadratique Q(sqrt(m)) : x^4 + 1 =
 *      (x^2 + sqrt(2) x + 1)(x^2 - sqrt(2) x + 1). Au-delà (résidus algébriques de degré
 *      plus grand), l'élément simple reste une intégrale non évaluée.
 *
 * Les logarithmes de même coefficient sont regroupés, ce qui redonne la forme minimale de
 * Rothstein-Trager quand les résidus sont rationnels : x/(x^2 - 4) s'intègre en
 * ln|x^2 - 4|/2 et non en ln|x - 2|/2 + ln|x + 2|/2.
 */

#include "Canonique.hpp"
#include "Journal.hpp"
#include "Primitives.hpp"

#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace symalgo {

namespace {

constexpr int DEGRE_FRACTION_MAX = 256;

const Polynome& polynomeUn() {
    static const Polynome* un = new Polynome({Nombre(1)});
    return *un;
}

bool egaux(const Polynome& a, const Polynome& b) {
    if (a.degre() != b.degre()) return false;
    for (int i = 0; i <= a.degre(); ++i) {
        if (a.coefficient(i) != b.coefficient(i)) return false;
    }
    return true;
}

Polynome quotient(const Polynome& a, const Polynome& b) {
    Polynome q, r;
    a.diviser(b, q, r);
    return q;
}

Polynome reste(const Polynome& a, const Polynome& b) {
    Polynome q, r;
    a.diviser(b, q, r);
    return r;
}

Polynome puissance(Polynome p, long long k) {
    Polynome r = polynomeUn();
    while (k > 0) {
        if (k & 1) r = r * p;
        k >>= 1;
        if (k > 0) p = p * p;
    }
    return r;
}

// a = c b pour une constante c (b non nul) ?
bool proportionnels(const Polynome& a, const Polynome& b, Nombre& c) {
    if (a.degre() != b.degre() || b.estNul()) return false;
    c = a.dominant() / b.dominant();
    return egaux(a, b * c);
}

// Primitive d'un polynôme (constante d'intégration nulle)
Polynome primitive(const Polynome& p) {
    std::vector<Nombre> c(static_cast<std::size_t>(p.degre() + 2), Nombre(0));
    for (int i = 0; i <= p.degre(); ++i) c[static_cast<std::size_t>(i + 1)] = p.coefficient(i) / Nombre(i + 1);
    return Polynome(std::move(c));
}

} // namespace

// ============================================================================
// Fraction
// ============================================================================

Fraction::Fraction() : m_d(polynomeUn()) {}

Fraction::Fraction(const Nombre& constante) : m_n(Polynome({constante})), m_d(polynomeUn()) {}

Fraction::Fraction(Polynome numerateur) : m_n(std::move(numerateur)), m_d(polynomeUn()) {}

Fraction::Fraction(Polynome numerateur, Polynome denominateur) : m_n(std::move(numerateur)), m_d(std::move(denominateur)) {
    reduire();
}

Fraction Fraction::variable() { return Fraction(Polynome({Nombre(0), Nombre(1)})); }

void Fraction::reduire() {
    if (m_d.estNul()) throw std::domain_error("Fraction : denominateur nul");
    if (m_n.estNul()) {
        m_d = polynomeUn();
        return;
    }
    if (m_d.degre() > 0) {
        const Polynome g = Polynome::pgcd(m_n, m_d);
        if (g.degre() > 0) {
            m_n = quotient(m_n, g);
            m_d = quotient(m_d, g);
        }
    }
    if (!m_d.dominant().estUn()) {
        const Nombre inverse = m_d.dominant().inverse();
        m_n = m_n * inverse;
        m_d = m_d * inverse;
    }
}

Fraction Fraction::operator+(const Fraction& b) const {
    if (estNulle()) return b;
    if (b.estNulle()) return *this;
    if (egaux(m_d, b.m_d)) {
        if (m_d.degre() == 0) return Fraction(m_n + b.m_n, m_d, DejaReduite{});
        return Fraction(m_n + b.m_n, m_d);
    }
    return Fraction(m_n * b.m_d + b.m_n * m_d, m_d * b.m_d);
}

Fraction Fraction::operator-() const { return Fraction(m_n * Nombre(-1), m_d, DejaReduite{}); }

Fraction Fraction::operator-(const Fraction& b) const { return *this + (-b); }

Fraction Fraction::operator*(const Fraction& b) const {
    if (estNulle() || b.estNulle()) return Fraction();
    if (m_d.degre() == 0 && b.m_d.degre() == 0) return Fraction(m_n * b.m_n, m_d, DejaReduite{});
    return Fraction(m_n * b.m_n, m_d * b.m_d);
}

Fraction Fraction::operator*(const Nombre& k) const {
    if (k.estZero()) return Fraction();
    return Fraction(m_n * k, m_d, DejaReduite{});
}

Fraction Fraction::operator/(const Fraction& b) const {
    if (b.estNulle()) throw std::domain_error("Fraction : division par zero");
    if (estNulle()) return Fraction();
    return Fraction(m_n * b.m_d, m_d * b.m_n);
}

bool Fraction::operator==(const Fraction& b) const { return egaux(m_n, b.m_n) && egaux(m_d, b.m_d); }

Fraction Fraction::puissance(long long k) const {
    if (k == 0) return Fraction(Nombre(1));
    if (k > 0) return Fraction(symalgo::puissance(m_n, k), symalgo::puissance(m_d, k), DejaReduite{});
    if (estNulle()) throw std::domain_error("Fraction : puissance negative de zero");
    // Inverse : seul le coefficient dominant du nouveau dénominateur est à ramener à 1
    const Nombre inverse = m_n.dominant().inverse();
    return Fraction(symalgo::puissance(m_d * inverse, -k), symalgo::puissance(m_n * inverse, -k), DejaReduite{});
}

Fraction Fraction::derivee() const {
    if (m_d.degre() == 0) return Fraction(m_n.derivee(), m_d, DejaReduite{});
    return Fraction(m_n.derivee() * m_d - m_n * m_d.derivee(), m_d * m_d);
}

Fraction Fraction::composer(const Fraction& g) const {
    if (estNulle()) return Fraction();
    const Polynome& p = g.m_n;
    const Polynome& q = g.m_d;
    // P(p/q) q^deg(P) = sum c_i p^i q^(deg - i), par le schéma de Horner
    auto homogene = [&](const Polynome& P) {
        Polynome r({P.dominant()}), qPuissance = polynomeUn();
        for (int i = P.degre() - 1; i >= 0; --i) {
            qPuissance = qPuissance * q;
            r = r * p + qPuissance * P.coefficient(i);
        }
        return r;
    };
    Polynome n = homogene(m_n), d = homogene(m_d);
    const int ecart = m_d.degre() - m_n.degre(); // f(g) = n q^deg(d) / (d q^deg(n))
    if (ecart > 0) n = n * symalgo::puissance(q, ecart);
    if (ecart < 0) d = d * symalgo::puissance(q, -ecart);
    return Fraction(std::move(n), std::move(d));
}

int Fraction::parite() const {
    // Fraction réduite : paire si n et d sont pairs, impaire si l'un est pair et l'autre impair
    auto pariteDe = [](const Polynome& p) {
        bool pair = true, impair = true;
        for (int i = 0; i <= p.degre(); ++i) {
            if (p.coefficient(i).estZero()) continue;
            (i % 2 == 0 ? impair : pair) = false;
        }
        return pair ? 1 : impair ? -1 : 0;
    };
    if (estNulle()) return 1;
    return pariteDe(m_n) * pariteDe(m_d);
}

Fraction Fraction::enPuissance(int k) const {
    auto reduit = [k](const Polynome& p) {
        std::vector<Nombre> c;
        for (int i = 0; i <= p.degre(); i += k) c.push_back(p.coefficient(i));
        return Polynome(std::move(c));
    };
    return Fraction(reduit(m_n), reduit(m_d), DejaReduite{});
}

// ============================================================================
// ElementQuadratique
// ============================================================================

ElementQuadratique ElementQuadratique::puissance(long long k) const {
    ElementQuadratique base = *this;
    if (k < 0) {
        // 1/(A + B y) = (A - B y)/(A^2 - B^2 y^2)
        const Fraction norme = a * a - b * b * *carre;
        if (norme.estNulle()) throw std::domain_error("ElementQuadratique : inverse de zero");
        base = {a / norme, -(b / norme), carre};
        k = -k;
    }
    ElementQuadratique r{Fraction(Nombre(1)), Fraction(), carre};
    while (k > 0) {
        if (k & 1) r = r * base;
        k >>= 1;
        if (k > 0) base = base * base;
    }
    return r;
}

ElementQuadratique ElementQuadratique::evaluer(const Fraction& f) const {
    auto horner = [this](const Polynome& p) {
        ElementQuadratique r{Fraction(), Fraction(), carre};
        for (int i = p.degre(); i >= 0; --i) {
            r = r * *this;
            r.a = r.a + Fraction(p.coefficient(i));
        }
        return r;
    };
    return horner(f.numerateur()) * horner(f.denominateur()).puissance(-1);
}

ExprPtr Fraction::versExpression(const ExprPtr& x) const {
    if (m_d.degre() == 0) return m_n.versExpression(x);
    Nombre contenu;
    const Polynome n = m_n.partiePrimitive(contenu);
    AccumulateurProduit p;
    p.multiplier(n.versExpression(x));
    for (const auto& [facteur, multiplicite] : m_d.sansCarre()) {
        Nombre c;
        const Polynome primitif = facteur.partiePrimitive(c); // facteur = c * primitif
        p.multiplier(ast_pow(primitif.versExpression(x), nombre(Nombre(-multiplicite))));
        contenu = contenu * c.puissanceEntiere(-multiplicite);
    }
    p.multiplier(nombre(contenu));
    return p.construire();
}

// ============================================================================
// Expression -> fraction
// ============================================================================

namespace {

bool tropGrande(const Fraction& f) {
    return f.numerateur().degre() > DEGRE_FRACTION_MAX || f.denominateur().degre() > DEGRE_FRACTION_MAX;
}

bool convertir(const ASTNode& e, const ASTNode* x, Fraction& f);

// base ^ exposant, exposant entier constant
bool convertirPuissance(const ASTNode& base, const ASTNode& exposant, const ASTNode* x, Fraction& f) {
    long long k;
    const Constante* c = comme<Constante>(&exposant);
    if (!c || !c->getNombre().versEntier(k) || k > DEGRE_FRACTION_MAX || k < -DEGRE_FRACTION_MAX) return false;
    Fraction b;
    if (!convertir(base, x, b)) return false;
    if (k < 0 && b.estNulle()) return false;
    const long long degre = std::max(b.numerateur().degre(), b.denominateur().degre());
    if (degre * (k < 0 ? -k : k) > DEGRE_FRACTION_MAX) return false;
    f = b.puissance(k);
    return true;
}

bool convertir(const ASTNode& e, const ASTNode* x, Fraction& f) {
    switch (e.type()) {
        case TypeNoeud::Constante: {
            const Nombre& c = static_cast<const Constante&>(e).getNombre();
            if (!c.estExact()) return false;
            f = Fraction(c);
            return true;
        }
        case TypeNoeud::Variable:
            if (&e != x) return false;
            f = Fraction::variable();
            return true;
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(e);
            if (!s.getConstante().estExact()) return false;
            Fraction total(s.getConstante());
            for (const Terme& t : s.getTermes()) {
                Fraction terme;
                if (!t.coefficient.estExact() || !convertir(*t.expression, x, terme)) return false;
                total = total + terme * t.coefficient;
                if (tropGrande(total)) return false;
            }
            f = std::move(total);
            return true;
        }
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(e);
            if (!p.getCoefficient().estExact()) return false;
            Fraction total(p.getCoefficient());
            for (const Facteur& facteur : p.getFacteurs()) {
                Fraction terme;
                if (!convertirPuissance(*facteur.base, *facteur.exposant, x, terme)) return false;
                total = total * terme;
                if (tropGrande(total)) return false;
            }
            f = std::move(total);
            return true;
        }
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return convertirPuissance(*p.getBase(), *p.getExposant(), x, f);
        }
        default: return false;
    }
}

} // namespace

bool versFraction(const ExprPtr& e, const ExprPtr& x, Fraction& f) {
    try {
        return convertir(*e, x.get(), f);
    } catch (const std::domain_error&) { // dénominateur identiquement nul
        return false;
    }
}

// ============================================================================
// Intégration
// ============================================================================

namespace {

/*
 * Nom : hermite
 * Description : Réduction d'Hermite, version linéaire de Mack : int a/d = g + int a2/dEtoile,
 *               dEtoile sans carré (d unitaire, a quelconque). À chaque tour, le
 *               dénominateur dMoins = pgcd(d, d') perd une puissance de chacun de ses
 *               facteurs : b et c sont tels que b (-dEtoile dMoins'/dMoins) + c dMoinsEtoile = a,
 *               et la partie b/dMoins passe dans g.
 */
void hermite(Polynome a, const Polynome& d, Fraction& g, Polynome& a2, Polynome& dEtoile) {
    Polynome dMoins = Polynome::pgcd(d, d.derivee());
    dEtoile = quotient(d, dMoins);
    // g est accumulée sur le dénominateur commun pgcd(d, d') : chaque dMoins divise le
    // précédent, et cofacteur = denominateur / dMoins. Une seule réduction, à la fin.
    const Polynome denominateur = dMoins;
    Polynome numerateur, cofacteur = polynomeUn();
    while (dMoins.degre() > 0) {
        const Polynome dMoinsDerivee = dMoins.derivee();
        const Polynome dMoins2 = Polynome::pgcd(dMoins, dMoinsDerivee);
        const Polynome dMoinsEtoile = quotient(dMoins, dMoins2);
        const Polynome p = quotient(dEtoile * dMoinsDerivee, dMoins) * Nombre(-1);
        Polynome u, v;
        Polynome::pgcdEtendu(p, dMoinsEtoile, u, v); // u p + v dMoinsEtoile = 1
        const Polynome b = reste(u * a, dMoinsEtoile);
        const Polynome c = quotient(a - b * p, dMoinsEtoile);
        a = c - b.derivee() * quotient(dEtoile, dMoinsEtoile);
        numerateur = numerateur + b * cofacteur;
        cofacteur = cofacteur * dMoinsEtoile;
        dMoins = dMoins2;
    }
    g = numerateur.estNul() ? Fraction() : Fraction(std::move(numerateur), denominateur);
    a2 = std::move(a);
}

// Vrai si p (degré >= 1) n'a aucune racine réelle
bool sansRacineReelle(const Polynome& p) {
    if (p.degre() % 2 == 1) return false;
    if (p.degre() == 2) {
        const Nombre& b = p.coefficient(1);
        return (b * b - Nombre(4) * p.coefficient(2) * p.coefficient(0)).signe() < 0;
    }
    Nombre borne(0); // borne de Cauchy : toutes les racines sont dans ]-borne, borne]
    const Nombre inverse = p.dominant().inverse();
    for (int i = 0; i < p.degre(); ++i) {
        const Nombre v = (p.coefficient(i) * inverse).valeurAbsolue();
        if (v > borne) borne = v;
    }
    borne = borne + Nombre(1);
    return p.nombreRacinesReelles(-borne, borne) == 0;
}

// Logarithmes c ln|argument| accumulés, regroupés par coefficient
struct Logarithmes {
    struct Groupe {
        Nombre coefficient;            // > 0
        Polynome dessus, dessous;      // c ln|dessus / dessous|
        bool positif;                  // dessus et dessous sans racine réelle
    };
    std::vector<Groupe> groupes;
    // Quand la variable est t = exp(u) : ln|t| = u, tenu à part des autres logarithmes
    bool variableExponentielle = false;
    Nombre coefficientVariable;
    bool quartiques = false; // les facteurs de degré 4 peuvent être scindés sur Q(sqrt(m))

    void ajouter(const Nombre& c, const Polynome& argument, bool positif) {
        if (c.estZero()) return;
        if (variableExponentielle && argument.degre() == 1 && argument.coefficient(0).estZero()) {
            coefficientVariable = coefficientVariable + c;
            return;
        }
        const Nombre module = c.valeurAbsolue();
        for (Groupe& g : groupes) {
            if (g.coefficient != module) continue;
            Polynome& cible = c.signe() > 0 ? g.dessus : g.dessous;
            cible = cible * argument;
            g.positif = g.positif && positif;
            return;
        }
        groupes.push_back({module, c.signe() > 0 ? argument : polynomeUn(), c.signe() > 0 ? polynomeUn() : argument, positif});
    }

    void versSomme(const ExprPtr& x, AccumulateurSomme& somme) const {
        if (!coefficientVariable.estZero()) somme.ajouter(static_cast<const Exponentielle&>(*x).m_argument, coefficientVariable);
        for (const Groupe& g : groupes) {
            // c ln|P/Q|, ou -c ln|Q| s'il n'y a que le dénominateur
            const bool inverse = g.dessus.degre() == 0;
            const ExprPtr argument = inverse ? g.dessous.versExpression(x)
                                             : g.dessus.versExpression(x) / g.dessous.versExpression(x);
            somme.ajouter(g.positif ? ast_ln(argument) : lnAbs(argument), inverse ? -g.coefficient : g.coefficient);
        }
    }
};

// ----------------------------------------------------------------------------
// Facteur irréductible de degré 4 : factorisation réelle sur un corps quadratique
// ----------------------------------------------------------------------------

// a + b sqrt(m), m > 0 rationnel non carré fixé par le contexte
struct Quad {
    Nombre a, b;
};

struct CorpsQuadratique {
    Nombre m;

    Quad plus(const Quad& u, const Quad& v) const { return {u.a + v.a, u.b + v.b}; }
    Quad moins(const Quad& u, const Quad& v) const { return {u.a - v.a, u.b - v.b}; }
    Quad fois(const Quad& u, const Quad& v) const { return {u.a * v.a + u.b * v.b * m, u.a * v.b + u.b * v.a}; }
    Quad fois(const Quad& u, const Nombre& k) const { return {u.a * k, u.b * k}; }
    Quad inverse(const Quad& u) const {
        const Nombre norme = (u.a * u.a - u.b * u.b * m).inverse();
        return {u.a * norme, -(u.b * norme)};
    }
    static Quad conjugue(const Quad& u) { return {u.a, -u.b}; }
    int signe(const Quad& u) const {
        const int sa = u.a.signe(), sb = u.b.signe();
        if (sb == 0 || sa == sb) return sa;
        if (sa == 0) return sb;
        return sa * (u.a * u.a - u.b * u.b * m).signe(); // signes opposés : le plus grand en module l'emporte
    }
    // a + b sqrt(m) en expression
    ExprPtr expression(const Quad& u) const { return nombre(u.a) + nombre(u.b) * ast_pow(nombre(m), frac(1, 2)); }
};

/*
 * Nom : factoriserQuartique
 * Description : Pour z^4 + P z^2 + Q z + R irréductible sur Q, cherche m > 0 rationnel et une
 *               factorisation réelle (z^2 + e1 z + e0)(z^2 + e1' z + e0'), e1 et e0 dans
 *               Q(sqrt(m)), le second facteur étant le conjugué du premier. Avec
 *               (z^2 + alpha z + beta)(z^2 - alpha z + delta) : beta + delta - alpha^2 = P,
 *               alpha (delta - beta) = Q, beta delta = R, et alpha^2 est racine de la
 *               résolvante cubique t^3 + 2P t^2 + (P^2 - 4R) t - Q^2. Renvoie faux si la
 *               factorisation réelle demande des nombres algébriques de degré > 2.
 */
bool factoriserQuartique(const Nombre& P, const Nombre& Q, const Nombre& R, Nombre& m, Quad& e1, Quad& e0) {
    const Nombre demi = Nombre::rationnel(1, 2);
    Nombre racine;
    if (Q.estZero()) {
        // alpha = 0 : beta et delta racines de t^2 - P t + R
        const Nombre D = P * P - Nombre(4) * R;
        if (D.signe() > 0 && !D.racineExacte(2, racine)) {
            m = D;
            e1 = {Nombre(0), Nombre(0)};
            e0 = {P * demi, demi};
            return true;
        }
        // beta = delta = rho, rho^2 = R, alpha^2 = 2 rho - P
        if (R.signe() > 0 && R.racineExacte(2, racine)) {
            for (const Nombre& rho : {racine, -racine}) {
                const Nombre carre = Nombre(2) * rho - P;
                Nombre inutile;
                if (carre.signe() <= 0 || carre.racineExacte(2, inutile)) continue;
                m = carre;
                e1 = {Nombre(0), Nombre(1)};
                e0 = {rho, Nombre(0)};
                return true;
            }
        }
        return false;
    }
    const Polynome resolvante({-(Q * Q), P * P - Nombre(4) * R, Nombre(2) * P, Nombre(1)});
    for (const RacineReelle& r : resolvante.racinesReelles()) {
        const Constante* valeur = r.exacte ? comme<Constante>(r.exacte) : nullptr;
        if (!valeur || valeur->getNombre().signe() <= 0) continue;
        m = valeur->getNombre(); // alpha = sqrt(m)
        if (m.racineExacte(2, racine)) continue; // le polynôme serait réductible
        e1 = {Nombre(0), Nombre(1)};
        e0 = {(P + m) * demi, -(Q / (Nombre(2) * m))}; // beta = (P + m - Q/alpha)/2
        return true;
    }
    return false;
}

/*
 * Nom : elementQuartique
 * Description : int a/p pour p irréductible de degré 4 (deg a < 4) quand p se factorise en
 *               deux trinômes réels conjugués q et q' sur Q(sqrt(m)) :
 *               a/p = (A z + B)/q + conjugué, z = x + p3/(4 p4), puis pour chaque trinôme
 *               (A/2) ln|q| + (B - A e1/2) int dz/q. Les deux logarithmes se regroupent en
 *               Re(A)/2 ln|p| + Im(A) sqrt(m)/2 ln|q/q'|.
 */
bool elementQuartique(const Polynome& a, const Polynome& p, const ExprPtr& x, Logarithmes& logs, AccumulateurSomme& somme) {
    const Nombre dominant = p.dominant();
    const Nombre h = p.coefficient(3) / (Nombre(4) * dominant); // z = x + h
    // p(z - h)/dominant = z^4 + P z^2 + Q z + R et n(z) = a(z - h)/dominant
    const Fraction decalage(Polynome({-h, Nombre(1)}));
    const Polynome pz = Fraction(p * dominant.inverse()).composer(decalage).numerateur();
    const Polynome n = Fraction(a * dominant.inverse()).composer(decalage).numerateur();
    CorpsQuadratique K;
    Quad e1, e0;
    if (!factoriserQuartique(pz.coefficient(2), pz.coefficient(1), pz.coefficient(0), K.m, e1, e0)) return false;

    // Dans K[z]/(q), q = z^2 + e1 z + e0 : éléments r1 z + r0
    struct Residu {
        Quad r1, r0;
    };
    const auto foisZ = [&](const Residu& r) { // z^2 = -e1 z - e0
        return Residu{K.moins(r.r0, K.fois(r.r1, e1)), K.fois(K.fois(r.r1, e0), Nombre(-1))};
    };
    const auto produit = [&](const Residu& u, const Residu& v) {
        const Residu haut = foisZ(Residu{K.fois(u.r1, v.r1), {Nombre(0), Nombre(0)}}); // u1 v1 z^2
        return Residu{K.plus(haut.r1, K.plus(K.fois(u.r1, v.r0), K.fois(u.r0, v.r1))), K.plus(haut.r0, K.fois(u.r0, v.r0))};
    };
    // n mod q (schéma de Horner)
    Residu reste{{Nombre(0), Nombre(0)}, {Nombre(0), Nombre(0)}};
    for (int i = n.degre(); i >= 0; --i) {
        reste = foisZ(reste);
        reste.r0.a = reste.r0.a + n.coefficient(i);
    }
    // Inverse de q' mod q : q' = q + (c1 z + c0), c = conj(e) - e
    const Quad c1 = K.moins(CorpsQuadratique::conjugue(e1), e1), c0 = K.moins(CorpsQuadratique::conjugue(e0), e0);
    // (c1 z + c0)(u z + v) = 1 : (c0 - c1 e1) u + c1 v = 0, -c1 e0 u + c0 v = 1
    const Quad m11 = K.moins(c0, K.fois(c1, e1)), m21 = K.fois(K.fois(c1, e0), Nombre(-1));
    const Quad determinant = K.plus(K.fois(m11, c0), K.fois(K.fois(c1, c1), e0));
    if (determinant.a.estZero() && determinant.b.estZero()) return false;
    const Quad inverseDet = K.inverse(determinant);
    (void)m21;
    const Residu inverse{K.fois(K.fois(c1, Nombre(-1)), inverseDet), K.fois(m11, inverseDet)};
    const Residu AB = produit(reste, inverse); // A z + B

    const ExprPtr racineM = ast_pow(nombre(K.m), frac(1, 2));
    const Nombre demi = Nombre::rationnel(1, 2);
    ExprPtr trinomes[2];
    bool positifs[2];
    for (int i = 0; i < 2; ++i) {
        const auto conj = [i](const Quad& u) { return i == 0 ? u : CorpsQuadratique::conjugue(u); };
        const Quad E1 = conj(e1), E0 = conj(e0), A = conj(AB.r1), B = conj(AB.r0);
        // q(x) = (x + h)^2 + E1 (x + h) + E0
        const Quad lineaire = K.plus(E1, {Nombre(2) * h, Nombre(0)});
        const Quad constant = K.plus(K.plus(K.fois(E1, h), E0), {h * h, Nombre(0)});
        AccumulateurSomme q;
        q.ajouter(ast_pow(x, 2.0), Nombre(1));
        q.ajouter(x, lineaire.a);
        q.ajouter(racineM * x, lineaire.b);
        q.ajouter(K.expression(constant), Nombre(1));
        trinomes[i] = q.construire();
        // discriminant/4 = E1^2/4 - E0
        const Quad quartDelta = K.moins(K.fois(K.fois(E1, E1), Nombre::rationnel(1, 4)), E0);
        const int signe = K.signe(quartDelta);
        positifs[i] = signe < 0;
        const Quad k = K.moins(B, K.fois(K.fois(A, E1), demi));
        if (k.a.estZero() && k.b.estZero()) continue;
        if (signe == 0) return false;
        const ExprPtr centre = x + K.expression(K.plus(K.fois(E1, demi), {h, Nombre(0)}));
        const ExprPtr s = ast_pow(K.expression(signe < 0 ? K.fois(quartDelta, Nombre(-1)) : quartDelta), frac(1, 2));
        if (signe < 0) {
            somme.ajouter(K.expression(k) * ast_atan(centre / s) / s, Nombre(1));
        } else {
            somme.ajouter(K.expression(k) * ast_ln(ast_abs((centre - s) / (centre + s))) / s, demi);
        }
    }
    // (A/2) ln|q| + (A'/2) ln|q'| = Re(A)/2 ln|p| + Im(A) sqrt(m)/2 ln|q/q'|
    logs.ajouter(AB.r1.a * demi, p, positifs[0] && positifs[1]);
    if (!AB.r1.b.estZero()) {
        const ExprPtr rapport = trinomes[0] / trinomes[1];
        somme.ajouter(racineM * (positifs[0] && positifs[1] ? ast_ln(rapport) : ast_ln(ast_abs(rapport))), AB.r1.b * demi);
    }
    return true;
}

/*
 * Nom : elementSimple
 * Description : int a/p pour un facteur p irréductible sur Q, primitif (deg a < deg p).
 *               Renvoie faux si la primitive demande des nombres algébriques de degré > 2.
 *                 degré 1, p = alpha x + beta : (a/alpha) ln|p| ;
 *                 degré 2, p = A x^2 + B x + C, a = m x + n : (m/2A) ln|p| + k int dx/p,
 *                   k = n - m B/(2A), et avec h = B/(2A), s = sqrt(|B^2 - 4AC|)/(2A) :
 *                   int dx/p = atan((x + h)/s)/(A s)                 (discriminant < 0)
 *                            = ln|(x + h - s)/(x + h + s)|/(2 A s)   (discriminant > 0) ;
 *                 degré >= 3 : c ln|p| si a = c p' ;
 *                 degré 4 : deux trinômes réels conjugués sur Q(sqrt(m)), s'ils existent.
 */
bool elementSimple(const Polynome& a, const Polynome& p, const ExprPtr& x, Logarithmes& logs, AccumulateurSomme& somme) {
    if (a.estNul()) return true;
    if (p.degre() == 1) {
        logs.ajouter(a.coefficient(0) / p.coefficient(1), p, false);
        return true;
    }
    if (p.degre() == 2) {
        const Nombre& A = p.coefficient(2);
        const Nombre& B = p.coefficient(1);
        const Nombre delta = B * B - Nombre(4) * A * p.coefficient(0);
        const Nombre deuxA = Nombre(2) * A;
        const Nombre& m = a.coefficient(1);
        logs.ajouter(m / deuxA, p, delta.signe() < 0);
        const Nombre k = a.coefficient(0) - m * B / deuxA;
        if (k.estZero()) return true;
        const ExprPtr centre = x + nombre(B / deuxA);
        const ExprPtr s = ast_pow(nombre(delta.valeurAbsolue() / (deuxA * deuxA)), frac(1, 2));
        if (delta.signe() < 0) {
            somme.ajouter(ast_atan(centre / s) / s, k / A);
        } else {
            somme.ajouter(ast_ln(ast_abs((centre - s) / (centre + s))) / s, k / deuxA);
        }
        return true;
    }
    Nombre c;
    if (proportionnels(a, p.derivee(), c)) {
        logs.ajouter(c, p, sansRacineReelle(p));
        return true;
    }
    return logs.quartiques && p.degre() == 4 && elementQuartique(a, p, x, logs, somme);
}

/*
 * Nom : partieLogarithmique
 * Description : int a/d pour d sans carré (deg a < deg d) : ajoute les logarithmes et
 *               arctangentes à somme, et à nonIntegre les éléments simples sans primitive
 *               à coefficients rationnels ou quadratiques.
 */
void partieLogarithmique(const Polynome& a, const Polynome& d, const ExprPtr& x, AccumulateurSomme& somme, Fraction& nonIntegre,
                         bool quartiques) {
    if (a.estNul()) return;
    Nombre contenu;
    const Polynome dPrimitif = d.partiePrimitive(contenu); // d = contenu * dPrimitif
    const Polynome numerateur = a * contenu.inverse();     // a/d = numerateur / dPrimitif
    Logarithmes logs;
    logs.variableExponentielle = x->type() == TypeNoeud::Exponentielle;
    logs.quartiques = quartiques;
    Nombre c;
    if (proportionnels(numerateur, dPrimitif.derivee(), c)) {
        // Dérivée logarithmique : c ln|d|, sans factoriser
        logs.ajouter(c, dPrimitif, sansRacineReelle(dPrimitif));
    } else {
        const auto facteurs = dPrimitif.facteursIrreductibles(c);
        AccumulateurSomme elements; // pour le journal seulement
        for (const auto& [facteur, multiplicite] : facteurs) {
            // numerateur/dPrimitif = sum a_j/p_j : a_j = numerateur (dPrimitif/p_j)^-1 mod p_j
            Polynome aj = numerateur;
            if (facteurs.size() > 1) {
                Polynome u, v;
                Polynome::pgcdEtendu(quotient(dPrimitif, facteur), facteur, u, v);
                aj = reste(numerateur * u, facteur);
                journal::noter([&](Journal&) { elements.ajouter(Fraction(aj, facteur).versExpression(x), Nombre(1)); });
            }
            if (!elementSimple(aj, facteur, x, logs, somme)) nonIntegre = nonIntegre + Fraction(aj, facteur);
        }
        if (facteurs.size() > 1) {
            journal::noter([&](Journal& j) {
                j.ajouter("integration.elements_simples", "Décomposition en éléments simples",
                          Fraction(numerateur, dPrimitif).versExpression(x), elements.construire());
            });
        }
    }
    logs.versSomme(x, somme);
}

} // namespace

namespace {

/*
 * Nom : substitutionMonomiale
 * Description : f(x) dx = h(u) du avec u = x^k, k >= 2 le plus grand possible : c'est le cas
 *               quand x f(x) est une fraction en x^k (x/(x^4 + 1), x^2/(x^6 + 1)). Renvoie
 *               faux s'il n'y a pas de tel k.
 */
bool substitutionMonomiale(const Fraction& f, Fraction& h, int& k) {
    const Fraction g = f * Fraction::variable();
    k = 0;
    for (const Polynome* p : {&g.numerateur(), &g.denominateur()}) {
        for (int i = 1; i <= p->degre(); ++i) {
            if (!p->coefficient(i).estZero()) k = std::gcd(k, i);
        }
    }
    if (k < 2) return false;
    auto enU = [k](const Polynome& p) {
        std::vector<Nombre> c;
        for (int i = 0; i <= p.degre(); i += k) c.push_back(p.coefficient(i));
        return Polynome(std::move(c));
    };
    // g(x) = G(u), dx/x = du/(k u)
    h = Fraction(enU(g.numerateur()), enU(g.denominateur()) * Polynome({Nombre(0), Nombre(k)}));
    return true;
}

} // namespace

namespace {

/*
 * Nom : integrerFraction
 * Description : Coeur de primitiveFraction : partie rationnelle (polynôme + réduction
 *               d'Hermite) dans rationnelle, logarithmes et arctangentes ajoutés à somme, ce
 *               qui n'a pas pu être intégré dans nonIntegre. quartiques autorise la
 *               factorisation réelle des facteurs irréductibles de degré 4.
 */
void integrerFraction(const Fraction& f, const ExprPtr& x, Fraction& rationnelle, AccumulateurSomme& somme, Fraction& nonIntegre,
                      bool quartiques) {
    Polynome q, r;
    f.numerateur().diviser(f.denominateur(), q, r);
    rationnelle = Fraction();
    if (!r.estNul()) {
        Polynome a2, dEtoile, q2, r2;
        hermite(r, f.denominateur(), rationnelle, a2, dEtoile);
        a2.diviser(dEtoile, q2, r2);
        q = q + q2;
        if (!rationnelle.estNulle()) {
            journal::noter([&](Journal& j) {
                j.ajouter("integration.hermite",
                          "Réduction d'Hermite : les pôles multiples donnent une fraction rationnelle, il reste un "
                          "dénominateur sans facteur carré",
                          Fraction(r, f.denominateur()).versExpression(x),
                          rationnelle.versExpression(x) + fabriquer<IntegraleNonEvaluee>(Fraction(a2, dEtoile).versExpression(x)));
            });
        }
        partieLogarithmique(r2, dEtoile, x, somme, nonIntegre, quartiques);
    }
    if (!q.estNul()) rationnelle = rationnelle + Fraction(primitive(q));
}

} // namespace

ExprPtr primitiveFraction(const Fraction& f, const ExprPtr& x, bool* complete) {
    AccumulateurSomme somme;
    Fraction rationnelle, nonIntegre;
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    // 1. Facteurs irréductibles de degré 1 et 2
    integrerFraction(f, x, rationnelle, somme, nonIntegre, false);
    if (!nonIntegre.estNulle()) {
        // 2. u = x^k abaisse le degré des facteurs irréductibles (x/(x^4 + 1) : atan(x^2)/2)
        Fraction h;
        int k;
        if (substitutionMonomiale(f, h, k)) {
            if (j) j->tronquer(etapes); // la tentative en x laisse la place au changement de variable
            const ExprPtr u = ast_pow(x, nombre(Nombre(k)));
            journal::noter([&](Journal& jr) {
                jr.ajouter("integration.substitution", "Changement de variable u = " + u->texte(), f.versExpression(x),
                           h.versExpression(var("u")));
            });
            bool reussie;
            ExprPtr resultat;
            {
                const journal::SousNiveau sous;
                resultat = primitiveFraction(h, u, &reussie);
            }
            if (reussie) {
                if (complete) *complete = true;
                return resultat;
            }
        }
        // 3. Facteurs de degré 4 scindés sur un corps quadratique réel
        if (j) j->tronquer(etapes);
        somme = AccumulateurSomme();
        nonIntegre = Fraction();
        integrerFraction(f, x, rationnelle, somme, nonIntegre, true);
        if (!nonIntegre.estNulle()) somme.ajouter(fabriquer<IntegraleNonEvaluee>(nonIntegre.versExpression(x)), Nombre(1));
    }
    if (complete) *complete = nonIntegre.estNulle();
    // Partie polynomiale et pôles multiples séparés : x^2/2 + 3*x - (8*x - 7)/(2*(x - 1)^2)
    Polynome q, r;
    rationnelle.numerateur().diviser(rationnelle.denominateur(), q, r);
    if (!q.estNul()) somme.ajouter(q.versExpression(x), Nombre(1));
    if (!r.estNul()) somme.ajouter(Fraction(r, rationnelle.denominateur()).versExpression(x), Nombre(1));
    return somme.construire();
}

bool primitiveFractionSeparee(const Fraction& f, const ExprPtr& x, Fraction& rationnelle, ExprPtr& transcendante) {
    for (const bool quartiques : {false, true}) {
        AccumulateurSomme somme;
        Fraction nonIntegre;
        Journal* const j = journal::actif();
        const std::size_t etapes = j ? j->etapes().size() : 0;
        integrerFraction(f, x, rationnelle, somme, nonIntegre, quartiques);
        if (nonIntegre.estNulle()) {
            transcendante = somme.construire();
            return true;
        }
        if (j) j->tronquer(etapes);
    }
    return false;
}

ExprPtr parFractionRationnelle(const ExprPtr& f) {
    const TypeNoeud type = f->type();
    if (type != TypeNoeud::Produit && type != TypeNoeud::Puissance) return nullptr;
    const ExprPtr x = variableDe(*f);
    Fraction fraction;
    // Un produit de polynômes est laissé au changement de variable (x^2 (x^3 + 1)^5), puis
    // au développement
    if (!x || !versFraction(f, x, fraction) || fraction.estPolynome()) return nullptr;
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    const ExprPtr resultat = primitiveFraction(fraction, x);
    // Rien n'a pu être intégré : aux méthodes suivantes d'essayer
    if (resultat->type() == TypeNoeud::IntegraleNonEvaluee) {
        if (j) j->tronquer(etapes);
        return nullptr;
    }
    return resultat;
}

ExprPtr parDeveloppement(const ExprPtr& f) {
    const TypeNoeud type = f->type();
    if (type != TypeNoeud::Produit && type != TypeNoeud::Puissance) return nullptr;
    const ExprPtr g = developper(f);
    if (g.get() == f.get() || g->type() != TypeNoeud::Somme) return nullptr;
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    journal::noter([&](Journal& jr) { jr.ajouter("integration.developpement", "Développement de l'intégrande", f, g); });
    ExprPtr G;
    {
        const journal::SousNiveau sous;
        G = g->integrer();
    }
    if (contientNonEvaluee(*G)) {
        if (j) j->tronquer(etapes);
        return nullptr;
    }
    return G;
}

} // namespace symalgo
