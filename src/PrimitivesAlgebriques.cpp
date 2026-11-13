/**
 * @file PrimitivesAlgebriques.cpp
 * @brief Intégration des fractions rationnelles de x et d'un radical.
 *
 * 1. Radical d'une fonction affine, R(x, (a x + b)^(p/q)) : t = (a x + b)^(1/q) donne
 *    x = (t^q - b)/a, dx = (q/a) t^(q-1) dt et une fraction rationnelle en t.
 *
 * 2. Racine carrée d'un trinôme, y = sqrt(Q), Q = a x^2 + b x + c. Toute fraction
 *    rationnelle de x et y a une forme unique A(x) + B(x) y (ElementQuadratique, y^2 = Q).
 *    int A est une intégrale de fraction rationnelle ; pour int B y = int (B Q)/y :
 *      - si N = B Q est un polynôme : int N/y = S y + lambda int dx/y, où le polynôme S
 *        (deg S = deg N - 1) et lambda se calculent du plus haut degré au plus bas en
 *        identifiant N = S' Q + S Q'/2 + lambda ;
 *      - sinon, substitution d'Euler :
 *          a > 0 : y = sqrt(a) y1, y1 = sqrt(x^2 + b x + c), t = x + y1,
 *                  x = (t^2 - c)/(2t + b), y1 = (t^2 + b t + c)/(2t + b) ;
 *          Q à racines rationnelles alpha, beta : t = y/(x - alpha),
 *                  x = (alpha t^2 - a beta)/(t^2 - a), y = t (x - alpha).
 *        La partie rationnelle de la primitive en t est récrite sous la forme A(x) + B(x) y
 *        en l'évaluant dans l'extension quadratique ; seuls les logarithmes et arctangentes
 *        gardent l'expression de t.
 */

#include "Canonique.hpp"
#include "Journal.hpp"
#include "Primitives.hpp"

#include <numeric>
#include <stdexcept>
#include <vector>

namespace symalgo {

namespace {

constexpr int DEGRE_MAX = 96;
constexpr long long ORDRE_MAX = 12; // ordre q du radical

// Les puissances d'exposant non entier de f : une seule base, ppcm des dénominateurs
struct Radicaux {
    const ASTNode* base = nullptr;
    long long ordre = 1;
};

bool chercherPuissance(const ASTNode& base, const ASTNode& exposant, Radicaux& r);

// Faux si f n'est pas une fraction rationnelle de x et des puissances fractionnaires d'une même base
bool chercher(const ASTNode& e, Radicaux& r) {
    if (!e.contientVariable()) return true;
    switch (e.type()) {
        case TypeNoeud::Variable: return true;
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) {
                if (!chercher(*t.expression, r)) return false;
            }
            return true;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                if (!chercherPuissance(*f.base, *f.exposant, r)) return false;
            }
            return true;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return chercherPuissance(*p.getBase(), *p.getExposant(), r);
        }
        default: return false;
    }
}

bool chercherPuissance(const ASTNode& base, const ASTNode& exposant, Radicaux& r) {
    if (!base.contientVariable()) return !exposant.contientVariable();
    const Constante* c = comme<Constante>(&exposant);
    long long num, den;
    if (!c || !c->getNombre().versFraction(num, den)) return false;
    if (den == 1) return chercher(base, r);
    if (r.base && r.base != &base) return false;
    r.base = &base;
    r.ordre = std::lcm(r.ordre, den);
    return r.ordre <= ORDRE_MAX;
}

// ----------------------------------------------------------------------------
// 1. Radical d'une fonction affine
// ----------------------------------------------------------------------------

struct ReglesAffines {
    const ASTNode* x;
    const ASTNode* base;
    long long ordre;
    Fraction xDeT; // (t^q - b)/a

    Feuille feuille(const ASTNode& e, ElementFraction& r) const {
        if (&e == x) {
            r.f = xDeT;
            return Feuille::Convertie;
        }
        if (&e == base) { // a x + b = t^q
            r.f = Fraction::variable().puissance(ordre);
            return Feuille::Convertie;
        }
        return e.type() == TypeNoeud::Variable ? Feuille::Refusee : Feuille::Descendre;
    }

    bool puissance(const ASTNode& b, const ASTNode& exposant, ElementFraction& r) const {
        long long num, den;
        if (&b != base || !static_cast<const Constante&>(exposant).getNombre().versFraction(num, den)) return false;
        const long long k = num * (ordre / den);
        if (k > DEGRE_MAX || k < -DEGRE_MAX) return false;
        r.f = Fraction::variable().puissance(k);
        return true;
    }
};

ExprPtr radicalAffine(const ExprPtr& f, const ExprPtr& x, const Radicaux& radicaux, const Polynome& base) {
    const Nombre& a = base.coefficient(1);
    const Nombre& b = base.coefficient(0);
    const long long q = radicaux.ordre;
    const Fraction t = Fraction::variable();
    ReglesAffines regles{x.get(), radicaux.base, q, (t.puissance(q) - Fraction(b)) * a.inverse()};
    ElementFraction h;
    try {
        if (!convertirRationnel(*f, h, regles, ElementFraction{Fraction(Nombre(1))}, DEGRE_MAX)) return nullptr;
    } catch (const std::domain_error&) {
        return nullptr;
    }
    const Fraction integrande = h.f * t.puissance(q - 1) * (Nombre(q) / a);
    const ExprPtr variable = ast_pow(radicaux.base->clone(), nombre(Nombre::rationnel(1, q)));
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    journal::noter([&](Journal& jr) {
        jr.ajouter("integration.substitution", "Changement de variable t = " + variable->texte(), f, integrande.versExpression(var("t")));
    });
    bool complete;
    ExprPtr F;
    {
        const journal::SousNiveau sous;
        F = primitiveFraction(integrande, variable, &complete);
    }
    if (!complete) {
        if (j) j->tronquer(etapes);
        return nullptr;
    }
    return F;
}

// ----------------------------------------------------------------------------
// 2. Racine carrée d'un trinôme
// ----------------------------------------------------------------------------

struct ReglesQuadratiques {
    const ASTNode* x;
    const ASTNode* base;
    const Fraction* carre; // Q

    Feuille feuille(const ASTNode& e, ElementQuadratique& r) const {
        if (&e == x) {
            r = {Fraction::variable(), Fraction(), carre};
            return Feuille::Convertie;
        }
        if (&e == base) {
            r = {*carre, Fraction(), carre};
            return Feuille::Convertie;
        }
        return e.type() == TypeNoeud::Variable ? Feuille::Refusee : Feuille::Descendre;
    }

    // Q^(k/2), k impair : Q^((k - 1)/2) y
    bool puissance(const ASTNode& b, const ASTNode& exposant, ElementQuadratique& r) const {
        long long num, den;
        if (&b != base || !static_cast<const Constante&>(exposant).getNombre().versFraction(num, den) || den != 2) return false;
        if (num > DEGRE_MAX || num < -DEGRE_MAX) return false;
        r = {Fraction(), carre->puissance((num - 1) / 2), carre};
        return true;
    }
};

/*
 * Nom : coefficientsIndetermines
 * Description : int N/y dx = S y + lambda int dx/y pour un polynôme N (y^2 = Q de degré 2) :
 *               le coefficient de x^j de S' Q + S Q'/2 est a j s_(j-1) + b (j + 1/2) s_j + c (j + 1) s_(j+1).
 */
void coefficientsIndetermines(const Polynome& N, const Polynome& Q, Polynome& S, Nombre& lambda) {
    const Nombre& a = Q.coefficient(2);
    const Nombre& b = Q.coefficient(1);
    const Nombre& c = Q.coefficient(0);
    const int n = N.degre();
    std::vector<Nombre> s(static_cast<std::size_t>(n + 2), Nombre(0));
    const Nombre demi = Nombre::rationnel(1, 2);
    for (int j = n; j >= 1; --j) {
        const std::size_t i = static_cast<std::size_t>(j);
        s[i - 1] = (N.coefficient(j) - b * (Nombre(j) + demi) * s[i] - c * Nombre(j + 1) * s[i + 1]) / (a * Nombre(j));
    }
    lambda = N.coefficient(0) - b * demi * s[0] - c * s[1];
    s.resize(static_cast<std::size_t>(n > 0 ? n : 0));
    S = Polynome(std::move(s));
}

// A(x) + B(x) racine, en expression
ExprPtr versExpression(const ElementQuadratique& e, const ExprPtr& x, const ExprPtr& racine) {
    AccumulateurSomme s;
    if (!e.a.estNulle()) s.ajouter(e.a.versExpression(x), Nombre(1));
    if (!e.b.estNulle()) s.ajouter(e.b.versExpression(x) * racine, Nombre(1));
    return s.construire();
}

/*
 * Nom : parEuler
 * Description : int B y dx par une substitution d'Euler rationnelle sur Q (voir l'en-tête du
 *               fichier). Pointeur nul si aucune ne s'applique (a < 0 et racines irrationnelles)
 *               ou si la fraction en t ne s'intègre pas complètement.
 */
ExprPtr parEuler(const Fraction& B, const Polynome& Q, const ExprPtr& x, const ExprPtr& y) {
    const Nombre& a = Q.coefficient(2);
    const Fraction t = Fraction::variable(), X = Fraction::variable();
    Fraction xDeT, yDeT, carre;         // x(t), la racine en fonction de t, son carré en x
    ElementQuadratique tDeX;            // t dans l'extension
    ExprPtr racine, facteur = un();     // racine utilisée, facteur constant sorti de l'intégrale
    if (a.signe() > 0) {
        const Polynome Q1 = Q * a.inverse(); // x^2 + b x + c
        const Fraction b(Q1.coefficient(1)), c(Q1.coefficient(0));
        const Fraction d = t * Nombre(2) + b;
        xDeT = (t * t - c) / d;
        yDeT = (t * t + b * t + c) / d;
        carre = Fraction(Q1);
        tDeX = {X, Fraction(Nombre(1)), &carre};
        facteur = ast_pow(nombre(a), frac(1, 2));
        racine = y / facteur;
    } else {
        Nombre racineDelta;
        const Nombre& b = Q.coefficient(1);
        const Nombre delta = b * b - Nombre(4) * a * Q.coefficient(0);
        if (delta.signe() <= 0 || !delta.racineExacte(2, racineDelta)) return nullptr;
        const Nombre alpha = (-b + racineDelta) / (Nombre(2) * a), beta = (-b - racineDelta) / (Nombre(2) * a);
        xDeT = (t * t * alpha - Fraction(a * beta)) / (t * t - Fraction(a));
        yDeT = t * (xDeT - Fraction(alpha));
        carre = Fraction(Q);
        tDeX = {Fraction(), Fraction(Nombre(1)) / (X - Fraction(alpha)), &carre};
        racine = y;
    }
    const Fraction h = B.composer(xDeT) * yDeT * xDeT.derivee();
    if (h.depasse(2 * DEGRE_MAX)) return nullptr;
    const ExprPtr variable = versExpression(tDeX, x, racine);
    Fraction rationnelle;
    ExprPtr transcendante;
    journal::noter([&](Journal& j) {
        j.ajouter("integration.euler", "Substitution d'Euler t = " + variable->texte(), B.versExpression(x) * y,
                  facteur * h.versExpression(var("t")));
    });
    {
        const journal::SousNiveau sous;
        if (!primitiveFractionSeparee(h, variable, rationnelle, transcendante)) return nullptr;
    }
    ElementQuadratique partie = tDeX.evaluer(rationnelle);
    // Constante d'intégration nulle : le terme constant de la partie polynomiale est retiré
    Polynome quotient, reste;
    partie.a.numerateur().diviser(partie.a.denominateur(), quotient, reste);
    if (!quotient.estNul()) partie.a = partie.a - Fraction(quotient.coefficient(0));
    return facteur * (versExpression(partie, x, racine) + transcendante);
}

ExprPtr radicalQuadratique(const ExprPtr& f, const ExprPtr& x, const Radicaux& radicaux, const Polynome& Q) {
    const Fraction carre(Q);
    ReglesQuadratiques regles{x.get(), radicaux.base, &carre};
    const ElementQuadratique unite{Fraction(Nombre(1)), Fraction(), &carre};
    ElementQuadratique e = unite;
    try {
        if (!convertirRationnel(*f, e, regles, unite, DEGRE_MAX)) return nullptr;
    } catch (const std::domain_error&) {
        return nullptr;
    }
    const ExprPtr y = ast_pow(radicaux.base->clone(), frac(1, 2));
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    auto echec = [&]() {
        if (j) j->tronquer(etapes);
        return ExprPtr();
    };
    const Nombre& a = Q.coefficient(2);
    const Nombre delta = Q.coefficient(1) * Q.coefficient(1) - Nombre(4) * a * Q.coefficient(0);
    if (a.signe() < 0 && delta.signe() <= 0) return nullptr; // jamais de racine réelle
    if (delta.estZero()) {
        // Q = a (x + h)^2 : y = sqrt(a) |x + h|, et sur chaque demi-droite
        // int A + B y = int A + sqrt(a) sign(x + h) int B (x + h)
        const Fraction centre = Fraction::variable() + Fraction(Q.coefficient(1) / (Nombre(2) * a));
        bool completeA, completeB;
        const ExprPtr FA = primitiveFraction(e.a, x, &completeA);
        const ExprPtr FB = primitiveFraction(e.b * centre, x, &completeB);
        if (!completeA || !completeB) return echec();
        return FA + ast_pow(nombre(a), frac(1, 2)) * ast_sign(centre.versExpression(x)) * FB;
    }
    AccumulateurSomme resultat;
    if (!e.a.estNulle()) {
        bool complete;
        resultat.ajouter(primitiveFraction(e.a, x, &complete), Nombre(1));
        if (!complete) return echec();
    }
    if (!e.b.estNulle()) {
        const Fraction N = e.b * carre; // B y = N / y
        if (N.estPolynome()) {
            Polynome S;
            Nombre lambda;
            coefficientsIndetermines(N.numerateur(), Q, S, lambda);
            if (!S.estNul()) resultat.ajouter(S.versExpression(x) * y, Nombre(1));
            if (!lambda.estZero()) {
                const ExprPtr J = tableEtendue(un() / y); // asin ou ln
                if (!J) return echec();
                resultat.ajouter(J, lambda);
            }
        } else {
            const ExprPtr F = parEuler(e.b, Q, x, y);
            if (!F) return echec();
            resultat.ajouter(F, Nombre(1));
        }
    }
    return resultat.construire();
}

} // namespace

ExprPtr parRadicaux(const ExprPtr& f) {
    const TypeNoeud type = f->type();
    if (type != TypeNoeud::Produit && type != TypeNoeud::Puissance) return nullptr;
    Radicaux radicaux;
    if (!chercher(*f, radicaux) || !radicaux.base) return nullptr;
    const ExprPtr x = variableDe(*f);
    Polynome base;
    ExprPtr v;
    if (!x || !Polynome::depuisExpression(radicaux.base->clone(), base, &v) || !base.estExact() || v.get() != x.get()) return nullptr;
    try {
        if (base.degre() == 1) return radicalAffine(f, x, radicaux, base);
        if (base.degre() == 2 && radicaux.ordre == 2) return radicalQuadratique(f, x, radicaux, base);
    } catch (const std::domain_error&) { // intégrande nulle part définie : 1/(y - sqrt(a) (x + h)) pour Q carré
    }
    return nullptr;
}

} // namespace symalgo
