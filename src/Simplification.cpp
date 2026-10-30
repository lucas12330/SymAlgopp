/**
 * @file Simplification.cpp
 * @brief Simplification avec options et réduction des fractions rationnelles par PGCD.
 *
 * La réduction rationnelle est appelée par Produit::calculerSimplification (Regles.cpp) sur
 * chaque produit reconstruit. simplifier() servant à l'intégration, aux EDO, au gradient...,
 * le cas courant (produit sans quotient de polynômes) est rejeté par un simple parcours des
 * facteurs, sans allocation ; les polynômes ne sont construits qu'ensuite.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"
#include "Polynome.hpp"
#include "Simplification.hpp"

#include <climits>
#include <cstddef>
#include <utility>
#include <vector>

namespace symalgo {

ExprPtr ASTNode::simplifier(const OptionsSimplification& options) const {
    if (options.rationnel) return simplifier(); // options par défaut : mémo du noeud
    // Sans règle facultative, il ne reste que la forme canonique, déjà imposée à la
    // construction : reconstruire le noeud à partir de ses enfants redonnerait le même
    // noeud (hash-consing). Une future règle facultative demandera ici un parcours avec
    // cache local à l'appel (comme derivee), pour ne pas polluer le mémo par défaut.
    return clone();
}

namespace {

// Exposant entier non nul
bool exposantEntier(const ExprPtr& exposant, long long& k) {
    const Constante* c = comme<Constante>(exposant);
    return c && c->getNombre().versEntier(k) && k != 0 && k != LLONG_MIN; // -k doit exister
}

// Facteur pouvant être un polynôme en x : base x ou somme contenant x, exposant entier
bool facteurCandidat(const Facteur& f, long long& k) {
    const TypeNoeud t = f.base->type();
    return (t == TypeNoeud::Variable || (t == TypeNoeud::Somme && f.base->contientVariable())) &&
           exposantEntier(f.exposant, k);
}

// Deux bases candidates peuvent-elles avoir un facteur commun ? x ne divise une somme
// polynomiale que si son terme constant est nul (cas courant des dérivées : x/(x + 1))
bool facteurCommunPossible(const ASTNode& a, const ASTNode& b) {
    if (a.type() == TypeNoeud::Variable) {
        return b.type() == TypeNoeud::Somme && static_cast<const Somme&>(b).getConstante().estZero();
    }
    if (b.type() == TypeNoeud::Variable) {
        return static_cast<const Somme&>(a).getConstante().estZero();
    }
    return true;
}

// Filtre sans allocation : une base du numérateur et une du dénominateur pouvant se simplifier
bool reductionPossible(const std::vector<Facteur>& facteurs) {
    for (const Facteur& d : facteurs) {
        long long kd;
        if (!facteurCandidat(d, kd) || kd > 0) continue;
        for (const Facteur& n : facteurs) {
            long long kn;
            if (facteurCandidat(n, kn) && kn > 0 && facteurCommunPossible(*n.base, *d.base)) return true;
        }
    }
    return false;
}

// Facteur polynomial P^exposant (exposant > 0) ; origine : expression d'origine tant que
// P n'a pas été divisé (elle est alors gardée telle quelle), nullptr sinon
struct FacteurPolynomial {
    Polynome p;
    long long exposant;
    ExprPtr origine;
};

// P <- P / g (division exacte)
void diviserPar(FacteurPolynomial& f, const Polynome& g) {
    Polynome quotient, reste;
    f.p.diviser(g, quotient, reste);
    f.p = std::move(quotient);
    f.origine = ExprPtr();
}

// Divise deux à deux numérateurs et dénominateurs par leurs PGCD jusqu'à ce qu'ils soient
// premiers entre eux : P^a / Q^b = (P/g)^a / (Q/g)^b * g^(a - b), avec g = pgcd(P, Q).
// Le degré total diminue à chaque division, d'où la terminaison. Renvoie vrai si un
// facteur commun a été trouvé.
bool reduire(std::vector<FacteurPolynomial>& numerateurs, std::vector<FacteurPolynomial>& denominateurs) {
    bool reduit = false, change = true;
    while (change) {
        change = false;
        for (std::size_t i = 0; i < numerateurs.size(); ++i) {
            for (std::size_t j = 0; j < denominateurs.size() && numerateurs[i].p.degre() > 0; ++j) {
                if (denominateurs[j].p.degre() < 1) continue;
                const Polynome g = Polynome::pgcd(numerateurs[i].p, denominateurs[j].p);
                if (g.degre() < 1) continue;
                diviserPar(numerateurs[i], g);
                diviserPar(denominateurs[j], g);
                // Après division par le PGCD, P/g et Q/g sont premiers entre eux ; g^(a - b)
                // est un nouveau facteur, à comparer à son tour à l'autre côté
                const long long a = numerateurs[i].exposant, b = denominateurs[j].exposant;
                if (a > b) numerateurs.push_back({g, a - b, ExprPtr()});
                if (b > a) denominateurs.push_back({g, b - a, ExprPtr()});
                reduit = change = true;
            }
        }
    }
    return reduit;
}

// Multiplie le résultat par f^signe*exposant : facteur d'origine intact, constante, ou
// partie primitive (contenu remonté dans le coefficient)
void multiplierPar(AccumulateurProduit& resultat, const FacteurPolynomial& f, long long signe, const ExprPtr& x) {
    const long long exposant = signe * f.exposant;
    if (f.origine) {
        resultat.multiplier(ast_pow(f.origine, nombre(Nombre(exposant))));
        return;
    }
    if (f.p.degre() < 1) {
        resultat.multiplier(nombre(f.p.coefficient(0).puissanceEntiere(exposant)));
        return;
    }
    Nombre contenu;
    const Polynome primitive = f.p.partiePrimitive(contenu);
    resultat.multiplier(nombre(contenu.puissanceEntiere(exposant)));
    resultat.multiplier(ast_pow(primitive.versExpression(x), nombre(Nombre(exposant))));
}

} // namespace

ExprPtr reduireQuotientRationnel(const ExprPtr& e) {
    // Rejet rapide, sans allocation
    const Produit* produit_ = comme<Produit>(e);
    if (!produit_ || !produit_->contientVariable() || produit_->plusieursVariables()) return e;
    const std::vector<Facteur>& facteurs = produit_->getFacteurs();
    if (!reductionPossible(facteurs)) return e;

    // Séparation : polynômes exacts du numérateur et du dénominateur, autres facteurs
    std::vector<FacteurPolynomial> numerateurs, denominateurs;
    AccumulateurProduit resultat;
    resultat.multiplier(nombre(produit_->getCoefficient()));
    ExprPtr x;
    for (const Facteur& f : facteurs) {
        long long k;
        Polynome p;
        ExprPtr variable;
        if (facteurCandidat(f, k) && Polynome::depuisExpression(f.base, p, &variable) && p.estExact() &&
            variable && p.degre() > 0) {
            x = variable;
            (k > 0 ? numerateurs : denominateurs).push_back({std::move(p), k > 0 ? k : -k, f.base});
        } else {
            resultat.multiplier(ast_pow(f.base, f.exposant)); // gardé tel quel
        }
    }
    if (numerateurs.empty() || denominateurs.empty() || !reduire(numerateurs, denominateurs)) return e;

    for (const FacteurPolynomial& f : numerateurs) multiplierPar(resultat, f, 1, x);
    for (const FacteurPolynomial& f : denominateurs) multiplierPar(resultat, f, -1, x);
    return resultat.construire();
}

} // namespace symalgo
