/**
 * @file Canonique.hpp
 * @brief Outils internes de construction canonique (non destinés aux utilisateurs).
 */

#pragma once

#include <climits>
#include <cmath>
#include <unordered_map>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

// Marqueur « exposant non entier » pour l'évaluation rapide des puissances
constexpr long long PAS_ENTIER = LLONG_MIN;

// Nombre d'éléments au-delà duquel les accumulateurs indexent leurs éléments par adresse
constexpr std::size_t SEUIL_INDEX = 16;

/*
 * Nom : AccumulateurSomme
 * Description : Construit une Somme canonique : aplatit les sommes, regroupe les termes
 *               identiques (par adresse, grâce au hash-consing) en additionnant leurs
 *               coefficients, supprime les termes nuls et trie le résultat.
 * Utilisation : AccumulateurSomme s; s.ajouter(a, Nombre(1)); s.ajouter(b, Nombre(-1));
 *               ExprPtr difference = s.construire();
 */
class AccumulateurSomme {
public:
    void ajouter(const ExprPtr& e, const Nombre& facteur); // ajoute facteur * e
    ExprPtr construire();

private:
    void ajouterTerme(const ExprPtr& terme, const Nombre& coefficient);

    Nombre m_constante;
    std::vector<Terme> m_termes;
    // Index des termes par adresse, construit seulement au-delà de SEUIL_INDEX termes
    // (en dessous, une recherche linéaire évite toute allocation)
    std::unordered_map<const ASTNode*, std::size_t> m_indices;
};

/*
 * Nom : AccumulateurProduit
 * Description : Construit un Produit canonique : aplatit les produits, additionne les
 *               exposants d'une même base, replie les constantes dans le coefficient et
 *               distribue un coefficient numérique sur une somme seule.
 */
class AccumulateurProduit {
public:
    void multiplier(const ExprPtr& e);
    ExprPtr construire();

private:
    void ajouterFacteur(const ExprPtr& base, const ExprPtr& exposant);

    Nombre m_coefficient{1};
    std::vector<Facteur> m_facteurs;
    std::unordered_map<const ASTNode*, std::size_t> m_indices; // idem AccumulateurSomme
};

/*
 * Nom : produitSansCoefficient
 * Description : Le produit de coefficient 1 formé des facteurs de p.
 */
ExprPtr produitSansCoefficient(const Produit& p);

/*
 * Nom : un
 * Description : La constante exacte 1 (partagée).
 */
const ExprPtr& un();

/*
 * Nom : signeReel
 * Description : Signe de x (-1, 0 ou 1 ; 0 aussi pour NaN), sans branchement.
 */
inline double signeReel(double x) { return static_cast<double>((x > 0.0) - (x < 0.0)); }

/*
 * Nom : coshReel, sinhReel, tanhReel
 * Description : Fonctions hyperboliques réelles construites sur une seule exponentielle
 *               (exp ou expm1), 7 à 9 fois plus rapides que std::cosh / std::sinh / std::tanh
 *               de la glibc pour une erreur d'au plus 2 ulp (mesurée dans test_evaluateur).
 *               Fonctions paires/impaires : on travaille sur |x| puis on rétablit le signe.
 *                 cosh(a) = (e + 1/e) / 2, e = exp(a)
 *                 sinh(a) = (t + t/(t+1)) / 2 avec t = expm1(a) si a < 1 (pas d'annulation
 *                           près de 0), sinon (e - 1/e) / 2
 *                 tanh(a) = t / (t + 2) avec t = expm1(2a) si a < 0,55, sinon 1 - 2/(exp(2a)+1)
 *               Au-delà de 709 on évite le débordement de exp(a) par (exp(a/2)/2) * exp(a/2) ;
 *               l'infini et NaN sont propagés comme std::.
 */
inline double coshReel(double x) {
    const double a = std::fabs(x);
    if (a <= 709.0) {
        const double e = std::exp(a);
        return 0.5 * (e + 1.0 / e);
    }
    if (a <= 710.4758600739439) {
        const double w = std::exp(0.5 * a);
        return (0.5 * w) * w;
    }
    return a != a ? a : HUGE_VAL;
}

inline double sinhReel(double x) {
    const double a = std::fabs(x);
    if (a < 3.725290298461914e-9) return x; // |x| < 2^-28 : sinh(x) = x (écart relatif < 2^-57)
    double r;
    if (a < 1.0) {
        const double t = std::expm1(a);
        r = 0.5 * (t + t / (t + 1.0));
    } else if (a <= 709.0) {
        const double e = std::exp(a);
        r = 0.5 * (e - 1.0 / e);
    } else if (a <= 710.4758600739439) {
        const double w = std::exp(0.5 * a);
        r = (0.5 * w) * w;
    } else {
        r = a != a ? a : HUGE_VAL;
    }
    return std::copysign(r, x);
}

inline double tanhReel(double x) {
    const double a = std::fabs(x);
    if (a < 3.725290298461914e-9) return x; // |x| < 2^-28 : tanh(x) = x (écart relatif < 2^-57)
    if (!(a < 22.0)) return a != a ? x : std::copysign(1.0, x); // tanh(22) arrondit à 1
    double r;
    if (a < 0.55) {
        const double t = std::expm1(2.0 * a);
        r = t / (t + 2.0);
    } else {
        r = 1.0 - 2.0 / (std::exp(2.0 * a) + 1.0);
    }
    return std::copysign(r, x);
}

/*
 * Nom : puissanceEntiereReelle
 * Description : b^k par exponentiation rapide (bien plus rapide que std::pow pour les
 *               petits exposants entiers, cas courant des polynômes).
 */
inline double puissanceEntiereReelle(double b, long long k) {
    switch (k) {
        case 1: return b;
        case 2: return b * b;
        case 3: return b * b * b;
        case -1: return 1.0 / b;
        case -2: return 1.0 / (b * b);
        default: break;
    }
    const bool negatif = k < 0;
    unsigned long long n = negatif ? 0ULL - static_cast<unsigned long long>(k) : static_cast<unsigned long long>(k);
    double r = 1.0;
    while (n) {
        if (n & 1ULL) r *= b;
        b *= b;
        n >>= 1;
    }
    return negatif ? 1.0 / r : r;
}

} // namespace symalgo
