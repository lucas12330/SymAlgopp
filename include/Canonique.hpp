/**
 * @file Canonique.hpp
 * @brief Outils internes de construction canonique (non destinés aux utilisateurs).
 */

#pragma once

#include <climits>
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
