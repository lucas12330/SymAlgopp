/**
 * @file Series.cpp
 * @brief Développements limités par arithmétique des séries de Taylor tronquées.
 *
 * Chaque noeud calcule directement les coefficients de Taylor de sa valeur (technique de
 * la différentiation automatique) : récurrences pour exp, ln, sin/cos, tan et u^p, en
 * O(n^2) par noeud au lieu de la croissance exponentielle des dérivées successives.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace symalgo {

namespace {

using Serie = std::vector<double>;

// Produit de Cauchy tronqué
Serie produitSeries(const Serie& u, const Serie& v) {
    Serie w(u.size(), 0.0);
    for (std::size_t k = 0; k < w.size(); ++k)
        for (std::size_t j = 0; j <= k; ++j) w[k] += u[j] * v[k - j];
    return w;
}

bool serieTaylor(const ASTNode& e, double a, int n, Serie& w);

/*
 * Nom : seriePuissance
 * Description : Série de u^p pour p constant : récurrence k u_0 w_k = sum (p j - (k - j)) u_j w_{k-j}
 *               si u(a) != 0, puissance entière positive sinon.
 */
bool seriePuissance(const Serie& u, double p, Serie& w) {
    const std::size_t taille = u.size();
    w.assign(taille, 0.0);
    if (u[0] != 0.0) {
        w[0] = std::pow(u[0], p);
        if (!std::isfinite(w[0])) return false;
        for (std::size_t k = 1; k < taille; ++k) {
            double somme = 0.0;
            for (std::size_t j = 1; j <= k; ++j) somme += (p * double(j) - double(k - j)) * u[j] * w[k - j];
            w[k] = somme / (double(k) * u[0]);
        }
        return true;
    }
    // u(a) = 0 : seule une puissance entière positive est développable
    if (p < 0.0 || std::floor(p) != p) return false;
    if (p > double(taille)) return true; // u^p = O((x-a)^p), nul jusqu'à l'ordre n
    w[0] = 1.0;
    for (int i = 0; i < static_cast<int>(p); ++i) w = produitSeries(w, u);
    return true;
}

// Série d'un facteur base^exposant
bool serieFacteur(const ExprPtr& base, const ExprPtr& exposant, double a, int n, Serie& w) {
    if (exposant->contientVariable()) {
        // u^v = exp(v ln u)
        return serieTaylor(*ast_exp(exposant * ast_ln(base)), a, n, w);
    }
    const Constante* c = comme<Constante>(exposant);
    if (!c) return false; // exposant symbolique (paramètre) ou constant non numérique
    Serie u;
    if (!serieTaylor(*base, a, n, u)) return false;
    return seriePuissance(u, c->getValeurConstante(), w);
}

/*
 * Nom : serieTaylor
 * Description : Coefficients de Taylor de e en a jusqu'à l'ordre n. Renvoie faux si un
 *               noeud n'est pas pris en charge (paramètre, noeud non évalué, singularité).
 */
bool serieTaylor(const ASTNode& e, double a, int n, Serie& w) {
    const std::size_t taille = static_cast<std::size_t>(n) + 1;
    w.assign(taille, 0.0);
    switch (e.type()) {
        case TypeNoeud::Constante: w[0] = e.getValeurConstante(); return true;
        case TypeNoeud::Pi: w[0] = 3.14159265358979323846; return true;
        case TypeNoeud::Variable:
            w[0] = a;
            if (n >= 1) w[1] = 1.0;
            return true;
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(e);
            w[0] = s.getConstante().versDouble();
            Serie t;
            for (const Terme& terme : s.getTermes()) {
                if (!serieTaylor(*terme.expression, a, n, t)) return false;
                const double c = terme.coefficient.versDouble();
                for (std::size_t k = 0; k < taille; ++k) w[k] += c * t[k];
            }
            return true;
        }
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(e);
            w[0] = p.getCoefficient().versDouble();
            Serie f;
            for (const Facteur& facteur : p.getFacteurs()) {
                if (!serieFacteur(facteur.base, facteur.exposant, a, n, f)) return false;
                w = produitSeries(w, f);
            }
            return true;
        }
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return serieFacteur(p.getBase(), p.getExposant(), a, n, w);
        }
        default: break;
    }
    const FonctionUnaire* f = comme<FonctionUnaire>(&e);
    if (!f) return false; // Parametre, noeuds non évalués
    Serie u;
    if (!serieTaylor(*f->m_argument, a, n, u)) return false;
    switch (e.type()) {
        case TypeNoeud::Exponentielle:
            // w' = u' w : k w_k = sum_{j=1..k} j u_j w_{k-j}
            w[0] = std::exp(u[0]);
            for (std::size_t k = 1; k < taille; ++k) {
                double somme = 0.0;
                for (std::size_t j = 1; j <= k; ++j) somme += double(j) * u[j] * w[k - j];
                w[k] = somme / double(k);
            }
            return true;
        case TypeNoeud::ArcSinus:
        case TypeNoeud::ArcCosinus:
        case TypeNoeud::ArcTangente: {
            // w' = u' s avec s = (1 - u^2)^(-1/2) (asin, acos au signe près) ou (1 + u^2)^(-1)
            const bool tangente = e.type() == TypeNoeud::ArcTangente;
            Serie v = produitSeries(u, u), s;
            for (double& vk : v) vk = tangente ? vk : -vk;
            v[0] += 1.0;
            if (!seriePuissance(v, tangente ? -1.0 : -0.5, s)) return false; // |u(a)| = 1 : non dérivable
            const double signe = e.type() == TypeNoeud::ArcCosinus ? -1.0 : 1.0;
            w[0] = tangente ? std::atan(u[0]) : (signe > 0 ? std::asin(u[0]) : std::acos(u[0]));
            if (!std::isfinite(w[0])) return false;
            for (std::size_t k = 1; k < taille; ++k) {
                double somme = 0.0;
                for (std::size_t j = 1; j <= k; ++j) somme += double(j) * u[j] * s[k - j];
                w[k] = signe * somme / double(k);
            }
            return true;
        }
        case TypeNoeud::Logarithme:
            // w' u = u' : k u_0 w_k = k u_k - sum_{j=1..k-1} j w_j u_{k-j}
            if (u[0] <= 0.0) return false;
            w[0] = std::log(u[0]);
            for (std::size_t k = 1; k < taille; ++k) {
                double somme = double(k) * u[k];
                for (std::size_t j = 1; j < k; ++j) somme -= double(j) * w[j] * u[k - j];
                w[k] = somme / (double(k) * u[0]);
            }
            return true;
        default: break;
    }
    // sin et cos se calculent ensemble : s' = u' c, c' = -u' s
    Serie sinus(taille, 0.0), cosinus(taille, 0.0);
    sinus[0] = std::sin(u[0]);
    cosinus[0] = std::cos(u[0]);
    for (std::size_t k = 1; k < taille; ++k) {
        double ss = 0.0, sc = 0.0;
        for (std::size_t j = 1; j <= k; ++j) {
            ss += double(j) * u[j] * cosinus[k - j];
            sc -= double(j) * u[j] * sinus[k - j];
        }
        sinus[k] = ss / double(k);
        cosinus[k] = sc / double(k);
    }
    switch (e.type()) {
        case TypeNoeud::Sinus: w = sinus; return true;
        case TypeNoeud::Cosinus: w = cosinus; return true;
        case TypeNoeud::Tangente:
            if (std::abs(cosinus[0]) < 1e-15) return false; // pôle de tan
            for (std::size_t k = 0; k < taille; ++k) {
                double somme = sinus[k];
                for (std::size_t j = 1; j <= k; ++j) somme -= cosinus[j] * w[k - j];
                w[k] = somme / cosinus[0];
            }
            return true;
        default: return false;
    }
}

// Coefficients de Taylor par dérivations symboliques successives (repli)
Serie serieParDerivation(const ASTNode& e, double a, int n) {
    Serie c(static_cast<std::size_t>(n) + 1, 0.0);
    ExprPtr deriv = e.clone();
    double factorielle = 1.0;
    for (int k = 0; k <= n; ++k) {
        if (k > 0) {
            deriv = deriv->derivee();
            factorielle *= k;
        }
        c[k] = deriv->eval(a) / factorielle;
    }
    return c;
}

} // namespace

ExprPtr ASTNode::DL(double a, int ordre) const {
    if (ordre < 0) throw std::invalid_argument("DL : ordre negatif");
    Serie c;
    if (!serieTaylor(*this, a, ordre, c)) c = serieParDerivation(*this, a, ordre);

    double echelle = 0.0;
    for (int k = 0; k <= ordre; ++k) {
        if (!std::isfinite(c[k])) {
            throw std::domain_error("DL : la fonction n'est pas developpable a l'ordre " + std::to_string(k) +
                                    " au point demande");
        }
        echelle = std::max(echelle, std::abs(c[k]));
    }
    // Les coefficients négligeables devant les autres (résidus d'arrondi, ex. cos(pi/2))
    // sont omis ; le seuil est relatif car 1/k! devient vite très petit
    const ExprPtr ecart = var("x") - a;
    AccumulateurSomme resultat;
    for (int k = 0; k <= ordre; ++k) {
        if (c[k] == 0.0 || std::abs(c[k]) <= 1e-15 * echelle) continue;
        resultat.ajouter(ast_pow(ecart, static_cast<double>(k)), Nombre::depuisDouble(c[k]));
    }
    return resultat.construire();
}

} // namespace symalgo
