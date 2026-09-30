/**
 * @file Polynome.cpp
 * @brief Développement des expressions et polynômes à coefficients exacts.
 */

#include "Polynome.hpp"
#include "Canonique.hpp"

#include <unordered_map>
#include <vector>

namespace symalgo {

namespace {

// Monôme d'une somme développée : coefficient * expression (un() pour la constante)
struct Monome {
    Nombre coefficient;
    ExprPtr expression;
};

std::vector<Monome> monomes(const ExprPtr& e) {
    switch (e->type()) {
        case TypeNoeud::Constante: return {{static_cast<const Constante&>(*e).getNombre(), un()}};
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(*e);
            std::vector<Monome> m;
            if (!s.getConstante().estZero()) m.push_back({s.getConstante(), un()});
            for (const Terme& t : s.getTermes()) m.push_back({t.coefficient, t.expression});
            return m;
        }
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(*e);
            return {{p.getCoefficient(), produitSansCoefficient(p)}};
        }
        default: return {{Nombre(1), e}};
    }
}

// (sum a_i) * (sum b_j) = sum a_i b_j
ExprPtr distribuer(const ExprPtr& a, const ExprPtr& b) {
    if (a->type() != TypeNoeud::Somme && b->type() != TypeNoeud::Somme) return a * b;
    AccumulateurSomme s;
    for (const Monome& ma : monomes(a)) {
        for (const Monome& mb : monomes(b)) s.ajouter(ma.expression * mb.expression, ma.coefficient * mb.coefficient);
    }
    return s.construire();
}

// s^n développé, par exponentiation rapide (s déjà développée, n >= 1)
ExprPtr puissanceDeveloppee(ExprPtr s, long long n) {
    ExprPtr resultat = un();
    while (n > 0) {
        if (n & 1) resultat = distribuer(resultat, s);
        n >>= 1;
        if (n > 0) s = distribuer(s, s);
    }
    return resultat;
}

class Developpeur {
public:
    ExprPtr developper(const ExprPtr& e) {
        if (!e->contientVariable() && e->type() == TypeNoeud::Constante) return e;
        const auto it = m_memo.find(e.get());
        if (it != m_memo.end()) return it->second;
        ExprPtr r = calculer(e);
        m_memo.emplace(e.get(), r);
        return r;
    }

private:
    ExprPtr facteur(const ExprPtr& base, const ExprPtr& exposant) {
        const ExprPtr b = developper(base);
        long long n;
        const Constante* c = comme<Constante>(exposant);
        if (c && c->getNombre().versEntier(n) && b->type() == TypeNoeud::Somme) {
            if (n > 0) return puissanceDeveloppee(b, n);
            if (n < 0) return ast_pow(puissanceDeveloppee(b, -n), nombre(Nombre(-1))); // dénominateur développé
        }
        return ast_pow(b, developper(exposant));
    }

    ExprPtr calculer(const ExprPtr& e) {
        switch (e->type()) {
            case TypeNoeud::Somme: {
                const Somme& s = static_cast<const Somme&>(*e);
                AccumulateurSomme acc;
                acc.ajouter(nombre(s.getConstante()), Nombre(1));
                for (const Terme& t : s.getTermes()) acc.ajouter(developper(t.expression), t.coefficient);
                return acc.construire();
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(*e);
                // Numérateurs distribués entre eux ; les dénominateurs restent un facteur
                ExprPtr numerateur = nombre(p.getCoefficient());
                AccumulateurProduit denominateur;
                for (const Facteur& f : p.getFacteurs()) {
                    const ExprPtr v = facteur(f.base, f.exposant);
                    const Constante* c = comme<Constante>(f.exposant);
                    if (c && c->getNombre().signe() < 0) {
                        denominateur.multiplier(v);
                    } else {
                        numerateur = distribuer(numerateur, v);
                    }
                }
                return distribuer(numerateur, denominateur.construire());
            }
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*e);
                return facteur(p.getBase(), p.getExposant());
            }
            case TypeNoeud::Sinus: return ast_sin(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::Cosinus: return ast_cos(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::Tangente: return ast_tan(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::Exponentielle: return ast_exp(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::Logarithme: return ast_ln(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            default: return e; // constantes, variable, paramètres, noeuds non évalués
        }
    }

    std::unordered_map<const ASTNode*, ExprPtr> m_memo;
};

} // namespace

ExprPtr developper(const ExprPtr& e) { return Developpeur().developper(e); }

} // namespace symalgo
