/**
 * @file Solveur.cpp
 * @brief Résolution d'équations : isolement de l'inconnue, polynômes certifiés, familles
 *        trigonométriques, recherche numérique sur un intervalle.
 */

#include "Solveur.hpp"
#include "Canonique.hpp"
#include "Evaluateur.hpp"
#include "Journal.hpp"
#include "Polynome.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace symalgo {

namespace {

constexpr double NAN_SOLUTION = std::numeric_limits<double>::quiet_NaN();
constexpr int PROFONDEUR_MAX = 40;

// Étapes du journal (journal::noter ne fait qu'un test de pointeur tant qu'aucun journal n'est
// actif : ni chaîne ni expression n'est construite hors des lambdas).
void etape(const char* regle, const char* description, const ExprPtr& gauche, const ExprPtr& droite) {
    journal::noter([&](Journal& j) { j.ajouterEquation(regle, description, gauche, droite); });
}

void note(const char* regle, const char* description) {
    journal::noter([&](Journal& j) { j.ajouter(regle, description); });
}

bool estZeroExact(const ExprPtr& e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().estZero();
}

// Valeur d'une expression constante (sans variable ni paramètre)
bool valeurNumerique(const ExprPtr& e, double& v) {
    if (e->contientVariable()) return false;
    try {
        v = e->eval(0.0);
    } catch (const std::logic_error&) {
        return false;
    }
    return !std::isnan(v);
}

// Parcourt les sous-expressions distinctes de e
template <class Visiteur>
void parcourir(const ExprPtr& e, Visiteur&& visiter, std::unordered_set<const ASTNode*>& vus) {
    if (!vus.insert(e.get()).second) return;
    visiter(e);
    switch (e->type()) {
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(*e).getTermes()) parcourir(t.expression, visiter, vus);
            break;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(*e).getFacteurs()) {
                parcourir(f.base, visiter, vus);
                parcourir(f.exposant, visiter, vus);
            }
            break;
        case TypeNoeud::Puissance:
            parcourir(static_cast<const Puissance&>(*e).getBase(), visiter, vus);
            parcourir(static_cast<const Puissance&>(*e).getExposant(), visiter, vus);
            break;
        default:
            if (const FonctionUnaire* f = comme<FonctionUnaire>(e)) parcourir(f->m_argument, visiter, vus);
    }
}

ExprPtr trouverVariable(const ExprPtr& e) {
    ExprPtr variable;
    std::unordered_set<const ASTNode*> vus;
    parcourir(e, [&](const ExprPtr& n) { if (!variable && n->type() == TypeNoeud::Variable) variable = n; }, vus);
    return variable;
}

// Vrai si l'expression ne contient que des constantes exactes
bool estExacte(const ExprPtr& e) {
    bool exacte = true;
    std::unordered_set<const ASTNode*> vus;
    parcourir(e, [&](const ExprPtr& n) {
        if (const Constante* c = comme<Constante>(n); c && !c->getNombre().estExact()) exacte = false;
        if (const Somme* s = comme<Somme>(n); s && !s->getConstante().estExact()) exacte = false;
        if (const Produit* p = comme<Produit>(n); p && !p->getCoefficient().estExact()) exacte = false;
    }, vus);
    return exacte;
}

/*
 * Nom : Resolveur
 * Description : Isole l'inconnue m_x dans L = R en inversant les opérations une à une.
 *               Les entiers libres des familles trigonométriques en cours sont empilés dans
 *               m_entiers ; toute forme non résolue rend le résultat incomplet.
 */
class Resolveur {
public:
    explicit Resolveur(ExprPtr x) : m_x(std::move(x)) {}

    Solutions resoudre(const ExprPtr& gauche, const ExprPtr& droite) {
        isoler(gauche, droite, 0);
        return m_resultat;
    }

private:
    void incomplet(const char* raison) {
        m_resultat.complet = false;
        journal::noter([&](Journal& j) {
            j.ajouter("resolution.incomplet", std::string("Forme non résolue, le résultat peut être incomplet : ") + raison);
        });
    }

    void ajouter(const ExprPtr& valeur, int multiplicite, bool exacte) {
        Solution s;
        s.valeur = valeur;
        s.multiplicite = multiplicite;
        s.exacte = exacte && estExacte(valeur);
        s.entiers = m_entiers;
        if (s.entiers.empty()) {
            double v;
            if (!valeurNumerique(valeur, v)) {
                incomplet("la valeur de la solution n'est pas un nombre");
                return;
            }
            s.approximation = v;
        } else {
            s.approximation = NAN_SOLUTION;
        }
        m_resultat.liste.push_back(s);
    }

    std::string nouvelEntier() {
        ++m_compteurEntiers;
        return m_compteurEntiers == 1 ? "k" : "k" + std::to_string(m_compteurEntiers);
    }

    // Polynôme en x : racines certifiées (suites de Sturm), liste complète
    bool parPolynome(const ExprPtr& D, bool journaliser) {
        Polynome p;
        ExprPtr v;
        const ExprPtr developpee = developper(D);
        if (!Polynome::depuisExpression(developpee, p, &v) || !v || v.get() != m_x.get()) return false;
        if (p.estNul()) {
            if (journaliser) note("resolution.identite", "Le polynôme est nul : l'équation est vraie pour tout x réel");
            m_resultat.toutReel = true;
            return true;
        }
        if (journaliser) journal::noter([&](Journal& j) {
            const int d = p.degre();
            j.ajouterEquation("resolution.polynome",
                              d == 1 ? "L'équation est du premier degré"
                              : d == 2 ? "L'équation est du second degré"
                                       : "L'équation est polynomiale de degré " + std::to_string(d),
                              developpee, nombre(Nombre(0)));
            if (d == 2) {
                const Nombre delta = p.coefficient(1) * p.coefficient(1) - Nombre(4) * p.coefficient(2) * p.coefficient(0);
                j.ajouterEquation("resolution.discriminant",
                                  "Discriminant Δ = b² − 4ac, avec a = " + p.coefficient(2).texte() +
                                      ", b = " + p.coefficient(1).texte() + ", c = " + p.coefficient(0).texte(),
                                  param("Delta"), nombre(delta));
                j.ajouter("resolution.discriminant_signe", delta.signe() > 0   ? "Δ > 0 : deux racines réelles distinctes"
                                                           : delta.signe() == 0 ? "Δ = 0 : une racine double"
                                                                                : "Δ < 0 : aucune racine réelle");
            } else if (d >= 3) {
                j.ajouter("resolution.sturm", "Racines réelles isolées et comptées par la suite de Sturm");
            }
        });
        for (const RacineReelle& r : p.racinesReelles()) {
            const ExprPtr valeur = r.exacte ? r.exacte : nombre(Nombre::reel(r.valeur));
            if (journaliser) journal::noter([&](Journal& j) {
                std::string texte = r.exacte ? "Racine exacte" : "Racine irrationnelle (valeur approchée certifiée)";
                if (r.multiplicite == 2) texte += ", double";
                if (r.multiplicite > 2) texte += ", de multiplicité " + std::to_string(r.multiplicite);
                j.ajouterEquation("resolution.racine", texte, m_x, valeur, !r.exacte);
            });
            ajouter(valeur, r.multiplicite, r.exacte != nullptr);
        }
        return true;
    }

    /*
     * Nom : parSubstitution
     * Description : Si x n'apparaît qu'à travers un seul « atome » u (sin(x), exp(2x), ln(x)...)
     *               et que l'équation est polynomiale en u, résout le polynôme en t = u puis
     *               u = t pour chaque racine (ex. sin(x)^2 - sin(x) = 0).
     */
    bool parSubstitution(const ExprPtr& D, int profondeur) {
        const ExprPtr d = developper(regrouperExponentielles(developper(D)));
        std::vector<ExprPtr> atomes;
        bool variableNue = false;
        std::unordered_set<const ASTNode*> vus;
        collecterAtomes(d, atomes, variableNue, vus);
        if (variableNue || atomes.size() != 1) return false;
        const ExprPtr u = atomes[0];
        const ExprPtr t = var("_t");
        const ExprPtr E = substituer(d, u, t);
        if (contient(E, m_x)) return false;
        Polynome p;
        ExprPtr v;
        // Degré 1 : c'est déjà u = constante, que l'isolement direct traite (sinon boucle infinie)
        if (!Polynome::depuisExpression(developper(E), p, &v) || !v || p.degre() < 2) return false;
        journal::noter([&](Journal& j) {
            if (m_exponentiellesRegroupees) {
                j.ajouter("resolution.regroupement_exponentielles",
                          "Les exponentielles ont des arguments proportionnels : on les écrit comme puissances d'une même exponentielle");
            }
            j.ajouterEquation("resolution.changement_variable",
                              "Changement de variable : x n'intervient que par une expression, que l'on note t", param("t"), u);
            j.ajouterEquation("resolution.polynome_en_t", "L'équation devient polynomiale en t", E, nombre(Nombre(0)));
        });
        Solutions racines;
        {
            journal::SousNiveau sous;
            racines = Resolveur(t).resoudre(E, nombre(Nombre(0)));
        }
        if (!racines.complet) incomplet("équation en t incomplète");
        if (racines.liste.empty()) note("resolution.aucune_racine", "Aucune racine réelle en t : l'équation n'a pas de solution");
        for (const Solution& r : racines.liste) {
            etape("resolution.retour_variable", "Retour à l'inconnue x pour cette racine en t", u, r.valeur);
            journal::SousNiveau sous;
            isoler(u, r.valeur, profondeur + 1);
        }
        return true;
    }

    /*
     * Nom : regrouperExponentielles
     * Description : Réécrit les exponentielles exp(c_i * g) d'arguments proportionnels comme
     *               puissances entières d'une même exponentielle : exp(2x) = exp(x)^2, ce qui
     *               rend exp(2x) - 3 exp(x) + 2 polynomiale en exp(x).
     */
    ExprPtr regrouperExponentielles(const ExprPtr& d) {
        m_exponentiellesRegroupees = false;
        std::vector<ExprPtr> exponentielles;
        std::unordered_set<const ASTNode*> vus;
        parcourir(d, [&](const ExprPtr& n) { if (n->type() == TypeNoeud::Exponentielle && contient(n, m_x)) exponentielles.push_back(n); }, vus);
        if (exponentielles.size() < 2) return d;
        // Argument = coefficient * noyau : les noyaux doivent coïncider
        ExprPtr noyau;
        std::vector<Nombre> coefficients;
        for (const ExprPtr& e : exponentielles) {
            const ExprPtr& arg = static_cast<const FonctionUnaire&>(*e).m_argument;
            Nombre c(1);
            ExprPtr g = arg;
            if (const Produit* p = comme<Produit>(arg)) {
                c = p->getCoefficient();
                g = produitSansCoefficient(*p);
            }
            if (!c.estExact() || (noyau && noyau.get() != g.get())) return d;
            noyau = g;
            coefficients.push_back(c);
        }
        // Pas commun : pgcd des coefficients rationnels
        Nombre numerateurs(0), denominateurs(1);
        for (const Nombre& c : coefficients) {
            numerateurs = Nombre::pgcd(numerateurs, c.numerateurNombre());
            const Nombre dc = c.denominateurNombre();
            denominateurs = denominateurs * dc / Nombre::pgcd(denominateurs, dc);
        }
        const Nombre pas = numerateurs / denominateurs;
        const ExprPtr base = ast_exp(nombre(pas) * noyau);
        ExprPtr resultat = d;
        for (std::size_t i = 0; i < exponentielles.size(); ++i) {
            resultat = substituer(resultat, exponentielles[i], ast_pow(base, nombre(coefficients[i] / pas)));
        }
        m_exponentiellesRegroupees = true;
        return resultat;
    }

    void collecterAtomes(const ExprPtr& e, std::vector<ExprPtr>& atomes, bool& variableNue,
                         std::unordered_set<const ASTNode*>& vus) {
        if (!contient(e, m_x) || !vus.insert(e.get()).second) return;
        auto ajouterAtome = [&](const ExprPtr& a) {
            if (std::none_of(atomes.begin(), atomes.end(), [&](const ExprPtr& b) { return b.get() == a.get(); })) atomes.push_back(a);
        };
        switch (e->type()) {
            case TypeNoeud::Variable: variableNue = true; break;
            case TypeNoeud::Somme:
                for (const Terme& t : static_cast<const Somme&>(*e).getTermes()) collecterAtomes(t.expression, atomes, variableNue, vus);
                break;
            case TypeNoeud::Produit:
                for (const Facteur& f : static_cast<const Produit&>(*e).getFacteurs()) {
                    collecterAtomes(ast_pow(f.base, f.exposant), atomes, variableNue, vus);
                }
                break;
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*e);
                long long k;
                const Constante* c = comme<Constante>(p.getExposant());
                if (c && c->getNombre().versEntier(k) && k > 0) {
                    collecterAtomes(p.getBase(), atomes, variableNue, vus); // puissance entière : polynomiale
                } else {
                    ajouterAtome(e);
                }
                break;
            }
            default: ajouterAtome(e); break;
        }
    }

    void isolerPuissance(const Puissance& p, const ExprPtr& R, int profondeur) {
        const ExprPtr& b = p.getBase();
        const ExprPtr& e = p.getExposant();
        if (!contient(e, m_x)) {
            const Constante* c = comme<Constante>(e);
            if (!c) {
                incomplet("exposant non numérique");
                return;
            }
            const Nombre& n = c->getNombre();
            if (estZeroExact(R)) {
                if (n.signe() > 0) {
                    etape("resolution.puissance_nulle", "Une puissance positive est nulle si sa base l'est", b, R);
                    isoler(b, R, profondeur + 1); // b^n = 0 => b = 0
                }
                return;
            }
            double v;
            if (!valeurNumerique(R, v)) {
                incomplet("second membre non numérique");
                return;
            }
            const ExprPtr inverse = nombre(n.inverse());
            long long k;
            if (n.versEntier(k) && k % 2 == 0) {
                if (v < 0.0) { // puissance paire négative : pas de solution réelle
                    note("resolution.hors_domaine", "Une puissance paire est positive : pas de solution réelle");
                    return;
                }
                const ExprPtr r = ast_pow(R, inverse);
                etape("resolution.racine_paire", "Puissance paire : on prend la racine, avec le signe +", b, r);
                isoler(b, r, profondeur + 1);
                if (v > 0.0) {
                    const ExprPtr oppose = -r;
                    etape("resolution.racine_paire", "Puissance paire : et avec le signe −", b, oppose);
                    isoler(b, oppose, profondeur + 1);
                }
            } else if (n.versEntier(k)) { // impaire
                const ExprPtr r = v >= 0.0 ? ast_pow(R, inverse) : -ast_pow(-R, inverse);
                etape("resolution.racine_impaire", "Puissance impaire : une seule racine réelle", b, r);
                isoler(b, r, profondeur + 1);
            } else if (v >= 0.0) { // exposant non entier : b >= 0
                const ExprPtr r = ast_pow(R, inverse);
                etape("resolution.racine_fractionnaire", "Exposant non entier : la base est positive, on élève à la puissance inverse", b, r);
                isoler(b, r, profondeur + 1);
            } else {
                note("resolution.hors_domaine", "Une puissance d'exposant non entier est positive : pas de solution réelle");
            }
            return;
        }
        if (!contient(b, m_x)) {
            // b^u = R avec b constante positive : u = ln(R) / ln(b)
            double base, v;
            if (!valeurNumerique(b, base) || !valeurNumerique(R, v) || base <= 0.0 || base == 1.0) {
                incomplet("base ou second membre non exploitable");
                return;
            }
            if (v > 0.0) {
                const ExprPtr suite = ast_ln(R) / ast_ln(b);
                etape("resolution.logarithme", "On applique ln des deux côtés (second membre strictement positif : domaine vérifié)", e, suite);
                isoler(e, suite, profondeur + 1);
            } else {
                note("resolution.hors_domaine", "Une puissance d'une base positive est strictement positive : pas de solution réelle");
            }
            return;
        }
        incomplet("l'inconnue figure dans la base et dans l'exposant");
    }

    void isolerTrigonometrique(TypeNoeud type, const ExprPtr& u, const ExprPtr& R, int profondeur) {
        const std::string k = nouvelEntier();
        const ExprPtr K = param(k);
        if (type == TypeNoeud::Tangente) {
            m_entiers.push_back(k);
            const ExprPtr suite = ast_atan(R) + K * pi(); // période pi
            etape("resolution.famille_trigonometrique", "tan est π-périodique : u = arctan(R) + kπ, k entier relatif", u, suite);
            isoler(u, suite, profondeur + 1);
            m_entiers.pop_back();
            return;
        }
        double v;
        if (!valeurNumerique(R, v)) {
            incomplet("second membre trigonométrique non numérique");
            return;
        }
        if (std::abs(v) > 1.0) { // pas de solution réelle
            note("resolution.hors_domaine", "sin et cos sont à valeurs dans [−1, 1] : pas de solution réelle");
            return;
        }
        const ExprPtr deuxKPi = cst(2.0) * K * pi();
        m_entiers.push_back(k);
        const bool sinus = type == TypeNoeud::Sinus;
        const ExprPtr A = sinus ? ast_asin(R) : ast_acos(R);
        const ExprPtr premiere = A + deuxKPi;
        etape("resolution.famille_trigonometrique",
              sinus ? "Famille trigonométrique : sin(u) = R donne u = arcsin(R) + 2kπ, k entier relatif"
                    : "Famille trigonométrique : cos(u) = R donne u = arccos(R) + 2kπ, k entier relatif",
              u, premiere);
        {
            journal::SousNiveau sous;
            isoler(u, premiere, profondeur + 1);
        }
        if (std::abs(v) != 1.0) {
            const ExprPtr seconde = sinus ? pi() - A + deuxKPi : -A + deuxKPi;
            etape("resolution.famille_trigonometrique",
                  sinus ? "Seconde famille : u = π − arcsin(R) + 2kπ" : "Seconde famille : u = −arccos(R) + 2kπ", u, seconde);
            journal::SousNiveau sous;
            isoler(u, seconde, profondeur + 1);
        }
        m_entiers.pop_back();
    }

    /*
     * Nom : isolerAbsoluOuHyperbolique
     * Description : Inverse abs, sign, sinh, cosh et tanh : abs(u) = R => u = +-R (R >= 0),
     *               sinh(u) = R => u = ln(R + (R^2+1)^(1/2)), cosh(u) = R => u = +-ln(R +
     *               (R^2-1)^(1/2)) (R >= 1), tanh(u) = R => u = ln((1+R)/(1-R))/2 (|R| < 1).
     *               sign(u) = 0 donne u = 0 ; sign(u) = +-1 est une inéquation (u > 0 ou
     *               u < 0), non représentable : le résultat est alors incomplet.
     */
    void isolerAbsoluOuHyperbolique(TypeNoeud type, const ExprPtr& u, const ExprPtr& R, int profondeur) {
        double v;
        if (!valeurNumerique(R, v)) {
            if (type != TypeNoeud::SinusHyperbolique) {
                incomplet("second membre non numérique");
                return;
            }
            v = 0.0; // sinh est bijective : aucune condition sur R
        }
        const ExprPtr deux = cst(2.0);
        switch (type) {
            case TypeNoeud::ValeurAbsolue:
                if (v < 0.0) {
                    note("resolution.hors_domaine", "Une valeur absolue est positive : pas de solution réelle");
                    return;
                }
                etape("resolution.valeur_absolue", "|u| = R (R ≥ 0) équivaut à u = R ou u = −R : première possibilité", u, R);
                isoler(u, R, profondeur + 1);
                if (v > 0.0) {
                    const ExprPtr oppose = -R;
                    etape("resolution.valeur_absolue", "Seconde possibilité", u, oppose);
                    isoler(u, oppose, profondeur + 1);
                }
                return;
            case TypeNoeud::Signe:
                if (v == 0.0) {
                    const ExprPtr zero = nombre(Nombre(0));
                    etape("resolution.signe", "sign(u) = 0 équivaut à u = 0", u, zero);
                    isoler(u, zero, profondeur + 1);
                } else if (std::abs(v) == 1.0) {
                    incomplet("sign(u) = ±1 est une inéquation (u > 0 ou u < 0)");
                } else {
                    note("resolution.hors_domaine", "sign ne prend que les valeurs −1, 0 et 1 : pas de solution");
                }
                return;
            case TypeNoeud::SinusHyperbolique: {
                const ExprPtr suite = ast_ln(R + ast_pow(ast_pow(R, deux) + 1.0, frac(1, 2)));
                etape("resolution.hyperbolique", "sinh est bijective : u = ln(R + √(R² + 1))", u, suite);
                isoler(u, suite, profondeur + 1);
                return;
            }
            case TypeNoeud::CosinusHyperbolique: {
                if (v < 1.0) {
                    note("resolution.hors_domaine", "cosh est supérieur ou égal à 1 : pas de solution réelle");
                    return;
                }
                const ExprPtr A = ast_ln(R + ast_pow(ast_pow(R, deux) - 1.0, frac(1, 2)));
                etape("resolution.hyperbolique", "cosh(u) = R (R ≥ 1) : u = ln(R + √(R² − 1))", u, A);
                isoler(u, A, profondeur + 1);
                if (v > 1.0) {
                    const ExprPtr oppose = -A;
                    etape("resolution.hyperbolique", "cosh est paire : et u = −ln(R + √(R² − 1))", u, oppose);
                    isoler(u, oppose, profondeur + 1);
                }
                return;
            }
            default: // TangenteHyperbolique
                if (std::abs(v) < 1.0) {
                    const ExprPtr suite = frac(1, 2) * ast_ln((1.0 + R) / (1.0 - R));
                    etape("resolution.hyperbolique", "tanh(u) = R (|R| < 1) : u = ln((1 + R)/(1 − R))/2", u, suite);
                    isoler(u, suite, profondeur + 1);
                } else {
                    note("resolution.hors_domaine", "tanh est à valeurs dans ]−1, 1[ : pas de solution réelle");
                }
                return;
        }
    }

    void isoler(ExprPtr L, ExprPtr R, int profondeur) {
        if (profondeur > PROFONDEUR_MAX) {
            incomplet("profondeur maximale atteinte");
            return;
        }
        if (contient(R, m_x)) { // l'inconnue d'un seul côté
            L = L - R;
            R = nombre(Nombre(0));
        }
        const ExprPtr developpee = developper(L - R);
        if (!contient(L, m_x) || !contient(developpee, m_x)) { // plus d'inconnue (ex. identité)
            const ExprPtr& D = developpee;
            double v;
            if (estZeroExact(D)) {
                if (m_entiers.empty()) m_resultat.toutReel = true;
            } else if (!valeurNumerique(D, v)) {
                incomplet("égalité symbolique indécidable"); // égalité symbolique indécidable
            }
            return;
        }
        // L = x : l'équation est déjà résolue (le polynôme de degré 1 n'apporte rien à expliquer)
        const bool directe = L.get() == m_x.get();
        if (directe) etape("resolution.isolee", "L'inconnue est isolée", L, R);
        if (m_entiers.empty() && parPolynome(L - R, !directe)) return;
        if (m_entiers.empty() && parSubstitution(L - R, profondeur)) return;

        switch (L->type()) {
            case TypeNoeud::Variable: ajouter(R->simplifier(), 1, true); return;
            case TypeNoeud::Somme: {
                // Les termes sans x passent à droite : c*t + reste = R  =>  t = (R - reste)/c
                const Somme& s = static_cast<const Somme&>(*L);
                AccumulateurSomme reste;
                reste.ajouter(nombre(s.getConstante()), Nombre(1));
                const Terme* seul = nullptr;
                int avecX = 0;
                for (const Terme& t : s.getTermes()) {
                    if (contient(t.expression, m_x)) {
                        seul = &t;
                        ++avecX;
                    } else {
                        reste.ajouter(t.expression, t.coefficient);
                    }
                }
                if (avecX != 1) {
                    incomplet("plusieurs termes dépendent de l'inconnue");
                    return;
                }
                const ExprPtr suite = (R - reste.construire()) / nombre(seul->coefficient);
                etape("resolution.isolement_somme", "On isole le terme qui contient x : on déplace les autres termes et on divise par son coefficient",
                      seul->expression, suite);
                isoler(seul->expression, suite, profondeur + 1);
                return;
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(*L);
                std::vector<ExprPtr> avecX;
                AccumulateurProduit autres;
                autres.multiplier(nombre(p.getCoefficient()));
                for (const Facteur& f : p.getFacteurs()) {
                    const ExprPtr facteur = ast_pow(f.base, f.exposant);
                    if (contient(facteur, m_x)) {
                        avecX.push_back(facteur);
                    } else {
                        autres.multiplier(facteur);
                    }
                }
                if (avecX.size() == 1) {
                    const ExprPtr suite = R / autres.construire();
                    etape("resolution.division", "On divise par les facteurs qui ne dépendent pas de x", avecX[0], suite);
                    isoler(avecX[0], suite, profondeur + 1);
                } else if (estZeroExact(R)) {
                    // Produit nul : un des facteurs est nul (un facteur d'exposant négatif ne l'est jamais)
                    etape("resolution.produit_nul", "Un produit est nul si l'un de ses facteurs l'est", L, R);
                    for (const ExprPtr& f : avecX) {
                        const Puissance* q = comme<Puissance>(f);
                        const Constante* e = q ? comme<Constante>(q->getExposant()) : nullptr;
                        if (!(e && e->getNombre().signe() < 0)) {
                            etape("resolution.facteur_nul", "Facteur à annuler", f, R);
                            journal::SousNiveau sous;
                            isoler(f, R, profondeur + 1);
                        } else {
                            note("resolution.facteur_denominateur", "Un facteur au dénominateur n'est jamais nul (il doit seulement être défini)");
                        }
                    }
                } else {
                    incomplet("produit de plusieurs facteurs en x égal à une constante non nulle");
                }
                return;
            }
            case TypeNoeud::Puissance: isolerPuissance(static_cast<const Puissance&>(*L), R, profondeur); return;
            case TypeNoeud::Exponentielle: {
                double v;
                if (!valeurNumerique(R, v)) {
                    incomplet("second membre de exp non numérique");
                } else if (v > 0.0) {
                    const ExprPtr& u = static_cast<const FonctionUnaire&>(*L).m_argument;
                    const ExprPtr suite = ast_ln(R);
                    etape("resolution.logarithme", "On applique ln des deux côtés (exp(u) = R exige R > 0 : domaine vérifié)", u, suite);
                    isoler(u, suite, profondeur + 1);
                } else {
                    note("resolution.hors_domaine", "exp est strictement positive : pas de solution réelle");
                }
                return;
            }
            case TypeNoeud::Logarithme: {
                const ExprPtr& u = static_cast<const FonctionUnaire&>(*L).m_argument;
                const ExprPtr suite = ast_exp(R);
                etape("resolution.exponentielle", "On applique exp des deux côtés : ln(u) = R équivaut à u = exp(R) (u > 0 vérifié à la fin)", u, suite);
                isoler(u, suite, profondeur + 1);
                return;
            }
            case TypeNoeud::ValeurAbsolue:
            case TypeNoeud::Signe:
            case TypeNoeud::SinusHyperbolique:
            case TypeNoeud::CosinusHyperbolique:
            case TypeNoeud::TangenteHyperbolique:
                isolerAbsoluOuHyperbolique(L->type(), static_cast<const FonctionUnaire&>(*L).m_argument, R, profondeur);
                return;
            case TypeNoeud::Sinus:
            case TypeNoeud::Cosinus:
            case TypeNoeud::Tangente:
                isolerTrigonometrique(L->type(), static_cast<const FonctionUnaire&>(*L).m_argument, R, profondeur);
                return;
            case TypeNoeud::ArcSinus:
            case TypeNoeud::ArcCosinus:
            case TypeNoeud::ArcTangente: {
                // Réciproques : l'image est un intervalle, hors duquel il n'y a pas de solution
                const double demiPi = 1.57079632679489661923;
                double v;
                if (!valeurNumerique(R, v)) {
                    incomplet("second membre d'une fonction réciproque non numérique");
                    return;
                }
                const ExprPtr& u = static_cast<const FonctionUnaire&>(*L).m_argument;
                const bool sinus = L->type() == TypeNoeud::ArcSinus, cosinus = L->type() == TypeNoeud::ArcCosinus;
                if ((sinus && std::abs(v) <= demiPi) || (cosinus && v >= 0.0 && v <= 2.0 * demiPi) ||
                    (!sinus && !cosinus && std::abs(v) < demiPi)) {
                    const ExprPtr suite = sinus ? ast_sin(R) : cosinus ? ast_cos(R) : ast_tan(R);
                    etape("resolution.reciproque", "On applique la fonction directe (le second membre est dans l'image de la réciproque : domaine vérifié)", u, suite);
                    isoler(u, suite, profondeur + 1);
                } else {
                    note("resolution.hors_domaine", "Le second membre est hors de l'image de la fonction réciproque : pas de solution");
                }
                return;
            }
            default: incomplet("forme non prise en charge"); return;
        }
    }

    ExprPtr m_x;
    Solutions m_resultat;
    std::vector<std::string> m_entiers;
    int m_compteurEntiers = 0;
    bool m_exponentiellesRegroupees = false; // dernier appel de regrouperExponentielles : réécriture faite
};

/*
 * Nom : nettoyer
 * Description : Supprime les doublons, élimine les solutions hors du domaine de l'équation
 *               (valeur non finie, ex. x = 0 dans x*ln(x) = 0) et trie par valeur.
 */
void nettoyer(Solutions& s, const ExprPtr& f, const ExprPtr& x) {
    const ExprPtr derivee = x ? f->derivee() : nullptr;
    std::vector<Solution> gardees;
    for (Solution& sol : s.liste) {
        if (!sol.estFamille() && x) {
            double fx, dfx = 0.0;
            try {
                fx = f->eval(sol.approximation);
                dfx = derivee->eval(sol.approximation);
            } catch (const std::logic_error&) {
                fx = 0.0; // paramètre symbolique : pas de vérification possible
            }
            const double tolerance = 1e-9 * (1.0 + std::abs(sol.approximation)) * (1.0 + (std::isfinite(dfx) ? std::abs(dfx) : 0.0));
            if (!std::isfinite(fx) || std::abs(fx) > tolerance) {
                journal::noter([&](Journal& j) {
                    j.ajouterEquation("resolution.solution_rejetee",
                                      std::isfinite(fx) ? "Solution rejetée : elle n'annule pas l'équation de départ (solution étrangère)"
                                                        : "Solution rejetée : hors du domaine de définition de l'équation (valeur non finie)",
                                      x, sol.valeur, !sol.exacte && !sol.estFamille());
                });
                continue;
            }
        }
        const bool doublon = std::any_of(gardees.begin(), gardees.end(), [&](const Solution& g) {
            if (g.valeur.get() == sol.valeur.get()) return true;
            return !g.estFamille() && !sol.estFamille() &&
                   std::abs(g.approximation - sol.approximation) <= 1e-12 * std::max(1.0, std::abs(g.approximation));
        });
        if (!doublon) {
            gardees.push_back(sol);
        } else {
            journal::noter([&](Journal& j) {
                j.ajouterEquation("resolution.doublon", "Solution déjà obtenue, ignorée", x, sol.valeur, !sol.exacte && !sol.estFamille());
            });
        }
    }
    std::stable_sort(gardees.begin(), gardees.end(), [](const Solution& a, const Solution& b) {
        if (a.estFamille() != b.estFamille()) return !a.estFamille();
        return !a.estFamille() && a.approximation < b.approximation;
    });
    s.liste = std::move(gardees);
}

// Méthode de Brent sur [a, b] avec f(a) et f(b) de signes opposés
template <class Fonction>
double brent(Fonction f, double a, double b, double fa, double fb) {
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a, fc = fa, d = b - a, e = d;
    for (int iteration = 0; iteration < 200; ++iteration) {
        if (fb == 0.0) return b;
        if ((fb > 0.0) == (fc > 0.0)) {
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        const double tolerance = 2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) + 1e-300;
        const double milieu = 0.5 * (c - b);
        if (std::abs(milieu) <= tolerance) return b;
        if (std::abs(e) >= tolerance && std::abs(fa) > std::abs(fb)) {
            // Interpolation (sécante ou quadratique inverse)
            double p, q;
            const double s = fb / fa;
            if (a == c) {
                p = 2.0 * milieu * s;
                q = 1.0 - s;
            } else {
                const double r = fb / fc, t = fa / fc;
                p = s * (2.0 * milieu * t * (t - r) - (b - a) * (r - 1.0));
                q = (t - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::abs(p);
            if (2.0 * p < std::min(3.0 * milieu * q - std::abs(tolerance * q), std::abs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = milieu;
                e = d;
            }
        } else {
            d = milieu;
            e = d;
        }
        a = b;
        fa = fb;
        b += std::abs(d) > tolerance ? d : (milieu > 0.0 ? tolerance : -tolerance);
        fb = f(b);
    }
    return b;
}

// Le solveur ne traite qu'une inconnue : les autres identifiants doivent être des paramètres
void exigerUneVariable(const ExprPtr& e) {
    if (e->plusieursVariables()) {
        throw std::invalid_argument("resoudre : l'expression '" + e->texte() +
                                    "' a plusieurs variables (une seule inconnue est admise)");
    }
}

} // namespace

namespace {

// Bilan écrit à la fin d'une résolution journalisée : une étape par solution retenue
void noterBilan(const Solutions& s, const ExprPtr& x) {
    journal::noter([&](Journal& j) {
        if (s.toutReel) {
            j.ajouter("resolution.bilan", "L'équation est vraie pour tout x réel");
        } else if (s.liste.empty()) {
            j.ajouter("resolution.bilan", s.complet ? "Aucune solution réelle" : "Aucune solution trouvée (résultat incomplet)");
        }
        for (const Solution& sol : s.liste) {
            std::string texte = sol.estFamille() ? "Famille de solutions (paramètres entiers relatifs)"
                                : sol.exacte     ? "Solution"
                                                 : "Solution (valeur approchée)";
            if (sol.multiplicite == 2) texte += ", racine double";
            if (sol.multiplicite > 2) texte += ", racine de multiplicité " + std::to_string(sol.multiplicite);
            j.ajouterEquation("resolution.solution", texte, x, sol.valeur, !sol.exacte && !sol.estFamille());
        }
        if (!s.complet && !s.liste.empty()) {
            j.ajouter("resolution.bilan", "Le résultat est incomplet : d'autres solutions peuvent exister");
        }
    });
}

} // namespace

Solutions resoudre(const ExprPtr& gauche, const ExprPtr& droite) {
    const ExprPtr f = gauche - droite;
    exigerUneVariable(f);
    const ExprPtr x = trouverVariable(f);
    etape("resolution.equation", "Équation à résoudre", gauche, droite);
    Solutions s;
    if (!x) {
        double v;
        const ExprPtr D = developper(f);
        if (estZeroExact(D)) {
            s.toutReel = true;
        } else if (!valeurNumerique(D, v)) {
            s.complet = false;
        }
        note("resolution.sans_inconnue", "Aucune inconnue : l'égalité est vérifiée ou non, sans calcul à faire");
        return s;
    }
    journal::noter([&](Journal& j) { // nombre(0) n'est construit que si un journal est actif
        if (!estZeroExact(droite))
            j.ajouterEquation("resolution.forme_nulle", "On se ramène à f(x) = 0 : tout dans le membre de gauche", f, nombre(Nombre(0)));
    });
    s = Resolveur(x).resoudre(gauche, droite);
    nettoyer(s, f, x);
    noterBilan(s, x);
    return s;
}

Solutions resoudre(const ExprPtr& gauche, const ExprPtr& droite, Journal& journal) {
    const journal::Activation activation(journal);
    return resoudre(gauche, droite);
}

Solutions resoudre(const ExprPtr& expression, Journal& journal) {
    const journal::Activation activation(journal);
    return resoudre(expression);
}

Solutions resoudre(const ExprPtr& expression) { return resoudre(expression, nombre(Nombre(0))); }

Solutions resoudreNumerique(const ExprPtr& f, double a, double b) {
    exigerUneVariable(f);
    Solutions s;
    s.complet = false;
    if (!(a < b)) return s;
    const ProgrammeEvaluation programme(f);
    if (!programme.estValide()) return s;
    const ProgrammeEvaluation programmeDerivee(f->derivee());
    auto valeur = [&](double t) { return programme.evaluer(t); };
    auto pente = [&](double t) { return programmeDerivee.evaluer(t); };

    constexpr std::size_t N = 4096;
    std::vector<double> xs(N + 1);
    for (std::size_t i = 0; i <= N; ++i) xs[i] = a + (b - a) * static_cast<double>(i) / static_cast<double>(N);
    const std::vector<double> ys = programme.evaluer(xs);
    double echelle = 1.0;
    for (double y : ys) {
        if (std::isfinite(y)) echelle = std::max(echelle, std::abs(y));
    }
    std::vector<double> racines;
    auto accepter = [&](double r) {
        const double fr = valeur(r);
        // Un changement de signe à travers un pôle (1/x en 0) n'est pas une racine
        if (std::isfinite(fr) && std::abs(fr) <= 1e-9 * echelle) racines.push_back(r);
    };
    for (std::size_t i = 0; i <= N; ++i) {
        if (ys[i] == 0.0) racines.push_back(xs[i]);
        if (i < N && std::isfinite(ys[i]) && std::isfinite(ys[i + 1]) && ys[i] * ys[i + 1] < 0.0) {
            accepter(brent(valeur, xs[i], xs[i + 1], ys[i], ys[i + 1]));
        }
    }
    // Racines de multiplicité paire (sans changement de signe) : points critiques où f s'annule
    if (programmeDerivee.estValide()) {
        const std::vector<double> dys = programmeDerivee.evaluer(xs);
        for (std::size_t i = 0; i < N; ++i) {
            if (std::isfinite(dys[i]) && std::isfinite(dys[i + 1]) && dys[i] * dys[i + 1] < 0.0) {
                const double c = brent(pente, xs[i], xs[i + 1], dys[i], dys[i + 1]);
                const double fc = valeur(c);
                if (std::isfinite(fc) && std::abs(fc) <= 1e-12 * echelle) racines.push_back(c);
            }
        }
    }
    std::sort(racines.begin(), racines.end());
    for (double r : racines) {
        if (!s.liste.empty() && std::abs(s.liste.back().approximation - r) <= 1e-9 * std::max(1.0, std::abs(r))) continue;
        Solution sol;
        sol.valeur = nombre(Nombre::reel(r));
        sol.approximation = r;
        sol.exacte = false;
        s.liste.push_back(sol);
    }
    return s;
}

Solutions resoudreSurIntervalle(const ExprPtr& expression, double a, double b) {
    exigerUneVariable(expression);
    journal::noter([&](Journal& j) {
        std::ostringstream os;
        os << "Résolution sur l'intervalle [" << a << ", " << b << "] : résolution exacte, puis dépliage des familles";
        j.ajouter("resolution.intervalle", os.str());
    });
    const Solutions exactes = resoudre(expression);
    Solutions s;
    s.complet = exactes.complet;
    s.toutReel = exactes.toutReel;
    constexpr long long MAX_ENTIERS = 100000;
    for (const Solution& sol : exactes.liste) {
        if (!sol.estFamille()) {
            if (sol.approximation >= a && sol.approximation <= b) s.liste.push_back(sol);
            continue;
        }
        if (sol.entiers.size() != 1) {
            s.complet = false;
            continue;
        }
        // Famille affine en k : v(k) = v0 + pas * k
        const ExprPtr k = param(sol.entiers[0]);
        auto enK = [&](long long n) { return substituer(sol.valeur, k, cst(static_cast<double>(n))); };
        double v0, v1, v2;
        if (!valeurNumerique(enK(0), v0) || !valeurNumerique(enK(1), v1) || !valeurNumerique(enK(2), v2) ||
            std::abs((v2 - v1) - (v1 - v0)) > 1e-9 * std::max(1.0, std::abs(v1 - v0))) {
            s.complet = false; // famille non affine : non dépliée
            note("resolution.famille_non_depliee", "Famille non affine en k : elle n'est pas dépliée sur l'intervalle");
            continue;
        }
        const double pas = v1 - v0;
        long long kMin = 0, kMax = 0;
        if (pas != 0.0) {
            const double k1 = (a - v0) / pas, k2 = (b - v0) / pas;
            kMin = static_cast<long long>(std::ceil(std::min(k1, k2) - 1e-9));
            kMax = static_cast<long long>(std::floor(std::max(k1, k2) + 1e-9));
        }
        if (kMax - kMin > MAX_ENTIERS) {
            s.complet = false;
            continue;
        }
        journal::noter([&](Journal& j) {
            j.ajouter("resolution.depliage", "Dépliage de la famille pour les entiers k de " + std::to_string(kMin) +
                                                 " à " + std::to_string(kMax) + " qui donnent une solution dans l'intervalle");
        });
        for (long long n = kMin; n <= kMax; ++n) {
            Solution depliee;
            depliee.valeur = enK(n);
            depliee.exacte = sol.exacte;
            if (!valeurNumerique(depliee.valeur, depliee.approximation)) continue;
            if (depliee.approximation >= a - 1e-12 && depliee.approximation <= b + 1e-12) s.liste.push_back(depliee);
        }
    }
    if (!s.complet && !s.toutReel) {
        note("resolution.numerique", "La résolution exacte est incomplète : recherche numérique des racines (échantillonnage et méthode de Brent)");
        for (const Solution& num : resoudreNumerique(expression, a, b).liste) s.liste.push_back(num);
    }
    const ExprPtr x = trouverVariable(expression);
    nettoyer(s, expression, x);
    if (x) noterBilan(s, x);
    return s;
}

Solutions resoudreSurIntervalle(const ExprPtr& expression, double a, double b, Journal& journal) {
    const journal::Activation activation(journal);
    return resoudreSurIntervalle(expression, a, b);
}

} // namespace symalgo
