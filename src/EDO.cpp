/**
 * @file EDO.cpp
 * @brief Résolution symbolique d'équations différentielles ordinaires non linéaires.
 *
 * Une équation du premier ordre affine en y' est mise sous la forme y' = f(x, y), ou
 * M dx + N dy = 0, puis les classes sont essayées de la moins chère à la plus générale :
 *
 *   quadrature        f sans y                       y = int f dx + C
 *   linéaire          f = -P(x) y + Q(x)             y = (int Q mu dx + C)/mu, mu = exp(int P)
 *   Bernoulli         f = a(x) y + b(x) y^n          linéaire en v = y^(1-n)
 *   séparable         f = g(x) h(y)                  int dy/h = int g dx + C
 *   affine            f = F(lambda x + y)            séparable en u = lambda x + y
 *   homogène          f(t x, t y) = f(x, y)          séparable en v = y/x (après translation si besoin)
 *   exacte            M_y = N_x                      Phi(x, y) = C ; sinon facteur intégrant mu(x),
 *                                                    mu(y), mu(x y) ou mu(x + y)
 *   Riccati           f = q0 + q1 y + q2 y^2         y = y1 + 1/v, y1 solution particulière
 *   Clairaut          y = x y' + g(y')               y = C x + g(C), et l'enveloppe
 *
 * Les classes polynomiales en y (linéaire, Bernoulli, Riccati) se lisent sur la décomposition
 * f = somme c_k(x) y^(e_k). Les autres tests de reconnaissance demandent de décider qu'une
 * expression est identiquement nulle : forme canonique, puis développement et
 * simplification, puis évaluation en des points irréguliers (estNul). Ce dernier étage n'est
 * qu'une forte présomption ; c'est pourquoi chaque famille obtenue est réinjectée dans
 * l'équation de départ (verifierEDO) et écartée si elle ne la vérifie pas.
 *
 * Second ordre affine en y'' : si y manque, l'équation est du premier ordre en p = y' ; si x
 * manque, p dp/dy = f(y, p) est du premier ordre en p(y), puis dy/p = dx.
 *
 * Les intégrales par rapport à une variable gèlent les autres (remplacées par des
 * paramètres) : elles héritent des limites de integrer() pour les coefficients symboliques.
 */

#include "EDO.hpp"

#include "Canonique.hpp"
#include "Lecture.hpp"
#include "Multivariable.hpp"
#include "Polynome.hpp"
#include "Primitives.hpp"
#include "Solveur.hpp"

#include <cctype>
#include <cmath>
#include <exception>
#include <utility>

namespace symalgo {

namespace {

// ============================================================================
// Outils : nullité, gel de variable, exponentielle simplifiée, isolement
// ============================================================================

bool estZero(const ExprPtr& e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().estZero();
}

ExprPtr zero() { return nombre(Nombre(0)); }

const std::string& nomDe(const ExprPtr& variable) { return static_cast<const Variable&>(*variable).getNom(); }

void collecterParametres(const ASTNode& e, std::vector<ExprPtr>& liste) {
    switch (e.type()) {
        case TypeNoeud::Parametre:
            for (const ExprPtr& p : liste) {
                if (p.get() == &e) return;
            }
            liste.push_back(e.clone());
            return;
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) collecterParametres(*t.expression, liste);
            return;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                collecterParametres(*f.base, liste);
                collecterParametres(*f.exposant, liste);
            }
            return;
        case TypeNoeud::Puissance:
            collecterParametres(*static_cast<const Puissance&>(e).getBase(), liste);
            collecterParametres(*static_cast<const Puissance&>(e).getExposant(), liste);
            return;
        default:
            if (const FonctionUnaire* f = comme<FonctionUnaire>(&e)) collecterParametres(*f->m_argument, liste);
    }
}

// Valeurs d'essai irrégulières, petites et grandes, pour tomber dans les domaines étroits
const double ESSAIS[] = {0.731, 1.317, 0.413, 1.873, 0.123, 0.577, 2.239, 0.291,
                         1.129, 0.067, 2.713, 0.867, 0.189, 1.621, 0.347, 3.141};
constexpr int NOMBRE_ESSAIS = 16;
// Constantes d'intégration : des deux signes et de tailles variées (domaines du type C < -exp(x))
const double CONSTANTES[] = {0.731, -1.317, 2.713, -0.413, 6.5, -5.3, 17.0, -23.0, 0.189, -0.067, 40.0, -61.0};
constexpr int NOMBRE_CONSTANTES = 12;
constexpr int NOMBRE_POINTS = 40;

enum class Verdict { Nul, NonNul, Indecis };

/*
 * Nom : nulNumerique
 * Description : Évalue e en des points irréguliers (variables et paramètres) : NonNul dès
 *               qu'une valeur dépasse 1e-7 fois la taille de ses termes, Nul si au moins 4
 *               points sont évaluables et tous nuls, Indecis sinon (hors domaine partout).
 *               Avec partout = faux, e est dite nulle si elle l'est en au moins 4 points,
 *               même si elle ne l'est pas ailleurs : c'est le cas d'une solution valable sur
 *               une partie du domaine seulement (y = (C exp(x/2) - 1)^2 pour y' = y + sqrt(y),
 *               tant que C exp(x/2) >= 1). Une coïncidence en 4 points irréguliers est exclue.
 */
Verdict nulNumerique(const ExprPtr& e, bool partout = true) {
    const std::vector<ExprPtr> vars = variables(e);
    std::vector<ExprPtr> parametres;
    collecterParametres(*e, parametres);
    std::vector<ExprPtr> termes; // pour l'échelle
    if (const Somme* s = comme<Somme>(e)) {
        for (const Terme& t : s->getTermes()) termes.push_back(t.expression * nombre(t.coefficient));
    }
    int valides = 0, nonNuls = 0;
    for (int k = 0; k < NOMBRE_POINTS; ++k) {
        Valeurs valeurs;
        for (std::size_t i = 0; i < vars.size(); ++i) {
            double v = ESSAIS[(static_cast<std::size_t>(k) * 3 + i * 5 + static_cast<std::size_t>(k) / NOMBRE_ESSAIS) % NOMBRE_ESSAIS];
            if (k % 5 == 4 && i == 0) v = -v; // quelques abscisses négatives
            valeurs[nomDe(vars[i])] = v;
        }
        ExprPtr f = e;
        std::vector<ExprPtr> fs = termes;
        for (std::size_t i = 0; i < parametres.size(); ++i) {
            const ExprPtr valeur = cst(CONSTANTES[(static_cast<std::size_t>(k) * 5 + i * 3 + 1) % NOMBRE_CONSTANTES]);
            f = substituer(f, parametres[i], valeur);
            for (ExprPtr& t : fs) t = substituer(t, parametres[i], valeur);
        }
        double valeur, echelle = 1.0;
        try {
            valeur = evaluer(f, valeurs);
            for (const ExprPtr& t : fs) {
                const double v = evaluer(t, valeurs);
                if (std::isfinite(v)) echelle += std::fabs(v);
            }
        } catch (const std::exception&) {
            continue;
        }
        if (!std::isfinite(valeur)) continue;
        if (std::fabs(valeur) > 1e-7 * echelle) {
            if (partout) return Verdict::NonNul;
            ++nonNuls;
            continue;
        }
        ++valides;
    }
    if (valides >= 4) return Verdict::Nul;
    return nonNuls > 0 ? Verdict::NonNul : Verdict::Indecis;
}

/*
 * Nom : estNul
 * Description : e est-elle identiquement nulle ? Certain si la forme canonique, ou le
 *               développement simplifié, est 0 ; sinon décidé numériquement (nulNumerique),
 *               ce qui n'est qu'une forte présomption : les solutions sont réinjectées.
 */
bool estNul(const ExprPtr& e, bool partout = true) {
    if (estZero(e)) return true;
    ExprPtr s;
    try {
        s = developper(e)->simplifier();
    } catch (const std::exception&) {
        return false;
    }
    if (estZero(s)) return true;
    if (s->type() == TypeNoeud::Constante) return false;
    return nulNumerique(s, partout) == Verdict::Nul;
}

/*
 * Nom : sansVariable
 * Description : Si e ne dépend pas de v (par sa forme, ou parce que sa dérivée par rapport à
 *               v est identiquement nulle), écrit dans reduite une expression égale à e où v
 *               n'apparaît plus, et renvoie vrai.
 */
bool sansVariable(const ExprPtr& e, const ExprPtr& v, ExprPtr& reduite) {
    if (!contient(e, v)) {
        reduite = e;
        return true;
    }
    const ExprPtr s = e->simplifier();
    if (!contient(s, v)) {
        reduite = s;
        return true;
    }
    if (!estNul(s->derivee(v))) return false;
    // e ne dépend pas de v : on y fixe v, en un point où elle est définie
    for (const ExprPtr& valeur : {frac(1), frac(2), frac(1, 2), frac(3)}) {
        const ExprPtr essai = substituer(s, v, valeur)->simplifier();
        if (contient(essai, v)) continue;
        if (estNul(essai - s)) {
            reduite = essai;
            return true;
        }
    }
    return false;
}

/*
 * Nom : integrerPar
 * Description : Primitive de e par rapport à v, les autres variables étant gelées (traitées
 *               comme des paramètres). Pointeur nul si la primitive n'est pas complète.
 */
ExprPtr integrerPar(const ExprPtr& e, const ExprPtr& v) {
    if (estZero(e)) return e;
    if (!contient(e, v)) return e * v;
    ExprPtr f = e;
    std::vector<std::pair<ExprPtr, ExprPtr>> gels;
    for (const ExprPtr& w : variables(e)) {
        if (w.get() == v.get()) continue;
        const ExprPtr gel = param("_" + nomDe(w));
        f = substituer(f, w, gel);
        gels.emplace_back(gel, w);
    }
    ExprPtr F = f->integrer();
    if (contientNonEvaluee(*F)) {
        const ExprPtr simple = f->simplifier();
        if (simple.get() == f.get()) return nullptr;
        F = simple->integrer();
        if (contientNonEvaluee(*F)) return nullptr;
    }
    for (const auto& [gel, w] : gels) F = substituer(F, gel, w);
    return F;
}

/*
 * Nom : exponentielle
 * Description : exp(argument), où les logarithmes sortent en puissances :
 *               exp(c ln(w) + r) = w^c exp(r). ln|w| perd sa valeur absolue quand c est
 *               entier : le signe est constant sur chaque intervalle et passe dans la
 *               constante d'intégration (ou ne change pas un facteur intégrant).
 */
ExprPtr exponentielle(const ExprPtr& argument) {
    AccumulateurProduit produit;
    AccumulateurSomme reste;
    const auto terme = [&](const ExprPtr& e, const Nombre& c) {
        if (e->type() != TypeNoeud::Logarithme) {
            reste.ajouter(e, c);
            return;
        }
        ExprPtr w = static_cast<const Logarithme&>(*e).m_argument;
        if (w->type() == TypeNoeud::ValeurAbsolue && c.estEntier()) w = static_cast<const ValeurAbsolue&>(*w).m_argument;
        produit.multiplier(ast_pow(w, nombre(c)));
    };
    if (const Somme* s = comme<Somme>(argument)) {
        reste.ajouter(nombre(s->getConstante()), Nombre(1));
        for (const Terme& t : s->getTermes()) terme(t.expression, t.coefficient);
    } else if (const Produit* p = comme<Produit>(argument);
               p && p->getFacteurs().size() == 1 && p->getFacteurs()[0].exposant.get() == un().get()) {
        terme(p->getFacteurs()[0].base, p->getCoefficient());
    } else {
        terme(argument, Nombre(1));
    }
    const ExprPtr r = reste.construire();
    if (!estZero(r)) produit.multiplier(ast_exp(r));
    return produit.construire();
}

// Exposant p/q avec p pair : u^(p/q) = d a deux solutions réelles opposées
bool exposantPair(const ExprPtr& exposant) {
    long long num, den;
    const Constante* c = comme<Constante>(exposant);
    return c && c->getNombre().versFraction(num, den) && num % 2 == 0;
}

/*
 * Nom : isoler
 * Description : Résout gauche = droite en v, quand v n'apparaît qu'à un endroit de gauche (et
 *               pas dans droite), en défaisant les opérations de l'extérieur vers
 *               l'intérieur. Une puissance paire ou une valeur absolue donne deux branches.
 *               Les fonctions périodiques ne donnent que leur branche principale. Renvoie
 *               faux si v ne peut pas être isolée ainsi.
 */
bool isoler(ExprPtr gauche, ExprPtr droite, const ExprPtr& v, std::vector<ExprPtr>& sorties) {
    for (int tour = 0; tour < 64; ++tour) {
        if (gauche.get() == v.get()) {
            sorties.push_back(droite);
            return true;
        }
        ExprPtr base, exposant; // puissance à défaire, le cas échéant
        switch (gauche->type()) {
            case TypeNoeud::Somme: {
                const Somme& s = static_cast<const Somme&>(*gauche);
                AccumulateurSomme sans;
                sans.ajouter(nombre(s.getConstante()), Nombre(1));
                const Terme* avec = nullptr;
                int nombreAvec = 0;
                for (const Terme& t : s.getTermes()) {
                    if (contient(t.expression, v)) {
                        avec = &t;
                        ++nombreAvec;
                    } else {
                        sans.ajouter(t.expression, t.coefficient);
                    }
                }
                if (nombreAvec != 1) return false;
                droite = (droite - sans.construire()) / nombre(avec->coefficient);
                gauche = avec->expression;
                continue;
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(*gauche);
                AccumulateurProduit sans;
                sans.multiplier(nombre(p.getCoefficient()));
                int nombreAvec = 0;
                for (const Facteur& f : p.getFacteurs()) {
                    if (contient(f.base, v) || contient(f.exposant, v)) {
                        base = f.base;
                        exposant = f.exposant;
                        ++nombreAvec;
                    } else {
                        sans.multiplier(ast_pow(f.base, f.exposant));
                    }
                }
                if (nombreAvec != 1) return false;
                droite = droite / sans.construire();
                break;
            }
            case TypeNoeud::Puissance:
                base = static_cast<const Puissance&>(*gauche).getBase();
                exposant = static_cast<const Puissance&>(*gauche).getExposant();
                break;
            case TypeNoeud::Logarithme:
                gauche = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                droite = exponentielle(droite);
                continue;
            case TypeNoeud::Exponentielle:
                gauche = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                droite = ast_ln(droite);
                continue;
            case TypeNoeud::ValeurAbsolue: {
                const ExprPtr& u = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                return isoler(u, droite, v, sorties) && isoler(u, -droite, v, sorties);
            }
            case TypeNoeud::Tangente:
            case TypeNoeud::ArcTangente:
            case TypeNoeud::Sinus:
            case TypeNoeud::ArcSinus:
            case TypeNoeud::Cosinus:
            case TypeNoeud::ArcCosinus: {
                static const std::pair<TypeNoeud, TypeNoeud> RECIPROQUES[] = {
                    {TypeNoeud::Tangente, TypeNoeud::ArcTangente}, {TypeNoeud::ArcTangente, TypeNoeud::Tangente},
                    {TypeNoeud::Sinus, TypeNoeud::ArcSinus},       {TypeNoeud::ArcSinus, TypeNoeud::Sinus},
                    {TypeNoeud::Cosinus, TypeNoeud::ArcCosinus},   {TypeNoeud::ArcCosinus, TypeNoeud::Cosinus}};
                for (const auto& [fonction, reciproque] : RECIPROQUES) {
                    if (fonction == gauche->type()) {
                        droite = appliquer(reciproque, droite);
                        break;
                    }
                }
                gauche = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                continue;
            }
            case TypeNoeud::TangenteHyperbolique: // tanh(u) = d : u = ln((1 + d)/(1 - d))/2
                gauche = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                droite = ast_ln((1.0 + droite) / (1.0 - droite)) / 2.0;
                continue;
            case TypeNoeud::SinusHyperbolique: // sinh(u) = d : u = ln(d + sqrt(d^2 + 1))
                gauche = static_cast<const FonctionUnaire&>(*gauche).m_argument;
                droite = ast_ln(droite + ast_pow(ast_pow(droite, 2.0) + 1.0, frac(1, 2)));
                continue;
            default: return false;
        }
        // base^exposant = droite
        if (!contient(exposant, v)) {
            const ExprPtr racine = ast_pow(droite, un() / exposant);
            if (exposantPair(exposant)) return isoler(base, racine, v, sorties) && isoler(base, -racine, v, sorties);
            gauche = base;
            droite = racine;
        } else if (!contient(base, v)) {
            gauche = exposant;
            droite = ast_ln(droite) / ast_ln(base);
        } else {
            return false;
        }
    }
    return false;
}

/*
 * Nom : expliciter
 * Description : Solutions en v de H(v) = G + C, C constante arbitraire (G sans v). Si H est
 *               un seul logarithme c ln|w|, la constante est absorbée : w = C exp(G/c), qui
 *               contient aussi w = 0. Renvoie faux si v ne s'isole pas.
 */
bool expliciter(const ExprPtr& H, const ExprPtr& G, const ExprPtr& v, const ExprPtr& C, std::vector<ExprPtr>& sorties) {
    ExprPtr logarithme = H;
    Nombre c(1);
    if (const Produit* p = comme<Produit>(H); p && p->getFacteurs().size() == 1 && p->getFacteurs()[0].exposant.get() == un().get()) {
        logarithme = p->getFacteurs()[0].base;
        c = p->getCoefficient();
    }
    if (logarithme->type() == TypeNoeud::Logarithme) {
        ExprPtr w = static_cast<const Logarithme&>(*logarithme).m_argument;
        if (w->type() == TypeNoeud::ValeurAbsolue) w = static_cast<const ValeurAbsolue&>(*w).m_argument;
        return isoler(w, C * exponentielle(G / nombre(c)), v, sorties);
    }
    return isoler(H, G + C, v, sorties);
}

// ============================================================================
// Vérification par réinjection
// ============================================================================

ExprPtr residu(const ExprPtr& equation, const InconnuesEDO& inc, const SolutionEDO& solution) {
    ExprPtr dy, d2y, y = inc.y;
    if (solution.y) {
        y = solution.y;
        dy = y->derivee(inc.x);
        d2y = dy->derivee(inc.x);
    } else {
        // relation R(x, y) = 0 : y' = -R_x/R_y, y'' = -(R_xx + 2 R_xy y' + R_yy y'^2)/R_y
        const ExprPtr& R = solution.relation;
        const ExprPtr Rx = R->derivee(inc.x), Ry = R->derivee(inc.y);
        dy = -Rx / Ry;
        d2y = -(Rx->derivee(inc.x) + 2.0 * Rx->derivee(inc.y) * dy + Ry->derivee(inc.y) * ast_pow(dy, 2.0)) / Ry;
    }
    // d2y et dy d'abord : leurs expressions contiennent encore la variable y (forme implicite)
    ExprPtr r = substituer(equation, inc.d2y, d2y);
    r = substituer(r, inc.dy, dy);
    if (solution.y) r = substituer(r, inc.y, y);
    return r;
}

// ============================================================================
// Décomposition f = somme c_k(x) y^(e_k)
// ============================================================================

struct TermeEnY {
    Nombre exposant;
    ExprPtr coefficient; // sans y
};

bool decomposer(const ExprPtr& f, const ExprPtr& y, std::vector<TermeEnY>& termes) {
    ExprPtr d;
    try {
        d = developper(f);
    } catch (const std::exception&) {
        return false;
    }
    const auto ajouter = [&](const ExprPtr& t, const Nombre& coefficient) {
        Nombre exposant(0);
        AccumulateurProduit c;
        c.multiplier(nombre(coefficient));
        const auto facteur = [&](const ExprPtr& base, const ExprPtr& e) {
            if (base.get() == y.get()) {
                const Constante* k = comme<Constante>(e);
                if (!k || !k->getNombre().estExact()) return false;
                exposant = exposant + k->getNombre();
                return true;
            }
            if (contient(base, y) || contient(e, y)) return false;
            c.multiplier(ast_pow(base, e));
            return true;
        };
        if (const Produit* p = comme<Produit>(t)) {
            c.multiplier(nombre(p->getCoefficient()));
            for (const Facteur& fa : p->getFacteurs()) {
                if (!facteur(fa.base, fa.exposant)) return false;
            }
        } else if (const Puissance* p = comme<Puissance>(t)) {
            if (!facteur(p->getBase(), p->getExposant())) return false;
        } else if (!facteur(t, un())) {
            return false;
        }
        const ExprPtr valeur = c.construire();
        for (TermeEnY& autre : termes) {
            if (autre.exposant == exposant) {
                autre.coefficient = autre.coefficient + valeur;
                return true;
            }
        }
        termes.push_back({exposant, valeur});
        return true;
    };
    if (const Somme* s = comme<Somme>(d)) {
        if (!s->getConstante().estZero() && !ajouter(un(), s->getConstante())) return false;
        for (const Terme& t : s->getTermes()) {
            if (!ajouter(t.expression, t.coefficient)) return false;
        }
    } else if (!ajouter(d, Nombre(1))) {
        return false;
    }
    for (std::size_t i = termes.size(); i-- > 0;) {
        if (estZero(termes[i].coefficient)) termes.erase(termes.begin() + static_cast<std::ptrdiff_t>(i));
    }
    return true;
}

ExprPtr coefficientDe(const std::vector<TermeEnY>& termes, const Nombre& exposant) {
    for (const TermeEnY& t : termes) {
        if (t.exposant == exposant) return t.coefficient;
    }
    return zero();
}

// ============================================================================
// Premier ordre
// ============================================================================

struct Probleme {
    InconnuesEDO inc;
    ExprPtr equation; // F(x, y, y') = 0, telle que donnée
    ExprPtr f;        // y' = f(x, y) (équation affine en y')
    ExprPtr M, N;     // M dx + N dy = 0
    int* constantes;  // compteur partagé des constantes C1, C2...

    ExprPtr constante() const { return param("C" + std::to_string(++*constantes)); }
};

ExprPtr variableAuxiliaire(const InconnuesEDO& inc, const char* nom, const char* secours) {
    return var(nomDe(inc.x) == nom || nomDe(inc.y) == nom ? secours : nom);
}

void ajouterExplicite(SolutionsEDO& r, ExprPtr y, bool singuliere = false) {
    // La forme simplifiée, si elle est plus courte : 1/(C1/x - ln(x)/x) s'écrit x/(C1 - ln(x))
    try {
        const ExprPtr simple = y->simplifier();
        if (simple->texte().size() < y->texte().size()) y = simple;
    } catch (const std::exception&) {
    }
    for (const SolutionEDO& s : r.liste) {
        if (s.y.get() == y.get()) return;
    }
    r.liste.push_back({y, ExprPtr(), singuliere});
}

/*
 * Nom : lineaire
 * Description : Solution générale de y' + P y = Q : (int Q mu dx + C)/mu, mu = exp(int P dx).
 *               Pointeur nul si une des deux primitives manque.
 */
ExprPtr lineaire(const ExprPtr& P, const ExprPtr& Q, const ExprPtr& x, const ExprPtr& C) {
    const ExprPtr I = integrerPar(P, x);
    if (!I) return nullptr;
    const ExprPtr mu = exponentielle(I);
    ExprPtr particuliere = zero();
    if (!estZero(Q)) {
        const ExprPtr J = integrerPar(Q * mu, x);
        if (!J) return nullptr;
        particuliere = developper(J / mu);
    }
    return particuliere + C * exponentielle(-I);
}

bool resoudreQuadrature(const Probleme& p, SolutionsEDO& r) {
    if (contient(p.f, p.inc.y)) return false;
    const ExprPtr F = integrerPar(p.f, p.inc.x);
    if (!F) return false;
    ajouterExplicite(r, F + p.constante());
    r.methode = "quadrature";
    return true;
}

bool resoudreLineaire(const Probleme& p, const std::vector<TermeEnY>& termes, SolutionsEDO& r) {
    for (const TermeEnY& t : termes) {
        if (!t.exposant.estZero() && !t.exposant.estUn()) return false;
    }
    const ExprPtr P = -coefficientDe(termes, Nombre(1)), Q = coefficientDe(termes, Nombre(0));
    const ExprPtr y = lineaire(P, Q, p.inc.x, p.constante());
    if (!y) return false;
    journal::noter([&](Journal& j) {
        j.ajouter("edo.lineaire", "Équation linéaire y' + P y = Q : y = (int Q mu dx + C)/mu avec mu = exp(int P dx)", P, Q);
    });
    ajouterExplicite(r, y);
    r.methode = "lineaire";
    return true;
}

bool resoudreBernoulli(const Probleme& p, const std::vector<TermeEnY>& termes, SolutionsEDO& r) {
    // f = a y + b y^n, n différent de 0 et de 1
    if (termes.empty() || termes.size() > 2) return false;
    const Nombre* n = nullptr;
    for (const TermeEnY& t : termes) {
        if (t.exposant.estUn()) continue;
        if (n || t.exposant.estZero()) return false;
        n = &t.exposant;
    }
    if (!n) return false;
    const Nombre k = Nombre(1) - *n; // v = y^k
    const ExprPtr a = coefficientDe(termes, Nombre(1)), b = coefficientDe(termes, *n);
    // v' = k a v + k b
    const ExprPtr v = lineaire(-(a * nombre(k)), b * nombre(k), p.inc.x, p.constante());
    if (!v) return false;
    journal::noter([&](Journal& j) {
        j.ajouter("edo.bernoulli", "Équation de Bernoulli : v = y^(" + k.texte() + ") vérifie une équation linéaire", p.f, v);
    });
    const ExprPtr y = ast_pow(v, nombre(k.inverse()));
    ajouterExplicite(r, y);
    long long num, den;
    if (k.versFraction(num, den) && num % 2 == 0) ajouterExplicite(r, -y); // y^k = v a deux racines réelles
    if (n->signe() > 0) ajouterExplicite(r, zero(), true);
    r.methode = "bernoulli";
    return true;
}

// Équilibres : zéros exacts de h(v), ajoutés comme solutions singulières par construire(v0)
template <class Construire>
void ajouterEquilibres(const ExprPtr& h, SolutionsEDO& r, const Construire& construire) {
    Solutions zeros;
    try {
        zeros = resoudre(h);
    } catch (const std::exception&) {
        return;
    }
    if (zeros.toutReel) return;
    for (const Solution& z : zeros.liste) {
        if (z.exacte && !z.estFamille()) ajouterExplicite(r, construire(z.valeur), true);
    }
}

/*
 * Nom : separer
 * Description : f = g(x) h(y) ? D'abord sur la forme (produit dont chaque facteur ne contient
 *               qu'une des deux variables), sinon avec g = f(x, y0) et h = f(x0, y)/f(x0, y0)
 *               en un point où f est finie et non nulle, l'identité f(x, y) f(x0, y0) =
 *               f(x, y0) f(x0, y) étant vérifiée par estNul.
 */
bool separer(const ExprPtr& f, const ExprPtr& x, const ExprPtr& y, ExprPtr& g, ExprPtr& h) {
    AccumulateurProduit enX, enY;
    bool mixte = false;
    const auto facteur = [&](const ExprPtr& e) {
        const bool aX = contient(e, x), aY = contient(e, y);
        if (aX && aY) mixte = true;
        (aY ? enY : enX).multiplier(e);
    };
    if (const Produit* p = comme<Produit>(f)) {
        enX.multiplier(nombre(p->getCoefficient()));
        for (const Facteur& fa : p->getFacteurs()) facteur(ast_pow(fa.base, fa.exposant));
    } else {
        facteur(f);
    }
    if (!mixte) {
        g = enX.construire();
        h = enY.construire();
        return true;
    }
    const ExprPtr points[][2] = {{frac(0), frac(0)}, {frac(1), frac(0)}, {frac(0), frac(1)},
                                 {frac(1), frac(1)}, {frac(2), frac(1)}, {frac(1), frac(2)}};
    for (const auto& point : points) {
        double valeur;
        try {
            valeur = evaluer(f, {{nomDe(x), point[0]->getValeurConstante()}, {nomDe(y), point[1]->getValeurConstante()}});
        } catch (const std::exception&) {
            return false; // paramètre ou noeud non évaluable
        }
        if (!std::isfinite(valeur) || std::fabs(valeur) < 1e-9) continue;
        const ExprPtr gx = substituer(f, y, point[1]), hy = substituer(f, x, point[0]), f00 = substituer(gx, x, point[0]);
        if (!estNul(f * f00 - gx * hy)) return false;
        g = gx;
        h = hy / f00;
        return true;
    }
    return false;
}

bool resoudreSeparable(const Probleme& p, SolutionsEDO& r) {
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    ExprPtr g, h;
    if (!separer(p.f, x, y, g, h)) return false;
    const ExprPtr H = integrerPar(un() / h, y), G = integrerPar(g, x);
    if (!H || !G) return false;
    journal::noter([&](Journal& j) {
        j.ajouter("edo.separable", "Variables séparables : int dy/h(y) = int g(x) dx + C", un() / h, g);
    });
    const ExprPtr C = p.constante();
    std::vector<ExprPtr> explicites;
    if (expliciter(H, G, y, C, explicites)) {
        for (const ExprPtr& e : explicites) ajouterExplicite(r, e);
    } else {
        r.liste.push_back({ExprPtr(), H - G - C, false});
    }
    ajouterEquilibres(h, r, [](const ExprPtr& y0) { return y0; });
    r.methode = "separable";
    return true;
}

bool resoudreAffine(const Probleme& p, SolutionsEDO& r) {
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    const ExprPtr fx = p.f->derivee(x), fy = p.f->derivee(y);
    if (estZero(fx) || estZero(fy)) return false;
    // f = F(lambda x + y) : f_x/f_y = lambda, constante
    ExprPtr sansX, lambda;
    if (!sansVariable(fx / fy, x, sansX) || !sansVariable(sansX, y, lambda)) return false;
    const ExprPtr u = variableAuxiliaire(p.inc, "u", "w");
    ExprPtr F;
    if (!sansVariable(substituer(p.f, y, u - lambda * x), x, F)) return false;
    // u' = lambda + F(u)
    const ExprPtr h = lambda + F;
    const ExprPtr H = integrerPar(un() / h, u);
    if (!H) return false;
    const ExprPtr combinaison = lambda * x + y;
    journal::noter([&](Journal& j) {
        j.ajouter("edo.affine", "y' = F(u) avec u = " + combinaison->texte() + " : u' = " + h->texte() + ", séparable", p.f, h);
    });
    const ExprPtr C = p.constante();
    std::vector<ExprPtr> explicites;
    if (expliciter(H, x, u, C, explicites)) {
        for (const ExprPtr& e : explicites) ajouterExplicite(r, e - lambda * x);
    } else {
        r.liste.push_back({ExprPtr(), substituer(H, u, combinaison) - x - C, false});
    }
    ajouterEquilibres(h, r, [&](const ExprPtr& u0) { return u0 - lambda * x; });
    r.methode = "affine";
    return true;
}

/*
 * Nom : homogene
 * Description : Résout y' = f(x, y) quand f(t x, t y) = f(x, y) : v = y/x vérifie
 *               x v' = F(v) - v. Les solutions sont exprimées en x et y.
 */
bool homogene(const Probleme& p, const ExprPtr& f, SolutionsEDO& r) {
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    // f(t x, t y) = f(x, y), y compris pour t < 0
    for (const ExprPtr& t : {frac(2), frac(-1), frac(1, 3)}) {
        const ExprPtr dilate = substituer(substituer(f, x, t * x), y, t * y);
        if (!estNul(dilate - f)) return false;
    }
    const ExprPtr v = variableAuxiliaire(p.inc, "v", "w");
    const ExprPtr F = substituer(substituer(f, x, un()), y, v); // f = F(y/x)
    const ExprPtr D = (F - v)->simplifier();                    // x v' = F(v) - v
    if (estZero(D)) return false;
    const ExprPtr H = integrerPar(un() / D, v);
    if (!H) return false;
    journal::noter([&](Journal& j) {
        j.ajouter("edo.homogene", "Équation homogène : v = " + (y / x)->texte() + " vérifie x v' = F(v) - v, séparable", f, D);
    });
    const ExprPtr C = p.constante();
    std::vector<ExprPtr> explicites;
    if (expliciter(H, ast_ln(ast_abs(x)), v, C, explicites)) {
        for (const ExprPtr& e : explicites) ajouterExplicite(r, developper(x * e));
    } else {
        r.liste.push_back({ExprPtr(), substituer(H, v, y / x) - ast_ln(ast_abs(x)) - C, false});
    }
    ajouterEquilibres(D, r, [&](const ExprPtr& v0) { return v0 * x; });
    return true;
}

// a x + b y + c à coefficients exacts
bool droite(const ASTNode& e, const ExprPtr& x, const ExprPtr& y, Nombre& a, Nombre& b, Nombre& c) {
    const Somme* s = comme<Somme>(&e);
    if (!s || !s->getConstante().estExact()) return false;
    a = b = Nombre(0);
    c = s->getConstante();
    for (const Terme& t : s->getTermes()) {
        if (!t.coefficient.estExact()) return false;
        if (t.expression.get() == x.get()) {
            a = t.coefficient;
        } else if (t.expression.get() == y.get()) {
            b = t.coefficient;
        } else {
            return false;
        }
    }
    return true;
}

void collecterDroites(const ASTNode& e, const ExprPtr& x, const ExprPtr& y, std::vector<std::vector<Nombre>>& droites) {
    Nombre a, b, c;
    if (droite(e, x, y, a, b, c)) {
        droites.push_back({a, b, c});
        return;
    }
    switch (e.type()) {
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) collecterDroites(*t.expression, x, y, droites);
            return;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                collecterDroites(*f.base, x, y, droites);
                collecterDroites(*f.exposant, x, y, droites);
            }
            return;
        case TypeNoeud::Puissance:
            collecterDroites(*static_cast<const Puissance&>(e).getBase(), x, y, droites);
            collecterDroites(*static_cast<const Puissance&>(e).getExposant(), x, y, droites);
            return;
        default:
            if (const FonctionUnaire* f = comme<FonctionUnaire>(&e)) collecterDroites(*f->m_argument, x, y, droites);
    }
}

/*
 * Nom : resoudreHomogene
 * Description : Équation homogène, ou qui le devient par la translation amenant à l'origine
 *               le point d'intersection de deux droites a x + b y + c présentes dans f
 *               (y' = F((a x + b y + c)/(d x + e y + f))).
 */
bool resoudreHomogene(const Probleme& p, SolutionsEDO& r) {
    r.methode = "homogene";
    if (homogene(p, p.f, r)) return true;
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    std::vector<std::vector<Nombre>> droites;
    collecterDroites(*p.f, x, y, droites);
    for (std::size_t i = 0; i < droites.size() && i < 4; ++i) {
        for (std::size_t k = i + 1; k < droites.size() && k < 4; ++k) {
            const std::vector<Nombre>& d1 = droites[i];
            const std::vector<Nombre>& d2 = droites[k];
            const Nombre determinant = d1[0] * d2[1] - d1[1] * d2[0];
            if (determinant.estZero() || (d1[2].estZero() && d2[2].estZero())) continue;
            const ExprPtr x0 = nombre((d1[1] * d2[2] - d1[2] * d2[1]) / determinant);
            const ExprPtr y0 = nombre((d1[2] * d2[0] - d1[0] * d2[2]) / determinant);
            // X = x - x0, Y = y - y0 : même équation en (X, Y), écrite avec les mêmes variables
            const ExprPtr translate = substituer(substituer(p.f, x, x + x0), y, y + y0);
            SolutionsEDO t;
            if (!homogene(p, translate, t)) continue;
            journal::noter([&](Journal& j) {
                j.ajouter("edo.translation", "Translation de l'origine en (" + x0->texte() + ", " + y0->texte() + ")", p.f, translate);
            });
            for (const SolutionEDO& sol : t.liste) {
                if (sol.y) {
                    ajouterExplicite(r, y0 + substituer(sol.y, x, x - x0), sol.singuliere);
                } else {
                    r.liste.push_back({ExprPtr(), substituer(substituer(sol.relation, x, x - x0), y, y - y0), sol.singuliere});
                }
            }
            return true;
        }
    }
    return false;
}

/*
 * Nom : potentiel
 * Description : Phi telle que dPhi = M dx + N dy si la forme est exacte (M_y = N_x) :
 *               Phi = int M dx + int (N - d/dy int M dx) dy. Pointeur nul sinon.
 */
ExprPtr potentiel(const ExprPtr& M, const ExprPtr& N, const ExprPtr& x, const ExprPtr& y) {
    if (!estNul(M->derivee(y) - N->derivee(x))) return nullptr;
    const ExprPtr F = integrerPar(M, x);
    if (!F) return nullptr;
    ExprPtr g;
    if (!sansVariable(N - F->derivee(y), x, g)) return nullptr;
    const ExprPtr G = integrerPar(g, y);
    return G ? F + G : nullptr;
}

bool resoudreExacte(const Probleme& p, SolutionsEDO& r) {
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    ExprPtr Phi = potentiel(p.M, p.N, x, y);
    r.methode = "exacte";
    if (!Phi) {
        // Facteur intégrant mu(x) : (M_y - N_x)/N sans y ; mu(y) : (N_x - M_y)/M sans x
        const ExprPtr ecart = p.M->derivee(y) - p.N->derivee(x);
        ExprPtr k, mu;
        if (sansVariable(ecart / p.N, y, k)) {
            if (const ExprPtr I = integrerPar(k, x)) mu = exponentielle(I);
        }
        if (!mu && sansVariable(-ecart / p.M, x, k)) {
            if (const ExprPtr I = integrerPar(k, y)) mu = exponentielle(I);
        }
        // mu(z), z = x y : mu'/mu = (N_x - M_y)/(x M - y N) ; z = x + y : (N_x - M_y)/(M - N)
        const ExprPtr z = variableAuxiliaire(p.inc, "z", "w");
        const std::pair<ExprPtr, ExprPtr> essais[] = {{x * p.M - y * p.N, x * y}, {p.M - p.N, x + y}};
        for (const auto& [diviseur, combinaison] : essais) {
            if (mu) break;
            if (estNul(diviseur)) continue;
            // y exprimé avec z : y = z/x ou y = z - x
            const ExprPtr yDeZ = combinaison->type() == TypeNoeud::Produit ? z / x : z - x;
            if (!sansVariable(substituer(-ecart / diviseur, y, yDeZ), x, k)) continue;
            if (const ExprPtr I = integrerPar(k, z)) mu = substituer(exponentielle(I), z, combinaison);
        }
        if (!mu) return false;
        Phi = potentiel((p.M * mu)->simplifier(), (p.N * mu)->simplifier(), x, y);
        if (!Phi) return false;
        journal::noter([&](Journal& j) { j.ajouter("edo.facteur_integrant", "Facteur intégrant", p.M * x + p.N * y, mu); });
        r.methode = "facteur_integrant";
    }
    journal::noter([&](Journal& j) {
        j.ajouter("edo.exacte", "Forme exacte M dx + N dy = dPhi : Phi(x, y) = C", p.M, Phi);
    });
    const ExprPtr C = p.constante();
    std::vector<ExprPtr> explicites;
    if (isoler(Phi, C, y, explicites)) {
        for (const ExprPtr& e : explicites) ajouterExplicite(r, e);
    } else {
        r.liste.push_back({ExprPtr(), Phi - C, false});
    }
    return true;
}

bool resoudreRiccati(const Probleme& p, const std::vector<TermeEnY>& termes, SolutionsEDO& r) {
    for (const TermeEnY& t : termes) {
        if (!t.exposant.estZero() && !t.exposant.estUn() && t.exposant != Nombre(2)) return false;
    }
    const ExprPtr q1 = coefficientDe(termes, Nombre(1)), q2 = coefficientDe(termes, Nombre(2));
    if (estZero(q2) || estZero(coefficientDe(termes, Nombre(0)))) return false;
    const ExprPtr& x = p.inc.x;
    // Solution particulière y1 = a x^k parmi de petites valeurs rationnelles
    static const int NUMERATEURS[] = {1, -1, 2, -2, 1, -1, 3, -3, 1, -1, 3, -3, 2, -2, 4, -4};
    static const int DENOMINATEURS[] = {1, 1, 1, 1, 2, 2, 1, 1, 3, 3, 2, 2, 3, 3, 1, 1};
    for (int k : {0, 1, -1, 2, -2}) {
        for (std::size_t i = 0; i < sizeof(NUMERATEURS) / sizeof(int); ++i) {
            const ExprPtr y1 = frac(NUMERATEURS[i], DENOMINATEURS[i]) * ast_pow(x, nombre(Nombre(k)));
            const ExprPtr ecart = y1->derivee(x) - substituer(p.f, p.inc.y, y1);
            double valeur;
            try {
                valeur = evaluer(ecart, {{nomDe(x), 1.317}});
            } catch (const std::exception&) {
                return false;
            }
            if (!(std::fabs(valeur) < 1e-9) || !estNul(ecart)) continue;
            // y = y1 + 1/v : v' + (q1 + 2 q2 y1) v = -q2
            const ExprPtr v = lineaire(q1 + 2.0 * q2 * y1, -q2, x, p.constante());
            if (!v) return false;
            journal::noter([&](Journal& j) {
                j.ajouter("edo.riccati", "Équation de Riccati, solution particulière y1 : y = y1 + 1/v, v linéaire", y1, v);
            });
            ajouterExplicite(r, y1 + un() / v);
            ajouterExplicite(r, y1, true);
            r.methode = "riccati";
            return true;
        }
    }
    return false;
}

/*
 * Nom : resoudreClairaut
 * Description : y = x p + g(p) (p = y') : droites y = C x + g(C), et leur enveloppe
 *               x = -g'(p), y = x p + g(p) quand p s'isole dans x = -g'(p).
 */
bool resoudreClairaut(const Probleme& p, SolutionsEDO& r) {
    const ExprPtr& x = p.inc.x;
    const ExprPtr& y = p.inc.y;
    const ExprPtr& dy = p.inc.dy;
    const ExprPtr Ey = p.equation->derivee(y);
    if (estZero(Ey) || contient(Ey, y)) return false; // équation affine en y
    const ExprPtr yDeP = (-substituer(p.equation, y, zero()) / Ey)->simplifier();
    if (!estNul(yDeP->derivee(x) - dy)) return false;
    ExprPtr g;
    if (!sansVariable(yDeP - x * dy, x, g)) return false;
    const ExprPtr C = p.constante();
    journal::noter([&](Journal& j) {
        j.ajouter("edo.clairaut", "Équation de Clairaut y = x y' + g(y') : droites y = C x + g(C) et leur enveloppe", yDeP, g);
    });
    ajouterExplicite(r, C * x + substituer(g, dy, C));
    std::vector<ExprPtr> pentes;
    if (isoler(g->derivee(dy), -x, dy, pentes)) {
        for (const ExprPtr& pente : pentes) ajouterExplicite(r, (x * pente + substituer(g, dy, pente))->simplifier(), true);
    }
    r.methode = "clairaut";
    return true;
}

// Ne garde que les solutions qui vérifient l'équation ; faux si la famille générale n'y est plus
bool filtrer(const ExprPtr& equation, const InconnuesEDO& inc, SolutionsEDO& r) {
    std::vector<SolutionEDO> gardees;
    bool generale = false;
    for (SolutionEDO& s : r.liste) {
        if (!verifierEDO(equation, inc, s)) continue;
        generale = generale || !s.singuliere;
        gardees.push_back(std::move(s));
    }
    r.liste = std::move(gardees);
    return generale;
}

/*
 * Nom : premierOrdre
 * Description : Résout equation(x, y, dy) = 0. constantes est le numéro de la dernière
 *               constante déjà utilisée (les nouvelles s'appellent C(constantes + 1)...).
 */
SolutionsEDO premierOrdre(const ExprPtr& equation, const InconnuesEDO& inc, int& constantes) {
    Probleme p{inc, equation, ExprPtr(), ExprPtr(), ExprPtr(), &constantes};
    const int depart = constantes;
    SolutionsEDO r;
    const auto essayer = [&](bool reussie) {
        if (reussie && filtrer(equation, inc, r)) return true;
        r = SolutionsEDO();
        constantes = depart;
        return false;
    };
    const ExprPtr A = equation->derivee(inc.dy);
    if (estZero(A)) return r; // pas d'équation différentielle
    if (contient(A, inc.dy)) {
        essayer(resoudreClairaut(p, r));
        return r;
    }
    // A y' + B = 0
    const ExprPtr B = substituer(equation, inc.dy, zero());
    p.f = (-B / A)->simplifier();
    p.M = B;
    p.N = A;
    if (!contient(A, inc.x) && !contient(A, inc.y)) {
        // y' = f donné sous forme résolue : M = -numérateur, N = dénominateur de f
        AccumulateurProduit denominateur;
        if (const Produit* produit = comme<Produit>(p.f)) {
            for (const Facteur& fa : produit->getFacteurs()) {
                const Constante* e = comme<Constante>(fa.exposant);
                if (e && e->getNombre().signe() < 0) denominateur.multiplier(ast_pow(fa.base, -fa.exposant));
            }
        } else if (const Puissance* puissance = comme<Puissance>(p.f)) {
            const Constante* e = comme<Constante>(puissance->getExposant());
            if (e && e->getNombre().signe() < 0) denominateur.multiplier(ast_pow(puissance->getBase(), -puissance->getExposant()));
        }
        p.N = denominateur.construire();
        p.M = -(p.f * p.N);
    }
    journal::noter([&](Journal& j) { j.ajouter("edo.forme", "Forme résolue y' = f(x, y)", equation, p.f); });

    std::vector<TermeEnY> termes;
    const bool polynomiale = decomposer(p.f, inc.y, termes);
    if (essayer(resoudreQuadrature(p, r))) return r;
    if (polynomiale && essayer(resoudreLineaire(p, termes, r))) return r;
    if (polynomiale && essayer(resoudreBernoulli(p, termes, r))) return r;
    if (essayer(resoudreSeparable(p, r))) return r;
    // Équation donnée sous la forme M + N y' = 0 : la forme exacte d'abord, elle garde M et N
    const bool formeDifferentielle = contient(A, inc.x) || contient(A, inc.y);
    if (formeDifferentielle && essayer(resoudreExacte(p, r))) return r;
    if (essayer(resoudreAffine(p, r))) return r;
    if (essayer(resoudreHomogene(p, r))) return r;
    if (!formeDifferentielle && essayer(resoudreExacte(p, r))) return r;
    if (polynomiale && essayer(resoudreRiccati(p, termes, r))) return r;
    return r;
}

// ============================================================================
// Second ordre réductible
// ============================================================================

SolutionsEDO secondOrdre(const ExprPtr& equation, const InconnuesEDO& inc) {
    SolutionsEDO r;
    const ExprPtr A = equation->derivee(inc.d2y);
    if (contient(A, inc.d2y)) return r; // non affine en y''
    const ExprPtr f = (-substituer(equation, inc.d2y, zero()) / A)->simplifier(); // y'' = f(x, y, y')
    int constantes = 0;
    const ExprPtr dp = var("d" + nomDe(inc.dy));
    if (!contient(f, inc.y)) {
        // y absente : équation du premier ordre en p(x) = y'
        const InconnuesEDO enP{inc.x, inc.dy, dp, var("d2" + nomDe(inc.dy))};
        const SolutionsEDO pentes = premierOrdre(dp - f, enP, constantes);
        journal::noter([&](Journal& j) {
            j.ajouter("edo.ordre2_sans_y", "y n'apparaît pas : équation du premier ordre en p = y', puis y = int p dx", equation, f);
        });
        const ExprPtr C = param("C" + std::to_string(++constantes));
        for (const SolutionEDO& s : pentes.liste) {
            if (!s.y) continue;
            const ExprPtr Y = integrerPar(s.y, inc.x);
            if (Y) r.liste.push_back({Y + C, ExprPtr(), s.singuliere});
        }
        r.methode = "ordre2_sans_y";
    } else if (!contient(f, inc.x)) {
        // x absente : y'' = p dp/dy, équation du premier ordre en p(y)
        const InconnuesEDO enP{inc.y, inc.dy, dp, var("d2" + nomDe(inc.dy))};
        const SolutionsEDO pentes = premierOrdre(inc.dy * dp - f, enP, constantes);
        journal::noter([&](Journal& j) {
            j.ajouter("edo.ordre2_sans_x", "x n'apparaît pas : y'' = p dp/dy, équation du premier ordre en p(y), puis int dy/p = x + C",
                      equation, f);
        });
        const ExprPtr C = param("C" + std::to_string(++constantes)); // la même pour les deux signes de p
        for (const SolutionEDO& s : pentes.liste) {
            if (!s.y) continue;
            if (estZero(s.y)) continue; // p = 0 : y constante, contenue dans les familles ou sans intérêt
            const ExprPtr H = integrerPar(un() / s.y, inc.y); // int dy/p(y) = x + C
            if (!H) continue;
            std::vector<ExprPtr> explicites;
            if (expliciter(H, inc.x, inc.y, C, explicites)) {
                for (const ExprPtr& e : explicites) r.liste.push_back({e, ExprPtr(), s.singuliere});
            } else {
                r.liste.push_back({ExprPtr(), H - inc.x - C, s.singuliere});
            }
        }
        r.methode = "ordre2_sans_x";
    }
    if (!filtrer(equation, inc, r)) r = SolutionsEDO();
    return r;
}

SolutionsEDO resoudreEquation(const ExprPtr& equation, const InconnuesEDO& inc) {
    if (contient(equation, inc.d2y)) return secondOrdre(equation, inc);
    int constantes = 0;
    return premierOrdre(equation, inc, constantes);
}

// y'' et y' deviennent les variables d2y et dy dans le texte
std::string remplacerPrimes(const std::string& texte, const std::string& fonction) {
    std::string r;
    for (std::size_t i = 0; i < texte.size();) {
        const bool debutDeMot = i == 0 || !(std::isalnum(static_cast<unsigned char>(texte[i - 1])) || texte[i - 1] == '_');
        if (debutDeMot && texte.compare(i, fonction.size(), fonction) == 0) {
            const std::size_t fin = i + fonction.size();
            if (texte.compare(fin, 2, "''") == 0) {
                r += "d2" + fonction;
                i = fin + 2;
                continue;
            }
            if (texte.compare(fin, 1, "'") == 0) {
                r += "d" + fonction;
                i = fin + 1;
                continue;
            }
        }
        r += texte[i++];
    }
    return r;
}

} // namespace

// ============================================================================
// API publique
// ============================================================================

std::string SolutionEDO::texte() const {
    if (y) return "y = " + y->texte();
    return relation ? relation->texte() + " = 0" : std::string();
}

bool verifierEDO(const ExprPtr& equation, const InconnuesEDO& inconnues, const SolutionEDO& solution) {
    if (!solution.y && !solution.relation) return false;
    try {
        return estNul(residu(equation, inconnues, solution), false);
    } catch (const std::exception&) {
        return false;
    }
}

SolutionsEDO resoudreEDO(const ExprPtr& equation, const InconnuesEDO& inconnues) {
    try {
        return resoudreEquation(equation, inconnues);
    } catch (const std::exception&) { // équation hors de portée : jamais d'exception pour un échec
        return SolutionsEDO();
    }
}

SolutionsEDO resoudreEDO(const ExprPtr& equation, const InconnuesEDO& inconnues, Journal& journal) {
    const journal::Activation activation(journal);
    return resoudreEDO(equation, inconnues);
}

SolutionsEDO resoudreEDO(const std::string& texte, const std::string& variable, const std::string& fonction) {
    const InconnuesEDO inc{var(variable), var(fonction), var("d" + fonction), var("d2" + fonction)};
    const EgaliteLue egalite =
        lireEquation(remplacerPrimes(texte, fonction), OptionsLecture{variable, {fonction, "d" + fonction, "d2" + fonction}});
    return resoudreEDO(egalite.gauche - egalite.droite, inc);
}

SolutionsEDO resoudreEDO(const std::string& texte, Journal& journal, const std::string& variable, const std::string& fonction) {
    const journal::Activation activation(journal);
    return resoudreEDO(texte, variable, fonction);
}

} // namespace symalgo
