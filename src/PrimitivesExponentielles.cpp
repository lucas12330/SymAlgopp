/**
 * @file PrimitivesExponentielles.cpp
 * @brief Intégration des expressions à exponentielles : classe exponentielle-polynôme,
 *        fractions rationnelles d'exponentielles, équation différentielle de Risch.
 *
 * 1. Classe exponentielle-polynôme (parExponentiellePolynome). Tout produit de polynômes en
 *    x, de exp, sin, cos, sinh, cosh d'arguments affines s'écrit comme une somme de monômes
 *
 *        c x^n H(W) T(V),   H dans {1, exp, cosh, sinh},  T dans {1, cos, sin},
 *
 *    W et V affines (les formules de produit en somme sont réelles : aucune exponentielle
 *    complexe n'apparaît dans les expressions). Avec w = W', b = V' et lambda = w + i b,
 *
 *        int x^n e^(lambda x) dx = e^(lambda x) sum_k (-1)^k n!/(n-k)! x^(n-k) / lambda^(k+1),
 *
 *    dont on prend la partie réelle ou imaginaire ; 1/lambda^(k+1) = (w - i b)^(k+1)/(w^2 + b^2)^(k+1)
 *    est calculé de proche en proche. Pour cosh et sinh, les deux exponentielles e^W et e^-W
 *    se regroupent : le résultat reste écrit en cosh et sinh.
 *
 * 2. Fractions rationnelles de exp(k u), sinh, cosh, tanh (parExponentielleRationnelle) :
 *    t = exp(u), dx = dt/(a t), puis intégration de la fraction rationnelle en t.
 *
 * 3. Équation de Risch (parEquationDeRisch) : int R exp(W) = S exp(W) avec S' + W' S = R.
 *    Les pôles de S sont ceux de R avec un ordre de moins, et son numérateur a un degré
 *    connu : S se trouve par un système triangulaire, ou n'existe pas (exp(x^2), exp(x)/x).
 */

#include "Canonique.hpp"
#include "Journal.hpp"
#include "Primitives.hpp"

#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace symalgo {

namespace {

constexpr int DEGRE_MAX = 64;              // degré en x d'un monôme
constexpr long long PUISSANCE_MAX = 32;    // exposant d'un sinus, d'un cosinus...
constexpr std::size_t MONOMES_MAX = 256;
constexpr int DEGRE_FRACTION_MAX = 128;

// ============================================================================
// 1. Classe exponentielle-polynôme
// ============================================================================

enum class Genre : std::uint8_t { Aucun, Exp, Pair, Impair }; // Pair : cosh ou cos ; Impair : sinh ou sin

struct Fonction {
    Genre genre = Genre::Aucun;
    ExprPtr argument, pente; // argument affine et sa dérivée
};

struct Monome {
    ExprPtr coefficient; // sans x
    int degre = 0;
    Fonction hyp, trig;
};

using Monomes = std::vector<Monome>;

// Signe conventionnel d'une expression sans x : exactement l'une de a et -a est « négative »
bool negatif(const ExprPtr& a) {
    switch (a->type()) {
        case TypeNoeud::Constante: return static_cast<const Constante&>(*a).getNombre().signe() < 0;
        case TypeNoeud::Produit: return static_cast<const Produit&>(*a).getCoefficient().signe() < 0;
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(*a);
            return s.getTermes().empty() ? s.getConstante().signe() < 0 : s.getTermes()[0].coefficient.signe() < 0;
        }
        default: return false;
    }
}

bool estZero(const ExprPtr& e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().estZero();
}

ExprPtr hyperbolique(Genre g, const ExprPtr& w) {
    return g == Genre::Exp ? ast_exp(w) : g == Genre::Pair ? ast_cosh(w) : ast_sinh(w);
}

ExprPtr trigonometrique(Genre g, const ExprPtr& v) { return g == Genre::Pair ? ast_cos(v) : ast_sin(v); }

/*
 * Nom : normaliser
 * Description : Forme unique d'un monôme : une fonction d'argument constant passe dans le
 *               coefficient, la pente d'un cosh, sinh, cos ou sin est rendue positive.
 *               Renvoie faux si le monôme est nul.
 */
bool normaliser(Monome& m) {
    for (Fonction* f : {&m.hyp, &m.trig}) {
        if (f->genre == Genre::Aucun) continue;
        if (!f->argument->contientVariable()) {
            m.coefficient = m.coefficient * (f == &m.hyp ? hyperbolique(f->genre, f->argument) : trigonometrique(f->genre, f->argument));
            *f = Fonction();
        } else if (f->genre != Genre::Exp && negatif(f->pente)) {
            f->argument = -f->argument;
            f->pente = -f->pente;
            if (f->genre == Genre::Impair) m.coefficient = -m.coefficient;
        }
    }
    return !estZero(m.coefficient);
}

bool memeFonction(const Fonction& a, const Fonction& b) {
    return a.genre == b.genre && a.argument.get() == b.argument.get();
}

void ajouter(Monomes& somme, Monome m) {
    if (!normaliser(m)) return;
    for (Monome& s : somme) {
        if (s.degre == m.degre && memeFonction(s.hyp, m.hyp) && memeFonction(s.trig, m.trig)) {
            s.coefficient = s.coefficient + m.coefficient;
            return;
        }
    }
    somme.push_back(std::move(m));
}

struct Morceau {
    Nombre facteur;
    Fonction fonction;
};

/*
 * Nom : produitFonctions
 * Description : a * b écrit comme somme d'au plus deux fonctions de la même famille
 *               (formules de produit en somme). epsilon = 1 pour la famille hyperbolique
 *               (sinh sinh = (cosh(+) - cosh(-))/2), -1 pour la famille trigonométrique
 *               (sin sin = (cos(-) - cos(+))/2).
 */
int produitFonctions(const Fonction& a, const Fonction& b, int epsilon, Morceau m[2]) {
    if (a.genre == Genre::Aucun || b.genre == Genre::Aucun) {
        m[0] = {Nombre(1), a.genre == Genre::Aucun ? b : a};
        return 1;
    }
    const Fonction plus{Genre::Aucun, a.argument + b.argument, a.pente + b.pente};
    if (a.genre == Genre::Exp && b.genre == Genre::Exp) {
        m[0] = {Nombre(1), {Genre::Exp, plus.argument, plus.pente}};
        return 1;
    }
    const Nombre demi = Nombre::rationnel(1, 2);
    if (a.genre == Genre::Exp || b.genre == Genre::Exp) {
        // e^E cosh(W) = (e^(E + W) + e^(E - W))/2, e^E sinh(W) = (e^(E + W) - e^(E - W))/2
        const Fonction& e = a.genre == Genre::Exp ? a : b;
        const Fonction& h = a.genre == Genre::Exp ? b : a;
        m[0] = {demi, {Genre::Exp, plus.argument, plus.pente}};
        m[1] = {h.genre == Genre::Pair ? demi : -demi, {Genre::Exp, e.argument - h.argument, e.pente - h.pente}};
        return 2;
    }
    const Fonction moins{Genre::Aucun, a.argument - b.argument, a.pente - b.pente};
    if (a.genre == b.genre) {
        // pair pair = (pair(+) + pair(-))/2 ; impair impair = epsilon (pair(+) - pair(-))/2
        const bool pair = a.genre == Genre::Pair;
        m[0] = {pair || epsilon > 0 ? demi : -demi, {Genre::Pair, plus.argument, plus.pente}};
        m[1] = {pair || epsilon < 0 ? demi : -demi, {Genre::Pair, moins.argument, moins.pente}};
        return 2;
    }
    // impair(a) pair(b) = (impair(a + b) + impair(a - b))/2 ; pair(a) impair(b) = (impair(a + b) - impair(a - b))/2
    m[0] = {demi, {Genre::Impair, plus.argument, plus.pente}};
    m[1] = {a.genre == Genre::Impair ? demi : -demi, {Genre::Impair, moins.argument, moins.pente}};
    return 2;
}

bool multiplier(Monomes& gauche, const Monomes& droite) {
    Monomes resultat;
    Morceau h[2], t[2];
    for (const Monome& a : gauche) {
        for (const Monome& b : droite) {
            const int degre = a.degre + b.degre;
            if (degre > DEGRE_MAX) return false;
            const ExprPtr coefficient = a.coefficient * b.coefficient;
            const int nh = produitFonctions(a.hyp, b.hyp, 1, h);
            const int nt = produitFonctions(a.trig, b.trig, -1, t);
            for (int i = 0; i < nh; ++i) {
                for (int j = 0; j < nt; ++j) {
                    ajouter(resultat, {coefficient * nombre(h[i].facteur * t[j].facteur), degre, h[i].fonction, t[j].fonction});
                }
            }
            if (resultat.size() > MONOMES_MAX) return false;
        }
    }
    gauche = std::move(resultat);
    return !gauche.empty();
}

bool multiplierPar(Monomes& gauche, Monome m) {
    const Monomes droite{std::move(m)};
    return multiplier(gauche, droite);
}

struct Lecture {
    const ASTNode* x;
    Monomes produit;
    bool transcendant = false;
    std::vector<std::pair<ExprPtr, ExprPtr>> puissances; // exp(u ln c) -> c^u, pour l'affichage
};

// produit *= base^exposant ; faux si le facteur n'est pas de la classe
bool lireFacteur(const ExprPtr& base, const ExprPtr& exposant, Lecture& l) {
    ExprPtr a;
    if (!base->contientVariable()) {
        if (!exposant->contientVariable()) return multiplierPar(l.produit, {ast_pow(base, exposant), 0, {}, {}});
        // c^(a x + b), c > 0 : exp((a x + b) ln c)
        if (!estPositifOuNul(base) || !coefficientLineaire(exposant, a)) return false;
        const ExprPtr lnc = ast_ln(base);
        const Fonction e{Genre::Exp, exposant * lnc, a * lnc};
        l.puissances.emplace_back(ast_exp(e.argument), ast_pow(base, exposant));
        l.transcendant = true;
        return multiplierPar(l.produit, {un(), 0, e, {}});
    }
    const Constante* c = comme<Constante>(exposant);
    if (!c || !c->getNombre().estExact()) return false;
    const TypeNoeud type = base->type();
    if (type == TypeNoeud::Exponentielle) { // exp(u)^r = exp(r u), r rationnel quelconque
        const ExprPtr& u = static_cast<const Exponentielle&>(*base).m_argument;
        if (!coefficientLineaire(u, a)) return false;
        l.transcendant = true;
        return multiplierPar(l.produit, {un(), 0, {Genre::Exp, u * exposant, a * exposant}, {}});
    }
    long long n;
    if (!c->getNombre().versEntier(n) || n < 1) return false;
    switch (type) {
        case TypeNoeud::Variable:
            return n <= DEGRE_MAX && multiplierPar(l.produit, {un(), static_cast<int>(n), {}, {}});
        case TypeNoeud::Sinus:
        case TypeNoeud::Cosinus:
        case TypeNoeud::SinusHyperbolique:
        case TypeNoeud::CosinusHyperbolique: {
            const ExprPtr& u = static_cast<const FonctionUnaire&>(*base).m_argument;
            if (n > PUISSANCE_MAX || !coefficientLineaire(u, a)) return false;
            const Genre genre = type == TypeNoeud::Cosinus || type == TypeNoeud::CosinusHyperbolique ? Genre::Pair : Genre::Impair;
            Monome m{un(), 0, {}, {}};
            (type == TypeNoeud::Sinus || type == TypeNoeud::Cosinus ? m.trig : m.hyp) = {genre, u, a};
            l.transcendant = true;
            const Monomes facteur{m};
            for (long long i = 0; i < n; ++i) {
                if (!multiplier(l.produit, facteur)) return false;
            }
            return true;
        }
        case TypeNoeud::Somme: { // polynôme en x à coefficients exacts
            Polynome p;
            ExprPtr v;
            if (!Polynome::depuisExpression(base, p, &v) || !p.estExact() || v.get() != l.x || p.degre() * n > DEGRE_MAX) return false;
            Polynome puissance({Nombre(1)});
            for (long long i = 0; i < n; ++i) puissance = puissance * p;
            Monomes facteur;
            for (int i = 0; i <= puissance.degre(); ++i) {
                if (!puissance.coefficient(i).estZero()) facteur.push_back({nombre(puissance.coefficient(i)), i, {}, {}});
            }
            return multiplier(l.produit, facteur);
        }
        default: return false;
    }
}

ExprPtr versExpression(const Monome& m, const ExprPtr& x) {
    ExprPtr e = m.coefficient * ast_pow(x, nombre(Nombre(m.degre)));
    if (m.hyp.genre != Genre::Aucun) e = e * hyperbolique(m.hyp.genre, m.hyp.argument);
    if (m.trig.genre != Genre::Aucun) e = e * trigonometrique(m.trig.genre, m.trig.argument);
    return e;
}

// Primitives accumulées, regroupées par exponentielle en facteur : (x - 1) exp(x)
struct Groupes {
    std::vector<std::pair<ExprPtr, AccumulateurSomme>> liste;

    AccumulateurSomme& pour(const ExprPtr& exponentielle) {
        for (auto& g : liste) {
            if (g.first.get() == exponentielle.get()) return g.second;
        }
        liste.emplace_back(exponentielle, AccumulateurSomme());
        return liste.back().second;
    }

    ExprPtr construire() {
        AccumulateurSomme total;
        for (auto& g : liste) {
            const ExprPtr s = g.second.construire();
            total.ajouter(g.first ? g.first * s : s, Nombre(1));
        }
        return total.construire();
    }
};

/*
 * Nom : integrerMonome
 * Description : Primitive de c x^n H(W) T(V) (voir l'en-tête du fichier). (p, q) parcourt
 *               les puissances (w - i b)^(k+1) ; pour cosh et sinh, le terme d'indice k porte
 *               sinh ou cosh selon la parité de k.
 */
void integrerMonome(const Monome& m, const ExprPtr& x, Groupes& groupes) {
    const int n = m.degre;
    const Genre hyp = m.hyp.genre, trig = m.trig.genre;
    if (hyp == Genre::Aucun && trig == Genre::Aucun) {
        groupes.pour(ExprPtr()).ajouter(m.coefficient * ast_pow(x, nombre(Nombre(n + 1))), Nombre::rationnel(1, n + 1));
        return;
    }
    const ExprPtr nul = nombre(Nombre(0));
    const ExprPtr& w = hyp == Genre::Aucun ? nul : m.hyp.pente;
    const ExprPtr& b = trig == Genre::Aucun ? nul : m.trig.pente;
    const ExprPtr module = ast_pow(w, 2.0) + ast_pow(b, 2.0);
    const bool exponentielle = hyp == Genre::Exp;
    AccumulateurSomme& cible = groupes.pour(exponentielle ? ast_exp(m.hyp.argument) : ExprPtr());
    ExprPtr ch, sh, cosV, sinV;
    if (hyp == Genre::Pair || hyp == Genre::Impair) {
        ch = ast_cosh(m.hyp.argument);
        sh = ast_sinh(m.hyp.argument);
    }
    if (trig != Genre::Aucun) {
        cosV = ast_cos(m.trig.argument);
        sinV = ast_sin(m.trig.argument);
    }
    ExprPtr p = w, q = -b;
    Nombre c(1); // (-1)^k n!/(n-k)!
    for (int k = 0; k <= n; ++k) {
        ExprPtr A = p, B = q;
        if (ch) {
            const bool sinhDevant = (k % 2 == 0) == (hyp == Genre::Pair);
            A = p * (sinhDevant ? sh : ch);
            B = q * (sinhDevant ? ch : sh);
        }
        const ExprPtr valeur = trig == Genre::Aucun ? A : trig == Genre::Pair ? A * cosV - B * sinV : A * sinV + B * cosV;
        cible.ajouter(m.coefficient * valeur * ast_pow(x, nombre(Nombre(n - k))) / ast_pow(module, nombre(Nombre(k + 1))), c);
        const ExprPtr suivant = p * w + q * b;
        q = q * w - p * b;
        p = suivant;
        c = c * Nombre(-(n - k));
    }
}

} // namespace

ExprPtr parExponentiellePolynome(const ExprPtr& f) {
    Lecture l;
    const ExprPtr x = variableDe(*f);
    if (!x) return nullptr;
    l.x = x.get();
    l.produit.push_back({un(), 0, {}, {}});
    if (const Produit* p = comme<Produit>(f)) {
        if (!p->getCoefficient().estUn()) l.produit[0].coefficient = nombre(p->getCoefficient());
        for (const Facteur& facteur : p->getFacteurs()) {
            if (!lireFacteur(facteur.base, facteur.exposant, l)) return nullptr;
        }
    } else if (const Puissance* p = comme<Puissance>(f)) {
        if (!lireFacteur(p->getBase(), p->getExposant(), l)) return nullptr;
    } else {
        return nullptr;
    }
    if (!l.transcendant) return nullptr;
    if (l.produit.size() > 1) {
        journal::noter([&](Journal& j) {
            AccumulateurSomme s;
            for (const Monome& m : l.produit) s.ajouter(versExpression(m, x), Nombre(1));
            j.ajouter("integration.linearisation", "Linéarisation : produits de fonctions écrits en somme", f, s.construire());
        });
    }
    Groupes groupes;
    for (const Monome& m : l.produit) integrerMonome(m, x, groupes);
    ExprPtr resultat = groupes.construire();
    for (const auto& [exponentielle, puissance] : l.puissances) resultat = substituer(resultat, exponentielle, puissance);
    return resultat;
}

// ============================================================================
// 2. Fractions rationnelles d'exponentielles et de fonctions hyperboliques
// ============================================================================

bool Multiples::ajouter(const ExprPtr& argument) {
    for (const ASTNode* a : arguments) {
        if (a == argument.get()) return true;
    }
    ExprPtr a;
    if (!coefficientLineaire(argument, a)) return false;
    if (arguments.empty()) m_premier = argument;
    arguments.push_back(argument.get());
    return arguments.size() <= 16;
}

bool Multiples::resoudre(long long multipleMax) {
    if (arguments.empty()) return false;
    std::vector<std::pair<long long, long long>> rapports;
    long long denominateur = 1;
    for (const ASTNode* a : arguments) {
        const ExprPtr rapport = a->clone() / m_premier;
        const Constante* r = comme<Constante>(rapport);
        long long num, den;
        if (!r || !r->getNombre().versFraction(num, den) || den > 64) return false;
        rapports.emplace_back(num, den);
        denominateur = std::lcm(denominateur, den);
        if (denominateur > 64) return false;
    }
    long long diviseur = 0;
    for (const auto& [num, den] : rapports) {
        multiples.push_back(num * (denominateur / den));
        diviseur = std::gcd(diviseur, multiples.back());
    }
    for (long long& k : multiples) {
        k /= diviseur;
        if (k > multipleMax || k < -multipleMax) return false;
    }
    u = m_premier * nombre(Nombre::rationnel(diviseur, denominateur));
    return coefficientLineaire(u, pente);
}

long long Multiples::multipleDe(const ASTNode* argument) const {
    for (std::size_t i = 0; i < arguments.size(); ++i) {
        if (arguments[i] == argument) return multiples[i];
    }
    return 0;
}

namespace {

struct ReglesExponentielles {
    const Multiples& m;

    Feuille feuille(const ASTNode& e, ElementFraction& r) const {
        const TypeNoeud type = e.type();
        if (type != TypeNoeud::Exponentielle && type != TypeNoeud::SinusHyperbolique && type != TypeNoeud::CosinusHyperbolique &&
            type != TypeNoeud::TangenteHyperbolique) {
            return Feuille::Descendre;
        }
        const long long k = m.multipleDe(static_cast<const FonctionUnaire&>(e).m_argument.get());
        const Fraction plus = Fraction::variable().puissance(k), moins = Fraction::variable().puissance(-k);
        const Nombre demi = Nombre::rationnel(1, 2);
        switch (type) {
            case TypeNoeud::Exponentielle: r.f = plus; break;
            case TypeNoeud::CosinusHyperbolique: r.f = (plus + moins) * demi; break;
            case TypeNoeud::SinusHyperbolique: r.f = (plus - moins) * demi; break;
            default: r.f = (plus - moins) / (plus + moins); break;
        }
        return Feuille::Convertie;
    }

    bool puissance(const ASTNode&, const ASTNode&, ElementFraction&) const { return false; }
};

} // namespace

ExprPtr parExponentielleRationnelle(const ExprPtr& f) {
    const TypeNoeud type = f->type();
    if (type != TypeNoeud::Produit && type != TypeNoeud::Puissance) return nullptr;
    Multiples m;
    const auto accepte = [](TypeNoeud t) {
        return t == TypeNoeud::Exponentielle || t == TypeNoeud::SinusHyperbolique || t == TypeNoeud::CosinusHyperbolique ||
               t == TypeNoeud::TangenteHyperbolique;
    };
    if (!collecter(*f, m, accepte) || !m.resoudre(DEGRE_FRACTION_MAX / 4)) return nullptr;
    ElementFraction h;
    const ReglesExponentielles regles{m};
    try {
        if (!convertirRationnel(*f, h, regles, ElementFraction{Fraction(Nombre(1))}, DEGRE_FRACTION_MAX)) return nullptr;
    } catch (const std::domain_error&) {
        return nullptr;
    }
    const ExprPtr t = ast_exp(m.u);
    const Fraction integrande = h.f / Fraction::variable(); // dx = dt/(a t)
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    journal::noter([&](Journal& jr) {
        jr.ajouter("integration.substitution", "Changement de variable t = " + t->texte() + ", dt = " + (m.pente * t)->texte() + " dx", f,
                   integrande.versExpression(var("t")) / m.pente);
    });
    bool complete;
    ExprPtr F;
    {
        const journal::SousNiveau sous;
        F = primitiveFraction(integrande, t, &complete);
    }
    if (!complete) {
        if (j) j->tronquer(etapes);
        return nullptr;
    }
    return F / m.pente;
}

// ============================================================================
// 3. Équation de Risch : S' + W' S = R
// ============================================================================

namespace {

bool memePolynome(const Polynome& a, const Polynome& b) { return (a - b).estNul(); }

/*
 * Nom : resoudreRisch
 * Description : Solution rationnelle S de S' + w S = g, w polynôme non nul. Avec g = n/d,
 *               D = pgcd(d, d') et E = d/D, S = T/D où le polynôme T vérifie
 *               d T' + (w d - E D') T = n D ; comme deg(w d - E D') >= deg d, le degré de T
 *               est deg(n D) - deg(w d - E D') et ses coefficients se calculent du plus haut
 *               au plus bas. Renvoie faux s'il n'y a pas de solution.
 */
bool resoudreRisch(const Polynome& w, const Fraction& g, Fraction& S) {
    if (g.estNulle()) {
        S = Fraction();
        return true;
    }
    const Polynome& n = g.numerateur();
    const Polynome& d = g.denominateur();
    const Polynome D = Polynome::pgcd(d, d.derivee());
    Polynome E, reste;
    d.diviser(D, E, reste);
    const Polynome b = w * d - E * D.derivee();
    Polynome c = n * D;
    const int degre = c.degre() - b.degre();
    if (degre < 0) return false;
    std::vector<Nombre> t(static_cast<std::size_t>(degre + 1), Nombre(0));
    for (int k = degre; k >= 0 && !c.estNul(); --k) {
        if (c.degre() > k + b.degre()) return false;
        if (c.degre() < k + b.degre()) continue;
        const Nombre tk = c.dominant() / b.dominant();
        t[static_cast<std::size_t>(k)] = tk;
        std::vector<Nombre> monome(static_cast<std::size_t>(k + 1), Nombre(0));
        monome.back() = tk;
        const Polynome m(std::move(monome));
        c = c - d * m.derivee() - b * m;
    }
    if (!c.estNul()) return false;
    S = Fraction(Polynome(std::move(t)), D);
    return true;
}

struct TermeExponentiel {
    Polynome W;
    ExprPtr exponentielle; // exp(W) tel qu'il sera écrit
    Fraction R;
};

// terme = coefficient * R(x) exp(W(x)) : l'ajoute à termes (regroupés par W)
bool lireTermeExponentiel(const ExprPtr& terme, const Nombre& coefficient, const ExprPtr& x, std::vector<TermeExponentiel>& termes) {
    AccumulateurProduit rationnel;
    AccumulateurSomme argument;
    bool exponentielle = false;
    auto facteur = [&](const ExprPtr& base, const ExprPtr& exposant) {
        if (base->type() != TypeNoeud::Exponentielle) {
            rationnel.multiplier(ast_pow(base, exposant));
            return true;
        }
        const Constante* k = comme<Constante>(exposant);
        if (!k || !k->getNombre().estExact()) return false;
        argument.ajouter(static_cast<const Exponentielle&>(*base).m_argument, k->getNombre());
        exponentielle = true;
        return true;
    };
    if (const Produit* p = comme<Produit>(terme)) {
        rationnel.multiplier(nombre(p->getCoefficient()));
        for (const Facteur& f : p->getFacteurs()) {
            if (!facteur(f.base, f.exposant)) return false;
        }
    } else if (const Puissance* p = comme<Puissance>(terme)) {
        if (!facteur(p->getBase(), p->getExposant())) return false;
    } else if (!facteur(terme, un())) {
        return false;
    }
    if (!exponentielle) return false;
    const ExprPtr W = argument.construire();
    TermeExponentiel t;
    ExprPtr v;
    if (!Polynome::depuisExpression(W, t.W, &v) || !t.W.estExact() || t.W.degre() < 1 || v.get() != x.get()) return false;
    if (!versFraction(rationnel.construire(), x, t.R)) return false;
    t.R = t.R * coefficient;
    for (TermeExponentiel& autre : termes) {
        if (memePolynome(autre.W, t.W)) {
            autre.R = autre.R + t.R;
            return true;
        }
    }
    t.exponentielle = ast_exp(W);
    termes.push_back(std::move(t));
    return true;
}

} // namespace

ExprPtr parEquationDeRisch(const ExprPtr& f) {
    const ExprPtr x = variableDe(*f);
    if (!x) return nullptr;
    std::vector<TermeExponentiel> termes;
    if (const Somme* s = comme<Somme>(f)) {
        if (!s->getConstante().estZero()) return nullptr;
        for (const Terme& t : s->getTermes()) {
            if (!lireTermeExponentiel(t.expression, t.coefficient, x, termes)) return nullptr;
        }
    } else if (!lireTermeExponentiel(f, Nombre(1), x, termes)) {
        return nullptr;
    }
    AccumulateurSomme resultat;
    for (const TermeExponentiel& t : termes) {
        Fraction S;
        if (!resoudreRisch(t.W.derivee(), t.R, S)) return nullptr;
        resultat.ajouter(S.versExpression(x) * t.exponentielle, Nombre(1));
    }
    return resultat.construire();
}

ExprPtr secondeChanceSomme(const ExprPtr& somme) { return parEquationDeRisch(somme); }

} // namespace symalgo
