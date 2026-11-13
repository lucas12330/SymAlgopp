/**
 * @file Primitives.cpp
 * @brief Intégration symbolique : point d'entrée ASTNode::integrer et suite de méthodes
 *        essayées dans l'ordre, de la moins chère à la plus chère (tableau METHODES).
 *
 *    1. règle du noeud (Regles.cpp, méthode protégée primitive()) : linéarité, facteur
 *       constant, primitives usuelles d'un argument affine F(ax+b)/a ;
 *    2. table étendue : 1/(a u^2 + b u + c) -> atan, 1/sqrt(...) -> asin ou ln, 1/cos^2 -> tan... ;
 *    3. polynômes en sin u et cos u impairs en l'un des deux (PrimitivesTrigonometriques.cpp) ;
 *    4. classe exponentielle-polynôme (PrimitivesExponentielles.cpp) ;
 *    5. fractions rationnelles : Hermite, éléments simples (PrimitivesRationnelles.cpp) ;
 *    6. fractions rationnelles de sin et cos : règles de Bioche, tan(u/2) ;
 *    7. fractions rationnelles d'exponentielles : t = exp(u) ;
 *    8. changement de variable « la dérivée divise » : u = g(x) pris parmi les
 *       sous-expressions de f ; si f/g' ne s'écrit qu'en fonction de u, on intègre en u ;
 *    9. intégration par parties : ln(v)^n, atan, asin, acos fois un cofacteur intégrable,
 *       puis u = ln(x) avec x = exp(u) (sin(ln(x))) ;
 *   10. radicaux : (a x + b)^(p/q), sqrt(a x^2 + b x + c) (PrimitivesAlgebriques.cpp) ;
 *   11. équation de Risch S' + W' S = R pour R(x) exp(W(x)) ;
 *   12. développement d'un produit en somme.
 *
 * Chaque méthode renvoie un pointeur nul si elle ne s'applique pas : elle rejette vite, sur
 * le type du noeud puis sur la forme. Une méthode qui en appelle d'autres (substitution,
 * parties) n'accepte qu'un résultat complet, sans intégrale non évaluée, et sa profondeur
 * d'imbrication est bornée.
 *
 * Coût : la règle du noeud est essayée d'abord, seule ; les autres méthodes ne servent
 * qu'aux intégrales qu'elle ne sait pas faire. Rien n'est ajouté à la construction
 * canonique, à la dérivée ni à l'évaluation.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"
#include "Journal.hpp"
#include "Polynome.hpp"
#include "Primitives.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace symalgo {

// ============================================================================
// Outils partagés avec Regles.cpp
// ============================================================================

bool coefficientLineaire(const ExprPtr& u, ExprPtr& a) {
    if (!u->contientVariable()) return false;
    const ExprPtr d = u->derivee();
    if (d->contientVariable()) return false;
    const Constante* c = comme<Constante>(d);
    if (c && (c->getNombre().estZero() || !c->getNombre().estFini())) return false;
    a = d;
    return true;
}

ExprPtr variableDe(const ASTNode& e) {
    const ASTNode* n = &e;
    while (n && n->contientVariable()) {
        const ASTNode* suivant = nullptr;
        switch (n->type()) {
            case TypeNoeud::Variable: return n->clone();
            case TypeNoeud::Somme:
                for (const Terme& t : static_cast<const Somme&>(*n).getTermes()) {
                    if (t.expression->contientVariable()) {
                        suivant = t.expression.get();
                        break;
                    }
                }
                break;
            case TypeNoeud::Produit:
                for (const Facteur& f : static_cast<const Produit&>(*n).getFacteurs()) {
                    if (f.base->contientVariable() || f.exposant->contientVariable()) {
                        suivant = f.base->contientVariable() ? f.base.get() : f.exposant.get();
                        break;
                    }
                }
                break;
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*n);
                suivant = p.getBase()->contientVariable() ? p.getBase().get() : p.getExposant().get();
                break;
            }
            case TypeNoeud::IntegraleNonEvaluee:
                suivant = static_cast<const IntegraleNonEvaluee&>(*n).getIntegrande().get();
                break;
            case TypeNoeud::LimiteNonEvaluee:
                suivant = static_cast<const LimiteNonEvaluee&>(*n).getExpression().get();
                break;
            default:
                if (const FonctionUnaire* f = comme<FonctionUnaire>(n)) suivant = f->m_argument.get();
        }
        n = suivant;
    }
    return ExprPtr();
}

namespace {

// Exposant entier pair (u^(2k) >= 0) ou racine d'ordre pair (u^(p/2k) >= 0 là où elle existe)
bool exposantPair(const ExprPtr& exposant) {
    long long k, den;
    const Constante* c = comme<Constante>(exposant);
    if (!c) return false;
    if (c->getNombre().versEntier(k)) return k % 2 == 0;
    return c->getNombre().versFraction(k, den) && den % 2 == 0;
}

// Trinôme a x^2 + b x + c (a > 0) sans racine réelle simple : positif ou nul
bool trinomePositif(const ExprPtr& u) {
    Polynome p;
    if (!Polynome::depuisExpression(u, p) || p.degre() != 2 || !p.estExact() || p.dominant().signe() <= 0) return false;
    const Nombre& b = p.coefficient(1);
    return (b * b - Nombre(4) * p.coefficient(2) * p.coefficient(0)).signe() <= 0;
}

} // namespace

bool estPositifOuNul(const ExprPtr& u) {
    switch (u->type()) {
        case TypeNoeud::Constante: return static_cast<const Constante&>(*u).getNombre().signe() > 0;
        case TypeNoeud::Pi:
        case TypeNoeud::Exponentielle:
        case TypeNoeud::ValeurAbsolue:
        case TypeNoeud::CosinusHyperbolique:
            return true;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(*u);
            return exposantPair(p.getExposant()) || estPositifOuNul(p.getBase());
        }
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(*u);
            if (p.getCoefficient().signe() <= 0) return false;
            for (const Facteur& f : p.getFacteurs()) {
                if (!exposantPair(f.exposant) && !estPositifOuNul(f.base)) return false;
            }
            return true;
        }
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(*u);
            bool termesPositifs = s.getConstante().signe() >= 0;
            for (const Terme& t : s.getTermes()) {
                if (!termesPositifs) break;
                termesPositifs = t.coefficient.signe() > 0 && estPositifOuNul(t.expression);
            }
            return termesPositifs || trinomePositif(u);
        }
        default: return false;
    }
}

ExprPtr lnAbs(const ExprPtr& u) { return ast_ln(estPositifOuNul(u) ? u : ast_abs(u)); }

bool estBaseExponentielle(const ExprPtr& c) {
    return !c->contientVariable() && (c->type() == TypeNoeud::Parametre || estPositifOuNul(c));
}

bool contientNonEvaluee(const ASTNode& e) {
    if (!e.contientVariable()) return false;
    switch (e.type()) {
        case TypeNoeud::IntegraleNonEvaluee: return true;
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) {
                if (contientNonEvaluee(*t.expression)) return true;
            }
            return false;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                if (contientNonEvaluee(*f.base) || contientNonEvaluee(*f.exposant)) return true;
            }
            return false;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return contientNonEvaluee(*p.getBase()) || contientNonEvaluee(*p.getExposant());
        }
        default: {
            const FonctionUnaire* f = comme<FonctionUnaire>(&e);
            return f && contientNonEvaluee(*f->m_argument);
        }
    }
}

namespace {

// ============================================================================
// Outils des méthodes
// ============================================================================

// f est-il l'intégrale non évaluée de lui-même (échec de la règle du noeud) ?
bool echecDe(const ExprPtr& r, const ASTNode& f) {
    const IntegraleNonEvaluee* i = comme<IntegraleNonEvaluee>(r);
    return i && i->getIntegrande().get() == &f;
}

/*
 * Nom : Remplacement
 * Description : Substitution cible -> remplacement (comme substituer()), qui reconnaît en
 *               plus les puissances de la cible écrites autrement :
 *                 cible exp(v) : exp(r v) = cible^r (r rationnel) — exp(2x) = exp(x)^2 ;
 *                 cible exp(v) : exp(r v + w) = cible^r exp(w) — exp(x + exp(x)) ;
 *                 cible x^n    : x^(k n) = cible^k (k entier) — x^4 = (x^2)^2 ;
 *                 cible ln(x)  : x = exp(cible) — sin(ln(x)) dx = sin(u) exp(u) du ;
 *               et réécrit abs(w) en w quand w est manifestement positif (retour à x :
 *               ln(abs(u)) avec u = x^2 + 1 devient ln(x^2 + 1)).
 */
class Remplacement {
public:
    Remplacement(const ExprPtr& cible, const ExprPtr& remplacement, bool inverserLn = false)
        : m_cible(cible), m_remplacement(remplacement) {
        if (const Exponentielle* e = comme<Exponentielle>(cible)) {
            m_argumentExp = e->m_argument;
            // v = k * terme : exp(k r terme + w) = cible^r exp(w)
            m_termeExp = m_argumentExp;
            m_coefficientExp = Nombre(1);
            if (const Produit* p = comme<Produit>(m_argumentExp)) {
                m_termeExp = produitSansCoefficient(*p);
                m_coefficientExp = p->getCoefficient();
            }
        }
        if (const Logarithme* l = comme<Logarithme>(cible); inverserLn && l && l->m_argument->type() == TypeNoeud::Variable) {
            m_variableLn = l->m_argument;
        }
        if (const Puissance* p = comme<Puissance>(cible); p && p->getBase()->type() == TypeNoeud::Variable) {
            if (const Constante* n = comme<Constante>(p->getExposant())) {
                m_variable = p->getBase();
                m_exposant = n->getNombre();
            }
        }
    }

    ExprPtr appliquerA(const ExprPtr& e) {
        if (e.get() == m_cible.get()) return m_remplacement;
        if (!e->contientVariable()) return e;
        const auto it = m_memo.find(e.get());
        if (it != m_memo.end()) return it->second;
        ExprPtr r;
        switch (e->type()) {
            case TypeNoeud::Somme: {
                const Somme& s = static_cast<const Somme&>(*e);
                AccumulateurSomme acc;
                acc.ajouter(nombre(s.getConstante()), Nombre(1));
                for (const Terme& t : s.getTermes()) acc.ajouter(appliquerA(t.expression), t.coefficient);
                r = acc.construire();
                break;
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(*e);
                AccumulateurProduit acc;
                acc.multiplier(nombre(p.getCoefficient()));
                for (const Facteur& f : p.getFacteurs()) acc.multiplier(puissance(f.base, f.exposant));
                r = acc.construire();
                break;
            }
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*e);
                r = puissance(p.getBase(), p.getExposant());
                break;
            }
            case TypeNoeud::Variable: {
                if (e.get() == m_variableLn.get()) { // cible ln(x) : x = exp(cible)
                    r = ast_exp(m_remplacement);
                    break;
                }
                // x = cible^(1/n) si 1/n est entier (cible x^(1/2) : x = cible^2)
                const Nombre k = m_variable ? m_exposant.inverse() : Nombre(0);
                r = e.get() == m_variable.get() && k.estEntier() ? ast_pow(m_remplacement, nombre(k)) : e;
                break;
            }
            case TypeNoeud::Exponentielle:
                r = exponentielle(static_cast<const Exponentielle&>(*e).m_argument);
                break;
            case TypeNoeud::ValeurAbsolue: {
                const ExprPtr a = appliquerA(static_cast<const ValeurAbsolue&>(*e).m_argument);
                r = estPositifOuNul(a) ? a : ast_abs(a);
                break;
            }
            case TypeNoeud::IntegraleNonEvaluee:
                r = fabriquer<IntegraleNonEvaluee>(appliquerA(static_cast<const IntegraleNonEvaluee&>(*e).getIntegrande()));
                break;
            default:
                if (const FonctionUnaire* f = comme<FonctionUnaire>(e)) {
                    r = appliquer(e->type(), appliquerA(f->m_argument));
                } else {
                    r = e;
                }
        }
        m_memo.emplace(e.get(), r);
        return r;
    }

private:
    // base^exposant, avec x^(k n) = cible^k quand la cible est x^n
    ExprPtr puissance(const ExprPtr& base, const ExprPtr& exposant) {
        if (m_variable && base.get() == m_variable.get()) {
            if (const Constante* m = comme<Constante>(exposant)) {
                const Nombre k = m->getNombre() / m_exposant;
                if (k.estEntier()) return ast_pow(m_remplacement, nombre(k));
            }
        }
        return ast_pow(appliquerA(base), appliquerA(exposant));
    }

    // exp(w), avec exp(r v) = cible^r quand la cible est exp(v)
    ExprPtr exponentielle(const ExprPtr& w) {
        if (m_argumentExp) {
            const ExprPtr rapport = w / m_argumentExp;
            const Constante* c = comme<Constante>(rapport);
            if (c && c->getNombre().estExact()) return ast_pow(m_remplacement, rapport);
            // exp(x + exp(x)) = cible * exp(exp(x)) : le terme proportionnel à v sort de la somme
            if (const Somme* s = comme<Somme>(w); s && m_termeExp->type() != TypeNoeud::Somme) {
                for (const Terme& t : s->getTermes()) {
                    if (t.expression.get() != m_termeExp.get() || !t.coefficient.estExact()) continue;
                    const ExprPtr reste = w - t.expression * nombre(t.coefficient);
                    return ast_pow(m_remplacement, nombre(t.coefficient / m_coefficientExp)) * ast_exp(appliquerA(reste));
                }
            }
        }
        const ExprPtr argument = appliquerA(w);
        // exp(r ln(z)) = z^r (retour à x après u = ln(x))
        if (const Produit* p = comme<Produit>(argument); p && p->getFacteurs().size() == 1) {
            const Facteur& f = p->getFacteurs()[0];
            if (f.exposant.get() == un().get() && f.base->type() == TypeNoeud::Logarithme) {
                return ast_pow(static_cast<const Logarithme&>(*f.base).m_argument, nombre(p->getCoefficient()));
            }
        }
        return ast_exp(argument);
    }

    ExprPtr m_cible, m_remplacement;
    ExprPtr m_argumentExp; // v si la cible est exp(v)
    ExprPtr m_termeExp;    // v sans son coefficient numérique
    Nombre m_coefficientExp;
    ExprPtr m_variableLn;  // x si la cible est ln(x)
    ExprPtr m_variable;    // x si la cible est x^n
    Nombre m_exposant;     // n
    std::unordered_map<const ASTNode*, ExprPtr> m_memo;
};

// ============================================================================
// Méthode 2 : table étendue
// ============================================================================

// Racine carrée exacte (radical canonique) d'un nombre positif
ExprPtr racine(const Nombre& n) { return ast_pow(nombre(n), frac(1, 2)); }

/*
 * Nom : trinome
 * Description : 1/(a x^2 + b x + c) et 1/sqrt(a x^2 + b x + c), coefficients rationnels, par
 *               la forme canonique a (u^2 + d/a) avec u = x + b/(2a), d = c - b^2/(4a) :
 *                 1/(a(u^2 + k^2))     -> atan(u/k) / (a k)
 *                 1/(a u^2)            -> -1/(a u)
 *                 1/sqrt(d - |a| u^2)  -> asin(u sqrt(|a|/d)) / sqrt(|a|)       (a < 0 < d)
 *                 1/sqrt(a(u^2 + d/a)) -> ln|u + sqrt(u^2 + d/a)| / sqrt(a)     (a > 0, d != 0)
 *               Les autres cas (racines réelles distinctes : fractions rationnelles) sont
 *               laissés aux méthodes suivantes.
 */
ExprPtr trinome(const ExprPtr& base, bool racineCarree) {
    Polynome p;
    ExprPtr x;
    if (!Polynome::depuisExpression(base, p, &x)) {
        // (2x + 1)^2 + 1 : développé seulement s'il est court
        const Somme* s = comme<Somme>(base);
        if (!s || s->getTermes().size() > 3) return nullptr;
        const ExprPtr d = developper(base);
        if (d.get() == base.get() || !Polynome::depuisExpression(d, p, &x)) return nullptr;
    }
    if (!x || p.degre() != 2 || !p.estExact()) return nullptr;
    const Nombre& a = p.coefficient(2);
    const Nombre& b = p.coefficient(1);
    const Nombre d = p.coefficient(0) - b * b / (Nombre(4) * a);
    const ExprPtr u = x + nombre(b / (Nombre(2) * a));
    if (!racineCarree) {
        if (d.estZero()) return -(nombre(a.inverse()) / u);
        const Nombre k2 = d / a;
        if (k2.signe() <= 0) return nullptr;
        const ExprPtr k = racine(k2);
        return ast_atan(u / k) / (nombre(a) * k);
    }
    if (a.signe() < 0 && d.signe() > 0) return ast_asin(u * racine(-a / d)) / racine(-a);
    if (a.signe() > 0 && !d.estZero()) {
        const ExprPtr w = u + ast_pow(base, frac(1, 2)) / racine(a); // u + sqrt(u^2 + d/a) > 0 si d > 0
        return (d.signe() > 0 ? ast_ln(w) : lnAbs(w)) / racine(a);
    }
    return nullptr;
}

} // namespace

ExprPtr tableEtendue(const ExprPtr& f) {
    const Puissance* p = comme<Puissance>(f);
    if (!p) return nullptr;
    const ExprPtr& base = p->getBase();
    const ExprPtr& exposant = p->getExposant();
    ExprPtr a;
    // c^(ax + b), c constante symbolique positive (pi, 2^(1/2)...) : c^u / (a ln c)
    if (!base->contientVariable()) {
        if (estBaseExponentielle(base) && coefficientLineaire(exposant, a)) return f / (a * ast_ln(base));
        return nullptr;
    }
    const Constante* c = comme<Constante>(exposant);
    if (!c) return nullptr;
    const Nombre& n = c->getNombre();
    if (n.estMoinsUn()) return trinome(base, false);
    if (n == Nombre::rationnel(-1, 2)) return trinome(base, true);
    const bool carre = n == Nombre(2);
    if (!carre && n != Nombre(-2)) return nullptr;
    // Carrés et inverses de carrés de fonctions trigonométriques et hyperboliques d'un
    // argument affine v = a x + b
    const FonctionUnaire* fn = comme<FonctionUnaire>(base);
    if (!fn || !coefficientLineaire(fn->m_argument, a)) return nullptr;
    const ExprPtr& v = fn->m_argument;
    const ExprPtr x = variableDe(*v);
    switch (base->type()) {
        case TypeNoeud::Cosinus: // 1/cos^2 v -> tan v
            return carre ? nullptr : ast_tan(v) / a;
        case TypeNoeud::Sinus: // 1/sin^2 v -> -cos v / sin v
            return carre ? nullptr : -(ast_cos(v) / ast_sin(v)) / a;
        case TypeNoeud::Tangente: // tan^2 v = 1/cos^2 v - 1
            return carre ? ast_tan(v) / a - x : nullptr;
        case TypeNoeud::CosinusHyperbolique: // 1/cosh^2 v -> tanh v
            return carre ? nullptr : ast_tanh(v) / a;
        case TypeNoeud::SinusHyperbolique: // 1/sinh^2 v -> -cosh v / sinh v
            return carre ? nullptr : -(ast_cosh(v) / ast_sinh(v)) / a;
        case TypeNoeud::TangenteHyperbolique: // tanh^2 v = 1 - 1/cosh^2 v
            return carre ? x - ast_tanh(v) / a : nullptr;
        default: return nullptr;
    }
}

namespace {

// ============================================================================
// Changement de variable « la dérivée divise »
// ============================================================================

constexpr std::size_t CANDIDATS_MAX = 8;   // changements de variable essayés par intégrale
constexpr std::size_t NOEUDS_PARCOURUS_MAX = 64;
constexpr int PROFONDEUR_SUBSTITUTION_MAX = 2; // substitutions imbriquées

int& profondeurSubstitution() {
    thread_local int profondeur = 0;
    return profondeur;
}

class Imbrication {
public:
    Imbrication() { ++profondeurSubstitution(); }
    ~Imbrication() { --profondeurSubstitution(); }
    Imbrication(const Imbrication&) = delete;
    Imbrication& operator=(const Imbrication&) = delete;
};

/*
 * Nom : Candidats
 * Description : Sous-expressions u de f candidates au changement de variable, dans l'ordre
 *               d'un parcours en profondeur, sans doublon (comparaison d'adresses grâce au
 *               hash-consing) : bases et exposants des puissances, fonctions et leurs
 *               arguments. Sont écartés x lui-même, les fonctions affines de x (déjà
 *               traitées par les règles des noeuds) et f. Nombre borné.
 */
class Candidats {
public:
    Candidats(const ExprPtr& f, const ExprPtr& x) : m_f(f.get()), m_x(x.get()) { parcourir(*f); }
    const std::vector<const ASTNode*>& liste() const { return m_liste; }

private:
    bool plein() const { return m_liste.size() >= CANDIDATS_MAX || m_vus.size() >= NOEUDS_PARCOURUS_MAX; }

    // x, a x ou a x + b
    bool affine(const ASTNode& u) const {
        const ASTNode* n = &u;
        if (const Somme* s = comme<Somme>(n); s && s->getTermes().size() == 1) n = s->getTermes()[0].expression.get();
        if (const Produit* p = comme<Produit>(n); p && p->getFacteurs().size() == 1 && p->getFacteurs()[0].exposant.get() == un().get()) {
            n = p->getFacteurs()[0].base.get();
        }
        return n == m_x;
    }

    void proposer(const ASTNode& u) {
        if (plein() || !u.contientVariable() || &u == m_f || affine(u)) return;
        for (const ASTNode* c : m_liste) {
            if (c == &u) return;
        }
        m_liste.push_back(&u);
    }

    void parcourir(const ASTNode& e) {
        if (plein() || !e.contientVariable() || &e == m_x) return;
        for (const ASTNode* v : m_vus) {
            if (v == &e) return;
        }
        m_vus.push_back(&e);
        switch (e.type()) {
            case TypeNoeud::Somme:
                for (const Terme& t : static_cast<const Somme&>(e).getTermes()) parcourir(*t.expression);
                break;
            case TypeNoeud::Produit:
                for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                    proposer(*f.base);
                    proposer(*f.exposant);
                }
                for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                    parcourir(*f.base);
                    parcourir(*f.exposant);
                }
                break;
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(e);
                proposer(*p.getBase());
                proposer(*p.getExposant());
                parcourir(*p.getBase());
                parcourir(*p.getExposant());
                break;
            }
            default:
                if (const FonctionUnaire* f = comme<FonctionUnaire>(&e)) {
                    proposer(e);
                    proposer(*f->m_argument);
                    parcourir(*f->m_argument);
                }
        }
    }

    const ASTNode* m_f;
    const ASTNode* m_x;
    std::vector<const ASTNode*> m_liste;
    std::vector<const ASTNode*> m_vus;
};

/*
 * Nom : substitution
 * Description : int f(x) dx = int g(u) du avec u = h(x) quand f = g(h(x)) h'(x) : pour chaque
 *               candidat u, q = f/u' (forme canonique, puis réduite par simplifier() si
 *               besoin) doit s'écrire sans x une fois u remplacé par une nouvelle variable.
 *               G = int g est calculée par integrer() (substitutions imbriquées bornées), puis
 *               u est remis à la place de la variable. Échec si G n'est que partiellement
 *               calculée. Avec inverserLn, seul u = ln(x) est essayé, et x est remplacé par
 *               exp(u) : sin(ln(x)) dx = sin(u) exp(u) du.
 */
ExprPtr substitution(const ExprPtr& f, bool inverserLn) {
    if (profondeurSubstitution() >= PROFONDEUR_SUBSTITUTION_MAX) return nullptr;
    const ExprPtr x = variableDe(*f);
    if (!x) return nullptr;
    // Nouvelle variable, distincte de x (nommée u pour le journal)
    const ExprPtr t = var(static_cast<const Variable&>(*x).getNom() == "u" ? "t" : "u");
    const std::uint16_t idT = t->identifiantVariable();
    auto sansX = [idT](const ExprPtr& g) { return !g->contientVariable() || g->identifiantVariable() == idT; };

    const Candidats candidats(f, x);
    for (const ASTNode* c : candidats.liste()) {
        // Seconde passe : seulement u = ln(x), avec x = exp(u)
        if (inverserLn && !(c->type() == TypeNoeud::Logarithme && static_cast<const Logarithme&>(*c).m_argument.get() == x.get())) continue;
        const ExprPtr u = c->clone();
        const ExprPtr du = u->derivee();
        const Constante* duConstante = comme<Constante>(du);
        if (duConstante && duConstante->getNombre().estZero()) continue; // sin^2 + cos^2...
        const ExprPtr q = f / du;
        ExprPtr g = Remplacement(u, t, inverserLn).appliquerA(q);
        if (!sansX(g)) {
            // (x + 1)/(2x + 2) : seule la réduction rationnelle fait disparaître x
            const ExprPtr qs = q->simplifier();
            if (qs.get() == q.get()) continue;
            g = Remplacement(u, t, inverserLn).appliquerA(qs);
            if (!sansX(g)) continue;
        }

        Journal* const j = journal::actif();
        const std::size_t etapes = j ? j->etapes().size() : 0;
        journal::noter([&](Journal& jr) {
            jr.ajouter("integration.substitution",
                       "Changement de variable " + t->texte() + " = " + u->texte() + ", d" + t->texte() + " = " +
                           du->texte() + " d" + x->texte(),
                       f, g);
        });
        ExprPtr G;
        {
            const Imbrication imbrication;
            const journal::SousNiveau sous;
            G = g->contientVariable() ? g->integrer() : g * t;
        }
        if (contientNonEvaluee(*G)) {
            if (j) j->tronquer(etapes);
            continue;
        }
        const ExprPtr resultat = Remplacement(t, u).appliquerA(G);
        journal::noter([&](Journal& jr) {
            jr.ajouter("integration.retour", "Retour à la variable " + x->texte() + " : " + t->texte() + " = " + u->texte(),
                       G, resultat);
        });
        return resultat;
    }
    return nullptr;
}

} // namespace

ExprPtr parSubstitution(const ExprPtr& f) { return substitution(f, false); }

ExprPtr parSubstitutionLogarithme(const ExprPtr& f) { return substitution(f, true); }

// ============================================================================
// Intégration par parties
// ============================================================================

namespace {

constexpr int PROFONDEUR_PARTIES_MAX = 6; // ln(x)^5 demande cinq intégrations par parties

int& profondeurParties() {
    thread_local int profondeur = 0;
    return profondeur;
}

// Fonctions que l'on dérive dans une intégration par parties : leur dérivée est algébrique
bool estLogarithmique(const ASTNode& e) {
    const TypeNoeud t = e.type();
    return (t == TypeNoeud::Logarithme || t == TypeNoeud::ArcTangente || t == TypeNoeud::ArcSinus || t == TypeNoeud::ArcCosinus) &&
           e.contientVariable();
}

bool entierPositif(const ExprPtr& exposant) {
    long long k;
    const Constante* c = comme<Constante>(exposant);
    return c && c->getNombre().versEntier(k) && k >= 1;
}

} // namespace

ExprPtr parParties(const ExprPtr& f) {
    if (profondeurParties() >= PROFONDEUR_PARTIES_MAX) return nullptr;
    // Le facteur L = ln(v)^n, atan(v), asin(v), acos(v)... à dériver
    ExprPtr L;
    if (estLogarithmique(*f)) {
        L = f;
    } else if (const Puissance* p = comme<Puissance>(f)) {
        if (estLogarithmique(*p->getBase()) && entierPositif(p->getExposant())) L = f;
    } else if (const Produit* p = comme<Produit>(f)) {
        for (const Facteur& facteur : p->getFacteurs()) {
            if (estLogarithmique(*facteur.base) && entierPositif(facteur.exposant)) {
                L = ast_pow(facteur.base, facteur.exposant);
                break;
            }
        }
    }
    if (!L) return nullptr;
    const ExprPtr x = variableDe(*f);
    const ExprPtr g = f / L;
    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    auto echec = [&]() {
        if (j) j->tronquer(etapes);
        return ExprPtr();
    };
    ++profondeurParties();
    struct Sortie {
        ~Sortie() { --profondeurParties(); }
    } sortie;
    ExprPtr G;
    {
        // Primitive du cofacteur : calcul annexe, sans étapes (l'étape « par parties » la donne)
        journal::actif() = nullptr;
        struct Retablir {
            Journal* journal;
            ~Retablir() { journal::actif() = journal; }
        } retablir{j};
        G = g->contientVariable() ? g->integrer() : g * x;
    }
    if (contientNonEvaluee(*G)) return echec();
    const ExprPtr reste = G * L->derivee();
    journal::noter([&](Journal& jr) {
        jr.ajouter("integration.parties",
                   "Intégration par parties : on dérive " + L->texte() + " et on intègre " + g->texte(), f,
                   G * L - fabriquer<IntegraleNonEvaluee>(reste));
    });
    ExprPtr I;
    {
        const journal::SousNiveau sous;
        I = reste->contientVariable() ? reste->integrer() : reste * x;
    }
    if (contientNonEvaluee(*I)) return echec();
    return G * L - I;
}

// ============================================================================
// Suite des méthodes
// ============================================================================

namespace {

struct Methode {
    ExprPtr (*appliquer)(const ExprPtr& f);
    const char* regle;
    const char* description;
};

// Du moins cher au plus cher ; chaque méthode rejette d'abord sur le type du noeud, puis sur
// la forme. Les méthodes sans règle notent elles-mêmes leurs étapes.
const Methode METHODES[] = {
    {tableEtendue, "integration.table", "Primitive usuelle"},
    {parTrigonometrieImpaire, "integration.rationnelle", "Primitive du polynôme en la nouvelle variable"},
    {parExponentiellePolynome, "integration.exponentielle_polynome",
     "Produit de polynômes, d'exponentielles, de sinus et de cosinus : primitive explicite de chaque terme"},
    {parFractionRationnelle, "integration.rationnelle",
     "Fraction rationnelle : division, réduction d'Hermite et éléments simples"},
    {parTrigonometrie, "integration.rationnelle", "Fraction rationnelle de la nouvelle variable"},
    {parExponentielleRationnelle, "integration.rationnelle", "Fraction rationnelle de la nouvelle variable"},
    {parSubstitution, nullptr, nullptr},
    {parParties, nullptr, nullptr},
    {parSubstitutionLogarithme, nullptr, nullptr},
    {parRadicaux, "integration.radicaux", "Fraction rationnelle de x et d'un radical"},
    {parEquationDeRisch, "integration.risch", "Primitive de la forme S(x) exp(W(x)) : S' + W' S = R"},
    {parDeveloppement, nullptr, nullptr},
};

} // namespace

// ============================================================================
// Point d'entrée
// ============================================================================

ExprPtr ASTNode::integrer() const {
    if (!contientVariable()) {
        const ExprPtr r = clone() * var("x");
        journal::noter([&](Journal& j) { j.ajouter("integration.constante", "Primitive d'une constante", clone(), r); });
        return r;
    }
    if (plusieursVariables()) return integraleNonEvaluee();
    // 1. Règle du noeud : linéarité (Somme, Produit) ou primitive usuelle d'un argument affine
    ExprPtr r = primitive();
    if (!echecDe(r, *this)) {
        if (m_type != TypeNoeud::Somme && m_type != TypeNoeud::Produit) {
            journal::noter([&](Journal& j) { j.ajouter("integration.table", "Primitive usuelle", clone(), r); });
        }
        return r;
    }
    // 2. Méthodes suivantes, dans l'ordre
    const ExprPtr f = clone();
    for (const Methode& m : METHODES) {
        if (ExprPtr resultat = m.appliquer(f)) {
            if (m.regle) journal::noter([&](Journal& j) { j.ajouter(m.regle, m.description, f, resultat); });
            return resultat;
        }
    }
    journal::noter([&](Journal& j) { j.ajouter("integration.non_evaluee", "Aucune méthode ne s'applique", f, r); });
    return r;
}

ExprPtr ASTNode::integrer(Journal& journal) const {
    // Pas de mémo de primitive à contourner : chaque appel refait toutes les étapes
    const journal::Activation activation(journal);
    return integrer();
}

} // namespace symalgo
