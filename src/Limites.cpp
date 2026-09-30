/**
 * @file Limites.cpp
 * @brief Calcul des limites en un point.
 *
 * Chaque noeud combine les limites de ses enfants : nombres (finis ou infinis),
 * expressions symboliques (paramètres) ou limites non déterminées. Les formes
 * indéterminées sont ramenées à un quotient N/D (dans la forme canonique, un quotient est
 * un produit à exposants négatifs) traité par la règle de L'Hôpital ; le signe d'un
 * infini c/0 vient du premier terme non nul du développement de Taylor de D.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"

#include <cmath>
#include <limits>

namespace symalgo {

namespace {

constexpr double INFINI = std::numeric_limits<double>::infinity();
constexpr double NAN_LIMITE = std::numeric_limits<double>::quiet_NaN();

// Profondeur maximale des réécritures récursives (L'Hôpital, 0·∞, 1^∞...)
constexpr int PROFONDEUR_MAX_REECRITURE = 8;
thread_local int g_profondeurReecriture = 0;

struct GardeProfondeur {
    GardeProfondeur() { ++g_profondeurReecriture; }
    ~GardeProfondeur() { --g_profondeurReecriture; }
    GardeProfondeur(const GardeProfondeur&) = delete;
    GardeProfondeur& operator=(const GardeProfondeur&) = delete;
};

bool reecriturePossible() { return g_profondeurReecriture < PROFONDEUR_MAX_REECRITURE; }

// Seuil en dessous duquel une limite numérique est considérée comme nulle
// (absorbe les résidus d'arrondi comme sin(pi) = 1.2e-16)
bool estNul(double v) { return std::abs(v) < 1e-12; }

enum class GenreLimite { Nombre, Symbolique, NonEvaluee };

/*
 * Nom : analyserLimite
 * Description : Classe le résultat d'un calcul de limite : nombre (écrit dans v),
 *               expression symbolique (paramètres) ou limite non déterminée. Une expression
 *               constante non numérique (ex. sin(2)) est évaluée.
 */
GenreLimite analyserLimite(const ExprPtr& limite, double& v) {
    if (limite->type() == TypeNoeud::LimiteNonEvaluee) return GenreLimite::NonEvaluee;
    try {
        v = limite->eval(0.0);
    } catch (const std::logic_error&) {
        return GenreLimite::Symbolique; // contient un paramètre
    }
    return std::isnan(v) ? GenreLimite::NonEvaluee : GenreLimite::Nombre;
}

// Constante v (entière exacte si possible), ou nullptr si v n'est pas un nombre
ExprPtr nombreOuNul(double v) { return std::isnan(v) ? nullptr : cst(v); }

bool estDeterminee(const ExprPtr& limite) {
    double v;
    return limite && analyserLimite(limite, v) != GenreLimite::NonEvaluee;
}

/*
 * Nom : premierTermeNonNul
 * Description : Ordre k et signe du premier terme non nul du développement de Taylor de f
 *               en a (f(a) étant nul). Au voisinage de a, f(x) ~ f^(k)(a)/k! (x-a)^k :
 *               k pair => même signe des deux côtés de a, k impair => changement de signe.
 */
bool premierTermeNonNul(const ExprPtr& f, double a, int& ordre, double& signe) {
    ExprPtr d = f;
    for (int k = 1; k <= PROFONDEUR_MAX_REECRITURE; ++k) {
        d = d->derivee();
        double v;
        try {
            v = d->eval(a);
        } catch (const std::logic_error&) {
            return false; // paramètre symbolique : signe inconnu
        }
        if (!std::isfinite(v)) return false;
        if (!estNul(v)) {
            ordre = k;
            signe = v > 0.0 ? 1.0 : -1.0;
            return true;
        }
    }
    return false;
}

/*
 * Nom : tendVersInfiniNonSigne
 * Description : Vrai si |f| tend vers +inf en a sans que f ait de limite (signes opposés de
 *               part et d'autre, ex. 1/x en 0) : f = c * b^(-k) avec b -> 0.
 */
bool tendVersInfiniNonSigne(const ExprPtr& f, double a) {
    // f = b^e (Puissance) ou c * b^e (Produit à un seul facteur)
    const ExprPtr* base = nullptr;
    const ExprPtr* exposant = nullptr;
    if (const Puissance* p = comme<Puissance>(f)) {
        base = &p->getBase();
        exposant = &p->getExposant();
    } else if (const Produit* q = comme<Produit>(f); q && q->getFacteurs().size() == 1) {
        base = &q->getFacteurs()[0].base;
        exposant = &q->getFacteurs()[0].exposant;
    } else {
        return false;
    }
    const Constante* e = comme<Constante>(*exposant);
    if (!e || e->getNombre().signe() >= 0) return false;
    double b;
    return analyserLimite((*base)->limite(a), b) == GenreLimite::Nombre && estNul(b);
}

/*
 * Nom : limiteQuotient
 * Description : Limite de N/D en a : L'Hôpital pour 0/0 et inf/inf (profondeur bornée),
 *               analyse de signe pour c/0. Renvoie nullptr si elle n'est pas déterminée.
 */
ExprPtr limiteQuotient(const ExprPtr& N, const ExprPtr& D, double a) {
    double n = 0.0, d = 0.0;
    const ExprPtr Ln = N->limite(a);
    const ExprPtr Ld = D->limite(a);
    GenreLimite tn = analyserLimite(Ln, n);
    GenreLimite td = analyserLimite(Ld, d);
    // Un infini sans signe (1/x en 0) reste un infini pour la règle de L'Hôpital
    if (tn == GenreLimite::NonEvaluee && tendVersInfiniNonSigne(N, a)) {
        tn = GenreLimite::Nombre;
        n = INFINI;
    }
    if (td == GenreLimite::NonEvaluee && tendVersInfiniNonSigne(D, a)) {
        td = GenreLimite::Nombre;
        d = INFINI;
    }
    if (tn == GenreLimite::NonEvaluee || td == GenreLimite::NonEvaluee) return nullptr;
    if (tn == GenreLimite::Symbolique || td == GenreLimite::Symbolique) return Ln / Ld;

    if ((estNul(n) && estNul(d)) || (std::isinf(n) && std::isinf(d))) {
        // Règle de L'Hôpital : lim N/D = lim N'/D'. Le quotient N'/D' est construit sous
        // forme canonique, ce qui simplifie les facteurs communs (x^-1 / -x^-2 = -x)
        if (!reecriturePossible()) return nullptr;
        GardeProfondeur garde;
        const ExprPtr r = (N->derivee() / D->derivee())->limite(a);
        return estDeterminee(r) ? r : nullptr;
    }
    if (estNul(d)) {
        // c / 0 : l'infini n'a un signe défini que si D garde le même signe des deux côtés
        int ordre;
        double signeD;
        if (premierTermeNonNul(D, a, ordre, signeD) && ordre % 2 == 0) {
            return cst((n > 0.0 ? 1.0 : -1.0) * signeD * INFINI);
        }
        return nullptr; // limites à gauche et à droite différentes
    }
    return nombreOuNul(n / d);
}

/*
 * Nom : limiteUnaire
 * Description : Limite de f(u) : f appliquée à la limite de l'argument u (nullptr si non
 *               déterminée).
 */
ExprPtr limiteUnaire(const ExprPtr& argument, double a, double (*fNum)(double), ExprPtr (*fSym)(const ExprPtr&)) {
    double u;
    const ExprPtr L = argument->limite(a);
    switch (analyserLimite(L, u)) {
        case GenreLimite::Nombre: return nombreOuNul(fNum(u));
        case GenreLimite::Symbolique: return fSym(L);
        case GenreLimite::NonEvaluee: break;
    }
    return nullptr;
}

bool exposantNumeriqueNegatif(const ExprPtr& e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().signe() < 0;
}

} // namespace

ExprPtr Constante::calculerLimite(double) const { return clone(); }
ExprPtr Variable::calculerLimite(double a) const { return cst(a); }
ExprPtr Parametre::calculerLimite(double) const { return clone(); }
ExprPtr Pi::calculerLimite(double) const { return clone(); }

ExprPtr Somme::calculerLimite(double a) const {
    double total = m_constante.versDouble();
    bool symbolique = false;
    AccumulateurSomme s;
    s.ajouter(nombre(m_constante), Nombre(1));
    for (std::size_t i = 0; i < m_termes.size(); ++i) {
        const ExprPtr L = m_termes[i].expression->limite(a);
        double v;
        switch (analyserLimite(L, v)) {
            case GenreLimite::NonEvaluee: return limiteNonEvaluee(a);
            case GenreLimite::Symbolique: symbolique = true; break;
            case GenreLimite::Nombre: total += m_termes[i].coefficient.versDouble() * v; break;
        }
        s.ajouter(L, m_termes[i].coefficient);
    }
    if (symbolique) return s.construire();
    if (ExprPtr r = nombreOuNul(total)) return r; // +inf + -inf : forme indéterminée non résolue
    return limiteNonEvaluee(a);
}

ExprPtr Produit::calculerLimite(double a) const {
    // Limites des facteurs b_i^e_i
    std::vector<double> valeurs;
    bool symbolique = false, indetermine = false, zero = false, infini = false;
    AccumulateurProduit p;
    p.multiplier(nombre(m_coefficient));
    for (const Facteur& f : m_facteurs) {
        const ExprPtr L = ast_pow(f.base, f.exposant)->limite(a);
        double v = 0.0;
        switch (analyserLimite(L, v)) {
            case GenreLimite::NonEvaluee: indetermine = true; break;
            case GenreLimite::Symbolique: symbolique = true; break;
            case GenreLimite::Nombre:
                zero = zero || estNul(v);
                infini = infini || std::isinf(v);
                break;
        }
        valeurs.push_back(v);
        p.multiplier(L);
    }
    if (!indetermine && !(zero && infini)) {
        if (symbolique) return p.construire();
        double produitValeurs = m_coefficient.versDouble();
        for (double v : valeurs) produitValeurs *= v;
        if (ExprPtr r = nombreOuNul(produitValeurs)) return r;
        return limiteNonEvaluee(a);
    }
    if (!reecriturePossible()) return limiteNonEvaluee(a);
    GardeProfondeur garde;
    // 1. Quotient naturel : exposants positifs au numérateur, négatifs au dénominateur
    //    (x * x^(-1) sous forme de facteurs séparés, sin(x) * x^(-1)...)
    AccumulateurProduit numerateur, denominateur;
    numerateur.multiplier(nombre(m_coefficient));
    bool aDenominateur = false;
    for (const Facteur& f : m_facteurs) {
        if (exposantNumeriqueNegatif(f.exposant)) {
            denominateur.multiplier(ast_pow(f.base, -f.exposant));
            aDenominateur = true;
        } else {
            numerateur.multiplier(ast_pow(f.base, f.exposant));
        }
    }
    if (aDenominateur) {
        if (ExprPtr r = limiteQuotient(numerateur.construire(), denominateur.construire(), a)) return r;
    }
    // 2. Forme 0 * inf : z * w = z / (1/w), sinon w / (1/z)
    if (zero && infini) {
        AccumulateurProduit nuls, infinis, inversesNuls, inversesInfinis;
        nuls.multiplier(nombre(m_coefficient));
        infinis.multiplier(nombre(m_coefficient));
        for (std::size_t i = 0; i < m_facteurs.size(); ++i) {
            const ExprPtr facteur = ast_pow(m_facteurs[i].base, m_facteurs[i].exposant);
            if (std::isinf(valeurs[i])) {
                infinis.multiplier(facteur);
                inversesInfinis.multiplier(ast_pow(facteur, nombre(Nombre(-1))));
            } else {
                nuls.multiplier(facteur);
                inversesNuls.multiplier(ast_pow(facteur, nombre(Nombre(-1))));
            }
        }
        if (ExprPtr r = limiteQuotient(nuls.construire(), inversesInfinis.construire(), a)) return r;
        if (ExprPtr r = limiteQuotient(infinis.construire(), inversesNuls.construire(), a)) return r;
    }
    return limiteNonEvaluee(a);
}

ExprPtr Puissance::calculerLimite(double a) const {
    // u^v = exp(v ln u) : lève les formes 1^inf, 0^0 et inf^0 quand l'exposant dépend de x
    auto parExponentielle = [&]() -> ExprPtr {
        if (!m_exposant->contientVariable() || !reecriturePossible()) return nullptr;
        GardeProfondeur garde;
        const ExprPtr r = ast_exp(m_exposant * ast_ln(m_base))->limite(a);
        return estDeterminee(r) ? r : nullptr;
    };
    // Exposant négatif : quotient 1 / b^(-e), dont le signe de l'infini dépend du côté
    if (exposantNumeriqueNegatif(m_exposant)) {
        if (!reecriturePossible()) return limiteNonEvaluee(a);
        GardeProfondeur garde;
        if (ExprPtr r = limiteQuotient(un(), ast_pow(m_base, -m_exposant), a)) return r;
        return limiteNonEvaluee(a);
    }
    double b = 0.0, e = 0.0;
    const ExprPtr Lb = m_base->limite(a);
    const ExprPtr Le = m_exposant->limite(a);
    const GenreLimite tb = analyserLimite(Lb, b);
    const GenreLimite te = analyserLimite(Le, e);
    if (tb == GenreLimite::NonEvaluee || te == GenreLimite::NonEvaluee) {
        if (ExprPtr r = parExponentielle()) return r;
        return limiteNonEvaluee(a);
    }
    if (tb == GenreLimite::Symbolique || te == GenreLimite::Symbolique) return ast_pow(Lb, Le);
    const bool indeterminee = (b == 1.0 && std::isinf(e)) || (estNul(b) && estNul(e)) || (std::isinf(b) && estNul(e));
    if (indeterminee) {
        if (ExprPtr r = parExponentielle()) return r;
        return limiteNonEvaluee(a);
    }
    if (estNul(b) && e < 0.0) {
        if (!reecriturePossible()) return limiteNonEvaluee(a);
        GardeProfondeur garde;
        if (ExprPtr r = limiteQuotient(un(), ast_pow(m_base, -m_exposant), a)) return r;
        return limiteNonEvaluee(a);
    }
    if (ExprPtr r = nombreOuNul(std::pow(b, e))) return r;
    return limiteNonEvaluee(a);
}

ExprPtr Sinus::calculerLimite(double a) const {
    // sin(+-inf) donne NaN : pas de limite
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::sin(u); }, ast_sin)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr Cosinus::calculerLimite(double a) const {
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::cos(u); }, ast_cos)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr Tangente::calculerLimite(double a) const {
    // Aux pôles (cos(u) = 0), tan tend vers +inf d'un côté et -inf de l'autre
    auto tanNum = [](double u) { return estNul(std::cos(u)) ? NAN_LIMITE : std::tan(u); };
    if (ExprPtr r = limiteUnaire(m_argument, a, tanNum, ast_tan)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr Exponentielle::calculerLimite(double a) const {
    // exp(+inf) = +inf, exp(-inf) = 0
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::exp(u); }, ast_exp)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr Logarithme::calculerLimite(double a) const {
    // ln(u) -> -inf quand u -> 0 ; hors du domaine (u < 0) la limite n'existe pas
    auto lnNum = [](double u) { return estNul(u) ? -INFINI : std::log(u); };
    if (ExprPtr r = limiteUnaire(m_argument, a, lnNum, ast_ln)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr ArcSinus::calculerLimite(double a) const {
    // Hors de [-1, 1], asin donne NaN : pas de limite réelle
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::asin(u); }, ast_asin)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr ArcCosinus::calculerLimite(double a) const {
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::acos(u); }, ast_acos)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr ArcTangente::calculerLimite(double a) const {
    // atan(+-inf) = +-pi/2
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::atan(u); }, ast_atan)) return r;
    return limiteNonEvaluee(a);
}

} // namespace symalgo
