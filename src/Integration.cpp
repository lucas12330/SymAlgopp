/**
 * @file Integration.cpp
 * @brief Intégrale définie : chemin exact (primitive certifiée sans singularité) et
 *        quadrature de Gauss-Kronrod G7-K15 adaptative sur l'évaluateur compilé.
 */

#include "Integration.hpp"
#include "Evaluateur.hpp"
#include "Multivariable.hpp"
#include "Solveur.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace symalgo {

namespace {

constexpr double NON_NOMBRE = std::numeric_limits<double>::quiet_NaN();
constexpr double INFINI = std::numeric_limits<double>::infinity();
constexpr double EPSILON = std::numeric_limits<double>::epsilon();

// ============================================================================
// Parcours des sous-expressions (dupliqué de Solveur.cpp, où il est privé)
// ============================================================================

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

bool contientNonEvalue(const ExprPtr& e) {
    bool trouve = false;
    std::unordered_set<const ASTNode*> vus;
    parcourir(e, [&](const ExprPtr& n) {
        if (n->type() == TypeNoeud::IntegraleNonEvaluee || n->type() == TypeNoeud::LimiteNonEvaluee) trouve = true;
    }, vus);
    return trouve;
}

// Valeur d'une expression sans variable ; faux si elle contient un paramètre ou un noeud non évalué
bool valeurConstante(const ExprPtr& e, double& v) {
    if (e->contientVariable()) return false;
    try {
        v = e->eval(0.0);
    } catch (const std::exception&) {
        return false;
    }
    return !std::isnan(v);
}

double evaluerSur(const ExprPtr& e, double x) {
    try {
        return e->eval(x);
    } catch (const std::exception&) {
        return NON_NOMBRE;
    }
}

// ============================================================================
// Détection de singularités sur [A, B]
// ============================================================================

enum class Condition { NonNul, Positif, NonNegatif };

/*
 * Nom : Analyseur
 * Description : Vérifie qu'une expression est définie et continue sur [A, B] : pour chaque
 *               dénominateur, ln, puissance non entière, tan, asin et acos, la condition
 *               correspondante est démontrée par les racines certifiées du solveur et un
 *               test de signe en des points intérieurs. Tout ce qui ne peut pas être
 *               démontré (solveur incomplet, paramètre, bornes inconnues) rend la
 *               vérification négative : l'appelant passe alors en numérique.
 *               En mode « ouvert », seul l'intérieur de ]A, B[ compte (une singularité
 *               intégrable aux bornes est tolérée) ; en mode « fermé », les bornes aussi.
 *               Les racines intérieures qui violent une condition sont mémorisées : ce sont
 *               les points où l'intervalle numérique doit être coupé.
 */
class Analyseur {
public:
    Analyseur(double a, double b, bool ferme, bool bornesConnues)
        : m_a(a), m_b(b), m_ferme(ferme), m_bornesConnues(bornesConnues) {
        m_epsilon = 1e-10 * std::max({1.0, std::abs(a), std::abs(b)});
    }

    bool sur(const ExprPtr& e) {
        m_ok = true;
        std::unordered_set<const ASTNode*> vus;
        parcourir(e, [this](const ExprPtr& n) { examiner(n); }, vus);
        return m_ok;
    }

    const std::vector<double>& coupures() const { return m_coupures; }

private:
    struct Racines {
        bool fiable = false;
        std::vector<double> valeurs; // triées, dans [A, B]
    };

    void examiner(const ExprPtr& n) {
        switch (n->type()) {
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*n);
                puissance(p.getBase(), p.getExposant());
                break;
            }
            case TypeNoeud::Produit:
                for (const Facteur& f : static_cast<const Produit&>(*n).getFacteurs()) puissance(f.base, f.exposant);
                break;
            case TypeNoeud::Logarithme:
                exiger(static_cast<const Logarithme&>(*n).m_argument, Condition::Positif);
                break;
            case TypeNoeud::Tangente:
                exiger(ast_cos(static_cast<const Tangente&>(*n).m_argument), Condition::NonNul);
                break;
            case TypeNoeud::ArcSinus:
            case TypeNoeud::ArcCosinus: {
                const ExprPtr& u = static_cast<const FonctionUnaire&>(*n).m_argument;
                exiger(1.0 - u * u, Condition::NonNegatif);
                break;
            }
            case TypeNoeud::IntegraleNonEvaluee:
            case TypeNoeud::LimiteNonEvaluee:
                m_ok = false;
                break;
            default:
                break;
        }
    }

    void puissance(const ExprPtr& base, const ExprPtr& exposant) {
        double q;
        if (!valeurConstante(exposant, q)) {
            // Exposant variable ou paramétré : la base doit être strictement positive
            exiger(base, Condition::Positif);
            return;
        }
        const bool entier = q == std::floor(q) && std::abs(q) < 1e15;
        if (entier) {
            if (q < 0.0) exiger(base, Condition::NonNul);
        } else {
            exiger(base, q > 0.0 ? Condition::NonNegatif : Condition::Positif);
        }
    }

    // Racines de u dans [A, B] (mémorisées), ou « non fiable »
    const Racines& racines(const ExprPtr& u) {
        auto [it, nouveau] = m_memo.try_emplace(u.get());
        if (!nouveau) return it->second;
        Racines& r = it->second;
        if (!m_bornesConnues || !(m_a < m_b)) return r;
        try {
            const Solutions s = resoudreSurIntervalle(u, m_a, m_b);
            if (!s.complet || s.toutReel) return r;
            for (const Solution& sol : s.liste) {
                if (!std::isfinite(sol.approximation)) return r;
                r.valeurs.push_back(sol.approximation);
            }
            std::sort(r.valeurs.begin(), r.valeurs.end());
            r.fiable = true;
        } catch (const std::exception&) {
            r.fiable = false;
        }
        return r;
    }

    bool interieur(double r) const { return r > m_a + m_epsilon && r < m_b - m_epsilon; }

    void exiger(const ExprPtr& u, Condition condition) {
        if (!verifier(u, condition)) m_ok = false;
    }

    bool verifier(const ExprPtr& u, Condition condition) {
        if (!u->contientVariable()) {
            double v;
            if (!valeurConstante(u, v)) return false;
            return condition == Condition::NonNul ? v != 0.0 : condition == Condition::Positif ? v > 0.0 : v >= 0.0;
        }
        if (u->plusieursVariables()) return false;
        const Racines& z = racines(u);
        if (!z.fiable) return false;

        // Racines qui comptent : toutes en mode fermé, les intérieures en mode ouvert ou pour u >= 0
        std::vector<double> dedans;
        bool surBorne = false;
        for (double r : z.valeurs) {
            if (interieur(r)) dedans.push_back(r);
            else surBorne = true;
        }

        bool ok;
        if (condition == Condition::NonNegatif) {
            // Signe >= 0 sur chaque morceau entre racines (une racine de contact est permise)
            ok = true;
            double debut = m_a;
            dedans.push_back(m_b);
            for (double fin : dedans) {
                if (fin > debut) {
                    const double v = evaluerSur(u, 0.5 * (debut + fin));
                    if (!(v >= 0.0)) ok = false;
                }
                debut = fin;
            }
            dedans.pop_back();
        } else {
            ok = dedans.empty() && !(m_ferme && surBorne);
            if (ok && condition == Condition::Positif) ok = evaluerSur(u, 0.5 * (m_a + m_b)) > 0.0;
        }
        if (!ok) m_coupures.insert(m_coupures.end(), dedans.begin(), dedans.end());
        return ok;
    }

    double m_a, m_b, m_epsilon;
    bool m_ferme, m_bornesConnues;
    bool m_ok = true;
    std::unordered_map<const ASTNode*, Racines> m_memo;
    std::vector<double> m_coupures;
};

// ============================================================================
// Quadrature de Gauss-Kronrod G7-K15 adaptative
// ============================================================================

// Noeuds de Kronrod sur [0, 1] (décroissants) et poids K15 / G7 (QUADPACK qk15)
constexpr double XGK[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                           0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                           0.586087235467691130294144838258730, 0.405845151377397166906606412076961,
                           0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
constexpr double WGK[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                           0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                           0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                           0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
constexpr double WG[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                          0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

constexpr int NOEUDS = 15;
constexpr int VAGUE_MAX = 4;                 // intervalles bissectés ensemble (un seul appel d'évaluation)
constexpr int LOT_MAX = 2 * VAGUE_MAX;       // enfants évalués ensemble
constexpr int PROFONDEUR_RESIDU = 46;        // au-delà, une évaluation non finie est un résidu d'arrondi
constexpr int STAGNATION_DIVERGENCE = 30;    // bissections consécutives sans décroissance : divergence
constexpr double RATIO_STAGNATION = 0.97;

// Domaine d'intégration réel, ramené à un intervalle fini par changement de variable
enum class Domaine {
    Fini,          // x = t
    DroiteInfinie, // [origine, +inf[ : x = origine + t/(1-t), t dans [0, 1[
    GaucheInfinie, // ]-inf, origine] : x = origine - t/(1-t), t dans [0, 1[
    Reels          // ]-inf, +inf[ : x = t/(1-t^2), t dans ]-1, 1[
};

struct Intervalle {
    double a = 0.0, b = 0.0;
    double integrale = 0.0;
    double erreur = 0.0;
    double plancher = 0.0; // erreur minimale atteignable (arrondis)
    int profondeur = 0;
    int stagnation = 0;
    bool nonFini = false;
    bool toutNonFini = false;
};

struct SortieQuadrature {
    double valeur = NON_NOMBRE;
    double erreur = INFINI;
    bool diverge = false;
};

class Quadrature {
public:
    Quadrature(const ProgrammeEvaluation& programme, Domaine domaine, double origine, const OptionsIntegration& options)
        : m_programme(programme), m_domaine(domaine), m_origine(origine), m_options(options) {
        m_tas.reserve(static_cast<std::size_t>(std::max(options.subdivisionsMax, 1)) + LOT_MAX + 2);
    }

    // tolAbs : tolérance absolue de ce morceau
    SortieQuadrature executer(double a, double b, double tolAbs) {
        m_tas.clear();
        Intervalle depart;
        depart.a = a;
        depart.b = b;
        evaluerLot(&depart, 1);
        if (depart.toutNonFini) return {};
        m_tas.push_back(depart);

        int bissections = 0;
        double residus = 0.0; // erreur des zones trop fines pour être évaluées (voir PROFONDEUR_RESIDU)
        const int maxBissections = std::max(m_options.subdivisionsMax, 1);
        Intervalle parents[VAGUE_MAX];
        Intervalle enfants[LOT_MAX];
        for (;;) {
            double integrale = 0.0, erreur = 0.0;
            bool nonFini = false;
            for (const Intervalle& iv : m_tas) {
                if (iv.nonFini) {
                    nonFini = true;
                } else {
                    integrale += iv.integrale;
                }
                erreur += iv.erreur;
            }
            const double tolerance = std::max(tolAbs, m_options.toleranceRelative * std::abs(integrale));
            if (!nonFini && erreur <= tolerance) return {integrale, erreur + residus, false};
            erreur += residus;
            if (bissections >= maxBissections) return finale(integrale, erreur, nonFini);

            const Intervalle& pire = m_tas.front();
            if (!pire.nonFini && pire.erreur <= pire.plancher * (1.0 + 1e-9)) return finale(integrale, erreur, nonFini);
            if (pire.b - pire.a <= 8.0 * EPSILON * std::max(std::abs(pire.a), std::abs(pire.b))) {
                return finale(integrale, erreur, nonFini); // plus de flottants entre les bornes
            }

            // Bissection des intervalles d'erreur supérieure à leur part de la tolérance
            const double part = tolerance / static_cast<double>(m_tas.size());
            int vague = 0;
            const int vagueMax = std::min(VAGUE_MAX, maxBissections - bissections);
            while (vague < vagueMax && !m_tas.empty()) {
                const Intervalle& haut = m_tas.front();
                if (vague > 0 && !haut.nonFini && haut.erreur <= part) break;
                std::pop_heap(m_tas.begin(), m_tas.end(), plusGrandeErreur);
                parents[vague] = m_tas.back();
                m_tas.pop_back();
                const double milieu = 0.5 * (parents[vague].a + parents[vague].b);
                enfants[2 * vague].a = parents[vague].a;
                enfants[2 * vague].b = milieu;
                enfants[2 * vague + 1].a = milieu;
                enfants[2 * vague + 1].b = parents[vague].b;
                ++vague;
            }
            evaluerLot(enfants, 2 * vague);
            bissections += vague;
            for (int i = 0; i < 2 * vague; ++i) {
                Intervalle& e = enfants[i];
                const Intervalle& p = parents[i / 2];
                e.profondeur = p.profondeur + 1;
                if (e.nonFini && e.profondeur > PROFONDEUR_RESIDU) {
                    // Valeur non finie en une zone infiniment fine (extrémité d'une
                    // singularité intégrable, aux limites des flottants) : contribution
                    // inconnue, majorée par celle du parent
                    residus += std::abs(p.integrale);
                    continue;
                }
                if (e.toutNonFini) return {};
                const double ie = std::abs(e.integrale), ip = std::abs(p.integrale);
                e.stagnation = (!e.nonFini && ie >= RATIO_STAGNATION * ip && ie > tolAbs) ? p.stagnation + 1 : 0;
                if (e.stagnation >= STAGNATION_DIVERGENCE) return {NON_NOMBRE, INFINI, true};
                m_tas.push_back(e);
                std::push_heap(m_tas.begin(), m_tas.end(), plusGrandeErreur);
            }
        }
    }

private:
    static bool plusGrandeErreur(const Intervalle& x, const Intervalle& y) { return x.erreur < y.erreur; }

    static SortieQuadrature finale(double integrale, double erreur, bool nonFini) {
        if (nonFini) return {};
        return {integrale, erreur, false};
    }

    // Évalue le K15 de n intervalles (a et b renseignés) en un seul appel groupé
    void evaluerLot(Intervalle* iv, int n) {
        const int total = NOEUDS * n;
        for (int i = 0; i < n; ++i) {
            const double c = 0.5 * (iv[i].a + iv[i].b), h = 0.5 * (iv[i].b - iv[i].a);
            double* t = m_t.data() + NOEUDS * i;
            for (int j = 0; j < 7; ++j) {
                t[j] = c - h * XGK[j];
                t[NOEUDS - 1 - j] = c + h * XGK[j];
            }
            t[7] = c;
        }
        const double* x = m_t.data();
        if (m_domaine != Domaine::Fini) {
            for (int k = 0; k < total; ++k) {
                const double t = m_t[k];
                if (m_domaine == Domaine::Reels) {
                    const double u = 1.0 - t * t;
                    m_x[k] = t / u;
                    m_jac[k] = (1.0 + t * t) / (u * u);
                } else {
                    const double u = 1.0 - t;
                    const double s = t / u;
                    m_x[k] = m_domaine == Domaine::DroiteInfinie ? m_origine + s : m_origine - s;
                    m_jac[k] = 1.0 / (u * u);
                }
            }
            x = m_x.data();
        }
        m_programme.evaluer(x, m_y.data(), static_cast<std::size_t>(total));
        if (m_domaine != Domaine::Fini) {
            for (int k = 0; k < total; ++k) m_y[k] = m_y[k] == 0.0 ? 0.0 : m_y[k] * m_jac[k];
        }
        for (int i = 0; i < n; ++i) kronrod(m_y.data() + NOEUDS * i, iv[i]);
    }

    static void kronrod(const double* v, Intervalle& iv) {
        int fini = 0;
        for (int k = 0; k < NOEUDS; ++k) fini += std::isfinite(v[k]) ? 1 : 0;
        if (fini < NOEUDS) {
            iv.nonFini = true;
            iv.toutNonFini = fini == 0;
            iv.integrale = 0.0;
            iv.erreur = INFINI;
            iv.plancher = 0.0;
            return;
        }
        iv.nonFini = false;
        iv.toutNonFini = false;
        const double h = 0.5 * (iv.b - iv.a);
        const double fc = v[7];
        double resg = WG[3] * fc, resk = WGK[7] * fc, resabs = WGK[7] * std::abs(fc);
        for (int j = 0; j < 3; ++j) { // noeuds de Gauss : XGK[1], XGK[3], XGK[5]
            const int m = 2 * j + 1;
            const double s = v[m] + v[NOEUDS - 1 - m];
            resg += WG[j] * s;
            resk += WGK[m] * s;
            resabs += WGK[m] * (std::abs(v[m]) + std::abs(v[NOEUDS - 1 - m]));
        }
        for (int j = 0; j < 4; ++j) { // noeuds propres à Kronrod : XGK[0], XGK[2], XGK[4], XGK[6]
            const int m = 2 * j;
            resk += WGK[m] * (v[m] + v[NOEUDS - 1 - m]);
            resabs += WGK[m] * (std::abs(v[m]) + std::abs(v[NOEUDS - 1 - m]));
        }
        const double moyenne = resk * 0.5;
        double resasc = WGK[7] * std::abs(fc - moyenne);
        for (int j = 0; j < 7; ++j) resasc += WGK[j] * (std::abs(v[j] - moyenne) + std::abs(v[NOEUDS - 1 - j] - moyenne));

        const double ah = std::abs(h);
        resabs *= ah;
        resasc *= ah;
        double erreur = std::abs((resk - resg) * h);
        if (resasc != 0.0 && erreur != 0.0) erreur = resasc * std::min(1.0, std::pow(200.0 * erreur / resasc, 1.5));
        iv.plancher = 50.0 * EPSILON * resabs;
        iv.erreur = std::max(iv.plancher, erreur);
        iv.integrale = resk * h;
    }

    const ProgrammeEvaluation& m_programme;
    Domaine m_domaine;
    double m_origine;
    OptionsIntegration m_options;
    std::vector<Intervalle> m_tas; // tas-max sur l'erreur
    std::array<double, NOEUDS * LOT_MAX> m_t{}, m_x{}, m_y{}, m_jac{};
};

// Intégrale numérique sur [A, B] (A < B, infinis permis), coupée aux points donnés
SortieQuadrature integrerNumerique(const ProgrammeEvaluation& programme, double A, double B,
                                   const std::vector<double>& coupures, const OptionsIntegration& options) {
    if (std::isinf(A) || std::isinf(B)) {
        Domaine domaine = Domaine::Reels;
        double origine = 0.0, a = -1.0, b = 1.0;
        if (std::isinf(A) && std::isinf(B)) {
            // x = t / (1 - t^2)
        } else {
            domaine = std::isinf(B) ? Domaine::DroiteInfinie : Domaine::GaucheInfinie;
            origine = std::isinf(B) ? A : B;
            a = 0.0;
            b = 1.0;
        }
        Quadrature q(programme, domaine, origine, options);
        return q.executer(a, b, options.toleranceAbsolue);
    }

    // Points de coupure intérieurs distincts
    std::vector<double> bornes;
    bornes.reserve(coupures.size() + 2);
    bornes.push_back(A);
    std::vector<double> tries(coupures);
    std::sort(tries.begin(), tries.end());
    for (double c : tries) {
        if (c > bornes.back() && c < B) bornes.push_back(c);
    }
    bornes.push_back(B);

    const double tolAbs = options.toleranceAbsolue / static_cast<double>(bornes.size() - 1);
    Quadrature q(programme, Domaine::Fini, 0.0, options);
    SortieQuadrature total{0.0, 0.0, false};
    for (std::size_t i = 0; i + 1 < bornes.size(); ++i) {
        const SortieQuadrature s = q.executer(bornes[i], bornes[i + 1], tolAbs);
        if (std::isnan(s.valeur)) return s;
        total.valeur += s.valeur;
        total.erreur += s.erreur;
    }
    return total;
}

// ============================================================================
// Chemin exact
// ============================================================================

// F(b) - F(a) quand F et f sont certifiées sans singularité sur [a, b] ; nul sinon
ExprPtr integraleExacte(const ExprPtr& f, const ExprPtr& a, const ExprPtr& b, double A, double B, bool bornesConnues,
                        bool fOk) {
    if (!fOk) return ExprPtr();
    const ExprPtr F = f->integrer();
    if (contientNonEvalue(F)) return ExprPtr();
    Analyseur ferme(A, B, true, bornesConnues);
    if (!ferme.sur(F)) return ExprPtr();

    try {
        const std::vector<ExprPtr> vars = variables(F);
        if (vars.size() > 1) return ExprPtr();
        ExprPtr resultat = vars.empty() ? cst(0.0) : (substituer(F, vars[0], b) - substituer(F, vars[0], a));
        resultat = resultat->simplifier();
        if (bornesConnues) {
            // Une valeur non finie (ln d'un nombre négatif...) invalide le résultat ; un
            // paramètre empêche seulement de l'évaluer
            double v = 0.0;
            try {
                v = resultat->eval(0.0);
            } catch (const std::logic_error&) {
            }
            if (!std::isfinite(v)) return ExprPtr();
        }
        return resultat;
    } catch (const std::exception&) {
        return ExprPtr();
    }
}

ResultatIntegrale oppose(ResultatIntegrale r) {
    if (r.exacte) r.exacte = (-r.exacte)->simplifier();
    r.valeur = -r.valeur;
    return r;
}

// a et b : bornes symboliques (nulles si seules des bornes double sans chemin exact sont données)
ResultatIntegrale calculer(const ExprPtr& f, ExprPtr a, ExprPtr b, double A, double B, bool bornesConnues,
                           const OptionsIntegration& options) {
    if (!f) throw std::invalid_argument("integrerDefinie : expression nulle");
    if (f->plusieursVariables()) {
        throw std::invalid_argument("integrerDefinie : l'integrale definie ne s'applique qu'aux expressions a une variable");
    }
    if ((a && a->contientVariable()) || (b && b->contientVariable())) {
        throw std::invalid_argument("integrerDefinie : les bornes doivent etre constantes");
    }
    if (bornesConnues && (std::isnan(A) || std::isnan(B))) throw std::invalid_argument("integrerDefinie : borne NaN");

    ResultatIntegrale r;
    if ((a && a == b) || (bornesConnues && !a && A == B)) {
        r.exacte = cst(0.0);
        r.valeur = 0.0;
        return r;
    }
    bool inverse = false;
    if (bornesConnues && A > B) {
        std::swap(A, B);
        std::swap(a, b);
        inverse = true;
    }

    if (bornesConnues && A == B) { // bornes distinctes mais de même valeur
        r.valeur = 0.0;
        return r;
    }

    // Singularités de f sur ]A, B[ (utiles au chemin exact et aux coupures numériques)
    const bool fini = bornesConnues && std::isfinite(A) && std::isfinite(B);
    Analyseur ouvert(A, B, false, fini);
    const bool fOk = ouvert.sur(f);

    if (options.exacteSiPossible && a && b && (fini || !bornesConnues)) {
        if (ExprPtr exacte = integraleExacte(f, a, b, A, B, fini, fOk)) {
            r.exacte = exacte;
            r.valeur = NON_NOMBRE;
            if (!valeurConstante(exacte, r.valeur)) r.valeur = NON_NOMBRE;
            r.erreurEstimee = 0.0;
            return inverse ? oppose(r) : r;
        }
    }

    if (!bornesConnues) {
        throw std::domain_error("integrerDefinie : bornes symboliques sans primitive exploitable, calcul numerique impossible");
    }
    const ProgrammeEvaluation programme(f);
    if (!programme.estValide()) {
        throw std::domain_error("integrerDefinie : integrale non calculable en exact et expression non evaluable "
                                "numeriquement (parametre ou noeud non evalue)");
    }
    const SortieQuadrature s = integrerNumerique(programme, A, B, ouvert.coupures(), options);
    r.valeur = inverse ? -s.valeur : s.valeur;
    r.erreurEstimee = s.erreur;
    return r;
}

} // namespace

/*
 * Nom : integrerDefinie
 * Description : Voir Integration.hpp.
 */
ResultatIntegrale integrerDefinie(const ExprPtr& f, const ExprPtr& a, const ExprPtr& b, const OptionsIntegration& options) {
    if (!a || !b) throw std::invalid_argument("integrerDefinie : borne nulle");
    double A = 0.0, B = 0.0;
    const bool connues = valeurConstante(a, A) && valeurConstante(b, B);
    if (!connues && (a->contientVariable() || b->contientVariable())) {
        throw std::invalid_argument("integrerDefinie : les bornes doivent etre constantes");
    }
    return calculer(f, a, b, A, B, connues, options);
}

ResultatIntegrale integrerDefinie(const ExprPtr& f, double a, double b, const OptionsIntegration& options) {
    // Chemin exact seulement pour des bornes entières (cst(0.5) serait un réel, pas 1/2)
    auto entier = [](double v) { return v == std::floor(v) && std::abs(v) < 9e15; };
    ExprPtr ea, eb;
    if (entier(a) && entier(b)) {
        ea = cst(a);
        eb = cst(b);
    }
    return calculer(f, ea, eb, a, b, true, options);
}

} // namespace symalgo
