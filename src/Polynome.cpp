/**
 * @file Polynome.cpp
 * @brief Développement des expressions et polynômes à coefficients exacts.
 */

#include "Polynome.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
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
            case TypeNoeud::ArcSinus: return ast_asin(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::ArcCosinus: return ast_acos(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            case TypeNoeud::ArcTangente: return ast_atan(developper(static_cast<const FonctionUnaire&>(*e).m_argument));
            default: return e; // constantes, variable, paramètres, noeuds non évalués
        }
    }

    std::unordered_map<const ASTNode*, ExprPtr> m_memo;
};

} // namespace

ExprPtr developper(const ExprPtr& e) { return Developpeur().developper(e); }

// ============================================================================
// Polynôme à coefficients exacts
// ============================================================================

Polynome::Polynome(std::vector<Nombre> coefficients) : m_coefficients(std::move(coefficients)) { normaliser(); }

void Polynome::normaliser() {
    while (!m_coefficients.empty() && m_coefficients.back().estZero()) m_coefficients.pop_back();
}

const Nombre& Polynome::coefficient(int i) const {
    static const Nombre ZERO(0);
    return i >= 0 && i <= degre() ? m_coefficients[static_cast<std::size_t>(i)] : ZERO;
}

namespace {

// Coefficient exact : les réels sont convertis en leur valeur exacte (rationnel dyadique)
bool coefficientExact(const Nombre& n, Nombre& exact, bool& toutExact) {
    if (n.estExact()) {
        exact = n;
        return true;
    }
    if (!n.estFini()) return false;
    exact = Nombre::exactDepuisDouble(n.versDouble());
    toutExact = false;
    return true;
}

// Monôme x^k (k entier positif) : écrit k et la variable
bool monomeEnX(const ExprPtr& e, long long& k, ExprPtr& variable) {
    if (e->type() == TypeNoeud::Variable) {
        k = 1;
        variable = e;
        return true;
    }
    if (const Puissance* p = comme<Puissance>(e)) {
        const Constante* c = comme<Constante>(p->getExposant());
        if (p->getBase()->type() == TypeNoeud::Variable && c && c->getNombre().versEntier(k) && k > 0) {
            variable = p->getBase();
            return true;
        }
    }
    return false;
}

} // namespace

bool Polynome::depuisExpression(const ExprPtr& e, Polynome& p, ExprPtr* variable) {
    std::vector<Nombre> c;
    ExprPtr x;
    bool toutExact = true;
    auto ajouter = [&](long long k, const Nombre& coef) {
        Nombre exact;
        if (!coefficientExact(coef, exact, toutExact)) return false;
        if (c.size() <= static_cast<std::size_t>(k)) c.resize(static_cast<std::size_t>(k) + 1, Nombre(0));
        c[static_cast<std::size_t>(k)] = c[static_cast<std::size_t>(k)] + exact;
        return true;
    };
    auto terme = [&](const ExprPtr& t, const Nombre& coef) {
        long long k;
        ExprPtr v;
        if (!monomeEnX(t, k, v)) return false;
        if (x && x.get() != v.get()) return false; // une seule variable
        x = v;
        return ajouter(k, coef);
    };
    switch (e->type()) {
        case TypeNoeud::Constante:
            if (!ajouter(0, static_cast<const Constante&>(*e).getNombre())) return false;
            break;
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(*e);
            if (!ajouter(0, s.getConstante())) return false;
            for (const Terme& t : s.getTermes()) {
                if (!terme(t.expression, t.coefficient)) return false;
            }
            break;
        }
        case TypeNoeud::Produit: {
            const Produit& pr = static_cast<const Produit&>(*e);
            if (pr.getFacteurs().size() != 1) return false;
            if (!terme(ast_pow(pr.getFacteurs()[0].base, pr.getFacteurs()[0].exposant), pr.getCoefficient())) return false;
            break;
        }
        default:
            if (!terme(e, Nombre(1))) return false;
    }
    p = Polynome(std::move(c));
    p.m_exact = toutExact;
    if (variable) *variable = x;
    return true;
}

Polynome Polynome::operator+(const Polynome& b) const {
    std::vector<Nombre> c(std::max(m_coefficients.size(), b.m_coefficients.size()), Nombre(0));
    for (std::size_t i = 0; i < c.size(); ++i) c[i] = coefficient(static_cast<int>(i)) + b.coefficient(static_cast<int>(i));
    return Polynome(std::move(c));
}

Polynome Polynome::operator-(const Polynome& b) const { return *this + b * Nombre(-1); }

Polynome Polynome::operator*(const Nombre& k) const {
    std::vector<Nombre> c = m_coefficients;
    for (Nombre& ci : c) ci = ci * k;
    return Polynome(std::move(c));
}

Polynome Polynome::operator*(const Polynome& b) const {
    if (estNul() || b.estNul()) return Polynome();
    std::vector<Nombre> c(m_coefficients.size() + b.m_coefficients.size() - 1, Nombre(0));
    for (std::size_t i = 0; i < m_coefficients.size(); ++i)
        for (std::size_t j = 0; j < b.m_coefficients.size(); ++j) c[i + j] = c[i + j] + m_coefficients[i] * b.m_coefficients[j];
    return Polynome(std::move(c));
}

void Polynome::diviser(const Polynome& d, Polynome& quotient, Polynome& reste) const {
    if (d.estNul()) throw std::domain_error("Polynome : division par le polynome nul");
    std::vector<Nombre> r = m_coefficients;
    const int n = degre(), m = d.degre();
    std::vector<Nombre> q(n >= m ? static_cast<std::size_t>(n - m + 1) : 0, Nombre(0));
    const Nombre inverseDominant = d.dominant().inverse();
    for (int k = n; k >= m; --k) {
        const Nombre facteur = r[static_cast<std::size_t>(k)] * inverseDominant;
        q[static_cast<std::size_t>(k - m)] = facteur;
        if (facteur.estZero()) continue;
        for (int j = 0; j <= m; ++j) {
            Nombre& rj = r[static_cast<std::size_t>(k - m + j)];
            rj = rj - facteur * d.m_coefficients[static_cast<std::size_t>(j)];
        }
    }
    quotient = Polynome(std::move(q));
    r.resize(static_cast<std::size_t>(std::max(0, std::min(n + 1, m))));
    reste = Polynome(std::move(r));
}

Polynome Polynome::derivee() const {
    std::vector<Nombre> c;
    for (std::size_t i = 1; i < m_coefficients.size(); ++i) c.push_back(m_coefficients[i] * Nombre(static_cast<long long>(i)));
    return Polynome(std::move(c));
}

Polynome Polynome::unitaire() const { return estNul() ? *this : *this * dominant().inverse(); }

Polynome Polynome::pgcd(Polynome a, Polynome b) {
    while (!b.estNul()) {
        Polynome q, r;
        a.diviser(b, q, r);
        a = std::move(b);
        b = std::move(r);
    }
    return a.unitaire();
}

std::vector<std::pair<Polynome, int>> Polynome::sansCarre() const {
    // Algorithme de Yun : f = prod f_i^i
    std::vector<std::pair<Polynome, int>> facteurs;
    if (degre() < 1) return facteurs;
    const Polynome f = unitaire();
    const Polynome fp = f.derivee();
    const Polynome a0 = pgcd(f, fp);
    Polynome q, r, b, c;
    f.diviser(a0, b, r);
    fp.diviser(a0, c, r);
    Polynome d = c - b.derivee();
    for (int i = 1; b.degre() > 0; ++i) {
        const Polynome a = pgcd(b, d);
        if (a.degre() > 0) facteurs.emplace_back(a, i);
        b.diviser(a, q, r);
        b = q;
        d.diviser(a, c, r);
        d = c - b.derivee();
    }
    return facteurs;
}

Nombre Polynome::evaluer(const Nombre& x) const {
    Nombre r(0);
    for (auto it = m_coefficients.rbegin(); it != m_coefficients.rend(); ++it) r = r * x + *it;
    return r;
}

double Polynome::evaluer(double x) const {
    double r = 0.0;
    for (auto it = m_coefficients.rbegin(); it != m_coefficients.rend(); ++it) r = r * x + it->versDouble();
    return r;
}

ExprPtr Polynome::versExpression(const ExprPtr& x) const {
    AccumulateurSomme s;
    for (std::size_t i = 0; i < m_coefficients.size(); ++i) {
        if (!m_coefficients[i].estZero()) s.ajouter(ast_pow(x, static_cast<double>(i)), m_coefficients[i]);
    }
    return s.construire();
}

// ============================================================================
// Racines réelles (suites de Sturm)
// ============================================================================

namespace {

std::vector<Polynome> suiteDeSturm(const Polynome& f) {
    std::vector<Polynome> suite = {f, f.derivee()};
    while (suite.back().degre() > 0) {
        Polynome q, r;
        suite[suite.size() - 2].diviser(suite.back(), q, r);
        if (r.estNul()) break;
        suite.push_back(r * Nombre(-1));
    }
    return suite;
}

int variations(const std::vector<Polynome>& suite, const Nombre& x) {
    int n = 0, precedent = 0;
    for (const Polynome& p : suite) {
        const int s = p.evaluer(x).signe();
        if (s == 0) continue;
        if (precedent != 0 && s != precedent) ++n;
        precedent = s;
    }
    return n;
}

// Borne de Cauchy : toutes les racines sont dans ]-B, B]
Nombre borneDeCauchy(const Polynome& f) {
    Nombre m(0);
    const Nombre inverse = f.dominant().inverse();
    for (int i = 0; i < f.degre(); ++i) {
        const Nombre v = (f.coefficient(i) * inverse).valeurAbsolue();
        if (v > m) m = v;
    }
    return m + Nombre(1);
}

// Diviseurs positifs de |n| (n non nul, |n| <= 10^12), sinon liste vide
std::vector<long long> diviseurs(const Nombre& n) {
    long long v;
    if (!n.versEntier(v) || v == 0) return {};
    v = v < 0 ? -v : v;
    if (v > 1000000000000LL) return {};
    std::vector<long long> petits, grands;
    for (long long d = 1; d * d <= v; ++d) {
        if (v % d == 0) {
            petits.push_back(d);
            if (d != v / d) grands.push_back(v / d);
        }
    }
    petits.insert(petits.end(), grands.rbegin(), grands.rend());
    return petits;
}

// Racine isolée d'un facteur sans carré, affinée jusqu'à la précision du double
struct RacineIsolee {
    Nombre gauche, droite; // seule racine dans ]gauche, droite]
    bool rationnelle = false;
    Nombre valeurExacte;
    double valeur = 0.0;
};

// Dénominateurs possibles des racines rationnelles : diviseurs du coefficient dominant
// de la forme entière primitive (théorème des racines rationnelles)
std::vector<long long> denominateursPossibles(const Polynome& f) {
    Nombre ppcm(1);
    for (int i = 0; i <= f.degre(); ++i) {
        const Nombre d = f.coefficient(i).denominateurNombre();
        ppcm = ppcm * d / Nombre::pgcd(ppcm, d);
    }
    return diviseurs(f.dominant() * ppcm);
}

// Cherche une racine rationnelle p/q dans ]a, b] près de la valeur approchée v
bool chercherRationnelle(const Polynome& f, const std::vector<long long>& denominateurs, const Nombre& a,
                         const Nombre& b, double v, Nombre& racine) {
    for (long long q : denominateurs) {
        const double p = std::round(v * static_cast<double>(q));
        if (std::abs(p) > 9.0e18) continue;
        const Nombre candidat = Nombre::rationnel(static_cast<long long>(p), q);
        if (a < candidat && candidat <= b && f.evaluer(candidat).estZero()) {
            racine = candidat;
            return true;
        }
    }
    return false;
}

bool intervalleEtroit(const Nombre& a, const Nombre& b, double largeurRelative) {
    const double da = a.versDouble(), db = b.versDouble();
    return db - da <= largeurRelative * std::max(1.0, std::max(std::abs(da), std::abs(db)));
}

/*
 * Nom : isolerRacines
 * Description : Isole les racines réelles d'un polynôme sans carré par dichotomie guidée
 *               par le théorème de Sturm (arithmétique exacte), puis les affine jusqu'à la
 *               précision du double. Les racines rationnelles sont reconnues dès que
 *               l'intervalle est assez étroit.
 */
std::vector<RacineIsolee> isolerRacines(const Polynome& f) {
    std::vector<RacineIsolee> racines;
    if (f.degre() < 1) return racines;
    const std::vector<Polynome> sturm = suiteDeSturm(f);
    const std::vector<long long> denominateurs = denominateursPossibles(f);
    const Nombre B = borneDeCauchy(f);
    const Nombre demi = Nombre::rationnel(1, 2);
    auto compter = [&](const Nombre& a, const Nombre& b) { return variations(sturm, a) - variations(sturm, b); };

    std::vector<std::pair<Nombre, Nombre>> pile = {{-B, B}};
    while (!pile.empty()) {
        auto [a, b] = pile.back();
        pile.pop_back();
        const int n = compter(a, b);
        if (n == 0) continue;
        if (n > 1) {
            const Nombre m = (a + b) * demi;
            pile.push_back({m, b});
            pile.push_back({a, m});
            continue;
        }
        RacineIsolee r;
        bool rationnelleTentee = false;
        // Affinage : une seule racine, simple, dans ]a, b]. Si f(a) != 0, le signe de f
        // suffit ; sinon (a est la racine de l'intervalle voisin) on compte par Sturm
        int signeA = f.evaluer(a).signe();
        for (int iteration = 0; iteration < 1200; ++iteration) { // 1200 > bits d'exposant + mantisse
            if (f.evaluer(b).estZero()) {
                r.rationnelle = true;
                r.valeurExacte = b;
                break;
            }
            if (!rationnelleTentee && intervalleEtroit(a, b, 1e-6)) {
                rationnelleTentee = true;
                if (chercherRationnelle(f, denominateurs, a, b, ((a + b) * demi).versDouble(), r.valeurExacte)) {
                    r.rationnelle = true;
                    break;
                }
            }
            // Arrêt quand les bornes sont des doubles voisins : précision maximale à toute échelle
            const double da = a.versDouble(), db = b.versDouble();
            if (std::nextafter(da, std::numeric_limits<double>::infinity()) >= db) break;
            const Nombre m = (a + b) * demi;
            bool racineAGauche;
            if (signeA != 0) {
                const int signeM = f.evaluer(m).signe();
                racineAGauche = signeM == 0 || signeM != signeA; // racine dans ]a, m]
            } else {
                racineAGauche = compter(a, m) == 1;
            }
            if (racineAGauche) {
                b = m;
            } else {
                a = m;
                signeA = f.evaluer(a).signe();
            }
        }
        if (!r.rationnelle) {
            r.rationnelle = chercherRationnelle(f, denominateurs, a, b, ((a + b) * demi).versDouble(), r.valeurExacte);
        }
        r.gauche = a;
        r.droite = b;
        if (r.rationnelle) {
            r.valeur = r.valeurExacte.versDouble();
        } else {
            // Le double le plus proche de la racine : celui des deux bornes dont f est le plus petit
            const double da = a.versDouble(), db = b.versDouble();
            r.valeur = f.evaluer(Nombre::exactDepuisDouble(da)).valeurAbsolue() <=
                               f.evaluer(Nombre::exactDepuisDouble(db)).valeurAbsolue()
                           ? da
                           : db;
        }
        racines.push_back(r);
    }
    return racines;
}

} // namespace

int Polynome::nombreRacinesReelles(const Nombre& a, const Nombre& b) const {
    if (degre() < 1) return 0;
    // Le théorème de Sturm compte les racines distinctes, même en présence de racines multiples
    const std::vector<Polynome> sturm = suiteDeSturm(*this);
    return variations(sturm, a) - variations(sturm, b);
}

std::vector<RacineReelle> Polynome::racinesReelles() const {
    std::vector<RacineReelle> resultat;
    for (const auto& [facteur, multiplicite] : sansCarre()) {
        std::vector<RacineIsolee> racines = isolerRacines(facteur);
        // Le reste sans les racines rationnelles : exact par radicaux s'il est de degré 2
        Polynome reste = facteur;
        for (const RacineIsolee& r : racines) {
            if (!r.rationnelle) continue;
            Polynome q, rr;
            reste.diviser(Polynome({-r.valeurExacte, Nombre(1)}), q, rr);
            reste = q;
        }
        std::vector<ExprPtr> quadratiques;
        if (reste.degre() == 2) {
            const Nombre& a = reste.coefficient(2);
            const Nombre& b = reste.coefficient(1);
            const Nombre& c = reste.coefficient(0);
            const Nombre delta = b * b - Nombre(4) * a * c;
            if (delta.signe() > 0) {
                const ExprPtr racineDelta = ast_pow(nombre(delta), frac(1, 2));
                const ExprPtr deuxA = nombre(Nombre(2) * a);
                // Ordre croissant pour correspondre aux racines isolées
                const ExprPtr r1 = (nombre(-b) - racineDelta) / deuxA;
                const ExprPtr r2 = (nombre(-b) + racineDelta) / deuxA;
                quadratiques = r1->eval(0.0) < r2->eval(0.0) ? std::vector<ExprPtr>{r1, r2} : std::vector<ExprPtr>{r2, r1};
            }
        }
        std::sort(racines.begin(), racines.end(), [](const RacineIsolee& x, const RacineIsolee& y) { return x.valeur < y.valeur; });
        std::size_t iq = 0;
        for (const RacineIsolee& r : racines) {
            RacineReelle rr;
            rr.valeur = r.valeur;
            rr.multiplicite = multiplicite;
            rr.gauche = r.gauche;
            rr.droite = r.droite;
            if (!m_exact) {
                // Coefficients approchés : valeur certifiée, mais pas de forme exacte
            } else if (r.rationnelle) {
                rr.exacte = nombre(r.valeurExacte);
            } else if (iq < quadratiques.size()) {
                rr.exacte = quadratiques[iq++];
            }
            resultat.push_back(rr);
        }
    }
    std::sort(resultat.begin(), resultat.end(), [](const RacineReelle& x, const RacineReelle& y) { return x.valeur < y.valeur; });
    return resultat;
}

// ============================================================================
// Factorisation
// ============================================================================

namespace {

// Forme entière primitive d'un polynôme : f = contenu * primitive
Polynome partiePrimitive(const Polynome& f, Nombre& contenu) {
    Nombre ppcm(1), pgcdNum(0);
    for (int i = 0; i <= f.degre(); ++i) {
        const Nombre d = f.coefficient(i).denominateurNombre();
        ppcm = ppcm * d / Nombre::pgcd(ppcm, d);
    }
    for (int i = 0; i <= f.degre(); ++i) pgcdNum = Nombre::pgcd(pgcdNum, (f.coefficient(i) * ppcm).numerateurNombre());
    contenu = pgcdNum / ppcm;
    if (f.dominant().signe() < 0) contenu = -contenu; // coefficient dominant positif
    return f * contenu.inverse();
}

} // namespace

ExprPtr factoriser(const ExprPtr& e) {
    const ExprPtr d = developper(e);
    Polynome p;
    ExprPtr x;
    if (!Polynome::depuisExpression(d, p, &x) || !x || p.degre() < 2 || !p.estExact()) return e;

    AccumulateurProduit produit_;
    produit_.multiplier(nombre(p.dominant())); // les facteurs de sansCarre sont unitaires
    for (const auto& [facteur, multiplicite] : p.sansCarre()) {
        std::vector<RacineIsolee> racines = isolerRacines(facteur);
        Polynome reste = facteur;
        for (const RacineIsolee& r : racines) {
            if (!r.rationnelle) continue;
            Polynome q, rr;
            reste.diviser(Polynome({-r.valeurExacte, Nombre(1)}), q, rr);
            reste = q;
            // x - p/q = (q x - p) / q
            const Nombre den = r.valeurExacte.denominateurNombre();
            const Nombre num = r.valeurExacte.numerateurNombre();
            const ExprPtr lineaire = nombre(den) * x - nombre(num);
            produit_.multiplier(ast_pow(lineaire, static_cast<double>(multiplicite)));
            produit_.multiplier(nombre(den.puissanceEntiere(-multiplicite)));
        }
        if (reste.degre() >= 1) {
            Nombre contenu;
            const Polynome primitive = partiePrimitive(reste, contenu);
            produit_.multiplier(ast_pow(primitive.versExpression(x), static_cast<double>(multiplicite)));
            produit_.multiplier(nombre(contenu.puissanceEntiere(multiplicite)));
        }
    }
    return produit_.construire();
}

} // namespace symalgo
