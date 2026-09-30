/**
 * @file Nombre.cpp
 * @brief Implémentation des nombres exacts (int64 avec repli sur GMP) et réels.
 */

#include "Nombre.hpp"

#include <gmpxx.h>

#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace symalgo {

struct Nombre::Grand {
    mpq_class q;
};

namespace {

constexpr std::int64_t INT64_MIN_ = std::numeric_limits<std::int64_t>::min();
constexpr double DEUX_PUISSANCE_53 = 9007199254740992.0;

mpz_class versMpz(std::int64_t v) {
    if constexpr (sizeof(long) >= sizeof(std::int64_t)) {
        return mpz_class(static_cast<long>(v));
    } else {
        return mpz_class(std::to_string(v));
    }
}

// Écrit z dans v s'il tient sur un int64 différent de INT64_MIN
bool depuisMpz(const mpz_class& z, std::int64_t& v) {
    if (mpz_sizeinbase(z.get_mpz_t(), 2) > 63) return false;
    if constexpr (sizeof(long) >= sizeof(std::int64_t)) {
        v = static_cast<std::int64_t>(z.get_si());
    } else {
        v = std::stoll(z.get_str());
    }
    return true;
}

// Réduit num/den (den != 0) avec den > 0 ; faux si une valeur vaut INT64_MIN
bool reduire(std::int64_t& num, std::int64_t& den) {
    if (num == INT64_MIN_ || den == INT64_MIN_) return false;
    if (den < 0) {
        num = -num;
        den = -den;
    }
    if (num == 0) {
        den = 1;
        return true;
    }
    const std::int64_t g = std::gcd(num, den);
    num /= g;
    den /= g;
    return true;
}

bool multiplier(std::int64_t a, std::int64_t b, std::int64_t& r) { return !__builtin_mul_overflow(a, b, &r); }
bool additionner(std::int64_t a, std::int64_t b, std::int64_t& r) { return !__builtin_add_overflow(a, b, &r); }

void combiner(std::size_t& graine, std::size_t valeur) {
    graine ^= valeur + 0x9e3779b97f4a7c15ULL + (graine << 6) + (graine >> 2);
}

} // namespace

// ============================================================================
// Construction, copie, destruction
// ============================================================================

Nombre::Nombre(long long entier) {
    if (entier == INT64_MIN_) {
        m_forme = Forme::Grand;
        m_grand = new Grand{mpq_class(versMpz(entier))};
    } else {
        m_petit = {static_cast<std::int64_t>(entier), 1};
    }
}

Nombre Nombre::rationnel(long long num, long long den) {
    if (den == 0) throw std::invalid_argument("Nombre : denominateur nul");
    std::int64_t n = num, d = den;
    if (reduire(n, d)) {
        Nombre r;
        r.m_petit = {n, d};
        return r;
    }
    Grand g{mpq_class(versMpz(num), versMpz(den))};
    g.q.canonicalize();
    return depuisGrand(std::move(g));
}

Nombre Nombre::reel(double v) {
    Nombre r;
    r.m_forme = Forme::Reel;
    r.m_reel = v;
    return r;
}

Nombre Nombre::depuisDouble(double v) {
    if (std::isfinite(v) && std::floor(v) == v && std::abs(v) <= DEUX_PUISSANCE_53) {
        return Nombre(static_cast<long long>(v));
    }
    return reel(v);
}

Nombre Nombre::depuisTexte(const std::string& texte) {
    Grand g;
    if (g.q.set_str(texte, 10) != 0) throw std::invalid_argument("Nombre : texte invalide « " + texte + " »");
    if (g.q.get_den() == 0) throw std::invalid_argument("Nombre : denominateur nul");
    g.q.canonicalize();
    return depuisGrand(std::move(g));
}

Nombre::Nombre(const Nombre& autre) : m_forme(autre.m_forme) {
    switch (m_forme) {
        case Forme::Petit: m_petit = autre.m_petit; break;
        case Forme::Reel: m_reel = autre.m_reel; break;
        case Forme::Grand: m_grand = new Grand(*autre.m_grand); break;
    }
}

Nombre::Nombre(Nombre&& autre) noexcept : m_forme(autre.m_forme) {
    switch (m_forme) {
        case Forme::Petit: m_petit = autre.m_petit; break;
        case Forme::Reel: m_reel = autre.m_reel; break;
        case Forme::Grand:
            m_grand = autre.m_grand;
            autre.m_forme = Forme::Petit;
            autre.m_petit = {0, 1};
            break;
    }
}

Nombre& Nombre::operator=(const Nombre& autre) {
    if (this != &autre) {
        Nombre copie(autre);
        *this = std::move(copie);
    }
    return *this;
}

Nombre& Nombre::operator=(Nombre&& autre) noexcept {
    if (this != &autre) {
        liberer();
        m_forme = autre.m_forme;
        switch (m_forme) {
            case Forme::Petit: m_petit = autre.m_petit; break;
            case Forme::Reel: m_reel = autre.m_reel; break;
            case Forme::Grand:
                m_grand = autre.m_grand;
                autre.m_forme = Forme::Petit;
                autre.m_petit = {0, 1};
                break;
        }
    }
    return *this;
}

Nombre::~Nombre() { liberer(); }

void Nombre::liberer() noexcept {
    if (m_forme == Forme::Grand) delete m_grand;
    m_forme = Forme::Petit;
    m_petit = {0, 1};
}

Nombre Nombre::depuisGrand(Grand&& g) {
    std::int64_t n, d;
    if (depuisMpz(g.q.get_num(), n) && depuisMpz(g.q.get_den(), d)) {
        Nombre r;
        r.m_petit = {n, d};
        return r;
    }
    Nombre r;
    r.m_forme = Forme::Grand;
    r.m_grand = new Grand(std::move(g));
    return r;
}

Nombre::Grand Nombre::versGrand() const {
    if (m_forme == Forme::Grand) return *m_grand;
    Grand g{mpq_class(versMpz(m_petit.num), versMpz(m_petit.den))};
    return g; // déjà réduit
}

// ============================================================================
// Propriétés
// ============================================================================

bool Nombre::estEntier() const {
    switch (m_forme) {
        case Forme::Petit: return m_petit.den == 1;
        case Forme::Grand: return m_grand->q.get_den() == 1;
        default: return false;
    }
}

bool Nombre::estZero() const {
    return (m_forme == Forme::Petit && m_petit.num == 0) || (m_forme == Forme::Reel && m_reel == 0.0);
}

bool Nombre::estUn() const {
    return (m_forme == Forme::Petit && m_petit.num == 1 && m_petit.den == 1) ||
           (m_forme == Forme::Reel && m_reel == 1.0);
}

bool Nombre::estMoinsUn() const {
    return (m_forme == Forme::Petit && m_petit.num == -1 && m_petit.den == 1) ||
           (m_forme == Forme::Reel && m_reel == -1.0);
}

int Nombre::signe() const {
    switch (m_forme) {
        case Forme::Petit: return (m_petit.num > 0) - (m_petit.num < 0);
        case Forme::Grand: return sgn(m_grand->q);
        default: return (m_reel > 0.0) - (m_reel < 0.0);
    }
}

bool Nombre::estFini() const { return m_forme != Forme::Reel || std::isfinite(m_reel); }

double Nombre::versDouble() const {
    switch (m_forme) {
        case Forme::Petit:
            return m_petit.den == 1 ? static_cast<double>(m_petit.num)
                                    : static_cast<double>(m_petit.num) / static_cast<double>(m_petit.den);
        case Forme::Grand: return m_grand->q.get_d();
        default: return m_reel;
    }
}

bool Nombre::versEntier(long long& n) const {
    if (m_forme != Forme::Petit || m_petit.den != 1) return false;
    n = m_petit.num;
    return true;
}

std::string Nombre::numerateur() const {
    if (m_forme == Forme::Reel) throw std::logic_error("Nombre::numerateur : nombre reel");
    return versGrand().q.get_num().get_str();
}

std::string Nombre::denominateur() const {
    if (m_forme == Forme::Reel) throw std::logic_error("Nombre::denominateur : nombre reel");
    return versGrand().q.get_den().get_str();
}

// ============================================================================
// Arithmétique
// ============================================================================

Nombre Nombre::operator+(const Nombre& b) const {
    if (m_forme == Forme::Reel || b.m_forme == Forme::Reel) return reel(versDouble() + b.versDouble());
    if (m_forme == Forme::Petit && b.m_forme == Forme::Petit) {
        // a/p + c/q = (a (q/g) + c (p/g)) / (p (q/g)) avec g = pgcd(p, q)
        const std::int64_t g = std::gcd(m_petit.den, b.m_petit.den);
        std::int64_t t1, t2, n, d;
        if (multiplier(m_petit.num, b.m_petit.den / g, t1) && multiplier(b.m_petit.num, m_petit.den / g, t2) &&
            additionner(t1, t2, n) && multiplier(m_petit.den, b.m_petit.den / g, d) && reduire(n, d)) {
            Nombre r;
            r.m_petit = {n, d};
            return r;
        }
    }
    Grand r{versGrand().q + b.versGrand().q};
    return depuisGrand(std::move(r));
}

Nombre Nombre::operator-() const {
    switch (m_forme) {
        case Forme::Petit: {
            Nombre r;
            r.m_petit = {-m_petit.num, m_petit.den}; // num != INT64_MIN par construction
            return r;
        }
        case Forme::Grand: return depuisGrand(Grand{-m_grand->q});
        default: return reel(-m_reel);
    }
}

Nombre Nombre::operator-(const Nombre& b) const { return *this + (-b); }

Nombre Nombre::operator*(const Nombre& b) const {
    if (m_forme == Forme::Reel || b.m_forme == Forme::Reel) return reel(versDouble() * b.versDouble());
    if (m_forme == Forme::Petit && b.m_forme == Forme::Petit) {
        // Réduction croisée : le résultat est directement irréductible
        const std::int64_t g1 = std::gcd(m_petit.num, b.m_petit.den);
        const std::int64_t g2 = std::gcd(b.m_petit.num, m_petit.den);
        const std::int64_t a = g1 ? m_petit.num / g1 : 0, d = g1 ? b.m_petit.den / g1 : b.m_petit.den;
        const std::int64_t c = g2 ? b.m_petit.num / g2 : 0, p = g2 ? m_petit.den / g2 : m_petit.den;
        std::int64_t n, q;
        if (multiplier(a, c, n) && multiplier(p, d, q) && reduire(n, q)) {
            Nombre r;
            r.m_petit = {n, q};
            return r;
        }
    }
    Grand r{versGrand().q * b.versGrand().q};
    return depuisGrand(std::move(r));
}

Nombre Nombre::inverse() const {
    switch (m_forme) {
        case Forme::Reel: return reel(1.0 / m_reel);
        case Forme::Petit: {
            if (m_petit.num == 0) throw std::domain_error("Nombre : division exacte par zero");
            std::int64_t n = m_petit.den, d = m_petit.num;
            reduire(n, d); // remet le signe au numérateur
            Nombre r;
            r.m_petit = {n, d};
            return r;
        }
        default: return depuisGrand(Grand{1 / m_grand->q});
    }
}

Nombre Nombre::operator/(const Nombre& b) const {
    if (m_forme == Forme::Reel || b.m_forme == Forme::Reel) return reel(versDouble() / b.versDouble());
    return *this * b.inverse();
}

Nombre Nombre::puissanceEntiere(long long e) const {
    if (m_forme == Forme::Reel) return reel(std::pow(m_reel, static_cast<double>(e)));
    if (e == std::numeric_limits<long long>::min()) throw std::domain_error("Nombre : exposant hors limites");
    if (e < 0) {
        if (estZero()) throw std::domain_error("Nombre : zero a une puissance negative");
        return inverse().puissanceEntiere(-e);
    }
    Nombre resultat(1), base = *this;
    while (e > 0) {
        if (e & 1) resultat = resultat * base;
        e >>= 1;
        if (e > 0) base = base * base;
    }
    return resultat;
}

bool Nombre::racineExacte(long long k, Nombre& racine) const {
    if (m_forme == Forme::Reel || k < 1 || signe() < 0) return false;
    if (k == 1) {
        racine = *this;
        return true;
    }
    const Grand g = versGrand();
    mpz_class rn, rd;
    const int exactN = mpz_root(rn.get_mpz_t(), g.q.get_num().get_mpz_t(), static_cast<unsigned long>(k));
    const int exactD = mpz_root(rd.get_mpz_t(), g.q.get_den().get_mpz_t(), static_cast<unsigned long>(k));
    if (!exactN || !exactD) return false;
    racine = depuisGrand(Grand{mpq_class(rn, rd)});
    return true;
}

// ============================================================================
// Comparaison, empreinte, affichage
// ============================================================================

bool Nombre::operator==(const Nombre& b) const {
    if (m_forme != b.m_forme) return false; // les petites valeurs sont toujours en forme Petit
    switch (m_forme) {
        case Forme::Petit: return m_petit.num == b.m_petit.num && m_petit.den == b.m_petit.den;
        case Forme::Grand: return m_grand->q == b.m_grand->q;
        default: return m_reel == b.m_reel;
    }
}

int Nombre::comparer(const Nombre& b) const {
    int c;
    if (m_forme == Forme::Reel || b.m_forme == Forme::Reel) {
        const double x = versDouble(), y = b.versDouble();
        const bool nx = std::isnan(x), ny = std::isnan(y);
        if (nx || ny) {
            c = nx && ny ? 0 : (nx ? 1 : -1); // NaN en dernier
        } else {
            c = (x > y) - (x < y);
        }
    } else if (m_forme == Forme::Petit && b.m_forme == Forme::Petit) {
        std::int64_t g, d;
        if (multiplier(m_petit.num, b.m_petit.den, g) && multiplier(b.m_petit.num, m_petit.den, d)) {
            c = (g > d) - (g < d);
        } else {
            c = cmp(versGrand().q, b.versGrand().q);
            c = (c > 0) - (c < 0);
        }
    } else {
        c = cmp(versGrand().q, b.versGrand().q);
        c = (c > 0) - (c < 0);
    }
    if (c != 0) return c;
    // À valeur égale : exact avant réel
    return static_cast<int>(m_forme == Forme::Reel) - static_cast<int>(b.m_forme == Forme::Reel);
}

std::size_t Nombre::hash() const {
    std::size_t h = static_cast<std::size_t>(m_forme);
    switch (m_forme) {
        case Forme::Petit:
            combiner(h, std::hash<std::int64_t>()(m_petit.num));
            combiner(h, std::hash<std::int64_t>()(m_petit.den));
            break;
        case Forme::Grand: combiner(h, std::hash<std::string>()(m_grand->q.get_str())); break;
        default: combiner(h, std::hash<double>()(m_reel == 0.0 ? 0.0 : m_reel)); break;
    }
    return h;
}

void Nombre::afficher(std::ostream& os) const {
    switch (m_forme) {
        case Forme::Petit:
            os << m_petit.num;
            if (m_petit.den != 1) os << '/' << m_petit.den;
            break;
        case Forme::Grand: os << m_grand->q.get_str(); break;
        default: os << m_reel; break;
    }
}

std::string Nombre::texte() const {
    std::ostringstream os;
    afficher(os);
    return os.str();
}

std::ostream& operator<<(std::ostream& os, const Nombre& n) {
    n.afficher(os);
    return os;
}

} // namespace symalgo
