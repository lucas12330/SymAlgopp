/**
 * @file Factorisation.cpp
 * @brief Euclide étendu et factorisation complète des polynômes sur Q.
 *
 * facteursIrreductibles : contenu et partie primitive, décomposition sans carré (Yun), puis,
 * pour chaque facteur sans carré, l'algorithme de Berlekamp-Zassenhaus :
 *   1. factorisation modulo un petit nombre premier p (algorithme de Berlekamp). Plusieurs
 *      premiers sont essayés : le degré d'un facteur sur Z est une somme de degrés de
 *      facteurs modulo chacun d'eux, ce qui suffit souvent à prouver l'irréductibilité ;
 *   2. remontée de Hensel linéaire jusqu'à p^k > 2 |lc| B, B borne de Mignotte des
 *      coefficients d'un facteur ;
 *   3. recombinaison : les produits de facteurs modulaires sont essayés par taille
 *      croissante ; un candidat est un facteur s'il divise exactement le polynôme sur Z.
 *
 * L'arithmétique modulaire tient sur 64 bits (produits sur 128 bits). Au-delà (degré,
 * coefficients ou borne trop grands), seuls les facteurs linéaires rationnels sont extraits
 * et le reste est rendu tel quel, signalé non certifié.
 */

#include "Polynome.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

namespace symalgo {

// ============================================================================
// Euclide étendu
// ============================================================================

Polynome Polynome::pgcdEtendu(const Polynome& a, const Polynome& b, Polynome& u, Polynome& v) {
    Polynome r0 = a, r1 = b;
    Polynome s0({Nombre(1)}), s1, t0, t1({Nombre(1)});
    while (!r1.estNul()) {
        Polynome q, r;
        r0.diviser(r1, q, r);
        Polynome s2 = s0 - q * s1;
        Polynome t2 = t0 - q * t1;
        r0 = std::move(r1);
        r1 = std::move(r);
        s0 = std::move(s1);
        s1 = std::move(s2);
        t0 = std::move(t1);
        t1 = std::move(t2);
    }
    if (r0.estNul()) {
        u = Polynome();
        v = Polynome();
        return r0;
    }
    const Nombre inverse = r0.dominant().inverse();
    u = s0 * inverse;
    v = t0 * inverse;
    return r0 * inverse;
}

namespace {

// ============================================================================
// Polynômes sur Z/mZ
// ============================================================================

using Entier = std::int64_t;
__extension__ typedef __int128 Large; // produits de deux résidus de 62 bits

// Coefficients dans [0, m[, c[i] devant x^i, sans zéro dominant (vide : polynôme nul)
using PolyMod = std::vector<Entier>;

constexpr int DEGRE_MAX = 62;                       // les ensembles de degrés tiennent sur 64 bits
constexpr std::size_t FACTEURS_MODULAIRES_MAX = 18; // au-delà, trop de recombinaisons
constexpr int PREMIERS_ESSAYES = 4;                 // premiers utilisables essayés
constexpr long double BORNE_MAX = 18014398509481984.0L; // 2^54 : p^k < 2^62 pour p < 2^8

const Entier PREMIERS[] = {3,  5,  7,  11, 13, 17, 19, 23, 29,  31,  37,  41,  43,  47,  53,
                           59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127};

inline Entier reduire(Entier a, Entier m) {
    a %= m;
    return a < 0 ? a + m : a;
}

inline Entier produitMod(Entier a, Entier b, Entier m) { return static_cast<Entier>(static_cast<Large>(a) * b % m); }

// Inverse de a modulo m (a et m premiers entre eux)
Entier inverseMod(Entier a, Entier m) {
    Entier r0 = m, r1 = reduire(a, m);
    Large s0 = 0, s1 = 1;
    while (r1 != 0) {
        const Entier q = r0 / r1;
        const Entier r2 = r0 - q * r1;
        const Large s2 = s0 - static_cast<Large>(q) * s1;
        r0 = r1;
        r1 = r2;
        s0 = s1;
        s1 = s2;
    }
    return reduire(static_cast<Entier>(s0 % m), m);
}

void rogner(PolyMod& a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

int degreDe(const PolyMod& a) { return static_cast<int>(a.size()) - 1; }

PolyMod soustraire(const PolyMod& a, const PolyMod& b, Entier m) {
    PolyMod c(std::max(a.size(), b.size()), 0);
    for (std::size_t i = 0; i < c.size(); ++i) {
        const Entier x = i < a.size() ? a[i] : 0, y = i < b.size() ? b[i] : 0;
        c[i] = x >= y ? x - y : x - y + m;
    }
    rogner(c);
    return c;
}

PolyMod multiplier(const PolyMod& a, const PolyMod& b, Entier m) {
    if (a.empty() || b.empty()) return {};
    PolyMod c(a.size() + b.size() - 1, 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;
        for (std::size_t j = 0; j < b.size(); ++j) {
            c[i + j] = static_cast<Entier>((static_cast<Large>(a[i]) * b[j] + c[i + j]) % m);
        }
    }
    rogner(c);
    return c;
}

// a = q b + r, deg r < deg b ; le coefficient dominant de b est inversible modulo m
void diviser(const PolyMod& a, const PolyMod& b, Entier m, PolyMod* q, PolyMod& r) {
    PolyMod reste = a;
    const int db = degreDe(b);
    const Entier inverse = b.back() == 1 ? 1 : inverseMod(b.back(), m);
    if (q) q->assign(static_cast<std::size_t>(std::max(0, degreDe(a) - db + 1)), 0);
    for (int k = degreDe(reste); k >= db; --k) {
        const Entier c = produitMod(reste[static_cast<std::size_t>(k)], inverse, m);
        if (q) (*q)[static_cast<std::size_t>(k - db)] = c;
        if (c == 0) continue;
        for (int j = 0; j <= db; ++j) {
            Entier& x = reste[static_cast<std::size_t>(k - db + j)];
            x = reduire(x - produitMod(c, b[static_cast<std::size_t>(j)], m), m);
        }
    }
    if (degreDe(reste) >= db) reste.resize(static_cast<std::size_t>(db));
    rogner(reste);
    if (q) rogner(*q);
    r = std::move(reste);
}

PolyMod rendreUnitaire(PolyMod a, Entier p) {
    if (a.empty() || a.back() == 1) return a;
    const Entier inverse = inverseMod(a.back(), p);
    for (Entier& c : a) c = produitMod(c, inverse, p);
    return a;
}

PolyMod pgcdMod(PolyMod a, PolyMod b, Entier p) {
    while (!b.empty()) {
        PolyMod r;
        diviser(a, b, p, nullptr, r);
        a = std::move(b);
        b = std::move(r);
    }
    return rendreUnitaire(std::move(a), p);
}

PolyMod deriveeMod(const PolyMod& a, Entier p) {
    PolyMod d;
    for (std::size_t i = 1; i < a.size(); ++i) d.push_back(produitMod(a[i], static_cast<Entier>(i) % p, p));
    rogner(d);
    return d;
}

// s a + t b = 1 modulo p (a et b premiers entre eux modulo p)
void bezoutMod(const PolyMod& a, const PolyMod& b, Entier p, PolyMod& s, PolyMod& t) {
    PolyMod r0 = a, r1 = b, s0 = {1}, s1, t0, t1 = {1};
    while (!r1.empty()) {
        PolyMod q, r;
        diviser(r0, r1, p, &q, r);
        PolyMod s2 = soustraire(s0, multiplier(q, s1, p), p);
        PolyMod t2 = soustraire(t0, multiplier(q, t1, p), p);
        r0 = std::move(r1);
        r1 = std::move(r);
        s0 = std::move(s1);
        s1 = std::move(s2);
        t0 = std::move(t1);
        t1 = std::move(t2);
    }
    const Entier inverse = inverseMod(r0[0], p);
    for (Entier& c : s0) c = produitMod(c, inverse, p);
    for (Entier& c : t0) c = produitMod(c, inverse, p);
    s = std::move(s0);
    t = std::move(t0);
}

// ============================================================================
// Berlekamp : facteurs irréductibles modulo un petit premier
// ============================================================================

/*
 * Nom : berlekamp
 * Description : Facteurs irréductibles unitaires de f (unitaire, sans carré) modulo p. Les
 *               polynômes v tels que v^p = v modulo f forment un espace de dimension r, le
 *               nombre de facteurs : c'est le noyau à gauche de Q - I, où la ligne i de Q
 *               est x^(i p) mod f. Chaque v sépare les facteurs : f = prod_s pgcd(f, v - s).
 */
std::vector<PolyMod> berlekamp(const PolyMod& f, Entier p) {
    const int n = degreDe(f);
    if (n <= 1) return {f};
    auto modF = [&](const PolyMod& a) {
        PolyMod r;
        diviser(a, f, p, nullptr, r);
        return r;
    };
    PolyMod xp = {1}, base = {0, 1}; // x^p mod f
    for (Entier e = p; e > 0; e >>= 1) {
        if (e & 1) xp = modF(multiplier(xp, base, p));
        if (e > 1) base = modF(multiplier(base, base, p));
    }
    // A = (Q - I) transposée : son noyau donne les v
    const std::size_t taille = static_cast<std::size_t>(n);
    std::vector<std::vector<Entier>> A(taille, std::vector<Entier>(taille, 0));
    PolyMod ligne = {1};
    for (std::size_t i = 0; i < taille; ++i) {
        for (std::size_t j = 0; j < ligne.size(); ++j) A[j][i] = ligne[j];
        A[i][i] = reduire(A[i][i] - 1, p);
        ligne = modF(multiplier(ligne, xp, p));
    }
    std::vector<int> lignePivot(taille, -1); // ligne du pivot de chaque colonne, ou -1
    std::size_t rang = 0;
    for (std::size_t col = 0; col < taille && rang < taille; ++col) {
        std::size_t pivot = rang;
        while (pivot < taille && A[pivot][col] == 0) ++pivot;
        if (pivot == taille) continue;
        std::swap(A[rang], A[pivot]);
        const Entier inverse = inverseMod(A[rang][col], p);
        for (Entier& c : A[rang]) c = produitMod(c, inverse, p);
        for (std::size_t l = 0; l < taille; ++l) {
            const Entier k = A[l][col];
            if (l == rang || k == 0) continue;
            for (std::size_t c = col; c < taille; ++c) A[l][c] = reduire(A[l][c] - produitMod(k, A[rang][c], p), p);
        }
        lignePivot[col] = static_cast<int>(rang++);
    }
    std::vector<PolyMod> noyau;
    for (std::size_t libre = 0; libre < taille; ++libre) {
        if (lignePivot[libre] >= 0) continue;
        PolyMod v(taille, 0);
        v[libre] = 1;
        for (std::size_t col = 0; col < taille; ++col) {
            if (lignePivot[col] >= 0) v[col] = reduire(-A[static_cast<std::size_t>(lignePivot[col])][libre], p);
        }
        rogner(v);
        noyau.push_back(std::move(v));
    }

    std::vector<PolyMod> facteurs = {f};
    for (const PolyMod& v : noyau) {
        if (facteurs.size() == noyau.size()) break;
        if (degreDe(v) < 1) continue; // v constant : ne sépare rien
        std::vector<PolyMod> suivants;
        for (const PolyMod& u : facteurs) {
            if (degreDe(u) <= 1) {
                suivants.push_back(u);
                continue;
            }
            PolyMod reste = u;
            for (Entier s = 0; s < p && degreDe(reste) >= 1; ++s) {
                PolyMod w = v;
                w[0] = reduire(w[0] - s, p);
                const PolyMod g = pgcdMod(reste, w, p);
                if (degreDe(g) < 1) continue;
                PolyMod q, r;
                diviser(reste, g, p, &q, r);
                reste = std::move(q);
                suivants.push_back(g);
            }
        }
        facteurs = std::move(suivants);
    }
    return facteurs;
}

// Bit d à 1 : un produit de facteurs de la liste peut être de degré d
std::uint64_t degresPossibles(const std::vector<PolyMod>& facteurs) {
    std::uint64_t bits = 1;
    for (const PolyMod& f : facteurs) bits |= bits << degreDe(f);
    return bits;
}

// ============================================================================
// Remontée de Hensel
// ============================================================================

/*
 * Nom : remonter
 * Description : F unitaire modulo M = p^k, F = prod a_i modulo p (a_i unitaires, premiers
 *               entre eux deux à deux) : renvoie les A_i unitaires, A_i = a_i modulo p et
 *               F = prod A_i modulo M. Remontée linéaire, un facteur après l'autre : pour
 *               F = A B modulo m, l'erreur e = (F - A B)/m se corrige par A += m alpha,
 *               B += m beta avec alpha b + beta a = e modulo p, deg alpha < deg a.
 */
std::vector<PolyMod> remonter(PolyMod F, const std::vector<PolyMod>& facteursP, Entier p, Entier M) {
    std::vector<PolyMod> leves;
    for (std::size_t i = 0; i + 1 < facteursP.size(); ++i) {
        const PolyMod& a = facteursP[i];
        PolyMod b = {1};
        for (std::size_t j = i + 1; j < facteursP.size(); ++j) b = multiplier(b, facteursP[j], p);
        PolyMod s, t;
        bezoutMod(a, b, p, s, t);
        PolyMod A = a, B = b;
        for (Entier m = p; m < M; m *= p) {
            const PolyMod E = soustraire(F, multiplier(A, B, M), M);
            PolyMod e(E.size());
            for (std::size_t k = 0; k < E.size(); ++k) e[k] = (E[k] / m) % p;
            rogner(e);
            if (e.empty()) continue;
            PolyMod alpha, beta, reste;
            diviser(multiplier(t, e, p), a, p, nullptr, alpha);
            diviser(soustraire(e, multiplier(alpha, b, p), p), a, p, &beta, reste);
            if (A.size() < alpha.size()) A.resize(alpha.size(), 0);
            if (B.size() < beta.size()) B.resize(beta.size(), 0);
            for (std::size_t k = 0; k < alpha.size(); ++k) A[k] = static_cast<Entier>((static_cast<Large>(m) * alpha[k] + A[k]) % M);
            for (std::size_t k = 0; k < beta.size(); ++k) B[k] = static_cast<Entier>((static_cast<Large>(m) * beta[k] + B[k]) % M);
        }
        leves.push_back(std::move(A));
        F = std::move(B);
    }
    leves.push_back(std::move(F));
    return leves;
}

// ============================================================================
// Zassenhaus
// ============================================================================

bool coefficientsEntiers(const Polynome& g, std::vector<Entier>& c) {
    c.resize(static_cast<std::size_t>(g.degre() + 1));
    for (int i = 0; i <= g.degre(); ++i) {
        long long v;
        if (!g.coefficient(i).versEntier(v)) return false;
        c[static_cast<std::size_t>(i)] = v;
    }
    return true;
}

PolyMod modulo(const std::vector<Entier>& c, Entier m) {
    PolyMod r(c.size());
    for (std::size_t i = 0; i < c.size(); ++i) r[i] = reduire(c[i], m);
    rogner(r);
    return r;
}

// Sous-ensemble suivant de même cardinal (ordre croissant des masques)
std::uint32_t sousEnsembleSuivant(std::uint32_t v) {
    const std::uint32_t c = v & (0u - v), r = v + c;
    return (((r ^ v) >> 2) / c) | r;
}

/*
 * Nom : zassenhaus
 * Description : Facteurs irréductibles de g (entier, primitif, sans carré, coefficient
 *               dominant positif, degré >= 3), ajoutés à facteurs. Renvoie faux sans rien
 *               ajouter si la factorisation sort du domaine traité (voir l'en-tête).
 */
bool zassenhaus(const Polynome& g, std::vector<Polynome>& facteurs) {
    const int n = g.degre();
    std::vector<Entier> c;
    if (n > DEGRE_MAX || !coefficientsEntiers(g, c)) return false;

    // Borne de Mignotte : un facteur h de g vérifie |h_i| <= 2^n ||g||_2, et l'on reconstruit lc * h
    long double norme = 0.0L;
    for (Entier ci : c) norme += static_cast<long double>(ci) * static_cast<long double>(ci);
    const long double borne = 2.0L * std::fabs(static_cast<long double>(c.back())) * std::pow(2.0L, n) * std::sqrt(norme);
    if (!(borne < BORNE_MAX)) return false;

    // 1. Modulo quelques petits premiers : degrés possibles et meilleure factorisation
    const std::uint64_t interieur = ((std::uint64_t(1) << n) - 1) & ~std::uint64_t(1); // degrés 1 .. n-1
    std::uint64_t degres = ~std::uint64_t(0);
    Entier p = 0;
    std::vector<PolyMod> modulaires;
    int essais = 0;
    for (Entier premier : PREMIERS) {
        if (c.back() % premier == 0) continue;
        const PolyMod gp = modulo(c, premier);
        if (degreDe(pgcdMod(gp, deriveeMod(gp, premier), premier)) > 0) continue; // carré modulo ce premier
        std::vector<PolyMod> fp = berlekamp(rendreUnitaire(gp, premier), premier);
        degres &= degresPossibles(fp);
        if ((degres & interieur) == 0) { // aucun degré de facteur propre possible
            facteurs.push_back(g);
            return true;
        }
        if (p == 0 || fp.size() < modulaires.size()) {
            p = premier;
            modulaires = std::move(fp);
        }
        if (++essais == PREMIERS_ESSAYES) break;
    }
    if (p == 0 || modulaires.size() > FACTEURS_MODULAIRES_MAX) return false;

    // 2. Remontée jusqu'à M = p^k > borne (M < 2^62)
    Entier M = p;
    while (static_cast<long double>(M) <= borne) M *= p;
    PolyMod F = modulo(c, M);
    const Entier inverseDominant = inverseMod(F.back(), M);
    for (Entier& ci : F) ci = produitMod(ci, inverseDominant, M);
    std::vector<PolyMod> pool = remonter(std::move(F), modulaires, p, M);

    // 3. Recombinaison par taille croissante
    Polynome reste = g;
    std::size_t taille = 1;
    while (2 * taille <= pool.size()) {
        long long dominant = 0, constant = 0;
        reste.dominant().versEntier(dominant);
        reste.coefficient(0).versEntier(constant);
        const Entier dominantMod = reduire(dominant, M);
        const std::uint32_t fin = std::uint32_t(1) << pool.size();
        bool trouve = false;
        for (std::uint32_t masque = (std::uint32_t(1) << taille) - 1; masque < fin; masque = sousEnsembleSuivant(masque)) {
            int degre = 0;
            Entier terme = dominantMod; // terme constant du candidat lc * prod G_i
            for (std::size_t i = 0; i < pool.size(); ++i) {
                if (!((masque >> i) & 1u)) continue;
                degre += degreDe(pool[i]);
                terme = produitMod(terme, pool[i][0], M);
            }
            if (!((degres >> degre) & 1u)) continue;
            // Test rapide : le terme constant du candidat divise lc * g(0)
            const Entier symetrique = terme > M / 2 ? terme - M : terme;
            const Large cible = static_cast<Large>(dominant) * constant;
            if (symetrique == 0 ? cible != 0 : cible % symetrique != 0) continue;

            PolyMod h = {dominantMod};
            for (std::size_t i = 0; i < pool.size(); ++i) {
                if ((masque >> i) & 1u) h = multiplier(h, pool[i], M);
            }
            std::vector<Nombre> coefficients;
            coefficients.reserve(h.size());
            for (Entier hi : h) coefficients.push_back(Nombre(hi > M / 2 ? hi - M : hi));
            Nombre contenu;
            const Polynome candidat = Polynome(std::move(coefficients)).partiePrimitive(contenu);
            Polynome q, r;
            reste.diviser(candidat, q, r);
            if (!r.estNul()) continue;
            facteurs.push_back(candidat);
            reste = std::move(q);
            for (std::size_t i = pool.size(); i-- > 0;) {
                if ((masque >> i) & 1u) pool.erase(pool.begin() + static_cast<std::ptrdiff_t>(i));
            }
            trouve = true;
            break;
        }
        if (!trouve) ++taille;
    }
    if (reste.degre() > 0) facteurs.push_back(std::move(reste));
    return true;
}

/*
 * Nom : factoriserSansCarre
 * Description : Facteurs irréductibles de g (entier, primitif, sans carré, coefficient
 *               dominant positif). Renvoie faux si le dernier facteur ajouté n'est pas
 *               certifié irréductible.
 */
bool factoriserSansCarre(const Polynome& g, std::vector<Polynome>& facteurs) {
    if (g.degre() < 1) return true;
    if (g.degre() == 1) {
        facteurs.push_back(g);
        return true;
    }
    if (g.degre() == 2) {
        // Réductible si et seulement si le discriminant est un carré : (2a x + b - r)(2a x + b + r) / (4a)
        const Nombre& a = g.coefficient(2);
        const Nombre& b = g.coefficient(1);
        const Nombre delta = b * b - Nombre(4) * a * g.coefficient(0);
        Nombre racine, contenu;
        if (delta.signe() > 0 && delta.racineExacte(2, racine)) {
            facteurs.push_back(Polynome({b - racine, Nombre(2) * a}).partiePrimitive(contenu));
            facteurs.push_back(Polynome({b + racine, Nombre(2) * a}).partiePrimitive(contenu));
        } else {
            facteurs.push_back(g);
        }
        return true;
    }
    if (zassenhaus(g, facteurs)) return true;
    // Hors du domaine de Zassenhaus : facteurs linéaires des racines rationnelles, reste non certifié
    Polynome reste = g;
    for (const RacineReelle& racine : g.racinesReelles()) {
        const Constante* valeur = racine.exacte ? comme<Constante>(racine.exacte) : nullptr;
        if (!valeur) continue;
        const Nombre& v = valeur->getNombre();
        const Polynome lineaire({-v.numerateurNombre(), v.denominateurNombre()});
        Polynome q, r;
        reste.diviser(lineaire, q, r);
        if (!r.estNul()) continue;
        Nombre contenu;
        reste = q.partiePrimitive(contenu);
        facteurs.push_back(lineaire);
    }
    if (reste.degre() < 1) return true;
    const bool certifie = reste.degre() == 1;
    facteurs.push_back(std::move(reste));
    return certifie;
}

// Ordre des facteurs : degré croissant, puis coefficients à partir du terme dominant
bool avant(const Polynome& a, const Polynome& b) {
    if (a.degre() != b.degre()) return a.degre() < b.degre();
    for (int i = a.degre(); i >= 0; --i) {
        const int c = a.coefficient(i).comparer(b.coefficient(i));
        if (c != 0) return c < 0;
    }
    return false;
}

} // namespace

std::vector<std::pair<Polynome, int>> Polynome::facteursIrreductibles(Nombre& contenu, bool* complet) const {
    std::vector<std::pair<Polynome, int>> resultat;
    if (complet) *complet = true;
    if (estNul()) {
        contenu = Nombre(0);
        return resultat;
    }
    const Polynome primitif = partiePrimitive(contenu);
    if (!m_exact) { // coefficients approchés : aucune factorisation exacte n'a de sens
        if (complet) *complet = false;
        if (primitif.degre() > 0) resultat.emplace_back(primitif, 1);
        return resultat;
    }
    for (const auto& [facteur, multiplicite] : primitif.sansCarre()) {
        Nombre c;
        std::vector<Polynome> irreductibles;
        if (!factoriserSansCarre(facteur.partiePrimitive(c), irreductibles) && complet) *complet = false;
        for (Polynome& f : irreductibles) resultat.emplace_back(std::move(f), multiplicite);
    }
    std::sort(resultat.begin(), resultat.end(), [](const auto& a, const auto& b) { return avant(a.first, b.first); });
    return resultat;
}

// ============================================================================
// Factorisation d'une expression
// ============================================================================

ExprPtr factoriser(const ExprPtr& e) {
    const ExprPtr d = developper(e);
    Polynome p;
    ExprPtr x;
    if (!Polynome::depuisExpression(d, p, &x) || !x || p.degre() < 2 || !p.estExact()) return e;
    Nombre contenu;
    AccumulateurProduit produit_;
    for (const auto& [facteur, multiplicite] : p.facteursIrreductibles(contenu)) {
        produit_.multiplier(ast_pow(facteur.versExpression(x), static_cast<double>(multiplicite)));
    }
    produit_.multiplier(nombre(contenu));
    return produit_.construire();
}

} // namespace symalgo
