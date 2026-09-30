/**
 * @file Nombre.hpp
 * @brief Nombre exact (rationnel de taille arbitraire) ou réel (double).
 *
 * Un Nombre exact est un rationnel irréductible p/q (q > 0). Tant qu'il tient sur
 * 64 bits, il est stocké directement (aucune allocation) ; au-delà, il bascule sur un
 * rationnel GMP (mpq_class), et revient au stockage direct dès que le résultat tient
 * de nouveau. Un Nombre réel est un double.
 *
 * Arithmétique : exact op exact reste exact ; dès qu'un opérande est réel, le résultat
 * est réel (comme dans GiNaC ou SymPy).
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>

// La dépendance à GMP reste confinée à Nombre.cpp (type opaque Nombre::Grand)

namespace symalgo {

class Nombre {
public:
    /*
     * Nom : Nombre
     * Description : Entier exact (0 par défaut).
     * Utilisation : Nombre n(42);
     */
    Nombre(long long entier = 0);

    /*
     * Nom : rationnel
     * Description : Rationnel exact num/den, réduit (lève std::invalid_argument si den = 0).
     * Utilisation : Nombre tiers = Nombre::rationnel(1, 3);
     */
    static Nombre rationnel(long long num, long long den);

    /*
     * Nom : reel
     * Description : Nombre réel (inexact) de valeur v.
     */
    static Nombre reel(double v);

    /*
     * Nom : depuisDouble
     * Description : Conversion d'un double saisi par l'utilisateur : une valeur entière
     *               (|v| <= 2^53) devient un entier exact, toute autre valeur un réel.
     * Utilisation : Nombre::depuisDouble(2.0) est l'entier 2, depuisDouble(0.5) le réel 0.5.
     */
    static Nombre depuisDouble(double v);

    /*
     * Nom : depuisTexte
     * Description : Rationnel exact de taille arbitraire, écrit « p » ou « p/q » en base 10.
     * Utilisation : Nombre n = Nombre::depuisTexte("123456789012345678901234567890/7");
     */
    static Nombre depuisTexte(const std::string& texte);

    Nombre(const Nombre& autre);
    Nombre(Nombre&& autre) noexcept;
    Nombre& operator=(const Nombre& autre);
    Nombre& operator=(Nombre&& autre) noexcept;
    ~Nombre();

    bool estExact() const { return m_forme != Forme::Reel; }
    bool estEntier() const;          // exact et de dénominateur 1
    bool estZero() const;
    bool estUn() const;              // vaut 1 (exact ou réel)
    bool estMoinsUn() const;
    int signe() const;               // -1, 0 ou 1 (0 aussi pour NaN)
    bool estFini() const;            // faux pour un réel infini ou NaN

    /*
     * Nom : versDouble
     * Description : Valeur approchée en double.
     */
    double versDouble() const;

    /*
     * Nom : versEntier
     * Description : Écrit la valeur si c'est un entier exact tenant sur 64 bits.
     */
    bool versEntier(long long& n) const;

    /*
     * Nom : numerateur / denominateur
     * Description : Composantes d'un nombre exact, en base 10 (taille arbitraire).
     */
    std::string numerateur() const;
    std::string denominateur() const;

    Nombre operator+(const Nombre& b) const;
    Nombre operator-(const Nombre& b) const;
    Nombre operator*(const Nombre& b) const;
    Nombre operator/(const Nombre& b) const; // lève std::domain_error pour une division exacte par 0
    Nombre operator-() const;
    Nombre inverse() const;

    /*
     * Nom : puissanceEntiere
     * Description : Puissance d'exposant entier (exacte pour un nombre exact).
     *               Lève std::domain_error pour 0 exact à une puissance négative.
     */
    Nombre puissanceEntiere(long long e) const;

    /*
     * Nom : racineExacte
     * Description : Si ce nombre exact positif est la puissance k-ième d'un rationnel,
     *               écrit cette racine et renvoie vrai (ex. 4/9 et k = 2 donnent 2/3).
     */
    bool racineExacte(long long k, Nombre& racine) const;

    /*
     * Nom : operator== / comparer
     * Description : Égalité stricte (un exact et un réel ne sont jamais égaux) et ordre
     *               total utilisé pour trier les expressions (par valeur, puis exact avant réel).
     */
    bool operator==(const Nombre& b) const;
    bool operator!=(const Nombre& b) const { return !(*this == b); }
    int comparer(const Nombre& b) const;

    std::size_t hash() const;

    /*
     * Nom : afficher
     * Description : Écrit « 3 », « -1/2 » (exact) ou la valeur du double (réel).
     */
    void afficher(std::ostream& os) const;
    std::string texte() const;

private:
    enum class Forme : std::uint8_t { Petit, Grand, Reel };

    struct Grand; // rationnel GMP (défini dans Nombre.cpp)

    Forme m_forme = Forme::Petit;
    union {
        struct {
            std::int64_t num;
            std::int64_t den;
        } m_petit;
        double m_reel;
        Grand* m_grand;
    };

    static Nombre depuisGrand(Grand&& g); // redescend en Petit si possible
    Grand versGrand() const;
    void liberer() noexcept;
};

std::ostream& operator<<(std::ostream& os, const Nombre& n);

} // namespace symalgo
