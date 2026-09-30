/**
 * @file test_nombre.cpp
 * @brief Tests unitaires de Nombre (rationnels exacts avec repli GMP, réels).
 */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "Nombre.hpp"
#include "test_framework.hpp"

using namespace symalgo;

TEST_CASE(construction_et_reduction) {
    CHECK_EQ(Nombre(42).texte(), std::string("42"));
    CHECK_EQ(Nombre::rationnel(2, 6).texte(), std::string("1/3"));
    CHECK_EQ(Nombre::rationnel(1, -2).texte(), std::string("-1/2"));
    CHECK_EQ(Nombre::rationnel(0, -5).texte(), std::string("0"));
    CHECK_EQ(Nombre::rationnel(-4, -2).texte(), std::string("2"));
    CHECK_THROWS(Nombre::rationnel(1, 0), std::invalid_argument);
    CHECK(Nombre::rationnel(6, 3).estEntier());
    CHECK(!Nombre::rationnel(1, 3).estEntier());
}

TEST_CASE(depuis_double) {
    CHECK(Nombre::depuisDouble(2.0).estExact());
    CHECK(Nombre::depuisDouble(-7.0).estEntier());
    CHECK(!Nombre::depuisDouble(0.5).estExact());
    CHECK(!Nombre::depuisDouble(1e300).estExact()); // hors de la plage exacte des doubles
    CHECK(!Nombre::depuisDouble(std::numeric_limits<double>::infinity()).estExact());
}

TEST_CASE(arithmetique_exacte) {
    const Nombre tiers = Nombre::rationnel(1, 3), sixieme = Nombre::rationnel(1, 6);
    CHECK_EQ((tiers + sixieme).texte(), std::string("1/2"));
    CHECK_EQ((tiers - sixieme).texte(), std::string("1/6"));
    CHECK_EQ((tiers * Nombre(3)).texte(), std::string("1"));
    CHECK_EQ((tiers / sixieme).texte(), std::string("2"));
    CHECK_EQ((-tiers).texte(), std::string("-1/3"));
    CHECK_EQ(Nombre::rationnel(-2, 3).inverse().texte(), std::string("-3/2"));
    CHECK_THROWS(Nombre(1) / Nombre(0), std::domain_error);
    CHECK_THROWS(Nombre(0).inverse(), std::domain_error);
}

TEST_CASE(debordement_vers_gmp_et_retour) {
    const Nombre grand(4000000000000000000LL); // 4e18, proche de la limite int64
    const Nombre somme = grand + grand;          // 8e18 : tient encore
    CHECK_EQ(somme.texte(), std::string("8000000000000000000"));
    const Nombre carre = grand * grand;          // 1,6e37 : bascule sur GMP
    CHECK_EQ(carre.texte(), std::string("16000000000000000000000000000000000000"));
    // Retour à la forme compacte dès que le résultat tient sur 64 bits
    const Nombre retour = carre / grand;
    CHECK(retour == grand);
    CHECK_EQ(retour.hash(), grand.hash());
    // 2^100 exact, puis 2^100 / 2^99 = 2
    const Nombre p100 = Nombre(2).puissanceEntiere(100);
    CHECK_EQ(p100.texte(), std::string("1267650600228229401496703205376"));
    CHECK(p100 / Nombre(2).puissanceEntiere(99) == Nombre(2));
    // Fractions avec grand dénominateur
    const Nombre f = Nombre(1) / Nombre(3).puissanceEntiere(50);
    CHECK_EQ(f.denominateur(), std::string("717897987691852588770249"));
    CHECK((f * Nombre(3).puissanceEntiere(50)).estUn());
}

TEST_CASE(texte_de_taille_arbitraire) {
    const Nombre n = Nombre::depuisTexte("123456789012345678901234567890/10");
    CHECK_EQ(n.texte(), std::string("12345678901234567890123456789"));
    CHECK_THROWS(Nombre::depuisTexte("abc"), std::invalid_argument);
}

TEST_CASE(melange_exact_reel) {
    const Nombre r = Nombre::rationnel(1, 2) + Nombre::reel(0.25);
    CHECK(!r.estExact());
    CHECK_NEAR(r.versDouble(), 0.75, 1e-15);
    CHECK(!(Nombre(2) == Nombre::reel(2.0))); // exact et réel ne sont jamais égaux
    CHECK(Nombre::reel(1.0).estUn());
}

TEST_CASE(puissances_et_racines) {
    CHECK_EQ(Nombre::rationnel(2, 3).puissanceEntiere(3).texte(), std::string("8/27"));
    CHECK_EQ(Nombre::rationnel(2, 3).puissanceEntiere(-2).texte(), std::string("9/4"));
    CHECK_EQ(Nombre(5).puissanceEntiere(0).texte(), std::string("1"));
    CHECK_THROWS(Nombre(0).puissanceEntiere(-1), std::domain_error);
    Nombre racine;
    CHECK(Nombre::rationnel(4, 9).racineExacte(2, racine));
    CHECK_EQ(racine.texte(), std::string("2/3"));
    CHECK(Nombre(27).racineExacte(3, racine));
    CHECK_EQ(racine.texte(), std::string("3"));
    CHECK(!Nombre(2).racineExacte(2, racine));   // sqrt(2) n'est pas rationnel
    CHECK(!Nombre(-4).racineExacte(2, racine));  // pas de racine réelle
}

TEST_CASE(ordre_total) {
    CHECK(Nombre::rationnel(1, 3).comparer(Nombre::rationnel(1, 2)) < 0);
    CHECK(Nombre(-1).comparer(Nombre(0)) < 0);
    CHECK(Nombre(2).comparer(Nombre::reel(2.0)) < 0); // exact avant réel à valeur égale
    CHECK(Nombre::reel(1.5).comparer(Nombre::rationnel(3, 2)) > 0);
    CHECK(Nombre(2).puissanceEntiere(100).comparer(Nombre(2).puissanceEntiere(99)) > 0);
    CHECK(Nombre::reel(std::nan("")).comparer(Nombre(1000)) > 0);
    CHECK_EQ(Nombre(3).comparer(Nombre(3)), 0);
}

TEST_CASE(proprietes) {
    CHECK(Nombre(0).estZero());
    CHECK(Nombre(1).estUn());
    CHECK(Nombre(-1).estMoinsUn());
    CHECK_EQ(Nombre::rationnel(-3, 4).signe(), -1);
    CHECK_EQ(Nombre(0).signe(), 0);
    long long n = 0;
    CHECK(Nombre(-12).versEntier(n));
    CHECK_EQ(n, -12LL);
    CHECK(!Nombre::rationnel(1, 2).versEntier(n));
    CHECK(!Nombre::reel(std::numeric_limits<double>::infinity()).estFini());
}

int main() { return test::executerTous(); }
