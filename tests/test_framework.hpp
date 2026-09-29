/**
 * @file test_framework.hpp
 * @brief Mini-framework de tests unitaires sans dépendance pour SymAlgo++.
 *
 * Utilisation :
 *   #include "test_framework.hpp"
 *   TEST_CASE(mon_test) {
 *       CHECK(1 + 1 == 2);
 *       CHECK_NEAR(0.1 + 0.2, 0.3, 1e-12);
 *   }
 *   int main() { return test::executerTous(); }
 *
 * Chaque échec affiche le fichier, la ligne et l'expression fautive ; le
 * programme renvoie un code non nul si au moins une vérification a échoué.
 */

#pragma once

#include <cmath>
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace test {

using CasDeTest = std::pair<const char*, std::function<void()>>;

inline std::vector<CasDeTest>& registre() {
    static std::vector<CasDeTest> cas;
    return cas;
}

inline int& echecsDuCasCourant() {
    static int echecs = 0;
    return echecs;
}

struct Enregistreur {
    Enregistreur(const char* nom, std::function<void()> fonction) {
        registre().emplace_back(nom, std::move(fonction));
    }
};

inline void signalerEchec(const char* fichier, int ligne, const std::string& message) {
    ++echecsDuCasCourant();
    std::cerr << "    ECHEC " << fichier << ":" << ligne << " : " << message << "\n";
}

/*
 * Nom : executerTous
 * Description : Exécute tous les cas enregistrés, affiche un bilan et renvoie
 *               0 si tout est vert, 1 sinon. Une exception inattendue fait
 *               échouer le cas sans interrompre les suivants.
 */
inline int executerTous() {
    int casEchoues = 0;
    for (const auto& [nom, fonction] : registre()) {
        echecsDuCasCourant() = 0;
        try {
            fonction();
        } catch (const std::exception& e) {
            signalerEchec(nom, 0, std::string("exception inattendue : ") + e.what());
        } catch (...) {
            signalerEchec(nom, 0, "exception inattendue de type inconnu");
        }
        const bool ok = echecsDuCasCourant() == 0;
        if (!ok) ++casEchoues;
        std::cout << (ok ? "[ OK ] " : "[FAIL] ") << nom << "\n";
    }
    std::cout << "\n" << (registre().size() - casEchoues) << "/" << registre().size()
              << " cas de test réussis\n";
    return casEchoues == 0 ? 0 : 1;
}

} // namespace test

#define TEST_CONCAT_IMPL(a, b) a##b
#define TEST_CONCAT(a, b) TEST_CONCAT_IMPL(a, b)

#define TEST_CASE(nom)                                                          \
    static void nom();                                                          \
    static const test::Enregistreur TEST_CONCAT(enregistreur_, nom)(#nom, nom); \
    static void nom()

#define CHECK(condition)                                                        \
    do {                                                                        \
        if (!(condition)) test::signalerEchec(__FILE__, __LINE__, #condition);  \
    } while (0)

#define CHECK_NEAR(obtenu, attendu, tolerance)                                  \
    do {                                                                        \
        const double test_o_ = (obtenu);                                        \
        const double test_a_ = (attendu);                                       \
        if (!(std::abs(test_o_ - test_a_) <= (tolerance))) {                    \
            std::ostringstream test_msg_;                                       \
            test_msg_.precision(17);                                            \
            test_msg_ << #obtenu << " = " << test_o_ << ", attendu " << test_a_ \
                      << " (tolérance " << (tolerance) << ")";                  \
            test::signalerEchec(__FILE__, __LINE__, test_msg_.str());           \
        }                                                                       \
    } while (0)

#define CHECK_EQ(obtenu, attendu)                                               \
    do {                                                                        \
        const auto test_o_ = (obtenu);                                          \
        const auto test_a_ = (attendu);                                         \
        if (!(test_o_ == test_a_)) {                                            \
            std::ostringstream test_msg_;                                       \
            test_msg_ << #obtenu << " = \"" << test_o_ << "\", attendu \""      \
                      << test_a_ << "\"";                                       \
            test::signalerEchec(__FILE__, __LINE__, test_msg_.str());           \
        }                                                                       \
    } while (0)

#define CHECK_THROWS(instruction, TypeException)                                \
    do {                                                                        \
        bool test_leve_ = false;                                                \
        try {                                                                   \
            instruction;                                                        \
        } catch (const TypeException&) {                                        \
            test_leve_ = true;                                                  \
        } catch (...) {                                                         \
        }                                                                       \
        if (!test_leve_)                                                        \
            test::signalerEchec(__FILE__, __LINE__,                             \
                                #instruction " n'a pas levé " #TypeException);  \
    } while (0)
