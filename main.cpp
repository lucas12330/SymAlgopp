/**
 * @file main.cpp
 * @brief Programme de test pour les fonctionnalités des Équations Différentielles
 */

#include <iostream>
#include "EquationDifferentielle.hpp"

int main() {
    std::cout << "===== DEMARRAGE DES TESTS : EQUATIONS DIFFERENTIELLES =====\n\n";

    // --- TEST 1 : Création simple (ex: Oscillateur Harmonique) ---
    std::cout << "Test 1 : Oscillateur Harmonique (y'' + 4y = 0)" << std::endl;
    EquationDifferentielle eqHarmonique;
    eqHarmonique.ajouterTerme(2, 1.0); // 1 * y''
    eqHarmonique.ajouterTerme(0, 4.0); // 4 * y
    std::cout << "Affichage attendu : 4*y + 1*y^(2) = 0" << std::endl;
    std::cout << "Affichage obtenu  : ";
    eqHarmonique.afficher();
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 2 : Équation plus complexe avec des nombres négatifs ---
    std::cout << "Test 2 : Équation plus complexe (3y''' - 2y' + 5y = 0)" << std::endl;
    EquationDifferentielle eqComplexe;
    eqComplexe.ajouterTerme(3, 3.0);  // 3 * y'''
    eqComplexe.ajouterTerme(1, -2.0); // -2 * y'
    eqComplexe.ajouterTerme(0, 5.0);  // 5 * y
    // N.B: Ton code actuel affiche "-2*y^1". Le test valide que l'ordre est bon !
    std::cout << "Affichage obtenu  : ";
    eqComplexe.afficher();
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 3 : Fusion de termes (utilisation du +=) ---
    std::cout << "Test 3 : Ajout multiple sur la même dérivée" << std::endl;
    EquationDifferentielle eqFusion;
    eqFusion.ajouterTerme(1, 2.0); // 2 * y'
    eqFusion.ajouterTerme(1, 3.5); // On ajoute 3.5 * y' -> devrait donner 5.5 * y'
    std::cout << "Affichage attendu : 5.5*y^1 = 0" << std::endl;
    std::cout << "Affichage obtenu  : ";
    eqFusion.afficher();
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 4 : Méthode eval() héritée ---
    std::cout << "Test 4 : Appel de la méthode virtual eval()" << std::endl;
    std::cout << "eval(42.0) = " << eqFusion.eval(42.0) << " (Attendu: 0)" << std::endl;
    std::cout << "--------------------------------------------------------\n\n";

    std::cout << "===== FIN DES TESTS =====\n";
    return 0;
}
