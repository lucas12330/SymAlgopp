#include <iostream>
#include <cmath>
#include "EquationDifferentielle.hpp"
#include "EquationClassique.hpp"

int main() {
    std::cout << "===== DEMARRAGE DES TESTS : EQUATIONS DIFFERENTIELLES =====\n\n";

    // --- TEST 1 : Oscillateur Harmonique (y'' + 4y = 0) ---
    std::cout << "Test 1 : Oscillateur Harmonique (y'' + 4y = 0)" << std::endl;
    EquationDifferentielle eqHarmonique;
    eqHarmonique.ajouterTerme(2, 1.0); // y''
    eqHarmonique.ajouterTerme(0, 4.0); // 4y
    std::cout << "Affichage attendu : y'' + 4*y = 0" << std::endl;
    std::cout << "Affichage obtenu  : ";
    eqHarmonique.afficher();
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 2 : Résolution Littérale (Oscillateur) ---
    std::cout << "Test 2 : Résolution Littérale de l'Oscillateur Harmonique" << std::endl;
    EquationClassique* eqSolHarmonique = eqHarmonique.resoudreLitteral();
    std::cout << "Solution exacte obtenue : ";
    eqSolHarmonique->afficher();
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 3 : Résolution Numérique (RK4) ---
    std::cout << "Test 3 : Évaluation Numérique (RK4)" << std::endl;
    // Conditions initiales pour cos(2x) : y(0) = 1, y'(0) = 0
    eqHarmonique.setConditionsInitiales({1.0, 0.0});
    
    double t = 3.14159265359 / 4.0; // t = pi/4
    // Solution exacte de y(t) = cos(2t), pour t = pi/4, y = cos(pi/2) = 0
    double y_num = eqHarmonique.eval(t);
    std::cout << "y(pi/4) via RK4  = " << y_num << std::endl;
    std::cout << "Attendu (RK4)    = 0.0 (approximatif)" << std::endl;
    std::cout << "--------------------------------------------------------\n\n";

    // --- TEST 4 : Équation plus complexe (3y''' - 2y' + 5y = 0) ---
    std::cout << "Test 4 : Équation complexe (3y''' - 2y' + 5y = 0)" << std::endl;
    EquationDifferentielle eqComplexe;
    eqComplexe.ajouterTerme(3, 3.0);  // 3 * y'''
    eqComplexe.ajouterTerme(1, -2.0); // -2 * y'
    eqComplexe.ajouterTerme(0, 5.0);  // 5 * y
    std::cout << "Affichage obtenu  : ";
    eqComplexe.afficher();
    
    EquationClassique* eqSolComplexe = eqComplexe.resoudreLitteral();
    std::cout << "Solution exacte : ";
    eqSolComplexe->afficher();
    std::cout << "--------------------------------------------------------\n\n";

    delete eqSolHarmonique;
    delete eqSolComplexe;

    std::cout << "===== FIN DES TESTS =====\n";
    return 0;
}
