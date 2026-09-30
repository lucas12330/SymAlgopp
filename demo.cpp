#include <iostream>
#include <iomanip>
#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"
#include "ASTNode.hpp"
#include "Lecture.hpp"
#include "Solveur.hpp"

using namespace symalgo;

// Exemple 1 : Utilisation des expressions symboliques classiques (Physique)
void exemple_physique() {
    std::cout << "--- 1. Exemple Physique : Energie Cinetique ---\n";
    // E_k(v) = 0.5 * m * v^2 (On prend m = 10 kg)
    auto V = var("v");
    auto expr = 0.5 * cst(10.0) * ast_pow(V, 2); 
    
    EquationClassique eq(expr);
    std::cout << "Equation de l'energie cinetique E(v) : ";
    eq.afficher();
    
    std::cout << "Evaluation pour v = 10 m/s : " << eq.eval(10.0) << " Joules\n";
    
    // La dérivée de l'énergie cinétique par rapport à la vitesse est la quantité de mouvement (p = m*v)
    const EquationClassique eq_derivee = eq.derivee();
    std::cout << "Derivee (Quantite de mouvement p = dE/dv) : ";
    eq_derivee.afficher();
    std::cout << "Evaluation pour v = 10 m/s : " << eq_derivee.eval(10.0) << " kg.m/s\n\n";
}

// Exemple 2 : Approximation symbolique (Développement limité)
void exemple_math_dl() {
    std::cout << "--- 2. Exemple Mathematique : Developpement Limite ---\n";
    // Fonction f(x) = exp(x) * sin(x)
    auto X = var("x");
    auto expr = ast_exp(X) * ast_sin(X);
    
    EquationClassique f(expr);
    std::cout << "Fonction f(x) : ";
    f.afficher();
    
    // On demande le DL en 0 à l'ordre 3
    const EquationClassique dl = f.DL(0.0, 3);
    std::cout << "Developpement limite en 0 a l'ordre 3 : ";
    dl.afficher();
    std::cout << "\n";
}

// Exemple 3 : Équations Différentielles (Ressort avec Frottement)
void exemple_edo() {
    std::cout << "--- 3. Exemple EDO : Ressort Amorti (m.x'' + c.x' + k.x = 0) ---\n";
    // Paramètres : masse m = 1, amortissement c = 0.2, raideur k = 4
    EquationDifferentielle eqRessort;
    eqRessort.ajouterTerme(2, 1.0);  // 1.0 * x''
    eqRessort.ajouterTerme(1, 0.2);  // 0.2 * x'
    eqRessort.ajouterTerme(0, 4.0);  // 4.0 * x
    
    std::cout << "Modele du systeme : ";
    eqRessort.afficher();
    
    // Conditions initiales : position x(0) = 1m, vitesse x'(0) = 0m/s
    eqRessort.setConditionsInitiales({1.0, 0.0});
    
    std::cout << "Position a t=2s (solution numerique RK4) : " << eqRessort.eval(2.0) << " m\n";
    
    // Résolution analytique littérale
    const EquationClassique eqGenerale = eqRessort.resoudreLitteral();
    std::cout << "Solution analytique generale x(t) = ";
    eqGenerale.afficher();

    // Solution exacte avec les conditions initiales, comparable a RK4
    const EquationClassique eqExacte = eqRessort.resoudreProblemeCauchy();
    std::cout << "Solution exacte avec x(0) = 1, x'(0) = 0 : ";
    eqExacte.afficher();
    std::cout << "Position a t=2s (solution exacte)        : " << eqExacte.eval(2.0) << " m\n";
}

// Exemple 4 : Équations saisies en texte, résolues exactement
void exemple_resolution() {
    std::cout << "\n--- 4. Exemple : Lecture et Resolution d'Equations ---\n";
    for (const char* texte : {"x^3 - x^2 - 2x + 2 = 0", "exp(2x) - 3exp(x) + 2 = 0", "sin(x) = 1/2", "2^x = 8"}) {
        const EquationClassique eq(texte);
        const Solutions s = eq.resoudre();
        std::cout << texte << "  :  ";
        for (std::size_t i = 0; i < s.liste.size(); ++i) {
            std::cout << (i ? " ; " : "") << "x = " << s.liste[i].valeur;
            if (s.liste[i].estFamille()) std::cout << " (k entier)";
        }
        std::cout << "\n";
    }
    // Sans forme exacte : solution numérique sur un intervalle
    const Solutions s = EquationClassique("exp(x) + x = 0").resoudre(-5.0, 5.0);
    std::cout << "exp(x) + x = 0  :  x ~ " << std::setprecision(15) << s.liste[0].approximation << "\n";
}

int main() {
    std::cout << "========================================\n";
    std::cout << "=== TEST EN CONDITION REELLE : SymAlgo++ \n";
    std::cout << "========================================\n\n";
    
    exemple_physique();
    exemple_math_dl();
    exemple_edo();
    exemple_resolution();
    
    return 0;
}
