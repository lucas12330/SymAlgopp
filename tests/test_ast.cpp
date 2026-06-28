#include <iostream>
#include "EquationClassique.hpp"
#include "ASTNode.hpp"

int main() {
    std::cout << "===== DEMARRAGE DES TESTS : AST / EQUATION CLASSIQUE =====\n\n";

    // 1. Construit: x^2 + x*5 + 6
    auto X = var("x");
    auto eq1_expr = ast_pow(X, 2) + (X * 5) + 6;
    EquationClassique eq1(eq1_expr);
    
    std::cout << "Equation 1               : "; eq1.afficher();
    std::cout << "Evaluation eq1(x=2)      : " << eq1.eval(2.0) << " (Attendu: 20)\n";
    
    EquationClassique* eq1_derive = eq1.derivee();
    std::cout << "Derivee eq1 (simplifiee) : "; eq1_derive->afficher();
    std::cout << "Eval Derivee(x=2)        : " << eq1_derive->eval(2.0) << " (Attendu: 9)\n\n";

    // 2. Construit: 2*sin(x) + 3*sin(x)
    auto expr2 = (cst(2) * ast_sin(X)) + (cst(3) * ast_sin(X));
    EquationClassique eq2(expr2);
    std::cout << "Equation 2 (avant simp)  : "; eq2.afficher();
    eq2.simplifier();
    std::cout << "Equation 2 (apres simp)  : "; eq2.afficher();
    
    // Test derivee sin -> cos
    EquationClassique* eq2_derive = eq2.derivee();
    std::cout << "Derivee eq2              : "; eq2_derive->afficher();

    // 3. Test Integrer
    EquationClassique* eq1_int = eq1.integrer();
    std::cout << "\nIntegrale eq1            : "; eq1_int->afficher();

    // 4. Test Limite (L'Hopital 0/0)
    auto N = ast_sin(X);
    auto D = X;
    EquationClassique eq3(N / D);
    std::cout << "\nEquation 3 (N/D)         : "; eq3.afficher();
    EquationClassique* eq3_lim = eq3.limite(0.0);
    std::cout << "Limite eq3 en x=0        : "; eq3_lim->afficher();

    // 5. Test DL (Taylor)
    auto eq4_expr = ast_pow(X, 3) + ast_pow(X, 2) + X + 1;
    EquationClassique eq4(eq4_expr);
    std::cout << "\nEquation 4               : "; eq4.afficher();
    EquationClassique* eq4_dl = eq4.DL(0.0, 2);
    std::cout << "DL de eq4 en x=0 (ord 2) : "; eq4_dl->afficher();
    
    // 6. Test Generer points
    auto points = eq4.genererPointsTrace(0.0, 1.0, 0.1);
    std::cout << "\nTracé eq4 (nb points)    : " << points.size() << "\n";

    delete eq1_derive;
    delete eq2_derive;
    delete eq1_int;
    delete eq3_lim;
    delete eq4_dl;

    std::cout << "\n===== FIN DES TESTS =====\n";
    return 0;
}
