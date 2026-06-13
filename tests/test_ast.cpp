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

    delete eq1_derive;
    delete eq2_derive;

    std::cout << "\n===== FIN DES TESTS =====\n";
    return 0;
}
