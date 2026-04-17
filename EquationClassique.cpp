#include "EquationClassique.hpp"
#include <iostream>

EquationClassique::EquationClassique(ExprPtr racine) : m_racine(racine) {}

EquationClassique::EquationClassique() : m_racine(cst(0.0)) {}

double EquationClassique::eval(double x) const {
    if (m_racine) return m_racine->eval(x);
    return 0.0;
}

EquationClassique* EquationClassique::derivee() {
    if (m_racine) {
        // Appelle la dérivation symbolique puis simplifie l'expression immédiatement
        ExprPtr deriveeFormelle = m_racine->derivee()->simplifier();
        return new EquationClassique(deriveeFormelle);
    }
    return new EquationClassique(cst(0.0));
}

void EquationClassique::simplifier() {
    if (m_racine) {
        m_racine = m_racine->simplifier();
    }
}

void EquationClassique::afficher() const {
    if (m_racine) {
        m_racine->afficher(std::cout);
        std::cout << " = 0" << std::endl;
    } else {
        std::cout << "0 = 0" << std::endl;
    }
}
