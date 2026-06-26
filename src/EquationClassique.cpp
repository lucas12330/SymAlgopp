#include "EquationClassique.hpp"
#include <iostream>

/*
 * Nom : EquationClassique
 * Description : Constructeur de l'équation avec une racine AST spécifique.
 * Utilisation : EquationClassique eq(mon_ast);
 */
EquationClassique::EquationClassique(ExprPtr racine) : m_racine(racine) {}

/*
 * Nom : EquationClassique
 * Description : Constructeur par défaut, initialise la racine à une constante 0.
 * Utilisation : EquationClassique eq;
 */
EquationClassique::EquationClassique() : m_racine(cst(0.0)) {}

/*
 * Nom : eval
 * Description : Calcule la valeur de l'équation en un point x.
 * Utilisation : double y = eq.eval(x);
 */
double EquationClassique::eval(double x) const {
    if (m_racine) return m_racine->eval(x);
    return 0.0;
}

/*
 * Nom : derivee
 * Description : Calcule la dérivée symbolique de l'équation et renvoie un nouveau pointeur d'EquationClassique.
 * Utilisation : EquationClassique* derivee_eq = eq.derivee();
 */
EquationClassique* EquationClassique::derivee() {
    if (m_racine) {
        // Appelle la dérivation symbolique puis simplifie l'expression immédiatement
        ExprPtr deriveeFormelle = m_racine->derivee()->simplifier();
        return new EquationClassique(deriveeFormelle);
    }
    return new EquationClassique(cst(0.0));
}

/*
 * Nom : simplifier
 * Description : Simplifie l'expression mathématique de l'équation (factorisation, etc.).
 * Utilisation : eq.simplifier();
 */
void EquationClassique::simplifier() {
    if (m_racine) {
        m_racine = m_racine->simplifier();
    }
}

/*
 * Nom : afficher
 * Description : Affiche l'équation symbolique dans le terminal suivie de " = 0".
 * Utilisation : eq.afficher();
 */
void EquationClassique::afficher() const {
    if (m_racine) {
        m_racine->afficher(std::cout);
        std::cout << " = 0" << std::endl;
    } else {
        std::cout << "0 = 0" << std::endl;
    }
}
