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

EquationClassique* EquationClassique::integrer() const {
    if (m_racine) {
        ExprPtr intFormelle = m_racine->integrer()->simplifier();
        return new EquationClassique(intFormelle);
    }
    return new EquationClassique(cst(0.0));
}

double EquationClassique::integraleDefinie(double a, double b) const {
    EquationClassique* primitive = this->integrer();
    double result = primitive->eval(b) - primitive->eval(a);
    delete primitive;
    return result;
}

EquationClassique* EquationClassique::limite(double x0) const {
    if (m_racine) {
        ExprPtr limFormelle = m_racine->limite(x0)->simplifier();
        return new EquationClassique(limFormelle);
    }
    return new EquationClassique(cst(0.0));
}

EquationClassique* EquationClassique::DL(double x0, int ordre) const {
    if (m_racine) {
        ExprPtr dlFormelle = m_racine->DL(x0, ordre);
        return new EquationClassique(dlFormelle);
    }
    return new EquationClassique(cst(0.0));
}

void EquationClassique::echantillonnageAdaptatif(double x1, double y1, double x2, double y2, std::vector<std::pair<double, double>>& pts, double tolerance, int depth) const {
    if (depth > 10) return; // Limite de récursion
    double xm = (x1 + x2) / 2.0;
    double ym = this->eval(xm);
    double y_interp = (y1 + y2) / 2.0;

    if (std::abs(ym - y_interp) > tolerance) {
        echantillonnageAdaptatif(x1, y1, xm, ym, pts, tolerance, depth + 1);
        pts.push_back({xm, ym});
        echantillonnageAdaptatif(xm, ym, x2, y2, pts, tolerance, depth + 1);
    }
}

std::vector<std::pair<double, double>> EquationClassique::genererPointsTrace(double xMin, double xMax, double tolerance) const {
    std::vector<std::pair<double, double>> points;
    if (xMin >= xMax) return points;

    // Découpage initial grossier (10 segments) pour éviter de rater les grandes variations
    int segments = 10;
    double pas = (xMax - xMin) / segments;
    
    double currX = xMin;
    double currY = this->eval(currX);
    points.push_back({currX, currY});

    for (int i = 1; i <= segments; ++i) {
        double nextX = xMin + i * pas;
        double nextY = this->eval(nextX);
        
        echantillonnageAdaptatif(currX, currY, nextX, nextY, points, tolerance, 0);
        points.push_back({nextX, nextY});
        
        currX = nextX;
        currY = nextY;
    }
    
    return points;
}
