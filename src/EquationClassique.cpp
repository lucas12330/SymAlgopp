#include "EquationClassique.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace symalgo {

/*
 * Nom : EquationClassique
 * Description : Constructeur de l'équation avec une racine AST spécifique.
 * Utilisation : EquationClassique eq(mon_ast);
 */
EquationClassique::EquationClassique(ExprPtr racine) : m_racine(std::move(racine)) {
    if (!m_racine) throw std::invalid_argument("EquationClassique : expression nulle");
}

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
double EquationClassique::eval(double x) const { return m_racine->eval(x); }

/*
 * Nom : derivee
 * Description : Calcule la dérivée symbolique simplifiée de l'équation.
 * Utilisation : EquationClassique d = eq.derivee();
 */
EquationClassique EquationClassique::derivee() const {
    return EquationClassique(m_racine->derivee()->simplifier());
}

std::unique_ptr<Equation> EquationClassique::deriveeGenerique() const {
    return std::make_unique<EquationClassique>(derivee());
}

/*
 * Nom : simplifier
 * Description : Simplifie l'expression mathématique de l'équation (factorisation, etc.).
 * Utilisation : eq.simplifier();
 */
void EquationClassique::simplifier() { m_racine = m_racine->simplifier(); }

/*
 * Nom : afficher
 * Description : Affiche l'équation symbolique dans le terminal suivie de " = 0".
 * Utilisation : eq.afficher();
 */
void EquationClassique::afficher() const {
    m_racine->afficher(std::cout);
    std::cout << " = 0" << std::endl;
}

/*
 * Nom : integrer
 * Description : Calcule une primitive symbolique simplifiée (IntegraleNonEvaluee si aucune
 *               règle ne s'applique).
 * Utilisation : EquationClassique F = eq.integrer();
 */
EquationClassique EquationClassique::integrer() const {
    return EquationClassique(m_racine->integrer()->simplifier());
}

/*
 * Nom : limite
 * Description : Calcule la limite en x0 (LimiteNonEvaluee si elle n'est pas déterminée).
 * Utilisation : EquationClassique l = eq.limite(0.0);
 */
EquationClassique EquationClassique::limite(double x0) const {
    return EquationClassique(m_racine->limite(x0)->simplifier());
}

/*
 * Nom : DL
 * Description : Développement limité de Taylor en x0 à l'ordre donné.
 * Utilisation : EquationClassique dl = eq.DL(0.0, 3);
 */
EquationClassique EquationClassique::DL(double x0, int ordre) const {
    return EquationClassique(m_racine->DL(x0, ordre));
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

} // namespace symalgo
