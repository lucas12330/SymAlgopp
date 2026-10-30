// Programme consommateur : vérifie qu'un projet extérieur peut utiliser la bibliothèque
// installée (en-tête unique, GMP, Eigen). Code de retour non nul au premier écart.
#include <symalgopp>

#include <cmath>
#include <iostream>

using namespace symalgo;

int main() {
    // fonctions de plusieurs variables : d/dy (x^2 y + sin(x y)) = x^2 + x cos(x y)
    const ExprPtr f = lire("x^2*y + sin(x*y)", {"x", {"y"}});
    const double v = evaluer(f->derivee("y"), {{"x", 2.0}, {"y", 0.5}});
    if (std::abs(v - (4.0 + 2.0 * std::cos(1.0))) > 1e-12) return 1;

    // nombres exacts (GMP) : la puissance de 2 dépasse 64 bits
    if (lire("2^100 + 1")->texte() != "1267650600228229401496703205377") return 2;

    // résolution exacte
    if (EquationClassique("x^2 = 2").resoudre().liste.size() != 2) return 3;

    // Eigen dans l'API publique
    EquationDifferentielle eq;
    eq.ajouterTerme(2, 1.0);
    eq.ajouterTerme(0, 4.0);
    const Eigen::MatrixXd A = eq.getMatriceCompagnon();
    if (A.rows() != 2 || A(1, 0) != -4.0) return 4;

    std::cout << "symalgopp : installation verifiee\n";
    return 0;
}
