/**
 * @file EquationClassique.hpp
 * @author Lucas Bezanilla (MODIFIED)
 * @brief Fichier d'en-tête pour la sous-classe EquationClassique.
 *
 * Cette classe a été restructurée pour opérer comme un wrapper (manipulateur)
 * dynamique autour d'un Abstract Syntax Tree (ASTNode), ce qui lui accorde
 * la résolution symbolique.
 */

#pragma once

#include "Equation.hpp"
#include "ASTNode.hpp"

class EquationClassique : public Equation {
public:
    /*
     * Nom : EquationClassique
     * Description : Construit l'équation autour de son arborescence.
     * Utilisation : EquationClassique eq(racine_ast);
     */
    explicit EquationClassique(ExprPtr racine);
    
    /*
     * Nom : EquationClassique
     * Description : Constructeur par défaut (initialise à zéro).
     * Utilisation : EquationClassique eq;
     */
    EquationClassique();
    
    /*
     * Nom : ~EquationClassique
     * Description : Destructeur par défaut.
     * Utilisation : Appelé automatiquement à la destruction de l'objet.
     */
    ~EquationClassique() override = default;

    /*
     * Nom : eval
     * Description : Surcharge de l'interface parente, évalue l'équation pour une valeur x donnée via l'arbre AST.
     * Utilisation : double resultat = eq.eval(valeur);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Génère formellement la dérivée et retourne un nouvel objet dérivé alloué.
     * Utilisation : EquationClassique* eq_derivee = eq.derivee();
     */
    EquationClassique* derivee() override;

    /*
     * Nom : simplifier
     * Description : Méthode propre pour forcer la factorisation et simplification de l'arbre.
     * Utilisation : eq.simplifier();
     */
    void simplifier();

    /*
     * Nom : afficher
     * Description : Affiche l'équation symbolique sur la sortie standard.
     * Utilisation : eq.afficher();
     */
    void afficher() const;

    /*
     * Nom : integrer
     * Description : Calcule l'intégrale formelle de l'équation.
     * Utilisation : EquationClassique* eq_int = eq.integrer();
     */
    EquationClassique* integrer() const;

    /*
     * Nom : limite
     * Description : Calcule la limite symbolique (avec L'Hôpital si besoin).
     * Utilisation : EquationClassique* eq_lim = eq.limite(x0);
     */
    EquationClassique* limite(double x0) const;

    /*
     * Nom : DL
     * Description : Calcule le développement limité de l'équation.
     * Utilisation : EquationClassique* eq_dl = eq.DL(x0, ordre);
     */
    EquationClassique* DL(double x0, int ordre) const;

    /*
     * Nom : genererPointsTrace
     * Description : Génère un tableau de points (x, y) optimisé par échantillonnage adaptatif.
     * Utilisation : auto pts = eq.genererPointsTrace(-10, 10, 0.05);
     */
    std::vector<std::pair<double, double>> genererPointsTrace(double xMin, double xMax, double tolerance = 1e-3) const;

private:
    ExprPtr m_racine;
    void echantillonnageAdaptatif(double x1, double y1, double x2, double y2, std::vector<std::pair<double, double>>& pts, double tolerance, int depth) const;
};
