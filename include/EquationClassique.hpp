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
    // Construit l'équation autour de son arborescence
    explicit EquationClassique(ExprPtr racine);
    
    // Constructeur par défaut (zéro)
    EquationClassique();
    
    ~EquationClassique() override = default;

    // Surcharge de l'interface parente (évaluation simple)
    double eval(double x) const override;
    
    // Génère formellement la dérivée (retourne un nouvel objet dérivé alloué)
    EquationClassique* derivee() override;

    // Méthode propre pour forcer la factorisation et simplification de l'arbre
    void simplifier();

    // Affiche l'équation symbolique
    void afficher() const;

private:
    ExprPtr m_racine;
};
