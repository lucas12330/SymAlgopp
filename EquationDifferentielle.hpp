/**
 * @file EquationDifferentielle.hpp
 * @author Lucas Bezanilla
 * @date 2026-04-14
 * @brief Fichier d'en-tête pour la sous-classe EquationDifferentielle.
 *
 * Projet SymAlgo++ :
 * Cette classe hérite de la classe de base Equation. Elle est dédiée
 * à la représentation et la résolution (numérique ou symbolique) des
 * équations différentielles.
 */

#pragma once

#include "Equation.hpp"

class EquationDifferentielle : public Equation {
public:
    EquationDifferentielle();
    ~EquationDifferentielle() override = default;

    /*
     * TODO: Surcharger les méthodes de l'interface et ajouter les attributs
     * propres aux EDO (ex: ordre de la dérivée, variables intégrées).
     */
};


