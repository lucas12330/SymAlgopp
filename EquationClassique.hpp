/**
 * @file EquationClassique.hpp
 * @author Lucas Bezanilla
 * @date 2026-04-14
 * @brief Fichier d'en-tête pour la sous-classe EquationClassique.
 *
 * Projet SymAlgo++ :
 * Cette classe hérite de la classe de base Equation. Elle a vocation à
 * représenter et manipuler des équations algébriques ou fonctions standards
 * (ex: polynômes, fonctions trigonométriques).
 * Elle implémente les méthodes héritées de façon spécifique à ces cas de
 * figure.
 */

#pragma once

#include "Equation.hpp"

class EquationClassique : public Equation {
public:
  EquationClassique();
  ~EquationClassique() override = default;

  /*
   * TODO: Implémenter et surcharger les méthodes virtuelles avec "override"
   * ex: double eval(double x) const override;
   */
};
