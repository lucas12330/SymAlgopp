/**
 * @file Equation.hpp
 * @author Lucas Bezanilla
 * @date 2026-04-14
 * @brief Fichier d'en-tête pour la classe virtuelle de base Equation.
 *
 * Projet SymAlgo++ :
 * Ce fichier définit la "super classe" abstraite de laquelle hériteront toutes
 * les autres formes d'équations (Classiques, Différentielles).
 * Elle a pour rôle de définir l'interface commune (virtuelle pure) pour des
 * actions comme la résolution, l'intégration, la dérivation et le tracé.
 */

#pragma once

#include <vector>

class Equation {
public:
  // Constructeur par défaut et Destructeur virtuel (crucial pour le
  // polymorphisme)
  Equation() = default;
  virtual ~Equation() = default;

  /*
   * TODO: Définir ici les méthodes virtuelles pures de l'interface, par exemple
   * : virtual double eval(double x) const = 0; virtual Equation* derivee()
   * const = 0;
   */
};
