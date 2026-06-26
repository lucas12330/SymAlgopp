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
  /*
   * Nom : Equation
   * Description : Constructeur par défaut de la classe Equation.
   * Utilisation : Appelé automatiquement lors de l'instanciation des classes dérivées.
   */
  Equation() = default;

  /*
   * Nom : ~Equation
   * Description : Destructeur virtuel par défaut de la classe Equation (crucial pour le polymorphisme).
   * Utilisation : Appelé automatiquement à la destruction d'un objet dérivé de Equation.
   */
  virtual ~Equation() = default;

  /*
   * Nom : eval
   * Description : Evalue l'équation pour une valeur x donnée. Méthode virtuelle pure.
   * Utilisation : Implémentée par les classes dérivées. S'utilise via double resultat = obj.eval(valeur);
   */
  virtual double eval(double x) const = 0;

  /*
   * Nom : derivee
   * Description : Calcule la dérivée de l'équation.
   * Utilisation : Implémentée par les classes dérivées. S'utilise via Equation* derivee_eq = obj.derivee();
   */
  virtual Equation *derivee() { return nullptr; }
};
