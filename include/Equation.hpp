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

#include <memory>

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
   * Nom : deriveeGenerique
   * Description : Version polymorphe de la dérivation, pour manipuler n'importe quelle
   *               équation via la classe de base. Chaque classe dérivée propose aussi
   *               derivee(), qui renvoie directement son propre type par valeur.
   * Utilisation : std::unique_ptr<Equation> d = equation.deriveeGenerique(); d->eval(x);
   */
  virtual std::unique_ptr<Equation> deriveeGenerique() const = 0;
};
