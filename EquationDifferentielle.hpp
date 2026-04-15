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
#include <map>

class EquationDifferentielle : public Equation {
public:
  EquationDifferentielle();
  ~EquationDifferentielle() override = default;
  double eval(double x) const override;
  
  // Fonction d'ajout des terme et coefficient.
  void ajouterTerme(unsigned int rang, double coeff);
  // Fonction simple d'affichage
  void afficher() const;

  /*
   * TODO: Surcharger les méthodes de l'interface et ajouter les attributs
   * propres aux EDO (ex: ordre de la dérivée, variables intégrées).
   */
private:
  std::map<unsigned int, double> m_terme;
};
