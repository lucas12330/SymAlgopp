/**
 * @file EquationDifferentielle.cpp
 * @author Lucas Bezanilla
 * @date 2026-04-15
 * @brief Fichier d'en-tête pour la sous-classe EquationDifferentielle.
 *
 * Projet SymAlgo++ :
 * Cette classe hérite de la classe de base Equation. Elle est dédiée
 * à la représentation et la résolution (numérique ou symbolique) des
 * équations différentielles.
 */

#include "EquationDifferentielle.hpp"

EquationDifferentielle::EquationDifferentielle() {}

void EquationDifferentielle::ajouterTerme(unsigned int rang, double coeff) {
  m_terme[rang] += coeff;
}

void EquationDifferentielle::afficher() const {
  for (auto const &[rang, coeff] : m_terme) {
    std::cout << coeff << "*y^" << rang << " + ";
  }
  std::cout << "0 = 0" << std::endl;
}
double EquationDifferentielle::eval(double x) const {
  return 0.0; // À implémenter plus tard avec une vraie logique
}