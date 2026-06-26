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
#include <iostream>

/*
 * Nom : EquationDifferentielle
 * Description : Constructeur par défaut de la classe.
 * Utilisation : EquationDifferentielle eq;
 */
EquationDifferentielle::EquationDifferentielle() {}

/*
 * Nom : ajouterTerme
 * Description : Ajoute ou modifie le coefficient d'un terme selon son rang (dérivée).
 * Utilisation : eq.ajouterTerme(1, 2.5); // Ajoute 2.5 * y'
 */
void EquationDifferentielle::ajouterTerme(unsigned int rang, double coeff) {
  m_terme[rang] += coeff;
}

/*
 * Nom : afficher
 * Description : Parcourt et affiche chaque terme de l'équation différentielle.
 * Utilisation : eq.afficher();
 */
void EquationDifferentielle::afficher() const {
  for (auto const &[rang, coeff] : m_terme) {
    std::cout << coeff << "*y^" << rang << " + ";
  }
  std::cout << "0 = 0" << std::endl;
}

/*
 * Nom : eval
 * Description : Evalue l'équation (implémentation minimale pour l'instant).
 * Utilisation : double val = eq.eval(x);
 */
double EquationDifferentielle::eval(double x) const {
  return 0.0; // À implémenter plus tard avec une vraie logique
}