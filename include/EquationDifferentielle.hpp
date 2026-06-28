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
#include <vector>
#include <Eigen/Dense>

class EquationClassique;

class EquationDifferentielle : public Equation {
public:
  /*
   * Nom : EquationDifferentielle
   * Description : Constructeur par défaut de l'équation différentielle.
   * Utilisation : EquationDifferentielle eq_diff;
   */
  EquationDifferentielle();

  /*
   * Nom : ~EquationDifferentielle
   * Description : Destructeur par défaut.
   * Utilisation : Appelé automatiquement à la destruction de l'objet.
   */
  ~EquationDifferentielle() override = default;

  /*
   * Nom : eval
   * Description : Evalue l'équation différentielle numériquement (Runge-Kutta 4) depuis x=0.
   * Utilisation : double resultat = eq_diff.eval(valeur);
   */
  double eval(double x) const override;
  
  /*
   * Nom : ajouterTerme
   * Description : Ajoute un terme à l'équation différentielle selon le rang de dérivation et son coefficient.
   * Utilisation : eq_diff.ajouterTerme(2, 5.0); // Ajoute 5.0 * y''
   */
  void ajouterTerme(unsigned int rang, double coeff);

  /*
   * Nom : afficher
   * Description : Affiche proprement l'équation différentielle.
   * Utilisation : eq_diff.afficher();
   */
  void afficher() const;

  /*
   * Nom : setConditionsInitiales
   * Description : Définit les conditions initiales pour la résolution numérique.
   * Utilisation : eq_diff.setConditionsInitiales({y0, y'0, ...});
   */
  void setConditionsInitiales(const std::vector<double>& ci);

  /*
   * Nom : getMatriceCompagnon
   * Description : Renvoie la matrice d'état A associée au système d'ordre 1 (Y' = A*Y).
   */
  Eigen::MatrixXd getMatriceCompagnon() const;

  /*
   * Nom : resoudreLitteral
   * Description : Résout analytiquement l'équation linéaire et retourne la solution formelle.
   */
  EquationClassique* resoudreLitteral() const;

private:
  std::map<unsigned int, double> m_terme;
  std::vector<double> m_conditions_initiales;
};
