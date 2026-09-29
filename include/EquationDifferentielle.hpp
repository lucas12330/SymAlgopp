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

#include "ASTNode.hpp"
#include "Equation.hpp"
#include "EquationClassique.hpp"

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <vector>

namespace symalgo {

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
   * Nom : derivee
   * Description : Équation dont la solution est y' : si y est solution, y' vérifie la même
   *               équation (linéaire homogène à coefficients constants) avec les conditions
   *               initiales (y'(0), ..., y^(n)(0)), la dernière étant tirée de l'équation.
   * Utilisation : EquationDifferentielle dy = eq.derivee(); double pente = dy.eval(x);
   */
  EquationDifferentielle derivee() const;

  /*
   * Nom : deriveeGenerique
   * Description : Dérivée renvoyée via l'interface polymorphe de Equation.
   * Utilisation : std::unique_ptr<Equation> d = eq.deriveeGenerique();
   */
  std::unique_ptr<Equation> deriveeGenerique() const override;
  
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
   * Description : Résout analytiquement l'équation linéaire et retourne la solution générale,
   *               combinaison linéaire des solutions de base avec des constantes symboliques
   *               C1..Cn (paramètres). Les racines multiples du polynôme caractéristique
   *               donnent les termes x^k e^(rx).
   * Utilisation : EquationClassique sol = eq.resoudreLitteral();
   */
  EquationClassique resoudreLitteral() const;

  /*
   * Nom : resoudreProblemeCauchy
   * Description : Solution exacte satisfaisant les conditions initiales en x = 0
   *               (y(0), y'(0), ... ; les conditions manquantes valent 0). Contrairement à
   *               resoudreLitteral, le résultat est directement évaluable.
   * Utilisation : eq.setConditionsInitiales({1.0, 0.0}); EquationClassique y = eq.resoudreProblemeCauchy();
   */
  EquationClassique resoudreProblemeCauchy() const;

private:
  /*
   * Nom : baseDeSolutions
   * Description : Fonctions de base de l'espace des solutions (une par unité d'ordre).
   */
  std::vector<ExprPtr> baseDeSolutions() const;

  std::map<unsigned int, double> m_terme;
  std::vector<double> m_conditions_initiales;
};

} // namespace symalgo
