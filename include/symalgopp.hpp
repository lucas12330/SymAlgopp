/**
 * @file symalgopp.hpp
 * @brief En-tête unique : toute l'API publique de SymAlgo++.
 *
 * Utilisation : #include <symalgopp>      (ou <symalgopp.hpp>)
 *
 * Expressions, nombres exacts, lecture depuis du texte, fonctions de plusieurs variables,
 * polynômes, solveur, évaluation compilée, équations classiques et différentielles. Canonique.hpp reste interne.
 * Les en-têtes peuvent aussi être inclus un par un, pour réduire le temps de compilation
 * (EquationDifferentielle.hpp inclut Eigen).
 */

#pragma once

#include "Nombre.hpp"
#include "ASTNode.hpp"
#include "Lecture.hpp"
#include "Multivariable.hpp"
#include "Polynome.hpp"
#include "Solveur.hpp"
#include "Evaluateur.hpp"
#include "Equation.hpp"
#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"
