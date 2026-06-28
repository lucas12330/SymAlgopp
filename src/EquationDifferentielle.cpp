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
#include "EquationClassique.hpp"
#include "ASTNode.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>

EquationDifferentielle::EquationDifferentielle() {}

void EquationDifferentielle::ajouterTerme(unsigned int rang, double coeff) {
  m_terme[rang] += coeff;
}

void EquationDifferentielle::setConditionsInitiales(const std::vector<double>& ci) {
    m_conditions_initiales = ci;
}

void EquationDifferentielle::afficher() const {
    bool first = true;
    for (auto it = m_terme.rbegin(); it != m_terme.rend(); ++it) {
        if (it->second == 0.0) continue;
        
        if (!first) {
            if (it->second > 0) std::cout << " + ";
            else std::cout << " - ";
        } else {
            if (it->second < 0) std::cout << "-";
        }
        
        double abs_coeff = std::abs(it->second);
        if (abs_coeff != 1.0 || it->first == 0) std::cout << abs_coeff;
        if (abs_coeff != 1.0 && it->first > 0) std::cout << "*";
        
        if (it->first == 0) {
            if (abs_coeff == 1.0) std::cout << "y";
            else std::cout << "*y";
        } else if (it->first == 1) {
            std::cout << "y'";
        } else if (it->first == 2) {
            std::cout << "y''";
        } else {
            std::cout << "y^(" << it->first << ")";
        }
        first = false;
    }
    if (first) std::cout << "0";
    std::cout << " = 0" << std::endl;
}

Eigen::MatrixXd EquationDifferentielle::getMatriceCompagnon() const {
    if (m_terme.empty()) return Eigen::MatrixXd::Zero(1, 1);
    unsigned int n = m_terme.rbegin()->first;
    if (n == 0) return Eigen::MatrixXd::Zero(1, 1);
    
    double an = m_terme.rbegin()->second;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
    for (unsigned int i = 0; i < n - 1; ++i) {
        A(i, i + 1) = 1.0;
    }
    for (unsigned int j = 0; j < n; ++j) {
        double aj = 0.0;
        if (m_terme.count(j)) aj = m_terme.at(j);
        A(n - 1, j) = -aj / an;
    }
    return A;
}

EquationClassique* EquationDifferentielle::resoudreLitteral() const {
    if (m_terme.empty()) return new EquationClassique(cst(0.0));
    unsigned int n = m_terme.rbegin()->first;
    if (n == 0) return new EquationClassique(cst(0.0));
    
    Eigen::MatrixXd A = getMatriceCompagnon();
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    
    auto X = var("x");
    ExprPtr solution = cst(0.0);
    auto eigenvalues = solver.eigenvalues();
    int var_idx = 1;
    
    for (int i = 0; i < eigenvalues.size(); ++i) {
        std::complex<double> lambda = eigenvalues[i];
        
        if (std::abs(lambda.imag()) < 1e-6) {
            std::string c_name = "C" + std::to_string(var_idx++);
            ExprPtr terme = var(c_name) * ast_exp(cst(lambda.real()) * X);
            solution = solution + terme;
        } else {
            if (lambda.imag() > 0) {
                std::string c_name1 = "C" + std::to_string(var_idx++);
                std::string c_name2 = "C" + std::to_string(var_idx++);
                ExprPtr exp_part = ast_exp(cst(lambda.real()) * X);
                ExprPtr trig_part = var(c_name1) * ast_cos(cst(lambda.imag()) * X) + var(c_name2) * ast_sin(cst(lambda.imag()) * X);
                
                if (std::abs(lambda.real()) < 1e-6) {
                    solution = solution + trig_part;
                } else {
                    solution = solution + (exp_part * trig_part);
                }
            }
        }
    }
    return new EquationClassique(solution->simplifier());
}

double EquationDifferentielle::eval(double x) const {
    unsigned int n = m_terme.empty() ? 0 : m_terme.rbegin()->first;
    if (n == 0) return 0.0;
    
    Eigen::VectorXd Y = Eigen::VectorXd::Zero(n);
    for (size_t i = 0; i < std::min((size_t)n, m_conditions_initiales.size()); ++i) {
        Y(i) = m_conditions_initiales[i];
    }
    
    if (std::abs(x) < 1e-9) return Y(0);
    
    Eigen::MatrixXd A = getMatriceCompagnon();
    int steps = std::max(100, (int)(std::abs(x) / 0.01));
    double h = x / steps;
    
    for (int i = 0; i < steps; ++i) {
        Eigen::VectorXd k1 = A * Y;
        Eigen::VectorXd k2 = A * (Y + 0.5 * h * k1);
        Eigen::VectorXd k3 = A * (Y + 0.5 * h * k2);
        Eigen::VectorXd k4 = A * (Y + h * k3);
        Y += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    
    return Y(0);
}