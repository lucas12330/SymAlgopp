/**
 * @file Evaluateur.hpp
 * @brief Évaluation compilée des expressions.
 *
 * Une expression est compilée en une suite linéaire d'instructions exécutée par un petit
 * interpréteur, au lieu d'être parcourue récursivement (un appel virtuel et un saut de
 * pointeur par noeud) :
 *  - les sous-expressions partagées (hash-consing) ne sont calculées qu'une fois ;
 *  - les registres sont réutilisés dès qu'une valeur n'est plus utile ;
 *  - l'évaluation par blocs applique chaque instruction à un bloc de points avant de passer
 *    à la suivante : les boucles simples (additions, produits, puissances entières) sont
 *    vectorisées par le compilateur.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

class ProgrammeEvaluation {
public:
    /*
     * Nom : ProgrammeEvaluation
     * Description : Compile l'expression. Une expression non évaluable (paramètre, noeud non
     *               évalué) donne un programme invalide, dont l'évaluation lève
     *               std::logic_error comme l'évaluation directe.
     * Utilisation : ProgrammeEvaluation p(expr); double y = p.evaluer(2.0);
     */
    explicit ProgrammeEvaluation(const ExprPtr& expression);

    /*
     * Nom : evaluer
     * Description : Valeur de l'expression en x.
     */
    double evaluer(double x) const;

    /*
     * Nom : evaluer (par blocs)
     * Description : ys[i] = f(xs[i]) pour i < n, par blocs vectorisés.
     * Utilisation : p.evaluer(xs.data(), ys.data(), xs.size());
     */
    void evaluer(const double* xs, double* ys, std::size_t n) const;
    std::vector<double> evaluer(const std::vector<double>& xs) const;

    bool estValide() const { return m_erreur.empty(); }
    std::size_t nombreInstructions() const { return m_instructions.size(); }
    std::size_t nombreRegistres() const { return m_nombreRegistres; }

    /*
     * Nom : gainPartage
     * Description : Rapport entre le nombre de noeuds que visiterait l'évaluation récursive
     *               de l'arbre (un noeud partagé est compté autant de fois qu'il est utilisé)
     *               et le nombre d'instructions du programme. Au-delà de 1,3, le programme
     *               compilé est plus rapide même pour une évaluation ponctuelle.
     */
    double gainPartage() const { return m_gainPartage; }

private:
    enum class Code : std::uint8_t {
        Constante,  // r = valeur
        Variable,   // r = x
        Axpy,       // r = a + valeur * b   (terme d'une somme)
        Echelle,    // r = valeur * a
        Produit,    // r = a * b
        PuissanceEntiere, // r = a^entier
        Puissance,  // r = a^b
        Sinus,
        Cosinus,
        Tangente,
        Exponentielle,
        Logarithme,
        ArcSinus,
        ArcCosinus,
        ArcTangente,
    };

    struct Instruction {
        Code code;
        std::uint32_t destination = 0;
        std::uint32_t a = 0;
        std::uint32_t b = 0;
        double valeur = 0.0;
        long long entier = 0;
    };

    void verifierValidite() const;

    std::vector<Instruction> m_instructions;
    std::uint32_t m_resultat = 0;
    std::uint32_t m_nombreRegistres = 0;
    std::string m_erreur; // raison de l'invalidité, vide si le programme est valide
    double m_gainPartage = 1.0;

    friend class CompilateurExpression;
};

} // namespace symalgo
