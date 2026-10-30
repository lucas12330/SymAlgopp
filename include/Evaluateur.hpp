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
     * Nom : ProgrammeEvaluation (plusieurs variables, plusieurs sorties)
     * Description : Compile ensemble plusieurs expressions de plusieurs variables : les
     *               sous-expressions communes (typiquement celles d'un gradient ou d'une
     *               hessienne) ne sont calculées qu'une fois. entrees[i] donne l'ordre des
     *               variables ; une variable de l'expression absente de la liste donne un
     *               programme invalide. Lève std::invalid_argument si une entrée n'est pas
     *               une variable ou apparaît deux fois.
     * Utilisation : ProgrammeEvaluation g(gradient(f), {var("x"), var("y")});
     *               std::vector<double> dfdxdy = g.evaluerEn({1.0, 2.0});
     */
    ProgrammeEvaluation(const std::vector<ExprPtr>& expressions, const std::vector<ExprPtr>& entrees);

    /*
     * Nom : evaluer
     * Description : Valeur de l'expression en x (programme à une variable et une sortie :
     *               sinon std::logic_error).
     */
    double evaluer(double x) const;

    /*
     * Nom : evaluerEn
     * Description : Valeurs des sorties au point donné (une valeur par entrée, dans l'ordre
     *               de la liste des entrées). Lève std::invalid_argument si le nombre de
     *               valeurs est incorrect.
     * Utilisation : p.evaluerEn(valeurs.data(), sorties.data());
     */
    void evaluerEn(const double* valeurs, double* sorties) const;
    std::vector<double> evaluerEn(const std::vector<double>& valeurs) const;

    /*
     * Nom : evaluer (par blocs)
     * Description : ys[i] = f(xs[i]) pour i < n, par blocs vectorisés.
     * Utilisation : p.evaluer(xs.data(), ys.data(), xs.size());
     */
    void evaluer(const double* xs, double* ys, std::size_t n) const;
    std::vector<double> evaluer(const std::vector<double>& xs) const;

    std::size_t nombreEntrees() const { return m_nombreEntrees; }
    std::size_t nombreSorties() const { return m_resultats.size(); }
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
        Variable,   // r = entrée numéro `entier`
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
    // Exécute le programme pour un point ; r doit contenir m_nombreRegistres valeurs
    template <class Entree>
    void executer(const Entree& entree, double* r) const;

    std::vector<Instruction> m_instructions;
    std::vector<std::uint32_t> m_resultats; // registre de chaque sortie
    std::uint32_t m_premierResultat = 0;    // m_resultats[0], lu par les évaluations à une sortie
    bool m_scalaire = true;                 // une entrée et une sortie : evaluer(x) est permis
    std::uint32_t m_nombreEntrees = 1;
    std::uint32_t m_nombreRegistres = 0;
    std::string m_erreur; // raison de l'invalidité, vide si le programme est valide
    double m_gainPartage = 1.0;

    friend class CompilateurExpression;
};

} // namespace symalgo
