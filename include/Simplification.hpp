/**
 * @file Simplification.hpp
 * @brief Options de simplification des expressions (ASTNode::simplifier).
 */

#pragma once

namespace symalgo {

/*
 * Nom : OptionsSimplification
 * Description : Règles facultatives appliquées par simplifier(options).
 *
 *   La forme canonique n'est PAS une option : elle est imposée à la construction de chaque
 *   noeud (x + x = 2*x, x*x^2 = x^3, x/x = 1, sin(pi/6) = 1/2...) et garantit qu'une
 *   expression n'existe qu'une fois en mémoire (hash-consing, égalité par adresse). Elle
 *   ne peut donc pas être désactivée ; aucune() ne garde qu'elle.
 *
 *   rationnel : un produit dont les facteurs à exposant entier sont des polynômes exacts en
 *     une seule variable est réduit par le PGCD de son numérateur et de son dénominateur :
 *     (x^2 - 1)/(x - 1) = x + 1, sin(x)*(x^2 - 1)/(x - 1) = sin(x)*(x + 1). Les facteurs
 *     modifiés sont écrits sous forme entière primitive (coefficient dominant positif, la
 *     constante remonte devant) ; les autres facteurs sont gardés tels quels. Comme
 *     cancel (SymPy) et normal (GiNaC), la réduction élargit le domaine : x + 1 est défini
 *     en x = 1, (x^2 - 1)/(x - 1) ne l'est pas. Sont laissés intacts : les expressions à
 *     plusieurs variables, les polynômes à coefficients réels (non exacts) ou symboliques
 *     (paramètres, pi...). Les sommes de quotients ne sont pas mises au même dénominateur.
 *
 * Utilisation : e->simplifier();                                    // options par défaut
 *               e->simplifier(OptionsSimplification::aucune());     // forme canonique seule
 *               OptionsSimplification o; o.rationnel = false; e->simplifier(o);
 */
struct OptionsSimplification {
    bool rationnel = true;

    /*
     * Nom : aucune
     * Description : Aucune règle facultative : l'expression est rendue telle quelle, sous
     *               la forme canonique qu'elle a depuis sa construction.
     */
    static OptionsSimplification aucune() {
        OptionsSimplification options;
        options.rationnel = false;
        return options;
    }
};

} // namespace symalgo
