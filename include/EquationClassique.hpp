/**
 * @file EquationClassique.hpp
 * @author Lucas Bezanilla (MODIFIED)
 * @brief Fichier d'en-tête pour la sous-classe EquationClassique.
 *
 * Cette classe a été restructurée pour opérer comme un wrapper (manipulateur)
 * dynamique autour d'un Abstract Syntax Tree (ASTNode), ce qui lui accorde
 * la résolution symbolique.
 */

#pragma once

#include "Equation.hpp"
#include "ASTNode.hpp"
#include "Evaluateur.hpp"
#include "Lecture.hpp"
#include "Solveur.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace symalgo {

class EquationClassique : public Equation {
public:
    /*
     * Nom : EquationClassique
     * Description : Construit l'équation autour de son arborescence (non nulle, sinon
     *               std::invalid_argument).
     * Utilisation : EquationClassique eq(racine_ast);
     */
    explicit EquationClassique(ExprPtr racine);
    explicit EquationClassique(std::nullptr_t) : EquationClassique(ExprPtr()) {} // lève, sans ambiguïté avec le texte

    /*
     * Nom : EquationClassique (texte)
     * Description : Équation lue depuis le texte (voir Lecture.hpp) : « gauche = droite »
     *               devient gauche - droite = 0. Lève ErreurLecture si le texte est mal formé.
     * Utilisation : EquationClassique eq("x^2 = 2");
     *               Solutions s = eq.resoudre();
     */
    explicit EquationClassique(const std::string& texte, const OptionsLecture& options = {});

    /*
     * Nom : EquationClassique
     * Description : Constructeur par défaut (initialise à zéro).
     * Utilisation : EquationClassique eq;
     */
    EquationClassique();

    /*
     * Nom : ~EquationClassique
     * Description : Destructeur par défaut.
     * Utilisation : Appelé automatiquement à la destruction de l'objet.
     */
    ~EquationClassique() override = default;

    /*
     * Nom : eval
     * Description : Surcharge de l'interface parente, évalue l'équation pour une valeur x donnée via l'arbre AST.
     * Utilisation : double resultat = eq.eval(valeur);
     */
    double eval(double x) const override;

    /*
     * Nom : eval (tableau)
     * Description : Évalue l'équation en chaque point, par l'évaluateur compilé et vectorisé
     *               (voir Evaluateur.hpp) : de loin le moyen le plus rapide pour de nombreux points.
     * Utilisation : std::vector<double> ys = eq.eval(xs);
     */
    std::vector<double> eval(const std::vector<double>& xs) const;

    /*
     * Nom : derivee
     * Description : Calcule la dérivée formelle (simplifiée) de l'équation.
     * Utilisation : EquationClassique d = eq.derivee();
     */
    EquationClassique derivee() const;

    /*
     * Nom : deriveeGenerique
     * Description : Dérivée renvoyée via l'interface polymorphe de Equation.
     * Utilisation : std::unique_ptr<Equation> d = eq.deriveeGenerique();
     */
    std::unique_ptr<Equation> deriveeGenerique() const override;

    /*
     * Nom : simplifier
     * Description : Méthode propre pour forcer la factorisation et simplification de l'arbre.
     * Utilisation : eq.simplifier();
     */
    void simplifier();

    /*
     * Nom : afficher
     * Description : Affiche l'équation symbolique sur la sortie standard.
     * Utilisation : eq.afficher();
     */
    void afficher() const;

    /*
     * Nom : integrer
     * Description : Calcule l'intégrale formelle de l'équation.
     * Utilisation : EquationClassique primitive = eq.integrer();
     */
    EquationClassique integrer() const;

    /*
     * Nom : limite
     * Description : Calcule la limite en x0 (avec L'Hôpital si besoin). Une limite qui n'a
     *               pas pu être déterminée est représentée par un noeud LimiteNonEvaluee.
     * Utilisation : EquationClassique l = eq.limite(x0);
     */
    EquationClassique limite(double x0) const;

    /*
     * Nom : DL
     * Description : Calcule le développement limité de l'équation.
     * Utilisation : EquationClassique dl = eq.DL(x0, ordre);
     */
    EquationClassique DL(double x0, int ordre) const;

    /*
     * Nom : resoudre
     * Description : Solutions réelles de l'équation (expression = 0) : exactes quand c'est
     *               possible, familles paramétrées par un entier pour les équations
     *               trigonométriques, avec un indicateur de complétude (voir Solveur.hpp).
     * Utilisation : Solutions s = eq.resoudre();
     */
    Solutions resoudre() const;

    /*
     * Nom : resoudre (intervalle)
     * Description : Solutions dans [a, b] : familles dépliées, recherche numérique en
     *               complément des formes non résolues exactement.
     * Utilisation : Solutions s = eq.resoudre(0.0, 10.0);
     */
    Solutions resoudre(double a, double b) const;

    /*
     * Nom : developper / factoriser
     * Description : Forme développée, et factorisation sur les rationnels (polynômes).
     * Utilisation : EquationClassique f = eq.factoriser();
     */
    EquationClassique developper() const;
    EquationClassique factoriser() const;

    /*
     * Nom : getExpression
     * Description : Donne accès à l'arbre (AST) de l'équation.
     * Utilisation : ExprPtr racine = eq.getExpression();
     */
    const ExprPtr& getExpression() const { return m_racine; }

    /*
     * Nom : genererPointsTrace
     * Description : Génère un tableau de points (x, y) optimisé par échantillonnage adaptatif.
     * Utilisation : auto pts = eq.genererPointsTrace(-10, 10, 0.05);
     */
    std::vector<std::pair<double, double>> genererPointsTrace(double xMin, double xMax, double tolerance = 1e-3) const;

private:
    // Nombre d'évaluations ponctuelles au-delà duquel l'expression est compilée : une
    // évaluation isolée ne paie pas le coût de la compilation
    static constexpr unsigned SEUIL_COMPILATION = 8;
    // Gain de partage minimal pour que l'évaluation ponctuelle compilée batte l'arbre
    // (mesuré : x3 pour des dérivées successives, gain 1,6 à 1,7 ; surcoût de 5 à 15 % pour
    // des expressions sans partage, gain 1 à 1,1)
    static constexpr double GAIN_PARTAGE_MIN = 1.3;

    const ProgrammeEvaluation& programme() const;

    ExprPtr m_racine;
    mutable std::shared_ptr<const ProgrammeEvaluation> m_programme; // compilé à la demande
    mutable unsigned m_evaluations = 0;
    mutable bool m_ponctuelCompile = false; // l'évaluation ponctuelle passe par le programme
    void echantillonnageAdaptatif(double x1, double y1, double x2, double y2, std::vector<std::pair<double, double>>& pts, double tolerance, int depth) const;
};

} // namespace symalgo
