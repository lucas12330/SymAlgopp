/**
 * @file Polynome.hpp
 * @brief Développement des expressions et polynômes à coefficients exacts.
 */

#pragma once

#include "ASTNode.hpp"

#include <utility>
#include <vector>

namespace symalgo {

/*
 * Nom : developper
 * Description : Développe l'expression : distribue les produits sur les sommes et
 *               développe les puissances entières positives de sommes, à tous les niveaux
 *               (y compris dans les arguments des fonctions). Le résultat est canonique :
 *               les termes semblables sont regroupés. Les dénominateurs sont développés
 *               mais restent au dénominateur.
 * Utilisation : developper(ast_pow(x + 1.0, 2.0))  // x^2 + 2*x + 1
 */
ExprPtr developper(const ExprPtr& e);

/*
 * Nom : RacineReelle
 * Description : Racine réelle d'un polynôme, certifiée : elle est la seule racine de
 *               l'intervalle ]gauche, droite]. Sa valeur exacte est donnée quand elle est
 *               connue (rationnelle, ou de degré 2 par radicaux), sinon approchée au double
 *               près.
 */
struct RacineReelle {
    ExprPtr exacte;     // nullptr si seule la valeur approchée est connue
    double valeur;      // approximation (au plus près en double)
    int multiplicite;
    Nombre gauche, droite;
};

/*
 * Nom : Polynome
 * Description : Polynôme d'une variable à coefficients exacts (Nombre), c[i] devant x^i.
 *               Les opérations sont exactes (rationnels de taille arbitraire).
 */
class Polynome {
public:
    Polynome() = default; // polynôme nul
    explicit Polynome(std::vector<Nombre> coefficients);

    /*
     * Nom : depuisExpression
     * Description : Reconnaît un polynôme dans une expression (développée au préalable).
     *               Les coefficients réels sont convertis en leur valeur exacte. Renvoie faux
     *               si l'expression n'est pas polynomiale à coefficients numériques. Si
     *               variable est fourni, il reçoit la variable rencontrée (nullptr si aucune).
     * Utilisation : Polynome p; if (Polynome::depuisExpression(e, p)) { ... }
     */
    static bool depuisExpression(const ExprPtr& e, Polynome& p, ExprPtr* variable = nullptr);

    int degre() const { return static_cast<int>(m_coefficients.size()) - 1; } // -1 : nul

    /*
     * Nom : estExact
     * Description : Faux si un coefficient provient d'un nombre réel (approché) : les racines
     *               restent certifiées numériquement, mais aucune forme exacte n'est donnée.
     */
    bool estExact() const { return m_exact; }
    bool estNul() const { return m_coefficients.empty(); }
    const Nombre& coefficient(int i) const;
    const Nombre& dominant() const { return m_coefficients.back(); }

    Polynome operator+(const Polynome& b) const;
    Polynome operator-(const Polynome& b) const;
    Polynome operator*(const Polynome& b) const;
    Polynome operator*(const Nombre& k) const;

    /*
     * Nom : diviser
     * Description : Division euclidienne exacte : *this = q * d + r, deg r < deg d.
     */
    void diviser(const Polynome& d, Polynome& quotient, Polynome& reste) const;

    Polynome derivee() const;
    Polynome unitaire() const; // coefficient dominant ramené à 1
    static Polynome pgcd(Polynome a, Polynome b);

    /*
     * Nom : sansCarre
     * Description : Décomposition sans facteur carré (algorithme de Yun) : le polynôme unitaire
     *               est le produit des f_i^i, f_i sans racine multiple. Renvoie les (f_i, i).
     */
    std::vector<std::pair<Polynome, int>> sansCarre() const;

    Nombre evaluer(const Nombre& x) const;
    double evaluer(double x) const;

    /*
     * Nom : nombreRacinesReelles
     * Description : Nombre exact de racines réelles distinctes dans ]a, b] (théorème de Sturm).
     */
    int nombreRacinesReelles(const Nombre& a, const Nombre& b) const;

    /*
     * Nom : racinesReelles
     * Description : Toutes les racines réelles, avec multiplicité, isolées par les suites de
     *               Sturm (liste certifiée complète), exactes quand c'est possible.
     */
    std::vector<RacineReelle> racinesReelles() const;

    /*
     * Nom : versExpression
     * Description : Expression sum c_i x^i dans la variable donnée.
     */
    ExprPtr versExpression(const ExprPtr& x) const;

private:
    void normaliser(); // retire les coefficients dominants nuls
    std::vector<Nombre> m_coefficients;
    bool m_exact = true;
};

/*
 * Nom : factoriser
 * Description : Factorise un polynôme sur les rationnels : facteurs linéaires issus des
 *               racines rationnelles (avec multiplicités) et facteurs restants sans carré,
 *               à coefficients entiers. Une expression non polynomiale est renvoyée
 *               inchangée. Limite de la forme canonique : un coefficient numérique devant
 *               un seul facteur somme est redistribué (2*(x - 1) s'écrit 2*x - 2).
 * Utilisation : factoriser(ast_pow(x, 3.0) - x)  // x*(x - 1)*(x + 1)
 */
ExprPtr factoriser(const ExprPtr& e);

} // namespace symalgo
