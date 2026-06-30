#pragma once
#include <memory>
#include <string>
#include <ostream>

class ASTNode;
using ExprPtr = std::shared_ptr<ASTNode>;

enum class NodeType {
    Constante, Variable, Addition, Soustraction, Multiplication, Division, Puissance,
    Sinus, Cosinus, tangeante, Exponentielle, Logarithme
};

// Classe abstraite de base pour tous les noeuds de l'arbre
class ASTNode : public std::enable_shared_from_this<ASTNode> {
public:
    /*
     * Nom : ~ASTNode
     * Description : Destructeur virtuel par défaut de la classe ASTNode.
     * Utilisation : Appelé automatiquement à la destruction d'un objet ASTNode.
     */
    virtual ~ASTNode() = default;
    
    /*
     * Nom : eval
     * Description : Evalue l'expression représentée par le noeud pour une valeur donnée.
     * Utilisation : double resultat = noeud->eval(x);
     */
    virtual double eval(double x) const = 0;
    
    /*
     * Nom : DL
     * Description : Calcule le développement limité symbolique de l'expression en a, à l'ordre donné.
     * Utilisation : ExprPtr dl = noeud->DL(a, ordre);
     */
    virtual ExprPtr DL(double a, int ordre) const;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée symbolique de l'expression.
     * Utilisation : ExprPtr d = noeud->derivee();
     */
    virtual ExprPtr derivee() const = 0;
    
    /*
     * Nom : simplifier
     * Description : Simplifie mathématiquement l'expression de l'arbre.
     * Utilisation : ExprPtr simp = noeud->simplifier();
     */
    virtual ExprPtr simplifier() const = 0;
    
    /*
     * Nom : afficher
     * Description : Affiche le contenu textuel de l'expression sur un flux.
     * Utilisation : noeud->afficher(std::cout);
     */
    virtual void afficher(std::ostream& os) const = 0;
    
    /*
     * Nom : clone
     * Description : Effectue une copie profonde de l'arbre courant.
     * Utilisation : ExprPtr copie = noeud->clone();
     */
    virtual ExprPtr clone() const = 0;

    /*
     * Nom : integrer
     * Description : Calcule la primitive symbolique de l'expression.
     * Utilisation : ExprPtr p = noeud->integrer();
     */
    virtual ExprPtr integrer() const = 0;
    
    /*
     * Nom : limite
     * Description : Calcule la limite symbolique quand x tend vers a.
     * Utilisation : ExprPtr l = noeud->limite(a);
     */
    virtual ExprPtr limite(double a) const = 0;

    
    /*
     * Nom : estEgal
     * Description : Compare la structure et le contenu mathématique du noeud avec un autre.
     * Utilisation : bool egal = noeud->estEgal(autre_noeud);
     */
    virtual bool estEgal(const ASTNode& autre) const = 0;
    
    /*
     * Nom : estConstante
     * Description : Indique si le noeud actuel représente une valeur constante mathématique.
     * Utilisation : bool cst = noeud->estConstante();
     */
    virtual bool estConstante() const { return false; }
    
    /*
     * Nom : getValeurConstante
     * Description : Renvoie la valeur numérique du noeud s'il s'agit d'une constante.
     * Utilisation : double val = noeud->getValeurConstante();
     */
    virtual double getValeurConstante() const { return 0.0; }
};

// --- Noeuds Terminaux ---

class Constante : public ASTNode {
    double m_valeur;
public:
    /*
     * Nom : Constante
     * Description : Constructeur initialisant la constante avec sa valeur numérique.
     * Utilisation : Constante c(5.0);
     */
    explicit Constante(double valeur);
    
    /*
     * Nom : eval
     * Description : Renvoie toujours la valeur de la constante, indépendamment de x.
     * Utilisation : double val = c.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : La dérivée d'une constante est toujours zéro.
     * Utilisation : ExprPtr d = c.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie la constante (renvoie une copie d'elle-même).
     * Utilisation : ExprPtr simp = c.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche la valeur de la constante.
     * Utilisation : c.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie exacte de cette constante.
     * Utilisation : ExprPtr copie = c.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est une constante de même valeur.
     * Utilisation : bool eq = c.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

    
    /*
     * Nom : estConstante
     * Description : Renvoie systématiquement vrai car ce noeud est une constante.
     * Utilisation : bool cst = c.estConstante();
     */
    bool estConstante() const override { return true; }
    
    /*
     * Nom : getValeurConstante
     * Description : Renvoie la valeur numérique exacte de cette constante.
     * Utilisation : double val = c.getValeurConstante();
     */
    double getValeurConstante() const override { return m_valeur; }
};

class Variable : public ASTNode {
    std::string m_nom;
public:
    /*
     * Nom : Variable
     * Description : Constructeur d'une variable mathématique avec un nom (par défaut "x").
     * Utilisation : Variable v("y");
     */
    explicit Variable(const std::string& nom = "x");
    
    /*
     * Nom : eval
     * Description : Renvoie la valeur x passée en paramètre (valeur de la variable).
     * Utilisation : double val = v.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : La dérivée de x par rapport à x est 1.
     * Utilisation : ExprPtr d = v.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Renvoie la variable elle-même (pas de simplification possible).
     * Utilisation : ExprPtr simp = v.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche le nom textuel de la variable.
     * Utilisation : v.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie identique de la variable.
     * Utilisation : ExprPtr copie = v.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si l'autre noeud est une variable avec le même nom.
     * Utilisation : bool eq = v.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

// --- Opérateurs Binaires ---

class OperateurBinaire : public ASTNode {
public:
    ExprPtr m_gauche;
    ExprPtr m_droite;
    
    /*
     * Nom : OperateurBinaire
     * Description : Constructeur de base pour tous les opérateurs prenant deux opérandes.
     * Utilisation : Appelé par les constructeurs des classes filles.
     */
    OperateurBinaire(ExprPtr gauche, ExprPtr droite);
};

class Addition : public OperateurBinaire {
public:
    /*
     * Nom : Addition
     * Description : Construit un noeud d'addition de deux expressions.
     * Utilisation : Addition add(expr1, expr2);
     */
    Addition(ExprPtr gauche, ExprPtr droite);
    
    /*
     * Nom : eval
     * Description : Evalue la somme des deux sous-arbres pour une valeur de x.
     * Utilisation : double val = add.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée d'une addition (u + v)' = u' + v'.
     * Utilisation : ExprPtr d = add.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie les zéros et regroupe les termes similaires de l'addition.
     * Utilisation : ExprPtr simp = add.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche l'addition au format (gauche + droite).
     * Utilisation : add.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie de l'addition et de ses deux enfants.
     * Utilisation : ExprPtr copie = add.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie l'égalité commutative (a+b = a+b ou a+b = b+a).
     * Utilisation : bool eq = add.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class Soustraction : public OperateurBinaire {
public:
    /*
     * Nom : Soustraction
     * Description : Construit un noeud de soustraction.
     * Utilisation : Soustraction sub(expr1, expr2);
     */
    Soustraction(ExprPtr gauche, ExprPtr droite);
    
    /*
     * Nom : eval
     * Description : Evalue la différence gauche - droite.
     * Utilisation : double val = sub.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée d'une soustraction (u - v)' = u' - v'.
     * Utilisation : ExprPtr d = sub.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie les zéros, constantes et annule (a - a).
     * Utilisation : ExprPtr simp = sub.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche la soustraction au format (gauche - droite).
     * Utilisation : sub.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie de la soustraction et de ses sous-arbres.
     * Utilisation : ExprPtr copie = sub.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie que le noeud actuel est une soustraction de termes identiques.
     * Utilisation : bool eq = sub.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class Multiplication : public OperateurBinaire {
public:
    /*
     * Nom : Multiplication
     * Description : Construit un noeud de multiplication de deux expressions.
     * Utilisation : Multiplication mul(expr1, expr2);
     */
    Multiplication(ExprPtr gauche, ExprPtr droite);
    
    /*
     * Nom : eval
     * Description : Evalue le produit des deux sous-arbres pour x.
     * Utilisation : double val = mul.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée d'un produit (uv)' = u'v + uv'.
     * Utilisation : ExprPtr d = mul.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie les multiplications par 0 ou 1 et associe les constantes.
     * Utilisation : ExprPtr simp = mul.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche le produit au format gauche * droite sans parenthèses.
     * Utilisation : mul.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie profonde de la multiplication.
     * Utilisation : ExprPtr copie = mul.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie l'égalité commutative pour la multiplication.
     * Utilisation : bool eq = mul.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class Division : public OperateurBinaire {
public:
    /*
     * Nom : Division
     * Description : Construit un noeud de division de deux expressions.
     * Utilisation : Division div(expr1, expr2);
     */
    Division(ExprPtr gauche, ExprPtr droite);
    
    /*
     * Nom : eval
     * Description : Evalue le quotient (gauche / droite).
     * Utilisation : double val = div.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée d'un quotient (u/v)' = (u'v - uv') / v^2.
     * Utilisation : ExprPtr d = div.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie les divisions par 1 et les termes identiques.
     * Utilisation : ExprPtr simp = div.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche la division au format (gauche / droite).
     * Utilisation : div.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie du noeud et de ses enfants.
     * Utilisation : ExprPtr copie = div.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est la même division exacte.
     * Utilisation : bool eq = div.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class Puissance : public OperateurBinaire {
public:
    /*
     * Nom : Puissance
     * Description : Construit un noeud de puissance (base^exposant).
     * Utilisation : Puissance p(base, exposant);
     */
    Puissance(ExprPtr base, ExprPtr exposant);
    
    /*
     * Nom : eval
     * Description : Evalue (gauche) élevé à la puissance (droite).
     * Utilisation : double val = p.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : Calcule la dérivée pour un exposant constant (u^n)' = n*u^{n-1}*u'.
     * Utilisation : ExprPtr d = p.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie les puissances 0 et 1.
     * Utilisation : ExprPtr simp = p.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche la puissance au format (base)^(exposant).
     * Utilisation : p.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Renvoie une copie de la puissance entière.
     * Utilisation : ExprPtr copie = p.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est une puissance avec mêmes base et exposant.
     * Utilisation : bool eq = p.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

// --- Fonctions Unaires ---

class FonctionUnaire : public ASTNode {
public:
    ExprPtr m_argument;
    
    /*
     * Nom : FonctionUnaire
     * Description : Constructeur de base pour les fonctions mathématiques à un paramètre (ex: sin, cos).
     * Utilisation : Appelé par les constructeurs de Sinus, Cosinus, etc.
     */
    explicit FonctionUnaire(ExprPtr arg);
};

class Sinus : public FonctionUnaire {
public:
    /*
     * Nom : Sinus
     * Description : Construit un noeud pour la fonction sinus.
     * Utilisation : Sinus s(expr);
     */
    explicit Sinus(ExprPtr arg);
    
    /*
     * Nom : eval
     * Description : Evalue sin(argument).
     * Utilisation : double val = s.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : La dérivée de sin(u) est cos(u)*u'.
     * Utilisation : ExprPtr d = s.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie en évaluant la constante si possible.
     * Utilisation : ExprPtr simp = s.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche sous la forme sin(argument).
     * Utilisation : s.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie du noeud sinus.
     * Utilisation : ExprPtr copie = s.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est un sinus avec le même argument.
     * Utilisation : bool eq = s.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class Cosinus : public FonctionUnaire {
public:
    /*
     * Nom : Cosinus
     * Description : Construit un noeud pour la fonction cosinus.
     * Utilisation : Cosinus c(expr);
     */
    explicit Cosinus(ExprPtr arg);
    
    /*
     * Nom : eval
     * Description : Evalue cos(argument).
     * Utilisation : double val = c.eval(x);
     */
    double eval(double x) const override;
    
    /*
     * Nom : derivee
     * Description : La dérivée de cos(u) est -sin(u)*u'.
     * Utilisation : ExprPtr d = c.derivee();
     */
    ExprPtr derivee() const override;
    
    /*
     * Nom : simplifier
     * Description : Simplifie en évaluant la constante si possible.
     * Utilisation : ExprPtr simp = c.simplifier();
     */
    ExprPtr simplifier() const override;
    
    /*
     * Nom : afficher
     * Description : Affiche sous la forme cos(argument).
     * Utilisation : c.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;
    
    /*
     * Nom : clone
     * Description : Crée une copie du noeud cosinus.
     * Utilisation : ExprPtr copie = c.clone();
     */
    ExprPtr clone() const override;
    
    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est un cosinus avec le même argument.
     * Utilisation : bool eq = c.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};

class tangeante : public FonctionUnaire {
    public:

    /*
     * Nom : tangeante
     * Description : Construit un noeud pour la fonction tangeante.
     * Utilisation : tangeante c(expr);
     */
    explicit tangeante(ExprPtr arg);

    /*
     * Nom : eval
     * Description : Evalue tan(argument).
     * Utilisation : double val = c.eval(x);
     */
    double eval(double x) const override;

    /*
     * Nom : derivee
     * Description : La dérivée de tan(u) est u'/(cos u )**2.
     * Utilisation : ExprPtr d = c.derivee();
     */
    ExprPtr derivee() const override;

    /*
     * Nom : simplifier
     * Description : Simplifie en évaluant la constante si possible.
     * Utilisation : ExprPtr simp = c.simplifier();
     */
    ExprPtr simplifier() const override;

    /*
     * Nom : afficher
     * Description : Affiche sous la forme tan(argument).
     * Utilisation : c.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;

    /*
     * Nom : clone
     * Description : Crée une copie du noeud tangeante.
     * Utilisation : ExprPtr copie = c.clone();
     */
    ExprPtr clone() const override;

    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est un cosinus avec le même argument.
     * Utilisation : bool eq = c.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;

};



class Logarithme : public FonctionUnaire {
public:
    explicit Logarithme(ExprPtr arg);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
    ExprPtr integrer() const override;
    ExprPtr limite(double a) const override;
};

// --- Fonctions Helpers & Surcharge d'Opérateurs ---

/*
 * Nom : cst
 * Description : Helper pour créer rapidement un noeud constante.
 * Utilisation : ExprPtr c = cst(5.0);
 */
ExprPtr cst(double valeur);

/*
 * Nom : var
 * Description : Helper pour créer rapidement un noeud variable.
 * Utilisation : ExprPtr v = var("x");
 */
ExprPtr var(const std::string& nom = "x");

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre deux ExprPtr.
 * Utilisation : ExprPtr res = expr1 + expr2;
 */
ExprPtr operator+(ExprPtr gauche, ExprPtr droite);

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre ExprPtr et double.
 * Utilisation : ExprPtr res = expr + 5.0;
 */
ExprPtr operator+(ExprPtr gauche, double droite);

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre double et ExprPtr.
 * Utilisation : ExprPtr res = 5.0 + expr;
 */
ExprPtr operator+(double gauche, ExprPtr droite);

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre deux ExprPtr.
 * Utilisation : ExprPtr res = expr1 - expr2;
 */
ExprPtr operator-(ExprPtr gauche, ExprPtr droite);

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre ExprPtr et double.
 * Utilisation : ExprPtr res = expr - 5.0;
 */
ExprPtr operator-(ExprPtr gauche, double droite);

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre double et ExprPtr.
 * Utilisation : ExprPtr res = 5.0 - expr;
 */
ExprPtr operator-(double gauche, ExprPtr droite);

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre deux ExprPtr.
 * Utilisation : ExprPtr res = expr1 * expr2;
 */
ExprPtr operator*(ExprPtr gauche, ExprPtr droite);

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre ExprPtr et double.
 * Utilisation : ExprPtr res = expr * 5.0;
 */
ExprPtr operator*(ExprPtr gauche, double droite);

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre double et ExprPtr.
 * Utilisation : ExprPtr res = 5.0 * expr;
 */
ExprPtr operator*(double gauche, ExprPtr droite);

/*
 * Nom : operator/
 * Description : Surcharge de la division entre deux ExprPtr.
 * Utilisation : ExprPtr res = expr1 / expr2;
 */
ExprPtr operator/(ExprPtr gauche, ExprPtr droite);

/*
 * Nom : operator/
 * Description : Surcharge de la division entre ExprPtr et double.
 * Utilisation : ExprPtr res = expr / 5.0;
 */
ExprPtr operator/(ExprPtr gauche, double droite);

/*
 * Nom : operator/
 * Description : Surcharge de la division entre double et ExprPtr.
 * Utilisation : ExprPtr res = 5.0 / expr;
 */
ExprPtr operator/(double gauche, ExprPtr droite);

/*
 * Nom : ast_pow
 * Description : Fonction pour créer une puissance avec deux expressions.
 * Utilisation : ExprPtr p = ast_pow(base_expr, exp_expr);
 */
ExprPtr ast_pow(ExprPtr base, ExprPtr exposant);

/*
 * Nom : ast_pow
 * Description : Fonction pour créer une puissance avec un exposant constant (double).
 * Utilisation : ExprPtr p = ast_pow(base_expr, 2.0);
 */
ExprPtr ast_pow(ExprPtr base, double exposant);

/*
 * Nom : ast_pow
 * Description : Fonction pour créer une puissance avec une base constante (double).
 * Utilisation : ExprPtr p = ast_pow(2.0, exp_expr);
 */
ExprPtr ast_pow(double base, ExprPtr exposant);

/*
 * Nom : ast_sin
 * Description : Fonction pour créer le noeud sinus d'une expression.
 * Utilisation : ExprPtr s = ast_sin(expr);
 */
ExprPtr ast_sin(ExprPtr arg);

/*
 * Nom : ast_cos
 * Description : Fonction pour créer le noeud cosinus d'une expression.
 * Utilisation : ExprPtr c = ast_cos(expr);
 */
ExprPtr ast_cos(ExprPtr arg);

/*
 * Nom : ast_tan
 * Description : Fonction pour créer le noeud tangeante d'une expression.
 * Utilisation : ExprPtr c = ast_tan(expr);
 */
 ExprPtr ast_tan(ExprPtr arg);




/*
 * Nom : ast_ln
 * Description : Fonction pour créer le noeud logarithme d'une expression.
 * Utilisation : ExprPtr l = ast_ln(expr);
 */
ExprPtr ast_ln(ExprPtr arg);
