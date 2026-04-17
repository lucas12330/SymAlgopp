#pragma once
#include <memory>
#include <string>
#include <ostream>

class ASTNode;
using ExprPtr = std::shared_ptr<ASTNode>;

// Classe abstraite de base pour tous les noeuds de l'arbre
class ASTNode : public std::enable_shared_from_this<ASTNode> {
public:
    virtual ~ASTNode() = default;
    
    virtual double eval(double x) const = 0;
    virtual ExprPtr derivee() const = 0;
    virtual ExprPtr simplifier() const = 0;
    virtual void afficher(std::ostream& os) const = 0;
    virtual ExprPtr clone() const = 0;
    virtual bool estEgal(const ASTNode& autre) const = 0;
    
    virtual bool estConstante() const { return false; }
    virtual double getValeurConstante() const { return 0.0; }
};

// --- Noeuds Terminaux ---

class Constante : public ASTNode {
    double m_valeur;
public:
    explicit Constante(double valeur);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
    bool estConstante() const override { return true; }
    double getValeurConstante() const override { return m_valeur; }
};

class Variable : public ASTNode {
    std::string m_nom;
public:
    explicit Variable(const std::string& nom = "x");
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

// --- Opérateurs Binaires ---

class OperateurBinaire : public ASTNode {
public:
    ExprPtr m_gauche;
    ExprPtr m_droite;
    OperateurBinaire(ExprPtr gauche, ExprPtr droite);
};

class Addition : public OperateurBinaire {
public:
    Addition(ExprPtr gauche, ExprPtr droite);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

class Soustraction : public OperateurBinaire {
public:
    Soustraction(ExprPtr gauche, ExprPtr droite);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

class Multiplication : public OperateurBinaire {
public:
    Multiplication(ExprPtr gauche, ExprPtr droite);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

class Division : public OperateurBinaire {
public:
    Division(ExprPtr gauche, ExprPtr droite);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

class Puissance : public OperateurBinaire {
public:
    Puissance(ExprPtr base, ExprPtr exposant);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

// --- Fonctions Unaires ---

class FonctionUnaire : public ASTNode {
public:
    ExprPtr m_argument;
    explicit FonctionUnaire(ExprPtr arg);
};

class Sinus : public FonctionUnaire {
public:
    explicit Sinus(ExprPtr arg);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};

class Cosinus : public FonctionUnaire {
public:
    explicit Cosinus(ExprPtr arg);
    double eval(double x) const override;
    ExprPtr derivee() const override;
    ExprPtr simplifier() const override;
    void afficher(std::ostream& os) const override;
    ExprPtr clone() const override;
    bool estEgal(const ASTNode& autre) const override;
};


// --- Fonctions Helpers & Surcharge d'Opérateurs ---
ExprPtr cst(double valeur);
ExprPtr var(const std::string& nom = "x");

ExprPtr operator+(ExprPtr gauche, ExprPtr droite);
ExprPtr operator+(ExprPtr gauche, double droite);
ExprPtr operator+(double gauche, ExprPtr droite);

ExprPtr operator-(ExprPtr gauche, ExprPtr droite);
ExprPtr operator-(ExprPtr gauche, double droite);
ExprPtr operator-(double gauche, ExprPtr droite);

ExprPtr operator*(ExprPtr gauche, ExprPtr droite);
ExprPtr operator*(ExprPtr gauche, double droite);
ExprPtr operator*(double gauche, ExprPtr droite);

ExprPtr operator/(ExprPtr gauche, ExprPtr droite);
ExprPtr operator/(ExprPtr gauche, double droite);
ExprPtr operator/(double gauche, ExprPtr droite);

ExprPtr ast_pow(ExprPtr base, ExprPtr exposant);
ExprPtr ast_pow(ExprPtr base, double exposant);
ExprPtr ast_pow(double base, ExprPtr exposant);

ExprPtr ast_sin(ExprPtr arg);
ExprPtr ast_cos(ExprPtr arg);
