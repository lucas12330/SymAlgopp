#include "ASTNode.hpp"
#include <cmath>
#include <iostream>

/*
 * Nom : extractCoeff
 * Description : Helper interne pour la factorisation de type a*U + b*U. Modifie coeff et u passés en référence.
 * Utilisation : double c; ExprPtr u; extractCoeff(expr, c, u);
 */
void extractCoeff(ExprPtr e, double& coeff, ExprPtr& u) {
    auto mult = std::dynamic_pointer_cast<Multiplication>(e);
    if (mult) {
        if (mult->m_gauche->estConstante()) {
            coeff = mult->m_gauche->getValeurConstante();
            u = mult->m_droite;
            return;
        }
        if (mult->m_droite->estConstante()) {
            coeff = mult->m_droite->getValeurConstante();
            u = mult->m_gauche;
            return;
        }
    }
    coeff = 1.0;
    u = e;
}

// ============== CONSTANTE ==================

/*
 * Nom : Constante
 * Description : Constructeur de constante qui stocke la valeur passée.
 * Utilisation : Constante c(3.14);
 */
Constante::Constante(double valeur) : m_valeur(valeur) {}

/*
 * Nom : eval
 * Description : Renvoie la valeur de la constante, ignore x.
 * Utilisation : double val = c.eval(x);
 */
double Constante::eval(double) const { return m_valeur; }

/*
 * Nom : derivee
 * Description : La dérivée d'une constante vaut 0. Renvoie un noeud Constante 0.
 * Utilisation : ExprPtr d = c.derivee();
 */
ExprPtr Constante::derivee() const { return cst(0.0); }

/*
 * Nom : simplifier
 * Description : Ne simplifie rien, renvoie une copie d'elle-même.
 * Utilisation : ExprPtr simp = c.simplifier();
 */
ExprPtr Constante::simplifier() const { return clone(); }

/*
 * Nom : afficher
 * Description : Ecrit la valeur numérique sur le flux.
 * Utilisation : c.afficher(std::cout);
 */
void Constante::afficher(std::ostream& os) const { os << m_valeur; }

/*
 * Nom : clone
 * Description : Renvoie un nouveau noeud constante avec la même valeur.
 * Utilisation : ExprPtr copie = c.clone();
 */
ExprPtr Constante::clone() const { return cst(m_valeur); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une constante et possède la même valeur aux erreurs flottantes près.
 * Utilisation : bool eq = c.estEgal(autre);
 */
bool Constante::estEgal(const ASTNode& autre) const {
    return autre.estConstante() && std::abs(autre.getValeurConstante() - m_valeur) < 1e-9;
}

// ============== VARIABLE ==================

/*
 * Nom : Variable
 * Description : Constructeur de variable qui l'initialise avec le nom donné.
 * Utilisation : Variable v("t");
 */
Variable::Variable(const std::string& nom) : m_nom(nom) {}

/*
 * Nom : eval
 * Description : Renvoie la valeur x, assumant que la variable est l'inconnue d'évaluation.
 * Utilisation : double val = v.eval(x);
 */
double Variable::eval(double x) const { return x; }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée de x qui est 1.
 * Utilisation : ExprPtr d = v.derivee();
 */
ExprPtr Variable::derivee() const { return cst(1.0); }

/*
 * Nom : simplifier
 * Description : Ne simplifie rien, renvoie une copie.
 * Utilisation : ExprPtr simp = v.simplifier();
 */
ExprPtr Variable::simplifier() const { return clone(); }

/*
 * Nom : afficher
 * Description : Imprime le nom de la variable.
 * Utilisation : v.afficher(std::cout);
 */
void Variable::afficher(std::ostream& os) const { os << m_nom; }

/*
 * Nom : clone
 * Description : Renvoie une copie de la variable avec le même nom.
 * Utilisation : ExprPtr copie = v.clone();
 */
ExprPtr Variable::clone() const { return var(m_nom); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une variable portant le même nom.
 * Utilisation : bool eq = v.estEgal(autre);
 */
bool Variable::estEgal(const ASTNode& autre) const {
    const Variable* v = dynamic_cast<const Variable*>(&autre);
    return v != nullptr && v->m_nom == m_nom;
}

// ============== OP Binaire Base ==================

/*
 * Nom : OperateurBinaire
 * Description : Initialise les sous-arbres gauche et droite d'un opérateur.
 * Utilisation : Appelé par les classes dérivées.
 */
OperateurBinaire::OperateurBinaire(ExprPtr gauche, ExprPtr droite) : m_gauche(gauche), m_droite(droite) {}

// ============== ADDITION ==================

/*
 * Nom : Addition
 * Description : Constructeur de l'addition.
 * Utilisation : Addition add(gauche, droite);
 */
Addition::Addition(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}

/*
 * Nom : eval
 * Description : Evalue chaque opérande puis les additionne.
 * Utilisation : double val = add.eval(x);
 */
double Addition::eval(double x) const { return m_gauche->eval(x) + m_droite->eval(x); }

/*
 * Nom : derivee
 * Description : Renvoie un noeud Addition des dérivées de chaque opérande.
 * Utilisation : ExprPtr d = add.derivee();
 */
ExprPtr Addition::derivee() const { return m_gauche->derivee() + m_droite->derivee(); }

/*
 * Nom : simplifier
 * Description : Simplifie les termes constants, supprime les zéros inutiles, et factorise a*U + b*U en (a+b)*U.
 * Utilisation : ExprPtr simp = add.simplifier();
 */
ExprPtr Addition::simplifier() const {
    auto g = m_gauche->simplifier();
    auto d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante()) return cst(g->getValeurConstante() + d->getValeurConstante());
    if (g->estConstante() && g->getValeurConstante() == 0.0) return d;
    if (d->estConstante() && d->getValeurConstante() == 0.0) return g;
    
    // Factorisation a*U + b*U = (a+b)*U
    double cG, cD;
    ExprPtr uG, uD;
    extractCoeff(g, cG, uG);
    extractCoeff(d, cD, uD);
    if (uG->estEgal(*uD)) {
        return cst(cG + cD) * uG;
    }
    
    return std::make_shared<Addition>(g, d);
}

/*
 * Nom : afficher
 * Description : Affiche (gauche + droite).
 * Utilisation : add.afficher(std::cout);
 */
void Addition::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " + "; m_droite->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie de l'addition et de ses enfants.
 * Utilisation : ExprPtr copie = add.clone();
 */
ExprPtr Addition::clone() const { return std::make_shared<Addition>(m_gauche->clone(), m_droite->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même addition (tient compte de la commutativité).
 * Utilisation : bool eq = add.estEgal(autre);
 */
bool Addition::estEgal(const ASTNode& autre) const {
    const Addition* a = dynamic_cast<const Addition*>(&autre);
    if (!a) return false;
    return (m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite))) ||
           (m_gauche->estEgal(*(a->m_droite)) && m_droite->estEgal(*(a->m_gauche)));
}

// ============== SOUSTRACTION ==================

/*
 * Nom : Soustraction
 * Description : Constructeur de la soustraction.
 * Utilisation : Soustraction sub(gauche, droite);
 */
Soustraction::Soustraction(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}

/*
 * Nom : eval
 * Description : Evalue chaque opérande puis effectue gauche - droite.
 * Utilisation : double val = sub.eval(x);
 */
double Soustraction::eval(double x) const { return m_gauche->eval(x) - m_droite->eval(x); }

/*
 * Nom : derivee
 * Description : Renvoie un noeud Soustraction des dérivées de chaque opérande.
 * Utilisation : ExprPtr d = sub.derivee();
 */
ExprPtr Soustraction::derivee() const { return m_gauche->derivee() - m_droite->derivee(); }

/*
 * Nom : simplifier
 * Description : Simplifie les termes constants, supprime -0 et factorise a*U - b*U en (a-b)*U.
 * Utilisation : ExprPtr simp = sub.simplifier();
 */
ExprPtr Soustraction::simplifier() const {
    auto g = m_gauche->simplifier();
    auto d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante()) return cst(g->getValeurConstante() - d->getValeurConstante());
    if (d->estConstante() && d->getValeurConstante() == 0.0) return g;
    
    double cG, cD;
    ExprPtr uG, uD;
    extractCoeff(g, cG, uG);
    extractCoeff(d, cD, uD);
    if (uG->estEgal(*uD)) {
        if (std::abs(cG - cD) < 1e-9) return cst(0.0);
        return cst(cG - cD) * uG;
    }
    
    return std::make_shared<Soustraction>(g, d);
}

/*
 * Nom : afficher
 * Description : Affiche (gauche - droite).
 * Utilisation : sub.afficher(std::cout);
 */
void Soustraction::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " - "; m_droite->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie de la soustraction et de ses enfants.
 * Utilisation : ExprPtr copie = sub.clone();
 */
ExprPtr Soustraction::clone() const { return std::make_shared<Soustraction>(m_gauche->clone(), m_droite->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même soustraction exacte.
 * Utilisation : bool eq = sub.estEgal(autre);
 */
bool Soustraction::estEgal(const ASTNode& autre) const {
    const Soustraction* a = dynamic_cast<const Soustraction*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== MULTIPLICATION ==================

/*
 * Nom : Multiplication
 * Description : Constructeur de la multiplication.
 * Utilisation : Multiplication mul(gauche, droite);
 */
Multiplication::Multiplication(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}

/*
 * Nom : eval
 * Description : Evalue chaque opérande puis effectue gauche * droite.
 * Utilisation : double val = mul.eval(x);
 */
double Multiplication::eval(double x) const { return m_gauche->eval(x) * m_droite->eval(x); }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée du produit : u'v + uv'.
 * Utilisation : ExprPtr d = mul.derivee();
 */
ExprPtr Multiplication::derivee() const {
    return (m_gauche->derivee() * m_droite) + (m_gauche * m_droite->derivee());
}

/*
 * Nom : simplifier
 * Description : Evalue les constantes, simplifie les multiplications par 0 ou 1. Place les constantes à gauche.
 * Utilisation : ExprPtr simp = mul.simplifier();
 */
ExprPtr Multiplication::simplifier() const {
    auto g = m_gauche->simplifier();
    auto d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante()) return cst(g->getValeurConstante() * d->getValeurConstante());
    if (g->estConstante()) {
        if (std::abs(g->getValeurConstante()) < 1e-9) return cst(0.0);
        if (std::abs(g->getValeurConstante() - 1.0) < 1e-9) return d;
    }
    if (d->estConstante()) {
        if (std::abs(d->getValeurConstante()) < 1e-9) return cst(0.0);
        if (std::abs(d->getValeurConstante() - 1.0) < 1e-9) return g;
        // Place toujours la constante à gauche
        return std::make_shared<Multiplication>(d, g); 
    }
    return std::make_shared<Multiplication>(g, d);
}

/*
 * Nom : afficher
 * Description : Affiche gauche * droite.
 * Utilisation : mul.afficher(std::cout);
 */
void Multiplication::afficher(std::ostream& os) const {
    m_gauche->afficher(os); os << " * "; m_droite->afficher(os);
}

/*
 * Nom : clone
 * Description : Renvoie une copie de la multiplication et de ses enfants.
 * Utilisation : ExprPtr copie = mul.clone();
 */
ExprPtr Multiplication::clone() const { return std::make_shared<Multiplication>(m_gauche->clone(), m_droite->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une multiplication identique (commutativité acceptée).
 * Utilisation : bool eq = mul.estEgal(autre);
 */
bool Multiplication::estEgal(const ASTNode& autre) const {
    const Multiplication* a = dynamic_cast<const Multiplication*>(&autre);
    if (!a) return false;
    return (m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite))) ||
           (m_gauche->estEgal(*(a->m_droite)) && m_droite->estEgal(*(a->m_gauche)));
}

// ============== DIVISION ==================

/*
 * Nom : Division
 * Description : Constructeur de la division.
 * Utilisation : Division div(gauche, droite);
 */
Division::Division(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}

/*
 * Nom : eval
 * Description : Evalue chaque opérande puis effectue gauche / droite.
 * Utilisation : double val = div.eval(x);
 */
double Division::eval(double x) const { return m_gauche->eval(x) / m_droite->eval(x); }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée du quotient : (u'v - uv') / v^2.
 * Utilisation : ExprPtr d = div.derivee();
 */
ExprPtr Division::derivee() const { // (u'v - uv') / v^2
    auto num = (m_gauche->derivee() * m_droite) - (m_gauche * m_droite->derivee());
    auto den = m_droite * m_droite;
    return num / den;
}

/*
 * Nom : simplifier
 * Description : Evalue les constantes, simplifie 0/u, u/1, et u/u.
 * Utilisation : ExprPtr simp = div.simplifier();
 */
ExprPtr Division::simplifier() const {
    auto g = m_gauche->simplifier();
    auto d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante() && d->getValeurConstante() != 0.0) 
        return cst(g->getValeurConstante() / d->getValeurConstante());
    if (g->estConstante() && g->getValeurConstante() == 0.0) return cst(0.0);
    if (d->estConstante() && d->getValeurConstante() == 1.0) return g;
    if (g->estEgal(*d)) return cst(1.0);
    return std::make_shared<Division>(g, d);
}

/*
 * Nom : afficher
 * Description : Affiche (gauche / droite).
 * Utilisation : div.afficher(std::cout);
 */
void Division::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " / "; m_droite->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie de la division et de ses enfants.
 * Utilisation : ExprPtr copie = div.clone();
 */
ExprPtr Division::clone() const { return std::make_shared<Division>(m_gauche->clone(), m_droite->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même division exacte.
 * Utilisation : bool eq = div.estEgal(autre);
 */
bool Division::estEgal(const ASTNode& autre) const {
    const Division* a = dynamic_cast<const Division*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== PUISSANCE ==================

/*
 * Nom : Puissance
 * Description : Constructeur de la puissance (base^exposant).
 * Utilisation : Puissance p(base, exposant);
 */
Puissance::Puissance(ExprPtr base, ExprPtr exposant) : OperateurBinaire(base, exposant) {}

/*
 * Nom : eval
 * Description : Evalue la base élevée à l'exposant avec pow().
 * Utilisation : double val = p.eval(x);
 */
double Puissance::eval(double x) const { return std::pow(m_gauche->eval(x), m_droite->eval(x)); }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée pour un exposant constant : n*u^{n-1}*u'. Renvoie 0 autrement.
 * Utilisation : ExprPtr d = p.derivee();
 */
ExprPtr Puissance::derivee() const {
    if (m_droite->estConstante()) { // (u^n)' = n*u^{n-1}*u'
        double n = m_droite->getValeurConstante();
        if (n == 0) return cst(0.0);
        return cst(n) * ast_pow(m_gauche, cst(n - 1)) * m_gauche->derivee();
    }
    // Cas general : (u^v)' = u^v * (v' * ln(u) + v * u' / u)
    auto ln_u = ast_ln(m_gauche);
    auto terme1 = m_droite->derivee() * ln_u;
    auto terme2 = m_droite * (m_gauche->derivee() / m_gauche);
    return ast_pow(m_gauche, m_droite) * (terme1 + terme2);
}

/*
 * Nom : simplifier
 * Description : Evalue les constantes, simplifie u^0 = 1, u^1 = u, 0^p = 0, 1^p = 1.
 * Utilisation : ExprPtr simp = p.simplifier();
 */
ExprPtr Puissance::simplifier() const {
    auto b = m_gauche->simplifier();
    auto p = m_droite->simplifier();
    if (b->estConstante() && p->estConstante()) return cst(std::pow(b->getValeurConstante(), p->getValeurConstante()));
    if (p->estConstante() && p->getValeurConstante() == 0.0) return cst(1.0);
    if (p->estConstante() && p->getValeurConstante() == 1.0) return b;
    if (b->estConstante() && b->getValeurConstante() == 0.0) return cst(0.0);
    if (b->estConstante() && b->getValeurConstante() == 1.0) return cst(1.0);
    return ast_pow(b, p);
}

/*
 * Nom : afficher
 * Description : Affiche (base)^(exposant).
 * Utilisation : p.afficher(std::cout);
 */
void Puissance::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << ")^("; m_droite->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie de la puissance et de ses enfants.
 * Utilisation : ExprPtr copie = p.clone();
 */
ExprPtr Puissance::clone() const { return std::make_shared<Puissance>(m_gauche->clone(), m_droite->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une puissance avec la même base et le même exposant.
 * Utilisation : bool eq = p.estEgal(autre);
 */
bool Puissance::estEgal(const ASTNode& autre) const {
    const Puissance* a = dynamic_cast<const Puissance*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== FONCTION UNAIRE BASE =================

/*
 * Nom : FonctionUnaire
 * Description : Constructeur de base des fonctions mathématiques unaires initialisant l'argument.
 * Utilisation : Appelé par les constructeurs des classes filles.
 */
FonctionUnaire::FonctionUnaire(ExprPtr arg) : m_argument(arg) {}

// ============== SINUS ==================

/*
 * Nom : Sinus
 * Description : Constructeur de la fonction sinus avec l'argument.
 * Utilisation : Sinus s(expr);
 */
Sinus::Sinus(ExprPtr arg) : FonctionUnaire(arg) {}

/*
 * Nom : eval
 * Description : Evalue sin(argument).
 * Utilisation : double val = s.eval(x);
 */
double Sinus::eval(double x) const { return std::sin(m_argument->eval(x)); }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée de sin(u) : cos(u)*u'.
 * Utilisation : ExprPtr d = s.derivee();
 */
ExprPtr Sinus::derivee() const { return ast_cos(m_argument) * m_argument->derivee(); }

/*
 * Nom : simplifier
 * Description : Evalue la constante si possible.
 * Utilisation : ExprPtr simp = s.simplifier();
 */
ExprPtr Sinus::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::sin(arg->getValeurConstante()));
    return ast_sin(arg);
}

/*
 * Nom : afficher
 * Description : Affiche sin(argument).
 * Utilisation : s.afficher(std::cout);
 */
void Sinus::afficher(std::ostream& os) const {
    os << "sin("; m_argument->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie du sinus et de son argument.
 * Utilisation : ExprPtr copie = s.clone();
 */
ExprPtr Sinus::clone() const { return std::make_shared<Sinus>(m_argument->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est un sinus avec le même argument exact.
 * Utilisation : bool eq = s.estEgal(autre);
 */
bool Sinus::estEgal(const ASTNode& autre) const {
    const Sinus* a = dynamic_cast<const Sinus*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

// ============== COSINUS ==================

/*
 * Nom : Cosinus
 * Description : Constructeur de la fonction cosinus avec l'argument.
 * Utilisation : Cosinus c(expr);
 */
Cosinus::Cosinus(ExprPtr arg) : FonctionUnaire(arg) {}

/*
 * Nom : eval
 * Description : Evalue cos(argument).
 * Utilisation : double val = c.eval(x);
 */
double Cosinus::eval(double x) const { return std::cos(m_argument->eval(x)); }

/*
 * Nom : derivee
 * Description : Renvoie la dérivée de cos(u) : -sin(u)*u'.
 * Utilisation : ExprPtr d = c.derivee();
 */
ExprPtr Cosinus::derivee() const { return (cst(-1.0) * ast_sin(m_argument)) * m_argument->derivee(); }

/*
 * Nom : simplifier
 * Description : Evalue la constante si possible.
 * Utilisation : ExprPtr simp = c.simplifier();
 */
ExprPtr Cosinus::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::cos(arg->getValeurConstante()));
    return ast_cos(arg);
}

/*
 * Nom : afficher
 * Description : Affiche cos(argument).
 * Utilisation : c.afficher(std::cout);
 */
void Cosinus::afficher(std::ostream& os) const {
    os << "cos("; m_argument->afficher(os); os << ")";
}

/*
 * Nom : clone
 * Description : Renvoie une copie du cosinus et de son argument.
 * Utilisation : ExprPtr copie = c.clone();
 */
ExprPtr Cosinus::clone() const { return std::make_shared<Cosinus>(m_argument->clone()); }

/*
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est un cosinus avec le même argument exact.
 * Utilisation : bool eq = c.estEgal(autre);
 */
bool Cosinus::estEgal(const ASTNode& autre) const {
    const Cosinus* a = dynamic_cast<const Cosinus*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}


// --- Surcharge d'Opérateurs =---

/*
 * Nom : cst
 * Description : Helper générant un noeud Constante à partir d'un double.
 * Utilisation : ExprPtr noeud = cst(3.14);
 */
ExprPtr cst(double valeur) { return std::make_shared<Constante>(valeur); }

/*
 * Nom : var
 * Description : Helper générant un noeud Variable à partir d'un nom de variable.
 * Utilisation : ExprPtr noeud = var("y");
 */
ExprPtr var(const std::string& nom) { return std::make_shared<Variable>(nom); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre deux ExprPtr générant un noeud Addition.
 * Utilisation : ExprPtr resultat = e1 + e2;
 */
ExprPtr operator+(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Addition>(gauche, droite); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre ExprPtr et double générant un noeud Addition.
 * Utilisation : ExprPtr resultat = e + 2.0;
 */
ExprPtr operator+(ExprPtr gauche, double droite) { return std::make_shared<Addition>(gauche, cst(droite)); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre double et ExprPtr générant un noeud Addition.
 * Utilisation : ExprPtr resultat = 2.0 + e;
 */
ExprPtr operator+(double gauche, ExprPtr droite) { return std::make_shared<Addition>(cst(gauche), droite); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre deux ExprPtr générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = e1 - e2;
 */
ExprPtr operator-(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Soustraction>(gauche, droite); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre ExprPtr et double générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = e - 2.0;
 */
ExprPtr operator-(ExprPtr gauche, double droite) { return std::make_shared<Soustraction>(gauche, cst(droite)); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre double et ExprPtr générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = 2.0 - e;
 */
ExprPtr operator-(double gauche, ExprPtr droite) { return std::make_shared<Soustraction>(cst(gauche), droite); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre deux ExprPtr générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = e1 * e2;
 */
ExprPtr operator*(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Multiplication>(gauche, droite); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre ExprPtr et double générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = e * 2.0;
 */
ExprPtr operator*(ExprPtr gauche, double droite) { return std::make_shared<Multiplication>(gauche, cst(droite)); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre double et ExprPtr générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = 2.0 * e;
 */
ExprPtr operator*(double gauche, ExprPtr droite) { return std::make_shared<Multiplication>(cst(gauche), droite); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre deux ExprPtr générant un noeud Division.
 * Utilisation : ExprPtr resultat = e1 / e2;
 */
ExprPtr operator/(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Division>(gauche, droite); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre ExprPtr et double générant un noeud Division.
 * Utilisation : ExprPtr resultat = e / 2.0;
 */
ExprPtr operator/(ExprPtr gauche, double droite) { return std::make_shared<Division>(gauche, cst(droite)); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre double et ExprPtr générant un noeud Division.
 * Utilisation : ExprPtr resultat = 2.0 / e;
 */
ExprPtr operator/(double gauche, ExprPtr droite) { return std::make_shared<Division>(cst(gauche), droite); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance entre deux ExprPtr.
 * Utilisation : ExprPtr resultat = ast_pow(base, exposant);
 */
ExprPtr ast_pow(ExprPtr base, ExprPtr exposant) { return std::make_shared<Puissance>(base, exposant); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance avec un exposant double constant.
 * Utilisation : ExprPtr resultat = ast_pow(base, 2.0);
 */
ExprPtr ast_pow(ExprPtr base, double exposant) { return std::make_shared<Puissance>(base, cst(exposant)); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance avec une base double constante.
 * Utilisation : ExprPtr resultat = ast_pow(2.0, exposant);
 */
ExprPtr ast_pow(double base, ExprPtr exposant) { return std::make_shared<Puissance>(cst(base), exposant); }

/*
 * Nom : ast_sin
 * Description : Helper générant un noeud Sinus.
 * Utilisation : ExprPtr resultat = ast_sin(e);
 */
ExprPtr ast_sin(ExprPtr arg) { return std::make_shared<Sinus>(arg); }

/*
 * Nom : ast_cos
 * Description : Helper générant un noeud Cosinus.
 * Utilisation : ExprPtr resultat = ast_cos(e);
 */
ExprPtr ast_cos(ExprPtr arg) { return std::make_shared<Cosinus>(arg); }

// ============================================================================
// IMPLÉMENTATION DES LIMITES, INTÉGRALES, DL ET LOGARITHME
// ============================================================================

#include <cmath>
#include <limits>

// --- ASTNode : Développement Limité ---
ExprPtr ASTNode::DL(double a, int ordre) const {
    ExprPtr result = cst(this->eval(a));
    ExprPtr deriv = this->clone();
    long long factorielle = 1;
    for (int i = 1; i <= ordre; ++i) {
        deriv = deriv->derivee();
        factorielle *= i;
        double coeff = deriv->eval(a) / factorielle;
        if (std::abs(coeff) > 1e-9) {
            ExprPtr terme = cst(coeff) * ast_pow(var("x") - a, i);
            result = result + terme;
        }
    }
    return result->simplifier();
}

// --- Constante ---
ExprPtr Constante::integrer() const {
    return cst(m_valeur) * var("x");
}
ExprPtr Constante::limite(double /*a*/) const {
    return cst(m_valeur);
}

// --- Variable ---
ExprPtr Variable::integrer() const {
    return cst(0.5) * ast_pow(var(m_nom), 2);
}
ExprPtr Variable::limite(double a) const {
    return cst(a);
}

// --- Addition ---
ExprPtr Addition::integrer() const {
    return m_gauche->integrer() + m_droite->integrer();
}
ExprPtr Addition::limite(double a) const {
    return m_gauche->limite(a) + m_droite->limite(a);
}

// --- Soustraction ---
ExprPtr Soustraction::integrer() const {
    return m_gauche->integrer() - m_droite->integrer();
}
ExprPtr Soustraction::limite(double a) const {
    return m_gauche->limite(a) - m_droite->limite(a);
}

// --- Multiplication ---
ExprPtr Multiplication::integrer() const {
    if (m_gauche->estConstante()) return m_gauche * m_droite->integrer();
    if (m_droite->estConstante()) return m_droite * m_gauche->integrer();
    // Intégration par parties non-gérée
    return cst(0.0);
}
ExprPtr Multiplication::limite(double a) const {
    return m_gauche->limite(a) * m_droite->limite(a);
}

// --- Division ---
ExprPtr Division::integrer() const {
    if (m_droite->estConstante()) return m_gauche->integrer() / m_droite;
    return cst(0.0);
}
ExprPtr Division::limite(double a) const {
    double n = m_gauche->eval(a);
    double d = m_droite->eval(a);
    if (std::abs(d) < 1e-9) {
        if (std::abs(n) < 1e-9) {
            // L'Hôpital : lim N/D = lim N'/D'
            return (m_gauche->derivee() / m_droite->derivee())->limite(a);
        } else {
            return cst(std::numeric_limits<double>::infinity());
        }
    }
    return m_gauche->limite(a) / m_droite->limite(a);
}

// --- Puissance ---
ExprPtr Puissance::integrer() const {
    if (m_droite->estConstante() && dynamic_cast<Variable*>(m_gauche.get())) {
        double n = m_droite->getValeurConstante();
        if (std::abs(n + 1.0) < 1e-9) {
            return ast_ln(m_gauche);
        }
        return ast_pow(m_gauche, cst(n + 1.0)) / cst(n + 1.0);
    }
    return cst(0.0);
}
ExprPtr Puissance::limite(double a) const {
    return ast_pow(m_gauche->limite(a), m_droite->limite(a));
}

// --- Sinus ---
ExprPtr Sinus::integrer() const {
    if (dynamic_cast<Variable*>(m_argument.get())) return cst(-1.0) * ast_cos(m_argument);
    return cst(0.0);
}
ExprPtr Sinus::limite(double a) const {
    return ast_sin(m_argument->limite(a));
}

// --- Cosinus ---
ExprPtr Cosinus::integrer() const {
    if (dynamic_cast<Variable*>(m_argument.get())) return ast_sin(m_argument);
    return cst(0.0);
}
ExprPtr Cosinus::limite(double a) const {
    return ast_cos(m_argument->limite(a));
}

// ============== EXPONENTIELLE ==================

Exponentielle::Exponentielle(ExprPtr arg) : FonctionUnaire(arg) {}

double Exponentielle::eval(double x) const { return std::exp(m_argument->eval(x)); }

ExprPtr Exponentielle::derivee() const {
    return ast_exp(m_argument) * m_argument->derivee();
}

ExprPtr Exponentielle::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::exp(arg->getValeurConstante()));
    if (auto ln_node = std::dynamic_pointer_cast<Logarithme>(arg)) {
        return ln_node->m_argument->simplifier();
    }
    return ast_exp(arg);
}

void Exponentielle::afficher(std::ostream& os) const {
    os << "exp("; m_argument->afficher(os); os << ")";
}

ExprPtr Exponentielle::clone() const { return ast_exp(m_argument->clone()); }

bool Exponentielle::estEgal(const ASTNode& autre) const {
    const Exponentielle* a = dynamic_cast<const Exponentielle*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

ExprPtr Exponentielle::integrer() const {
    if (dynamic_cast<Variable*>(m_argument.get())) return ast_exp(m_argument);
    double cG; ExprPtr uG;
    extractCoeff(m_argument, cG, uG);
    if (dynamic_cast<Variable*>(uG.get()) && cG != 0.0) {
        return cst(1.0 / cG) * ast_exp(m_argument);
    }
    return cst(0.0);
}

ExprPtr Exponentielle::limite(double a) const {
    return ast_exp(m_argument->limite(a));
}

ExprPtr ast_exp(ExprPtr arg) {
    return std::make_shared<Exponentielle>(arg);
}

// ============== LOGARITHME ==================

Logarithme::Logarithme(ExprPtr arg) : FonctionUnaire(arg) {}

double Logarithme::eval(double x) const { return std::log(m_argument->eval(x)); }

ExprPtr Logarithme::derivee() const {
    return m_argument->derivee() / m_argument;
}

ExprPtr Logarithme::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::log(arg->getValeurConstante()));
    return ast_ln(arg);
}

void Logarithme::afficher(std::ostream& os) const {
    os << "ln("; m_argument->afficher(os); os << ")";
}

ExprPtr Logarithme::clone() const { return ast_ln(m_argument->clone()); }

bool Logarithme::estEgal(const ASTNode& autre) const {
    const Logarithme* a = dynamic_cast<const Logarithme*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

ExprPtr Logarithme::integrer() const {
    if (dynamic_cast<Variable*>(m_argument.get())) {
        return m_argument * ast_ln(m_argument) - m_argument;
    }
    return cst(0.0);
}

ExprPtr Logarithme::limite(double a) const {
    double v = m_argument->eval(a);
    if (v <= 0.0) return cst(-std::numeric_limits<double>::infinity());
    return ast_ln(m_argument->limite(a));
}

ExprPtr ast_ln(ExprPtr arg) {
    return std::make_shared<Logarithme>(arg);
}
