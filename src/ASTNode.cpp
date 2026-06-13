#include "ASTNode.hpp"
#include <cmath>
#include <iostream>

// Helper interne pour la factorisation de type a*U + b*U
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
Constante::Constante(double valeur) : m_valeur(valeur) {}
double Constante::eval(double) const { return m_valeur; }
ExprPtr Constante::derivee() const { return cst(0.0); }
ExprPtr Constante::simplifier() const { return clone(); }
void Constante::afficher(std::ostream& os) const { os << m_valeur; }
ExprPtr Constante::clone() const { return cst(m_valeur); }
bool Constante::estEgal(const ASTNode& autre) const {
    return autre.estConstante() && std::abs(autre.getValeurConstante() - m_valeur) < 1e-9;
}

// ============== VARIABLE ==================
Variable::Variable(const std::string& nom) : m_nom(nom) {}
double Variable::eval(double x) const { return x; }
ExprPtr Variable::derivee() const { return cst(1.0); }
ExprPtr Variable::simplifier() const { return clone(); }
void Variable::afficher(std::ostream& os) const { os << m_nom; }
ExprPtr Variable::clone() const { return var(m_nom); }
bool Variable::estEgal(const ASTNode& autre) const {
    const Variable* v = dynamic_cast<const Variable*>(&autre);
    return v != nullptr && v->m_nom == m_nom;
}

// ============== OP Binaire Base ==================
OperateurBinaire::OperateurBinaire(ExprPtr gauche, ExprPtr droite) : m_gauche(gauche), m_droite(droite) {}

// ============== ADDITION ==================
Addition::Addition(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}
double Addition::eval(double x) const { return m_gauche->eval(x) + m_droite->eval(x); }
ExprPtr Addition::derivee() const { return m_gauche->derivee() + m_droite->derivee(); }
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
void Addition::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " + "; m_droite->afficher(os); os << ")";
}
ExprPtr Addition::clone() const { return std::make_shared<Addition>(m_gauche->clone(), m_droite->clone()); }
bool Addition::estEgal(const ASTNode& autre) const {
    const Addition* a = dynamic_cast<const Addition*>(&autre);
    if (!a) return false;
    return (m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite))) ||
           (m_gauche->estEgal(*(a->m_droite)) && m_droite->estEgal(*(a->m_gauche)));
}

// ============== SOUSTRACTION ==================
Soustraction::Soustraction(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}
double Soustraction::eval(double x) const { return m_gauche->eval(x) - m_droite->eval(x); }
ExprPtr Soustraction::derivee() const { return m_gauche->derivee() - m_droite->derivee(); }
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
void Soustraction::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " - "; m_droite->afficher(os); os << ")";
}
ExprPtr Soustraction::clone() const { return std::make_shared<Soustraction>(m_gauche->clone(), m_droite->clone()); }
bool Soustraction::estEgal(const ASTNode& autre) const {
    const Soustraction* a = dynamic_cast<const Soustraction*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== MULTIPLICATION ==================
Multiplication::Multiplication(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}
double Multiplication::eval(double x) const { return m_gauche->eval(x) * m_droite->eval(x); }
ExprPtr Multiplication::derivee() const {
    return (m_gauche->derivee() * m_droite) + (m_gauche * m_droite->derivee());
}
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
void Multiplication::afficher(std::ostream& os) const {
    m_gauche->afficher(os); os << " * "; m_droite->afficher(os);
}
ExprPtr Multiplication::clone() const { return std::make_shared<Multiplication>(m_gauche->clone(), m_droite->clone()); }
bool Multiplication::estEgal(const ASTNode& autre) const {
    const Multiplication* a = dynamic_cast<const Multiplication*>(&autre);
    if (!a) return false;
    return (m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite))) ||
           (m_gauche->estEgal(*(a->m_droite)) && m_droite->estEgal(*(a->m_gauche)));
}

// ============== DIVISION ==================
Division::Division(ExprPtr gauche, ExprPtr droite) : OperateurBinaire(gauche, droite) {}
double Division::eval(double x) const { return m_gauche->eval(x) / m_droite->eval(x); }
ExprPtr Division::derivee() const { // (u'v - uv') / v^2
    auto num = (m_gauche->derivee() * m_droite) - (m_gauche * m_droite->derivee());
    auto den = m_droite * m_droite;
    return num / den;
}
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
void Division::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << " / "; m_droite->afficher(os); os << ")";
}
ExprPtr Division::clone() const { return std::make_shared<Division>(m_gauche->clone(), m_droite->clone()); }
bool Division::estEgal(const ASTNode& autre) const {
    const Division* a = dynamic_cast<const Division*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== PUISSANCE ==================
Puissance::Puissance(ExprPtr base, ExprPtr exposant) : OperateurBinaire(base, exposant) {}
double Puissance::eval(double x) const { return std::pow(m_gauche->eval(x), m_droite->eval(x)); }
ExprPtr Puissance::derivee() const {
    if (m_droite->estConstante()) { // (u^n)' = n*u^{n-1}*u'
        double n = m_droite->getValeurConstante();
        if (n == 0) return cst(0.0);
        return cst(n) * ast_pow(m_gauche, cst(n - 1)) * m_gauche->derivee();
    }
    return cst(0.0); 
}
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
void Puissance::afficher(std::ostream& os) const {
    os << "("; m_gauche->afficher(os); os << ")^("; m_droite->afficher(os); os << ")";
}
ExprPtr Puissance::clone() const { return std::make_shared<Puissance>(m_gauche->clone(), m_droite->clone()); }
bool Puissance::estEgal(const ASTNode& autre) const {
    const Puissance* a = dynamic_cast<const Puissance*>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== FONCTION UNAIRE BASE =================
FonctionUnaire::FonctionUnaire(ExprPtr arg) : m_argument(arg) {}

// ============== SINUS ==================
Sinus::Sinus(ExprPtr arg) : FonctionUnaire(arg) {}
double Sinus::eval(double x) const { return std::sin(m_argument->eval(x)); }
ExprPtr Sinus::derivee() const { return ast_cos(m_argument) * m_argument->derivee(); }
ExprPtr Sinus::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::sin(arg->getValeurConstante()));
    return ast_sin(arg);
}
void Sinus::afficher(std::ostream& os) const {
    os << "sin("; m_argument->afficher(os); os << ")";
}
ExprPtr Sinus::clone() const { return std::make_shared<Sinus>(m_argument->clone()); }
bool Sinus::estEgal(const ASTNode& autre) const {
    const Sinus* a = dynamic_cast<const Sinus*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

// ============== COSINUS ==================
Cosinus::Cosinus(ExprPtr arg) : FonctionUnaire(arg) {}
double Cosinus::eval(double x) const { return std::cos(m_argument->eval(x)); }
ExprPtr Cosinus::derivee() const { return (cst(-1.0) * ast_sin(m_argument)) * m_argument->derivee(); }
ExprPtr Cosinus::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::cos(arg->getValeurConstante()));
    return ast_cos(arg);
}
void Cosinus::afficher(std::ostream& os) const {
    os << "cos("; m_argument->afficher(os); os << ")";
}
ExprPtr Cosinus::clone() const { return std::make_shared<Cosinus>(m_argument->clone()); }
bool Cosinus::estEgal(const ASTNode& autre) const {
    const Cosinus* a = dynamic_cast<const Cosinus*>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}


// --- Surcharge d'Opérateurs =---

ExprPtr cst(double valeur) { return std::make_shared<Constante>(valeur); }
ExprPtr var(const std::string& nom) { return std::make_shared<Variable>(nom); }

ExprPtr operator+(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Addition>(gauche, droite); }
ExprPtr operator+(ExprPtr gauche, double droite) { return std::make_shared<Addition>(gauche, cst(droite)); }
ExprPtr operator+(double gauche, ExprPtr droite) { return std::make_shared<Addition>(cst(gauche), droite); }

ExprPtr operator-(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Soustraction>(gauche, droite); }
ExprPtr operator-(ExprPtr gauche, double droite) { return std::make_shared<Soustraction>(gauche, cst(droite)); }
ExprPtr operator-(double gauche, ExprPtr droite) { return std::make_shared<Soustraction>(cst(gauche), droite); }

ExprPtr operator*(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Multiplication>(gauche, droite); }
ExprPtr operator*(ExprPtr gauche, double droite) { return std::make_shared<Multiplication>(gauche, cst(droite)); }
ExprPtr operator*(double gauche, ExprPtr droite) { return std::make_shared<Multiplication>(cst(gauche), droite); }

ExprPtr operator/(ExprPtr gauche, ExprPtr droite) { return std::make_shared<Division>(gauche, droite); }
ExprPtr operator/(ExprPtr gauche, double droite) { return std::make_shared<Division>(gauche, cst(droite)); }
ExprPtr operator/(double gauche, ExprPtr droite) { return std::make_shared<Division>(cst(gauche), droite); }

ExprPtr ast_pow(ExprPtr base, ExprPtr exposant) { return std::make_shared<Puissance>(base, exposant); }
ExprPtr ast_pow(ExprPtr base, double exposant) { return std::make_shared<Puissance>(base, cst(exposant)); }
ExprPtr ast_pow(double base, ExprPtr exposant) { return std::make_shared<Puissance>(cst(base), exposant); }

ExprPtr ast_sin(ExprPtr arg) { return std::make_shared<Sinus>(arg); }
ExprPtr ast_cos(ExprPtr arg) { return std::make_shared<Cosinus>(arg); }
