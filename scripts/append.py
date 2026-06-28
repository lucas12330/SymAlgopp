code = """
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
ExprPtr Constante::limite(double a) const {
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
"""

with open("src/ASTNode.cpp", "a") as f:
    f.write(code)
print("Code appended successfully!")
