/**
 * @file Regles.cpp
 * @brief Règles propres à chaque noeud : évaluation, dérivation, simplification et
 *        primitives.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"

#include <cmath>
#include <stdexcept>

namespace symalgo {

namespace {

/*
 * Nom : coefficientLineaire
 * Description : Si u est une fonction affine de x (u = a*x + b avec a non nul, a pouvant
 *               être symbolique comme pi ou C1), renvoie vrai et écrit a : la dérivée de u
 *               ne dépend alors pas de x. Sert à la règle d'intégration F(ax+b)/a.
 */
bool coefficientLineaire(const ExprPtr& u, ExprPtr& a) {
    if (!u->contientVariable()) return false;
    const ExprPtr d = u->derivee();
    if (d->contientVariable()) return false;
    const Constante* c = comme<Constante>(d);
    if (c && (c->getNombre().estZero() || !c->getNombre().estFini())) return false;
    a = d;
    return true;
}

// e / a
ExprPtr diviserPar(const ExprPtr& e, const ExprPtr& a) { return e / a; }

// Dérivée logarithmique d'un facteur b^e : e' ln(b) + e b'/b
ExprPtr deriveeLogarithmique(const ExprPtr& base, const ExprPtr& exposant, CacheDerivees& cache) {
    const ExprPtr db = base->derivee(cache);
    const ExprPtr de = exposant->derivee(cache);
    const bool dbNul = comme<Constante>(db) && comme<Constante>(db)->getNombre().estZero();
    const bool deNul = comme<Constante>(de) && comme<Constante>(de)->getNombre().estZero();
    if (dbNul && deNul) return nullptr;
    AccumulateurSomme s;
    if (!dbNul) s.ajouter(produit({exposant, db, ast_pow(base, nombre(Nombre(-1)))}), Nombre(1));
    if (!deNul) s.ajouter(de * ast_ln(base), Nombre(1));
    return s.construire();
}

} // namespace

// ============================================================================
// Constante, Variable, Parametre
// ============================================================================

double Constante::eval(double) const { return m_approx; }
ExprPtr Constante::calculerDerivee(CacheDerivees&) const { return nombre(Nombre(0)); }
ExprPtr Constante::calculerSimplification() const { return clone(); }
ExprPtr Constante::primitive() const { return clone() * var("x"); }

double Variable::eval(double x) const { return x; }
ExprPtr Variable::calculerDerivee(CacheDerivees&) const { return un(); }
ExprPtr Variable::calculerSimplification() const { return clone(); }
ExprPtr Variable::primitive() const { return frac(1, 2) * ast_pow(clone(), 2.0); }

double Parametre::eval(double) const {
    throw std::logic_error("Impossible d'evaluer le parametre symbolique '" + m_nom +
                           "' : il n'a pas de valeur numerique");
}
ExprPtr Parametre::calculerDerivee(CacheDerivees&) const { return nombre(Nombre(0)); }
ExprPtr Parametre::calculerSimplification() const { return clone(); }
ExprPtr Parametre::primitive() const { return clone() * var("x"); }

double Pi::eval(double) const { return 3.14159265358979323846; }
ExprPtr Pi::calculerDerivee(CacheDerivees&) const { return nombre(Nombre(0)); }
ExprPtr Pi::calculerSimplification() const { return clone(); }
ExprPtr Pi::primitive() const { return clone() * var("x"); }

// ============================================================================
// Somme
// ============================================================================

double Somme::eval(double x) const {
    double s = m_constante.versDouble();
    for (const Terme& t : m_termes) s += t.coefficient.versDouble() * t.expression->eval(x);
    return s;
}

ExprPtr Somme::calculerDerivee(CacheDerivees& cache) const {
    AccumulateurSomme s;
    for (const Terme& t : m_termes) s.ajouter(t.expression->derivee(cache), t.coefficient);
    return s.construire();
}

ExprPtr Somme::calculerSimplification() const {
    AccumulateurSomme s;
    s.ajouter(nombre(m_constante), Nombre(1));
    for (const Terme& t : m_termes) s.ajouter(t.expression->simplifier(), t.coefficient);
    return s.construire();
}

ExprPtr Somme::primitive() const {
    // Linéarité ; une somme dont un terme n'est pas intégrable garde le reste calculé
    AccumulateurSomme s;
    if (!m_constante.estZero()) s.ajouter(var("x"), m_constante);
    for (const Terme& t : m_termes) s.ajouter(t.expression->integrer(), t.coefficient);
    return s.construire();
}

// ============================================================================
// Produit
// ============================================================================

double Produit::eval(double x) const {
    double r = m_coefficient.versDouble();
    for (const Facteur& f : m_facteurs) {
        const double b = f.base->eval(x);
        long long k;
        const Constante* e = comme<Constante>(f.exposant);
        r *= e && e->getNombre().versEntier(k) ? puissanceEntiereReelle(b, k) : std::pow(b, f.exposant->eval(x));
    }
    return r;
}

ExprPtr Produit::calculerDerivee(CacheDerivees& cache) const {
    // (prod b_i^e_i)' = P * sum_i (e_i' ln b_i + e_i b_i' / b_i) : la multiplication
    // canonique fusionne ensuite b_i^e_i * b_i^(-1) en b_i^(e_i - 1)
    AccumulateurSomme s;
    const ExprPtr soi = clone();
    for (const Facteur& f : m_facteurs) {
        if (const ExprPtr dlog = deriveeLogarithmique(f.base, f.exposant, cache)) s.ajouter(soi * dlog, Nombre(1));
    }
    return s.construire();
}

ExprPtr Produit::calculerSimplification() const {
    AccumulateurProduit p;
    p.multiplier(nombre(m_coefficient));
    for (const Facteur& f : m_facteurs) p.multiplier(ast_pow(f.base->simplifier(), f.exposant->simplifier()));
    return p.construire();
}

ExprPtr Produit::primitive() const {
    // Les facteurs indépendants de x sortent de l'intégrale
    AccumulateurProduit constant, variable;
    constant.multiplier(nombre(m_coefficient));
    int nombreVariables = 0;
    for (const Facteur& f : m_facteurs) {
        const ExprPtr facteur = ast_pow(f.base, f.exposant);
        if (facteur->contientVariable()) {
            variable.multiplier(facteur);
            ++nombreVariables;
        } else {
            constant.multiplier(facteur);
        }
    }
    if (nombreVariables != 1) return integraleNonEvaluee(); // intégration par parties non gérée
    return constant.construire() * variable.construire()->integrer();
}

// ============================================================================
// Puissance
// ============================================================================

double Puissance::eval(double x) const {
    const double b = m_base->eval(x);
    return m_exposantEntier != PAS_ENTIER ? puissanceEntiereReelle(b, m_exposantEntier) : std::pow(b, m_exposant->eval(x));
}

ExprPtr Puissance::calculerDerivee(CacheDerivees& cache) const {
    const ExprPtr dlog = deriveeLogarithmique(m_base, m_exposant, cache);
    return dlog ? clone() * dlog : nombre(Nombre(0));
}

ExprPtr Puissance::calculerSimplification() const {
    return ast_pow(m_base->simplifier(), m_exposant->simplifier());
}

ExprPtr Puissance::primitive() const {
    ExprPtr a;
    // (a*x + b)^n, n constant : (ax+b)^(n+1) / (a(n+1)), ou ln(ax+b)/a si n = -1
    if (!m_exposant->contientVariable() && coefficientLineaire(m_base, a)) {
        if (const Constante* n = comme<Constante>(m_exposant)) {
            const Nombre& nv = n->getNombre();
            if (nv.estMoinsUn()) return diviserPar(ast_ln(m_base), a);
            const Nombre suivant = nv + Nombre(1);
            return diviserPar(ast_pow(m_base, nombre(suivant)), a * nombre(suivant));
        }
    }
    // c^(a*x + b), c constante positive différente de 1 : c^u / (a ln c)
    if (!m_base->contientVariable() && coefficientLineaire(m_exposant, a)) {
        const Constante* c = comme<Constante>(m_base);
        if (c && c->getNombre().signe() > 0 && !c->getNombre().estUn()) {
            return clone() / (a * ast_ln(m_base));
        }
    }
    return integraleNonEvaluee();
}

// ============================================================================
// Fonctions unaires
// ============================================================================

double Sinus::eval(double x) const { return std::sin(m_argument->eval(x)); }
double Cosinus::eval(double x) const { return std::cos(m_argument->eval(x)); }
double Tangente::eval(double x) const { return std::tan(m_argument->eval(x)); }
double Exponentielle::eval(double x) const { return std::exp(m_argument->eval(x)); }
double Logarithme::eval(double x) const { return std::log(m_argument->eval(x)); }

ExprPtr Sinus::calculerDerivee(CacheDerivees& cache) const { return ast_cos(m_argument) * m_argument->derivee(cache); }
ExprPtr Cosinus::calculerDerivee(CacheDerivees& cache) const { return -(ast_sin(m_argument) * m_argument->derivee(cache)); }
ExprPtr Tangente::calculerDerivee(CacheDerivees& cache) const {
    return (1.0 + ast_pow(clone(), 2.0)) * m_argument->derivee(cache); // 1 + tan^2
}
ExprPtr Exponentielle::calculerDerivee(CacheDerivees& cache) const { return clone() * m_argument->derivee(cache); }
ExprPtr Logarithme::calculerDerivee(CacheDerivees& cache) const { return m_argument->derivee(cache) / m_argument; }

ExprPtr Sinus::calculerSimplification() const { return ast_sin(m_argument->simplifier()); }
ExprPtr Cosinus::calculerSimplification() const { return ast_cos(m_argument->simplifier()); }
ExprPtr Tangente::calculerSimplification() const { return ast_tan(m_argument->simplifier()); }
ExprPtr Exponentielle::calculerSimplification() const { return ast_exp(m_argument->simplifier()); }
ExprPtr Logarithme::calculerSimplification() const { return ast_ln(m_argument->simplifier()); }

// Primitives par substitution linéaire : int f(ax+b) dx = F(ax+b) / a
ExprPtr Sinus::primitive() const {
    ExprPtr a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(-ast_cos(m_argument), a);
    return integraleNonEvaluee();
}

ExprPtr Cosinus::primitive() const {
    ExprPtr a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(ast_sin(m_argument), a);
    return integraleNonEvaluee();
}

ExprPtr Tangente::primitive() const {
    ExprPtr a;
    // -ln(cos(u)), valable là où cos(u) > 0
    if (coefficientLineaire(m_argument, a)) return diviserPar(-ast_ln(ast_cos(m_argument)), a);
    return integraleNonEvaluee();
}

ExprPtr Exponentielle::primitive() const {
    ExprPtr a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(clone(), a);
    return integraleNonEvaluee();
}

ExprPtr Logarithme::primitive() const {
    ExprPtr a;
    // u ln(u) - u
    if (coefficientLineaire(m_argument, a)) return diviserPar(m_argument * clone() - m_argument, a);
    return integraleNonEvaluee();
}

// ============================================================================
// Fonctions réciproques
// ============================================================================

double ArcSinus::eval(double x) const { return std::asin(m_argument->eval(x)); }
double ArcCosinus::eval(double x) const { return std::acos(m_argument->eval(x)); }
double ArcTangente::eval(double x) const { return std::atan(m_argument->eval(x)); }

namespace {

// (1 - u^2)^(-1/2)
ExprPtr inverseRacineUnMoinsCarre(const ExprPtr& u) { return ast_pow(1.0 - ast_pow(u, 2.0), frac(-1, 2)); }

} // namespace

ExprPtr ArcSinus::calculerDerivee(CacheDerivees& cache) const {
    return m_argument->derivee(cache) * inverseRacineUnMoinsCarre(m_argument);
}
ExprPtr ArcCosinus::calculerDerivee(CacheDerivees& cache) const {
    return -(m_argument->derivee(cache) * inverseRacineUnMoinsCarre(m_argument));
}
ExprPtr ArcTangente::calculerDerivee(CacheDerivees& cache) const {
    return m_argument->derivee(cache) / (1.0 + ast_pow(m_argument, 2.0));
}

ExprPtr ArcSinus::calculerSimplification() const { return ast_asin(m_argument->simplifier()); }
ExprPtr ArcCosinus::calculerSimplification() const { return ast_acos(m_argument->simplifier()); }
ExprPtr ArcTangente::calculerSimplification() const { return ast_atan(m_argument->simplifier()); }

// Primitives par intégration par parties, pour un argument affine u = a*x + b
ExprPtr ArcSinus::primitive() const {
    ExprPtr a;
    if (!coefficientLineaire(m_argument, a)) return integraleNonEvaluee();
    return diviserPar(m_argument * clone() + ast_pow(1.0 - ast_pow(m_argument, 2.0), frac(1, 2)), a);
}
ExprPtr ArcCosinus::primitive() const {
    ExprPtr a;
    if (!coefficientLineaire(m_argument, a)) return integraleNonEvaluee();
    return diviserPar(m_argument * clone() - ast_pow(1.0 - ast_pow(m_argument, 2.0), frac(1, 2)), a);
}
ExprPtr ArcTangente::primitive() const {
    ExprPtr a;
    if (!coefficientLineaire(m_argument, a)) return integraleNonEvaluee();
    return diviserPar(m_argument * clone() - frac(1, 2) * ast_ln(1.0 + ast_pow(m_argument, 2.0)), a);
}

// ============================================================================
// Noeuds non évalués
// ============================================================================

double IntegraleNonEvaluee::eval(double) const {
    throw std::logic_error("Impossible d'evaluer une primitive non calculee symboliquement");
}
// Théorème fondamental de l'analyse : (int f)' = f
ExprPtr IntegraleNonEvaluee::calculerDerivee(CacheDerivees&) const { return m_integrande; }
ExprPtr IntegraleNonEvaluee::calculerSimplification() const {
    return fabriquer<IntegraleNonEvaluee>(m_integrande->simplifier());
}
ExprPtr IntegraleNonEvaluee::primitive() const { return integraleNonEvaluee(); }
ExprPtr IntegraleNonEvaluee::calculerLimite(double a) const { return limiteNonEvaluee(a); }

double LimiteNonEvaluee::eval(double) const {
    throw std::logic_error("Impossible d'evaluer une limite non determinee");
}
ExprPtr LimiteNonEvaluee::calculerDerivee(CacheDerivees&) const { return nombre(Nombre(0)); }
ExprPtr LimiteNonEvaluee::calculerSimplification() const { return clone(); }
ExprPtr LimiteNonEvaluee::primitive() const { return clone() * var("x"); }
ExprPtr LimiteNonEvaluee::calculerLimite(double) const { return clone(); }

} // namespace symalgo
