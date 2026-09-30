#include "ASTNode.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric> // Pour std::gcd
#include <stdexcept>
#include <string>
#include <vector>

namespace symalgo {

namespace {

/*
 * Nom : extraireCoefficient
 * Description : Décompose e en coeff * u (coeff constante, 1 par défaut), pour la
 *               factorisation a*U + b*U = (a+b)*U. Le coefficient reste un noeud afin
 *               de préserver l'exactitude des fractions.
 * Utilisation : ExprPtr c, u; extraireCoefficient(expr, c, u);
 */
void extraireCoefficient(const ExprPtr& e, ExprPtr& coeff, ExprPtr& u) {
    if (const Multiplication* m = comme<Multiplication>(e.get())) {
        if (m->m_gauche->estConstante()) {
            coeff = m->m_gauche;
            u = m->m_droite;
            return;
        }
        if (m->m_droite->estConstante()) {
            coeff = m->m_droite;
            u = m->m_gauche;
            return;
        }
    }
    coeff = cst(1.0);
    u = e;
}

// Vrai si e est une constante valant exactement v
bool estValeur(const ExprPtr& e, double v) {
    return e->estConstante() && e->getValeurConstante() == v;
}

bool estFraction(const ExprPtr& e) { return comme<Fraction>(e.get()) != nullptr; }

/*
 * Nom : commeRationnel
 * Description : Écrit n/d si la constante e est un rationnel représentable exactement :
 *               une Fraction, ou une Constante entière (|v| <= 2^53).
 */
bool commeRationnel(const ExprPtr& e, int64_t& n, int64_t& d) {
    if (const Fraction* f = comme<Fraction>(e.get())) {
        n = f->getNum();
        d = f->getDen();
        return true;
    }
    if (const Constante* c = comme<Constante>(e.get())) {
        const double v = c->getValeurConstante();
        if (std::floor(v) == v && std::abs(v) <= 9007199254740992.0) {
            n = static_cast<int64_t>(v);
            d = 1;
            return true;
        }
    }
    return false;
}

/*
 * Nom : plierConstantes
 * Description : Calcule g op d pour deux constantes (op parmi + - * /). Le calcul est exact
 *               (résultat Fraction) dès qu'une Fraction est en jeu et que les deux opérandes
 *               sont rationnels, sauf débordement de int64 ; sinon il se fait en double.
 * Utilisation : ExprPtr r = plierConstantes(frac(1, 3), frac(1, 6), '+'); // (1/2)
 */
ExprPtr plierConstantes(const ExprPtr& g, const ExprPtr& d, char op) {
    int64_t n1, d1, n2, d2;
    if ((estFraction(g) || estFraction(d)) && commeRationnel(g, n1, d1) && commeRationnel(d, n2, d2)) {
        int64_t n = 0, q = 1, t1 = 0, t2 = 0;
        bool debordement = false;
        switch (op) {
            case '+':
            case '-':
                debordement = __builtin_mul_overflow(n1, d2, &t1) || __builtin_mul_overflow(n2, d1, &t2) ||
                              __builtin_mul_overflow(d1, d2, &q) ||
                              (op == '+' ? __builtin_add_overflow(t1, t2, &n) : __builtin_sub_overflow(t1, t2, &n));
                break;
            case '*':
                debordement = __builtin_mul_overflow(n1, n2, &n) || __builtin_mul_overflow(d1, d2, &q);
                break;
            default: // '/'
                debordement = n2 == 0 || __builtin_mul_overflow(n1, d2, &n) || __builtin_mul_overflow(d1, n2, &q);
                break;
        }
        if (!debordement) return frac(n, q);
    }
    const double a = g->getValeurConstante();
    const double b = d->getValeurConstante();
    switch (op) {
        case '+': return cst(a + b);
        case '-': return cst(a - b);
        case '*': return cst(a * b);
        default: return cst(a / b);
    }
}

// Écrit e = base^exposant (u seul vaut u^1)
void baseEtExposant(const ExprPtr& e, ExprPtr& base, ExprPtr& exposant) {
    if (const Puissance* p = comme<Puissance>(e.get())) {
        base = p->m_gauche;
        exposant = p->m_droite;
        return;
    }
    base = e;
    exposant = cst(1.0);
}

ExprPtr combinerPuissance(const ExprPtr& b, const ExprPtr& p);

/*
 * Nom : combinerProduit
 * Description : Simplifie g * d, g et d étant déjà simplifiés : repli des constantes (qui
 *               sont regroupées à gauche : c1 * (c2 * u) = (c1 c2) * u), éléments neutre et
 *               absorbant, et u^a * u^b = u^(a+b).
 */
ExprPtr combinerProduit(ExprPtr g, ExprPtr d) {
    if (d->estConstante() && !g->estConstante()) std::swap(g, d);
    if (g->estConstante()) {
        if (d->estConstante()) return plierConstantes(g, d, '*');
        if (estValeur(g, 0.0)) return cst(0.0);
        if (estValeur(g, 1.0)) return d;
        const Multiplication* m = comme<Multiplication>(d.get());
        if (m && m->m_gauche->estConstante()) {
            return combinerProduit(plierConstantes(g, m->m_gauche, '*'), m->m_droite);
        }
        return fabriquer<Multiplication>(g, d);
    }
    // Les constantes remontent à gauche : u * (c * v) = c * (u * v)
    if (const Multiplication* m = comme<Multiplication>(d.get()); m && m->m_gauche->estConstante()) {
        return combinerProduit(m->m_gauche, combinerProduit(g, m->m_droite));
    }
    if (const Multiplication* m = comme<Multiplication>(g.get()); m && m->m_gauche->estConstante()) {
        return combinerProduit(m->m_gauche, combinerProduit(m->m_droite, d));
    }
    ExprPtr bg, eg, bd, ed;
    baseEtExposant(g, bg, eg);
    baseEtExposant(d, bd, ed);
    if (eg->estConstante() && ed->estConstante() && bg->estEgal(*bd)) {
        return combinerPuissance(bg, plierConstantes(eg, ed, '+'));
    }
    return fabriquer<Multiplication>(g, d);
}

/*
 * Nom : combinerPuissance
 * Description : Simplifie b^p, b et p étant déjà simplifiés.
 */
ExprPtr combinerPuissance(const ExprPtr& b, const ExprPtr& p) {
    if (b->estConstante() && p->estConstante()) return cst(std::pow(b->getValeurConstante(), p->getValeurConstante()));
    if (estValeur(p, 0.0)) return cst(1.0);
    if (estValeur(p, 1.0)) return b;
    if (estValeur(b, 0.0)) return cst(0.0);
    if (estValeur(b, 1.0)) return cst(1.0);
    return ast_pow(b, p);
}

} // namespace

/*
 * Nom : coefficientLineaire
 * Description : Si u est une fonction affine de x (u = a*x + b avec a non nul), renvoie
 *               vrai et écrit a. Détection par la dérivée : u' doit se simplifier en une
 *               constante numérique. Sert à la règle d'intégration F(ax+b)/a.
 * Utilisation : double a; if (coefficientLineaire(u, a)) { ... }
 */
static bool coefficientLineaire(const ExprPtr& u, double& a) {
    if (!u->contientVariable()) return false;
    const ExprPtr d = u->derivee()->simplifier();
    if (!d->estConstante()) return false;
    a = d->getValeurConstante();
    return a != 0.0 && std::isfinite(a);
}

/*
 * Nom : diviserPar
 * Description : Renvoie e / a sous la forme (1/a) * e, ou e inchangée si a vaut 1.
 * Utilisation : ExprPtr r = diviserPar(ast_sin(u), 2.0);
 */
static ExprPtr diviserPar(const ExprPtr& e, double a) {
    if (a == 1.0) return e;
    return cst(1.0 / a) * e;
}

// ============== OUTILS POUR LES LIMITES ==================

namespace {

constexpr double INFINI = std::numeric_limits<double>::infinity();
constexpr double NAN_LIMITE = std::numeric_limits<double>::quiet_NaN();

// Profondeur maximale des réécritures récursives (L'Hôpital, 0·∞, 1^∞...)
constexpr int PROFONDEUR_MAX_REECRITURE = 8;
thread_local int g_profondeurReecriture = 0;

struct GardeProfondeur {
    GardeProfondeur() { ++g_profondeurReecriture; }
    ~GardeProfondeur() { --g_profondeurReecriture; }
    GardeProfondeur(const GardeProfondeur&) = delete;
    GardeProfondeur& operator=(const GardeProfondeur&) = delete;
};

bool reecriturePossible() { return g_profondeurReecriture < PROFONDEUR_MAX_REECRITURE; }

// Seuil en dessous duquel une limite numérique est considérée comme nulle
// (absorbe les résidus d'arrondi comme sin(pi) = 1.2e-16)
bool estNul(double v) { return std::abs(v) < 1e-12; }

enum class GenreLimite { Nombre, Symbolique, NonEvaluee };

/*
 * Nom : analyserLimite
 * Description : Classe le résultat d'un calcul de limite : nombre (fini ou infini, écrit
 *               dans v), expression symbolique (paramètres) ou limite non déterminée.
 */
GenreLimite analyserLimite(const ExprPtr& limite, double& v) {
    if (comme<LimiteNonEvaluee>(limite.get())) return GenreLimite::NonEvaluee;
    const ExprPtr s = limite->simplifier();
    if (comme<LimiteNonEvaluee>(s.get())) return GenreLimite::NonEvaluee;
    if (s->estConstante()) {
        v = s->getValeurConstante();
        return std::isnan(v) ? GenreLimite::NonEvaluee : GenreLimite::Nombre;
    }
    return GenreLimite::Symbolique;
}

// Constante v, ou nullptr si v n'est pas un nombre (forme indéterminée)
ExprPtr nombreOuNul(double v) { return std::isnan(v) ? nullptr : cst(v); }

// Vrai si le résultat d'une réécriture est une limite déterminée
bool estDeterminee(const ExprPtr& limite) {
    double v;
    return analyserLimite(limite, v) != GenreLimite::NonEvaluee;
}

/*
 * Nom : premierTermeNonNul
 * Description : Ordre k et signe du premier terme non nul du développement de Taylor de f
 *               en a (f(a) étant nul). Au voisinage de a, f(x) ~ f^(k)(a)/k! (x-a)^k :
 *               k pair => f garde le même signe des deux côtés de a, k impair => il change.
 */
bool premierTermeNonNul(const ExprPtr& f, double a, int& ordre, double& signe) {
    ExprPtr d = f;
    for (int k = 1; k <= PROFONDEUR_MAX_REECRITURE; ++k) {
        d = d->derivee()->simplifier();
        double v;
        try {
            v = d->eval(a);
        } catch (const std::logic_error&) {
            return false; // paramètre symbolique : signe inconnu
        }
        if (!std::isfinite(v)) return false;
        if (!estNul(v)) {
            ordre = k;
            signe = v > 0.0 ? 1.0 : -1.0;
            return true;
        }
    }
    return false;
}

/*
 * Nom : limiteUnaire
 * Description : Limite de f(u) : f appliquée à la limite de l'argument u. Renvoie nullptr
 *               si elle n'est pas déterminée.
 */
ExprPtr limiteUnaire(const ExprPtr& argument, double a,
                     double (*fNum)(double), ExprPtr (*fSym)(ExprPtr)) {
    double u;
    const ExprPtr L = argument->limite(a);
    switch (analyserLimite(L, u)) {
        case GenreLimite::Nombre: return nombreOuNul(fNum(u));
        case GenreLimite::Symbolique: return fSym(L);
        case GenreLimite::NonEvaluee: break;
    }
    return nullptr;
}

/*
 * Nom : commeQuotient
 * Description : Si e est un quotient (Division, ou puissance d'exposant constant négatif),
 *               renvoie vrai et écrit son numérateur et son dénominateur.
 */
bool commeQuotient(const ExprPtr& e, ExprPtr& num, ExprPtr& den) {
    if (const Division* d = comme<Division>(e.get())) {
        num = d->m_gauche;
        den = d->m_droite;
        return true;
    }
    const Puissance* p = comme<Puissance>(e.get());
    if (p && !p->m_droite->contientVariable()) {
        const ExprPtr n = p->m_droite->simplifier();
        if (n->estConstante() && n->getValeurConstante() < 0.0) {
            num = cst(1.0);
            den = ast_pow(p->m_gauche, cst(-n->getValeurConstante()));
            return true;
        }
    }
    return false;
}

/*
 * Nom : limiteParQuotientUnique
 * Description : Réécrit g op d (produit ou quotient) comme une seule fraction
 *               (ex. x * (1/x) -> x / x) puis en calcule la limite. Utile quand un
 *               facteur n'a pas de limite propre alors que l'ensemble en a une.
 *               Renvoie nullptr si la réécriture est impossible ou n'aboutit pas.
 */
ExprPtr limiteParQuotientUnique(const ExprPtr& g, const ExprPtr& d, bool estDivision, double a) {
    ExprPtr ng = g, dg = cst(1.0), nd = d, dd = cst(1.0);
    const bool qg = commeQuotient(g, ng, dg);
    const bool qd = commeQuotient(d, nd, dd);
    if ((!qg && !qd) || !reecriturePossible()) return nullptr;
    GardeProfondeur garde;
    // (ng/dg) * (nd/dd) = (ng nd) / (dg dd) ; (ng/dg) / (nd/dd) = (ng dd) / (dg nd)
    const ExprPtr num = estDivision ? ng * dd : ng * nd;
    const ExprPtr den = estDivision ? dg * nd : dg * dd;
    const ExprPtr r = (num->simplifier() / den->simplifier())->limite(a);
    return estDeterminee(r) ? r : nullptr;
}

} // namespace

// ============== ASTNODE ==================

/*
 * Nom : integrer
 * Description : Point d'entrée de l'intégration : une expression indépendante de x
 *               s'intègre en c*x, sinon on applique les règles du noeud.
 * Utilisation : ExprPtr p = noeud->integrer();
 */
ExprPtr ASTNode::integrer() const {
    if (!contientVariable()) return clone() * var("x");
    return primitive();
}

/*
 * Nom : limite
 * Description : Point d'entrée du calcul de limite en a.
 * Utilisation : ExprPtr l = noeud->limite(a);
 */
ExprPtr ASTNode::limite(double a) const { return calculerLimite(a); }

ExprPtr ASTNode::integraleNonEvaluee() const {
    return fabriquer<IntegraleNonEvaluee>(clone());
}

ExprPtr ASTNode::limiteNonEvaluee(double a) const {
    return fabriquer<LimiteNonEvaluee>(clone(), a);
}

// ============== CONSTANTE ==================

/*
 * Nom : Constante
 * Description : Constructeur de constante qui stocke la valeur passée.
 * Utilisation : Constante c(3.14);
 */
Constante::Constante(CleFabrique, double valeur) : ASTNode(TypeNoeud::Constante), m_valeur(valeur) {}

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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une constante et possède la même valeur aux erreurs flottantes près.
 * Utilisation : bool eq = c.estEgal(autre);
 */
bool Constante::estEgal(const ASTNode& autre) const {
    return autre.estConstante() && std::abs(autre.getValeurConstante() - m_valeur) < 1e-9;
}

// ============== FRACTION ==================

Fraction::Fraction(CleFabrique, int64_t num, int64_t den) : ASTNode(TypeNoeud::Fraction) {
    if (den == 0) {
        throw std::invalid_argument("Denominateur nul dans une Fraction");
    }
    int64_t g = std::gcd(num, den);
    m_num = num / g;
    m_den = den / g;
    if (m_den < 0) {
        m_num = -m_num;
        m_den = -m_den;
    }
    m_valeur_eval = static_cast<double>(m_num) / static_cast<double>(m_den);
}

double Fraction::eval(double /*x*/) const { return m_valeur_eval; }

ExprPtr Fraction::derivee() const { return cst(0.0); }

ExprPtr Fraction::simplifier() const { return clone(); }

void Fraction::afficher(std::ostream& os) const {
    if (m_den == 1) os << m_num;
    else os << "(" << m_num << "/" << m_den << ")";
}

bool Fraction::estEgal(const ASTNode& autre) const {
    if (const Fraction* f = comme<Fraction>(&autre)) {
        return m_num == f->m_num && m_den == f->m_den;
    }
    if (autre.estConstante()) {
        return std::abs(autre.getValeurConstante() - m_valeur_eval) < 1e-9;
    }
    return false;
}

ExprPtr Fraction::primitive() const { return frac(m_num, m_den) * var("x"); }
ExprPtr Fraction::calculerLimite(double /*a*/) const { return clone(); }

// ============== VARIABLE ==================

/*
 * Nom : Variable
 * Description : Constructeur de variable qui l'initialise avec le nom donné.
 * Utilisation : Variable v("t");
 */
Variable::Variable(CleFabrique, const std::string& nom) : ASTNode(TypeNoeud::Variable), m_nom(nom) {}

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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une variable portant le même nom.
 * Utilisation : bool eq = v.estEgal(autre);
 */
bool Variable::estEgal(const ASTNode& autre) const {
    const Variable* v = comme<Variable>(&autre);
    return v != nullptr && v->m_nom == m_nom;
}

// ============== PARAMETRE ==================

Parametre::Parametre(CleFabrique, const std::string& nom) : ASTNode(TypeNoeud::Parametre), m_nom(nom) {}

double Parametre::eval(double) const {
    throw std::logic_error("Impossible d'evaluer le parametre symbolique '" + m_nom +
                           "' : il n'a pas de valeur numerique");
}

ExprPtr Parametre::derivee() const { return cst(0.0); }

ExprPtr Parametre::simplifier() const { return clone(); }

void Parametre::afficher(std::ostream& os) const { os << m_nom; }

bool Parametre::estEgal(const ASTNode& autre) const {
    const Parametre* p = comme<Parametre>(&autre);
    return p != nullptr && p->m_nom == m_nom;
}

ExprPtr Parametre::primitive() const { return clone() * var("x"); }

ExprPtr Parametre::calculerLimite(double /*a*/) const { return clone(); }

// ============== OP Binaire Base ==================

/*
 * Nom : OperateurBinaire
 * Description : Initialise les sous-arbres gauche et droite d'un opérateur.
 * Utilisation : Appelé par les classes dérivées.
 */
OperateurBinaire::OperateurBinaire(TypeNoeud type, ExprPtr gauche, ExprPtr droite)
    : ASTNode(type), m_gauche(std::move(gauche)), m_droite(std::move(droite)) {}

// ============== ADDITION ==================

/*
 * Nom : Addition
 * Description : Constructeur de l'addition.
 * Utilisation : Addition add(gauche, droite);
 */
Addition::Addition(CleFabrique, ExprPtr gauche, ExprPtr droite) : OperateurBinaire(TypeNoeud::Addition, gauche, droite) {}

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
    const ExprPtr g = m_gauche->simplifier();
    const ExprPtr d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante()) return plierConstantes(g, d, '+');
    if (estValeur(g, 0.0)) return d;
    if (estValeur(d, 0.0)) return g;

    // Factorisation a*U + b*U = (a+b)*U
    ExprPtr cG, cD, uG, uD;
    extraireCoefficient(g, cG, uG);
    extraireCoefficient(d, cD, uD);
    if (uG->estEgal(*uD)) return combinerProduit(plierConstantes(cG, cD, '+'), uG);

    return fabriquer<Addition>(g, d);
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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même addition (tient compte de la commutativité).
 * Utilisation : bool eq = add.estEgal(autre);
 */
bool Addition::estEgal(const ASTNode& autre) const {
    const Addition* a = comme<Addition>(&autre);
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
Soustraction::Soustraction(CleFabrique, ExprPtr gauche, ExprPtr droite) : OperateurBinaire(TypeNoeud::Soustraction, gauche, droite) {}

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
    const ExprPtr g = m_gauche->simplifier();
    const ExprPtr d = m_droite->simplifier();
    if (g->estConstante() && d->estConstante()) return plierConstantes(g, d, '-');
    if (estValeur(d, 0.0)) return g;
    if (estValeur(g, 0.0)) return combinerProduit(cst(-1.0), d);

    // Factorisation a*U - b*U = (a-b)*U
    ExprPtr cG, cD, uG, uD;
    extraireCoefficient(g, cG, uG);
    extraireCoefficient(d, cD, uD);
    if (uG->estEgal(*uD)) return combinerProduit(plierConstantes(cG, cD, '-'), uG);

    return fabriquer<Soustraction>(g, d);
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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même soustraction exacte.
 * Utilisation : bool eq = sub.estEgal(autre);
 */
bool Soustraction::estEgal(const ASTNode& autre) const {
    const Soustraction* a = comme<Soustraction>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== MULTIPLICATION ==================

/*
 * Nom : Multiplication
 * Description : Constructeur de la multiplication.
 * Utilisation : Multiplication mul(gauche, droite);
 */
Multiplication::Multiplication(CleFabrique, ExprPtr gauche, ExprPtr droite) : OperateurBinaire(TypeNoeud::Multiplication, gauche, droite) {}

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
    return combinerProduit(m_gauche->simplifier(), m_droite->simplifier());
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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une multiplication identique (commutativité acceptée).
 * Utilisation : bool eq = mul.estEgal(autre);
 */
bool Multiplication::estEgal(const ASTNode& autre) const {
    const Multiplication* a = comme<Multiplication>(&autre);
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
Division::Division(CleFabrique, ExprPtr gauche, ExprPtr droite) : OperateurBinaire(TypeNoeud::Division, gauche, droite) {}

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
    const ExprPtr g = m_gauche->simplifier();
    const ExprPtr d = m_droite->simplifier();
    if (d->estConstante()) {
        if (estValeur(d, 0.0)) return fabriquer<Division>(g, d); // division par zéro conservée
        if (g->estConstante()) return plierConstantes(g, d, '/');
        if (estValeur(d, 1.0)) return g;
        // (c * u) / k = (c/k) * u
        const Multiplication* m = comme<Multiplication>(g.get());
        if (m && m->m_gauche->estConstante()) {
            return combinerProduit(plierConstantes(m->m_gauche, d, '/'), m->m_droite);
        }
        return fabriquer<Division>(g, d);
    }
    if (estValeur(g, 0.0)) return cst(0.0);
    if (g->estEgal(*d)) return cst(1.0);
    return fabriquer<Division>(g, d);
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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est la même division exacte.
 * Utilisation : bool eq = div.estEgal(autre);
 */
bool Division::estEgal(const ASTNode& autre) const {
    const Division* a = comme<Division>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== PUISSANCE ==================

/*
 * Nom : Puissance
 * Description : Constructeur de la puissance (base^exposant).
 * Utilisation : Puissance p(base, exposant);
 */
Puissance::Puissance(CleFabrique, ExprPtr base, ExprPtr exposant) : OperateurBinaire(TypeNoeud::Puissance, base, exposant) {}

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
    if (!m_droite->contientVariable()) { // exposant symbolique (paramètre) : même règle
        return m_droite * ast_pow(m_gauche, m_droite - 1.0) * m_gauche->derivee();
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
    return combinerPuissance(m_gauche->simplifier(), m_droite->simplifier());
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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est une puissance avec la même base et le même exposant.
 * Utilisation : bool eq = p.estEgal(autre);
 */
bool Puissance::estEgal(const ASTNode& autre) const {
    const Puissance* a = comme<Puissance>(&autre);
    return a && m_gauche->estEgal(*(a->m_gauche)) && m_droite->estEgal(*(a->m_droite));
}

// ============== FONCTION UNAIRE BASE =================

/*
 * Nom : FonctionUnaire
 * Description : Constructeur de base des fonctions mathématiques unaires initialisant l'argument.
 * Utilisation : Appelé par les constructeurs des classes filles.
 */
FonctionUnaire::FonctionUnaire(TypeNoeud type, ExprPtr arg) : ASTNode(type), m_argument(std::move(arg)) {}

// ============== SINUS ==================

/*
 * Nom : Sinus
 * Description : Constructeur de la fonction sinus avec l'argument.
 * Utilisation : Sinus s(expr);
 */
Sinus::Sinus(CleFabrique, ExprPtr arg) : FonctionUnaire(TypeNoeud::Sinus, arg) {}

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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est un sinus avec le même argument exact.
 * Utilisation : bool eq = s.estEgal(autre);
 */
bool Sinus::estEgal(const ASTNode& autre) const {
    const Sinus* a = comme<Sinus>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

// ============== COSINUS ==================

/*
 * Nom : Cosinus
 * Description : Constructeur de la fonction cosinus avec l'argument.
 * Utilisation : Cosinus c(expr);
 */
Cosinus::Cosinus(CleFabrique, ExprPtr arg) : FonctionUnaire(TypeNoeud::Cosinus, arg) {}

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
 * Nom : estEgal
 * Description : Vérifie si un autre noeud est un cosinus avec le même argument exact.
 * Utilisation : bool eq = c.estEgal(autre);
 */
bool Cosinus::estEgal(const ASTNode& autre) const {
    const Cosinus* a = comme<Cosinus>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

// ============== TANGENTE ==================

Tangente::Tangente(CleFabrique, ExprPtr arg) : FonctionUnaire(TypeNoeud::Tangente, arg) {}

double Tangente::eval(double x) const { return std::tan(m_argument->eval(x)); }

ExprPtr Tangente::derivee() const {
    // Dérivée de tan(u) = (1 + tan^2(u)) * u'
    auto tan_u = ast_tan(m_argument);
    return (cst(1.0) + ast_pow(tan_u, 2.0)) * m_argument->derivee();
}

ExprPtr Tangente::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::tan(arg->getValeurConstante()));
    return ast_tan(arg);
}

void Tangente::afficher(std::ostream& os) const {
    os << "tan("; m_argument->afficher(os); os << ")";
}

bool Tangente::estEgal(const ASTNode& autre) const {
    const Tangente* a = comme<Tangente>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

ExprPtr Tangente::primitive() const {
    double a;
    // Primitive de tan(u) : -ln(cos(u)), valable là où cos(u) > 0
    if (coefficientLineaire(m_argument, a)) return diviserPar(cst(-1.0) * ast_ln(ast_cos(m_argument)), a);
    return integraleNonEvaluee();
}

ExprPtr Tangente::calculerLimite(double a) const {
    // Aux pôles (cos(u) = 0), tan tend vers +inf d'un côté et -inf de l'autre
    auto tanNum = [](double u) { return estNul(std::cos(u)) ? NAN_LIMITE : std::tan(u); };
    if (ExprPtr r = limiteUnaire(m_argument, a, tanNum, ast_tan)) return r;
    return limiteNonEvaluee(a);
}

// --- Surcharge d'Opérateurs =---

/*
 * Nom : cst
 * Description : Helper générant un noeud Constante à partir d'un double.
 * Utilisation : ExprPtr noeud = cst(3.14);
 */
ExprPtr cst(double valeur) { return fabriquer<Constante>(valeur); }
ExprPtr frac(int64_t num, int64_t den) { return fabriquer<Fraction>(num, den); }

/*
 * Nom : var
 * Description : Helper générant un noeud Variable à partir d'un nom de variable.
 * Utilisation : ExprPtr noeud = var("y");
 */
ExprPtr var(const std::string& nom) { return fabriquer<Variable>(nom); }

/*
 * Nom : param
 * Description : Helper générant un noeud Parametre (constante symbolique).
 * Utilisation : ExprPtr noeud = param("C1");
 */
ExprPtr param(const std::string& nom) { return fabriquer<Parametre>(nom); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre deux ExprPtr générant un noeud Addition.
 * Utilisation : ExprPtr resultat = e1 + e2;
 */
ExprPtr operator+(ExprPtr gauche, ExprPtr droite) { return fabriquer<Addition>(gauche, droite); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre ExprPtr et double générant un noeud Addition.
 * Utilisation : ExprPtr resultat = e + 2.0;
 */
ExprPtr operator+(ExprPtr gauche, double droite) { return fabriquer<Addition>(gauche, cst(droite)); }

/*
 * Nom : operator+
 * Description : Surcharge de l'addition entre double et ExprPtr générant un noeud Addition.
 * Utilisation : ExprPtr resultat = 2.0 + e;
 */
ExprPtr operator+(double gauche, ExprPtr droite) { return fabriquer<Addition>(cst(gauche), droite); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre deux ExprPtr générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = e1 - e2;
 */
ExprPtr operator-(ExprPtr gauche, ExprPtr droite) { return fabriquer<Soustraction>(gauche, droite); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre ExprPtr et double générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = e - 2.0;
 */
ExprPtr operator-(ExprPtr gauche, double droite) { return fabriquer<Soustraction>(gauche, cst(droite)); }

/*
 * Nom : operator-
 * Description : Surcharge de la soustraction entre double et ExprPtr générant un noeud Soustraction.
 * Utilisation : ExprPtr resultat = 2.0 - e;
 */
ExprPtr operator-(double gauche, ExprPtr droite) { return fabriquer<Soustraction>(cst(gauche), droite); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre deux ExprPtr générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = e1 * e2;
 */
ExprPtr operator*(ExprPtr gauche, ExprPtr droite) { return fabriquer<Multiplication>(gauche, droite); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre ExprPtr et double générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = e * 2.0;
 */
ExprPtr operator*(ExprPtr gauche, double droite) { return fabriquer<Multiplication>(gauche, cst(droite)); }

/*
 * Nom : operator*
 * Description : Surcharge de la multiplication entre double et ExprPtr générant un noeud Multiplication.
 * Utilisation : ExprPtr resultat = 2.0 * e;
 */
ExprPtr operator*(double gauche, ExprPtr droite) { return fabriquer<Multiplication>(cst(gauche), droite); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre deux ExprPtr générant un noeud Division.
 * Utilisation : ExprPtr resultat = e1 / e2;
 */
ExprPtr operator/(ExprPtr gauche, ExprPtr droite) { return fabriquer<Division>(gauche, droite); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre ExprPtr et double générant un noeud Division.
 * Utilisation : ExprPtr resultat = e / 2.0;
 */
ExprPtr operator/(ExprPtr gauche, double droite) { return fabriquer<Division>(gauche, cst(droite)); }

/*
 * Nom : operator/
 * Description : Surcharge de la division entre double et ExprPtr générant un noeud Division.
 * Utilisation : ExprPtr resultat = 2.0 / e;
 */
ExprPtr operator/(double gauche, ExprPtr droite) { return fabriquer<Division>(cst(gauche), droite); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance entre deux ExprPtr.
 * Utilisation : ExprPtr resultat = ast_pow(base, exposant);
 */
ExprPtr ast_pow(ExprPtr base, ExprPtr exposant) { return fabriquer<Puissance>(base, exposant); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance avec un exposant double constant.
 * Utilisation : ExprPtr resultat = ast_pow(base, 2.0);
 */
ExprPtr ast_pow(ExprPtr base, double exposant) { return fabriquer<Puissance>(base, cst(exposant)); }

/*
 * Nom : ast_pow
 * Description : Helper générant un noeud Puissance avec une base double constante.
 * Utilisation : ExprPtr resultat = ast_pow(2.0, exposant);
 */
ExprPtr ast_pow(double base, ExprPtr exposant) { return fabriquer<Puissance>(cst(base), exposant); }

/*
 * Nom : ast_sin
 * Description : Helper générant un noeud Sinus.
 * Utilisation : ExprPtr resultat = ast_sin(e);
 */
ExprPtr ast_sin(ExprPtr arg) { return fabriquer<Sinus>(arg); }

/*
 * Nom : ast_cos
 * Description : Helper générant un noeud Cosinus.
 * Utilisation : ExprPtr resultat = ast_cos(e);
 */
ExprPtr ast_cos(ExprPtr arg) { return fabriquer<Cosinus>(arg); }

/*
 * Nom : ast_tan
 * Description : Helper générant un noeud Tangente.
 * Utilisation : ExprPtr resultat = ast_tan(e);
 */
ExprPtr ast_tan(ExprPtr arg) { return fabriquer<Tangente>(arg); }

// ============================================================================
// IMPLÉMENTATION DES LIMITES, INTÉGRALES, DL ET LOGARITHME
// ============================================================================

// --- ASTNode : Développement Limité ---

namespace {

using Serie = std::vector<double>;

// Produit de Cauchy tronqué
Serie produitSeries(const Serie& u, const Serie& v) {
    Serie w(u.size(), 0.0);
    for (size_t k = 0; k < w.size(); ++k)
        for (size_t j = 0; j <= k; ++j) w[k] += u[j] * v[k - j];
    return w;
}

/*
 * Nom : serieTaylor
 * Description : Coefficients de Taylor de e en a jusqu'à l'ordre n (arithmétique des séries
 *               tronquées, comme en différentiation automatique) : coût O(n^2) par noeud au
 *               lieu de la croissance exponentielle des dérivées symboliques successives.
 *               Renvoie faux si un noeud n'est pas pris en charge (paramètre, noeud non
 *               évalué, point singulier) ; l'appelant se replie alors sur les dérivées.
 */
bool serieTaylor(const ASTNode& e, double a, int n, Serie& w) {
    const size_t taille = static_cast<size_t>(n) + 1;
    w.assign(taille, 0.0);
    if (e.estConstante()) {
        w[0] = e.getValeurConstante();
        return true;
    }
    if (comme<Variable>(&e)) {
        w[0] = a;
        if (n >= 1) w[1] = 1.0;
        return true;
    }
    if (const OperateurBinaire* op = comme<OperateurBinaire>(&e)) {
        Serie u, v;
        if (!serieTaylor(*op->m_gauche, a, n, u) || !serieTaylor(*op->m_droite, a, n, v)) return false;
        if (comme<Addition>(op)) {
            for (size_t k = 0; k < taille; ++k) w[k] = u[k] + v[k];
        } else if (comme<Soustraction>(op)) {
            for (size_t k = 0; k < taille; ++k) w[k] = u[k] - v[k];
        } else if (comme<Multiplication>(op)) {
            w = produitSeries(u, v);
        } else if (comme<Division>(op)) {
            // q = u / v : q_k = (u_k - sum_{j>=1} v_j q_{k-j}) / v_0
            if (v[0] == 0.0) return false;
            for (size_t k = 0; k < taille; ++k) {
                double somme = u[k];
                for (size_t j = 1; j <= k; ++j) somme -= v[j] * w[k - j];
                w[k] = somme / v[0];
            }
        } else if (comme<Puissance>(op)) {
            if (op->m_droite->contientVariable()) {
                // u^v = exp(v ln u), défini pour u(a) > 0
                if (u[0] <= 0.0) return false;
                const ExprPtr forme = ast_exp(op->m_droite * ast_ln(op->m_gauche));
                return serieTaylor(*forme, a, n, w);
            }
            const double p = v[0];
            if (u[0] != 0.0) {
                // w = u^p : k u_0 w_k = sum_{j=1..k} (p j - (k - j)) u_j w_{k-j}
                w[0] = std::pow(u[0], p);
                if (!std::isfinite(w[0])) return false;
                for (size_t k = 1; k < taille; ++k) {
                    double somme = 0.0;
                    for (size_t j = 1; j <= k; ++j) {
                        somme += (p * double(j) - double(k - j)) * u[j] * w[k - j];
                    }
                    w[k] = somme / (double(k) * u[0]);
                }
            } else {
                // u(a) = 0 : seule une puissance entière positive est développable
                if (p < 0.0 || std::floor(p) != p) return false;
                w.assign(taille, 0.0);
                w[0] = 1.0;
                if (p > n) {
                    w[0] = 0.0; // u^p = O((x-a)^p), nul jusqu'à l'ordre n
                } else {
                    for (int i = 0; i < static_cast<int>(p); ++i) w = produitSeries(w, u);
                }
            }
        } else {
            return false;
        }
        return true;
    }
    if (const FonctionUnaire* f = comme<FonctionUnaire>(&e)) {
        Serie u;
        if (!serieTaylor(*f->m_argument, a, n, u)) return false;
        if (comme<Exponentielle>(f)) {
            // w' = u' w : k w_k = sum_{j=1..k} j u_j w_{k-j}
            w[0] = std::exp(u[0]);
            for (size_t k = 1; k < taille; ++k) {
                double somme = 0.0;
                for (size_t j = 1; j <= k; ++j) somme += double(j) * u[j] * w[k - j];
                w[k] = somme / double(k);
            }
            return true;
        }
        if (comme<Logarithme>(f)) {
            // w' u = u' : k u_0 w_k = k u_k - sum_{j=1..k-1} j w_j u_{k-j}
            if (u[0] <= 0.0) return false;
            w[0] = std::log(u[0]);
            for (size_t k = 1; k < taille; ++k) {
                double somme = double(k) * u[k];
                for (size_t j = 1; j < k; ++j) somme -= double(j) * w[j] * u[k - j];
                w[k] = somme / (double(k) * u[0]);
            }
            return true;
        }
        // sin et cos se calculent ensemble : s' = u' c, c' = -u' s
        Serie sinus(taille, 0.0), cosinus(taille, 0.0);
        sinus[0] = std::sin(u[0]);
        cosinus[0] = std::cos(u[0]);
        for (size_t k = 1; k < taille; ++k) {
            double ss = 0.0, sc = 0.0;
            for (size_t j = 1; j <= k; ++j) {
                ss += double(j) * u[j] * cosinus[k - j];
                sc -= double(j) * u[j] * sinus[k - j];
            }
            sinus[k] = ss / double(k);
            cosinus[k] = sc / double(k);
        }
        if (comme<Sinus>(f)) {
            w = sinus;
        } else if (comme<Cosinus>(f)) {
            w = cosinus;
        } else if (comme<Tangente>(f)) {
            if (std::abs(cosinus[0]) < 1e-15) return false; // pôle de tan
            for (size_t k = 0; k < taille; ++k) {
                double somme = sinus[k];
                for (size_t j = 1; j <= k; ++j) somme -= cosinus[j] * w[k - j];
                w[k] = somme / cosinus[0];
            }
        } else {
            return false;
        }
        return true;
    }
    return false; // Parametre, noeuds non évalués...
}

// Coefficients de Taylor par dérivations symboliques successives (repli)
Serie serieParDerivation(const ASTNode& e, double a, int n) {
    Serie c(static_cast<size_t>(n) + 1, 0.0);
    ExprPtr deriv = e.simplifier();
    double factorielle = 1.0;
    for (int k = 0; k <= n; ++k) {
        if (k > 0) {
            deriv = deriv->derivee()->simplifier();
            factorielle *= k;
        }
        c[k] = deriv->eval(a) / factorielle;
    }
    return c;
}

} // namespace

/*
 * Nom : DL
 * Description : Développement limité (Taylor) en a à l'ordre donné : sum c_k (x - a)^k.
 *               Lève std::domain_error si la fonction n'est pas développable en a.
 * Utilisation : ExprPtr dl = expr->DL(0.0, 5);
 */
ExprPtr ASTNode::DL(double a, int ordre) const {
    if (ordre < 0) throw std::invalid_argument("DL : ordre negatif");
    Serie c;
    if (!serieTaylor(*this, a, ordre, c)) c = serieParDerivation(*this, a, ordre);

    double echelle = 0.0;
    for (int k = 0; k <= ordre; ++k) {
        if (!std::isfinite(c[k])) {
            throw std::domain_error("DL : la fonction n'est pas developpable a l'ordre " +
                                    std::to_string(k) + " au point demande");
        }
        echelle = std::max(echelle, std::abs(c[k]));
    }
    // Les coefficients négligeables devant les autres (résidus d'arrondi, ex. cos(pi/2))
    // sont omis ; le seuil est relatif car 1/k! devient vite très petit
    const ExprPtr ecart = (var("x") - a)->simplifier();
    ExprPtr resultat = nullptr;
    for (int k = 0; k <= ordre; ++k) {
        if (std::abs(c[k]) <= 1e-15 * echelle || c[k] == 0.0) continue;
        const ExprPtr terme = k == 0 ? cst(c[k]) : cst(c[k]) * ast_pow(ecart, k);
        resultat = resultat ? resultat + terme : terme;
    }
    return resultat ? resultat->simplifier() : cst(0.0);
}

// --- Constante ---
ExprPtr Constante::primitive() const {
    return cst(m_valeur) * var("x");
}
ExprPtr Constante::calculerLimite(double /*a*/) const {
    return cst(m_valeur);
}

// --- Variable ---
ExprPtr Variable::primitive() const {
    return cst(0.5) * ast_pow(var(m_nom), 2);
}
ExprPtr Variable::calculerLimite(double a) const {
    return cst(a);
}

// --- Addition ---
ExprPtr Addition::primitive() const {
    return m_gauche->integrer() + m_droite->integrer();
}
ExprPtr Addition::calculerLimite(double a) const {
    double g = 0.0, d = 0.0;
    const ExprPtr Lg = m_gauche->limite(a);
    const ExprPtr Ld = m_droite->limite(a);
    const GenreLimite tg = analyserLimite(Lg, g);
    const GenreLimite td = analyserLimite(Ld, d);
    if (tg == GenreLimite::NonEvaluee || td == GenreLimite::NonEvaluee) return limiteNonEvaluee(a);
    if (tg == GenreLimite::Symbolique || td == GenreLimite::Symbolique) return Lg + Ld;
    // +inf + -inf donne NaN : forme indéterminée non résolue
    if (ExprPtr r = nombreOuNul(g + d)) return r;
    return limiteNonEvaluee(a);
}

// --- Soustraction ---
ExprPtr Soustraction::primitive() const {
    return m_gauche->integrer() - m_droite->integrer();
}
ExprPtr Soustraction::calculerLimite(double a) const {
    double g = 0.0, d = 0.0;
    const ExprPtr Lg = m_gauche->limite(a);
    const ExprPtr Ld = m_droite->limite(a);
    const GenreLimite tg = analyserLimite(Lg, g);
    const GenreLimite td = analyserLimite(Ld, d);
    if (tg == GenreLimite::NonEvaluee || td == GenreLimite::NonEvaluee) return limiteNonEvaluee(a);
    if (tg == GenreLimite::Symbolique || td == GenreLimite::Symbolique) return Lg - Ld;
    if (ExprPtr r = nombreOuNul(g - d)) return r;
    return limiteNonEvaluee(a);
}

// --- Multiplication ---
ExprPtr Multiplication::primitive() const {
    // Linéarité : un facteur indépendant de x sort de l'intégrale
    if (!m_gauche->contientVariable()) return m_gauche * m_droite->integrer();
    if (!m_droite->contientVariable()) return m_droite * m_gauche->integrer();
    if (m_gauche->estEgal(*m_droite)) return ast_pow(m_gauche, 2.0)->integrer();
    // Intégration par parties non gérée
    return integraleNonEvaluee();
}
ExprPtr Multiplication::calculerLimite(double a) const {
    double g = 0.0, d = 0.0;
    const ExprPtr Lg = m_gauche->limite(a);
    const ExprPtr Ld = m_droite->limite(a);
    const GenreLimite tg = analyserLimite(Lg, g);
    const GenreLimite td = analyserLimite(Ld, d);
    if (tg == GenreLimite::NonEvaluee || td == GenreLimite::NonEvaluee) {
        if (ExprPtr r = limiteParQuotientUnique(m_gauche, m_droite, false, a)) return r;
        return limiteNonEvaluee(a);
    }
    if (tg == GenreLimite::Symbolique || td == GenreLimite::Symbolique) return Lg * Ld;

    const bool zeroFoisInfini = (estNul(g) && std::isinf(d)) || (std::isinf(g) && estNul(d));
    if (!zeroFoisInfini) {
        if (ExprPtr r = nombreOuNul(g * d)) return r;
        return limiteNonEvaluee(a);
    }
    // Forme 0 * inf : d'abord sous forme de fraction unique, puis z * w = z / (1/w)
    // (forme 0/0), sinon w / (1/z) (forme inf/inf)
    if (ExprPtr r = limiteParQuotientUnique(m_gauche, m_droite, false, a)) return r;
    if (!reecriturePossible()) return limiteNonEvaluee(a);
    GardeProfondeur garde;
    const ExprPtr z = estNul(g) ? m_gauche : m_droite;
    const ExprPtr w = estNul(g) ? m_droite : m_gauche;
    ExprPtr r = (z / (cst(1.0) / w))->limite(a);
    if (estDeterminee(r)) return r;
    r = (w / (cst(1.0) / z))->limite(a);
    if (estDeterminee(r)) return r;
    return limiteNonEvaluee(a);
}

// --- Division ---
ExprPtr Division::primitive() const {
    if (!m_droite->contientVariable()) return m_gauche->integrer() / m_droite;
    if (!m_gauche->contientVariable()) {
        // c / v^n = c * v^(-n), puis règle des puissances
        const Puissance* p = comme<Puissance>(m_droite.get());
        if (p && !p->m_droite->contientVariable()) {
            return m_gauche * ast_pow(p->m_gauche, cst(-1.0) * p->m_droite)->integrer();
        }
        return m_gauche * ast_pow(m_droite, -1.0)->integrer();
    }
    return integraleNonEvaluee();
}
ExprPtr Division::calculerLimite(double a) const {
    double n = 0.0, d = 0.0;
    const ExprPtr Ln = m_gauche->limite(a);
    const ExprPtr Ld = m_droite->limite(a);
    const GenreLimite tn = analyserLimite(Ln, n);
    const GenreLimite td = analyserLimite(Ld, d);
    if (tn == GenreLimite::NonEvaluee || td == GenreLimite::NonEvaluee) {
        if (ExprPtr r = limiteParQuotientUnique(m_gauche, m_droite, true, a)) return r;
        return limiteNonEvaluee(a);
    }
    if (tn == GenreLimite::Symbolique || td == GenreLimite::Symbolique) return Ln / Ld;

    const bool zeroSurZero = estNul(n) && estNul(d);
    const bool infiniSurInfini = std::isinf(n) && std::isinf(d);
    if (zeroSurZero || infiniSurInfini) {
        // Règle de L'Hôpital : lim N/D = lim N'/D'
        if (!reecriturePossible()) return limiteNonEvaluee(a);
        GardeProfondeur garde;
        const ExprPtr r = (m_gauche->derivee()->simplifier() / m_droite->derivee()->simplifier())->limite(a);
        return estDeterminee(r) ? r : limiteNonEvaluee(a);
    }
    if (estNul(d)) {
        // c / 0 : l'infini n'a un signe défini que si D garde le même signe des deux côtés
        int ordre;
        double signeD;
        if (premierTermeNonNul(m_droite, a, ordre, signeD) && ordre % 2 == 0) {
            return cst((n > 0.0 ? 1.0 : -1.0) * signeD * INFINI);
        }
        return limiteNonEvaluee(a); // limites à gauche et à droite différentes
    }
    if (ExprPtr r = nombreOuNul(n / d)) return r;
    return limiteNonEvaluee(a);
}

// --- Puissance ---
ExprPtr Puissance::primitive() const {
    double a;
    // (a*x + b)^n avec n constant : (ax+b)^(n+1) / (a(n+1)), ou ln(ax+b) / a si n = -1
    if (!m_droite->contientVariable() && coefficientLineaire(m_gauche, a)) {
        const ExprPtr n = m_droite->simplifier();
        if (n->estConstante()) {
            const double nv = n->getValeurConstante();
            if (std::abs(nv + 1.0) < 1e-12) return diviserPar(ast_ln(m_gauche), a);
            return ast_pow(m_gauche, cst(nv + 1.0)) / cst(a * (nv + 1.0));
        }
    }
    // b^(a*x + c) avec b constant strictement positif et différent de 1 : b^u / (a ln b)
    if (!m_gauche->contientVariable() && coefficientLineaire(m_droite, a)) {
        const ExprPtr b = m_gauche->simplifier();
        if (b->estConstante() && b->getValeurConstante() > 0.0 && b->getValeurConstante() != 1.0) {
            return ast_pow(m_gauche, m_droite) / cst(a * std::log(b->getValeurConstante()));
        }
    }
    return integraleNonEvaluee();
}
ExprPtr Puissance::calculerLimite(double a) const {
    double b = 0.0, e = 0.0;
    const ExprPtr Lb = m_gauche->limite(a);
    const ExprPtr Le = m_droite->limite(a);
    const GenreLimite tb = analyserLimite(Lb, b);
    const GenreLimite te = analyserLimite(Le, e);

    // u^v = exp(v ln u) : lève les formes 1^inf, 0^0 et inf^0 quand l'exposant dépend de x
    auto parExponentielle = [&]() -> ExprPtr {
        if (!m_droite->contientVariable() || !reecriturePossible()) return nullptr;
        GardeProfondeur garde;
        const ExprPtr r = ast_exp(m_droite * ast_ln(m_gauche))->limite(a);
        return estDeterminee(r) ? r : nullptr;
    };

    if (tb == GenreLimite::NonEvaluee || te == GenreLimite::NonEvaluee) {
        if (ExprPtr r = parExponentielle()) return r;
        return limiteNonEvaluee(a);
    }
    if (tb == GenreLimite::Symbolique || te == GenreLimite::Symbolique) return ast_pow(Lb, Le);

    const bool indeterminee = (b == 1.0 && std::isinf(e)) || (estNul(b) && estNul(e)) ||
                              (std::isinf(b) && estNul(e));
    if (indeterminee && m_droite->contientVariable()) {
        if (ExprPtr r = parExponentielle()) return r;
        return limiteNonEvaluee(a);
    }
    if (estNul(b) && e < 0.0) {
        // 0^(-n) = 1 / 0^n : le signe de l'infini dépend du côté
        if (!reecriturePossible()) return limiteNonEvaluee(a);
        GardeProfondeur garde;
        const ExprPtr r = (cst(1.0) / ast_pow(m_gauche, (cst(-1.0) * m_droite)->simplifier()))->limite(a);
        return estDeterminee(r) ? r : limiteNonEvaluee(a);
    }
    if (ExprPtr r = nombreOuNul(std::pow(b, e))) return r;
    return limiteNonEvaluee(a);
}

// --- Sinus ---
ExprPtr Sinus::primitive() const {
    double a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(cst(-1.0) * ast_cos(m_argument), a);
    return integraleNonEvaluee();
}
ExprPtr Sinus::calculerLimite(double a) const {
    // sin(+-inf) donne NaN : pas de limite
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::sin(u); }, ast_sin)) return r;
    return limiteNonEvaluee(a);
}

// --- Cosinus ---
ExprPtr Cosinus::primitive() const {
    double a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(ast_sin(m_argument), a);
    return integraleNonEvaluee();
}
ExprPtr Cosinus::calculerLimite(double a) const {
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::cos(u); }, ast_cos)) return r;
    return limiteNonEvaluee(a);
}

// ============== EXPONENTIELLE ==================

Exponentielle::Exponentielle(CleFabrique, ExprPtr arg) : FonctionUnaire(TypeNoeud::Exponentielle, arg) {}

double Exponentielle::eval(double x) const { return std::exp(m_argument->eval(x)); }

ExprPtr Exponentielle::derivee() const {
    return ast_exp(m_argument) * m_argument->derivee();
}

ExprPtr Exponentielle::simplifier() const {
    auto arg = m_argument->simplifier();
    if (arg->estConstante()) return cst(std::exp(arg->getValeurConstante()));
    if (const Logarithme* ln_node = comme<Logarithme>(arg)) {
        return ln_node->m_argument->simplifier();
    }
    return ast_exp(arg);
}

void Exponentielle::afficher(std::ostream& os) const {
    os << "exp("; m_argument->afficher(os); os << ")";
}

bool Exponentielle::estEgal(const ASTNode& autre) const {
    const Exponentielle* a = comme<Exponentielle>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

ExprPtr Exponentielle::primitive() const {
    double a;
    if (coefficientLineaire(m_argument, a)) return diviserPar(ast_exp(m_argument), a);
    return integraleNonEvaluee();
}

ExprPtr Exponentielle::calculerLimite(double a) const {
    // exp(+inf) = +inf, exp(-inf) = 0
    if (ExprPtr r = limiteUnaire(m_argument, a, [](double u) { return std::exp(u); }, ast_exp)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr ast_exp(ExprPtr arg) {
    return fabriquer<Exponentielle>(arg);
}

// ============== LOGARITHME ==================

Logarithme::Logarithme(CleFabrique, ExprPtr arg) : FonctionUnaire(TypeNoeud::Logarithme, arg) {}

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

bool Logarithme::estEgal(const ASTNode& autre) const {
    const Logarithme* a = comme<Logarithme>(&autre);
    return a && m_argument->estEgal(*(a->m_argument));
}

ExprPtr Logarithme::primitive() const {
    double a;
    // Primitive de ln(u) : u ln(u) - u
    if (coefficientLineaire(m_argument, a)) {
        return diviserPar(m_argument * ast_ln(m_argument) - m_argument, a);
    }
    return integraleNonEvaluee();
}

ExprPtr Logarithme::calculerLimite(double a) const {
    // ln(u) -> -inf quand u -> 0 ; hors du domaine (u < 0) la limite n'existe pas
    auto lnNum = [](double u) { return estNul(u) ? -INFINI : std::log(u); };
    if (ExprPtr r = limiteUnaire(m_argument, a, lnNum, ast_ln)) return r;
    return limiteNonEvaluee(a);
}

ExprPtr ast_ln(ExprPtr arg) {
    return fabriquer<Logarithme>(arg);
}

// ============== INTEGRALE NON EVALUEE ==================

IntegraleNonEvaluee::IntegraleNonEvaluee(CleFabrique, ExprPtr integrande) : ASTNode(TypeNoeud::IntegraleNonEvaluee), m_integrande(std::move(integrande)) {}

double IntegraleNonEvaluee::eval(double) const {
    throw std::logic_error("Impossible d'evaluer une primitive non calculee symboliquement");
}

// Théorème fondamental de l'analyse : (∫f)' = f
ExprPtr IntegraleNonEvaluee::derivee() const { return m_integrande; }

ExprPtr IntegraleNonEvaluee::simplifier() const {
    return fabriquer<IntegraleNonEvaluee>(m_integrande->simplifier());
}

void IntegraleNonEvaluee::afficher(std::ostream& os) const {
    os << "integrale("; m_integrande->afficher(os); os << ")";
}

bool IntegraleNonEvaluee::estEgal(const ASTNode& autre) const {
    const IntegraleNonEvaluee* i = comme<IntegraleNonEvaluee>(&autre);
    return i && m_integrande->estEgal(*(i->m_integrande));
}

ExprPtr IntegraleNonEvaluee::primitive() const { return integraleNonEvaluee(); }

ExprPtr IntegraleNonEvaluee::calculerLimite(double a) const { return limiteNonEvaluee(a); }

// ============== LIMITE NON EVALUEE ==================

LimiteNonEvaluee::LimiteNonEvaluee(CleFabrique, ExprPtr expression, double point)
    : ASTNode(TypeNoeud::LimiteNonEvaluee), m_expression(std::move(expression)), m_point(point) {}

double LimiteNonEvaluee::eval(double) const {
    throw std::logic_error("Impossible d'evaluer une limite non determinee");
}

ExprPtr LimiteNonEvaluee::derivee() const { return cst(0.0); }

ExprPtr LimiteNonEvaluee::simplifier() const { return clone(); }

void LimiteNonEvaluee::afficher(std::ostream& os) const {
    os << "lim(x->" << m_point << ", "; m_expression->afficher(os); os << ")";
}

bool LimiteNonEvaluee::estEgal(const ASTNode& autre) const {
    const LimiteNonEvaluee* l = comme<LimiteNonEvaluee>(&autre);
    return l && l->m_point == m_point && m_expression->estEgal(*(l->m_expression));
}

ExprPtr LimiteNonEvaluee::primitive() const { return clone() * var("x"); }

ExprPtr LimiteNonEvaluee::calculerLimite(double) const { return clone(); }

} // namespace symalgo
