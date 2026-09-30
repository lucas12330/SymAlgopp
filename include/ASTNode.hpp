/**
 * @file ASTNode.hpp
 * @brief Arbre de syntaxe abstraite (AST) des expressions symboliques de SymAlgo++.
 *
 * Les expressions sont des graphes immuables de noeuds partagés :
 *  - chaque noeud est géré par un ExprPtr (compteur de références intrusif, Ref.hpp) ;
 *  - chaque expression n'existe qu'une fois en mémoire (hash-consing) : deux expressions
 *    identiques sont le même pointeur ;
 *  - les expressions sont construites directement sous forme canonique : sommes et
 *    produits n-aires, termes triés, coefficients rationnels exacts (Nombre). Ainsi
 *    x + x donne 2*x, x * x^2 donne x^3 et 2*(x + 1) donne 2*x + 2.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Nombre.hpp"
#include "Ref.hpp"

namespace symalgo {

class ASTNode;
class TableNoeuds;

/*
 * ExprPtr : pointeur partagé vers un noeud immuable de l'AST.
 */
using ExprPtr = Ref<ASTNode>;

/*
 * Nom : CacheDerivees
 * Description : Dérivées déjà calculées pendant un appel à derivee(), par noeud.
 */
using CacheDerivees = std::unordered_map<const ASTNode*, ExprPtr>;

/*
 * Nom : TypeNoeud
 * Description : Type concret d'un noeud, stocké sur un octet (aiguillage par switch et
 *               conversions typées sans RTTI, voir comme<T>).
 */
enum class TypeNoeud : std::uint8_t {
    Constante,
    Variable,
    Parametre,
    Pi,
    Somme,
    Produit,
    Puissance,
    // Fonctions unaires
    Sinus,
    Cosinus,
    Tangente,
    Exponentielle,
    Logarithme,
    ArcSinus,
    ArcCosinus,
    ArcTangente,
    // Noeuds non évalués
    IntegraleNonEvaluee,
    LimiteNonEvaluee,
};

/*
 * Nom : Terme
 * Description : Terme d'une Somme : coefficient * expression.
 */
struct Terme {
    ExprPtr expression;
    Nombre coefficient;
};

/*
 * Nom : Facteur
 * Description : Facteur d'un Produit : base ^ exposant.
 */
struct Facteur {
    ExprPtr base;
    ExprPtr exposant;
};

/*
 * Nom : Signature
 * Description : Description d'un noeud (type, valeurs, enfants) sans le construire. Sert au
 *               hash-consing : on cherche d'abord un noeud identique existant et l'on
 *               n'alloue un nouveau noeud que s'il n'existe pas encore.
 */
struct Signature {
    TypeNoeud type;
    const ASTNode* enfants[2] = {nullptr, nullptr};
    const Nombre* nombre = nullptr;                 // Constante, constante d'une Somme, coefficient d'un Produit
    double reel = 0.0;                              // point d'une LimiteNonEvaluee
    const std::string* nom = nullptr;               // Variable, Parametre
    const std::vector<Terme>* termes = nullptr;     // Somme
    const std::vector<Facteur>* facteurs = nullptr; // Produit

    static Signature unaire(TypeNoeud type, const ASTNode* a) {
        Signature s{type};
        s.enfants[0] = a;
        return s;
    }
};

std::size_t hashSignature(const Signature& s);
const ASTNode* chercherNoeud(const Signature& s, std::size_t hash);
ExprPtr enregistrerNoeud(ASTNode* nouveau, std::size_t hash);

/*
 * Nom : nombreNoeudsVivants
 * Description : Nombre de noeuds distincts actuellement en mémoire (diagnostic, tests).
 */
std::size_t nombreNoeudsVivants();

/*
 * Nom : CleFabrique
 * Description : Clé de construction (idiome « passkey ») : seuls les helpers de fabrication
 *               peuvent créer des noeuds, toujours gérés par un ExprPtr (jamais sur la pile).
 */
class CleFabrique {
    CleFabrique() = default;
    template <class T, class... Args>
    friend ExprPtr fabriquer(Args&&... args);
};

/*
 * Nom : fabriquer
 * Description : Crée un noeud de type T, ou renvoie le noeud identique déjà existant
 *               (hash-consing). Réservé à l'implémentation : les arguments doivent déjà être
 *               sous forme canonique (utiliser les helpers et opérateurs publics).
 */
template <class T, class... Args>
ExprPtr fabriquer(Args&&... args) {
    const Signature signature = T::signature(static_cast<const Args&>(args)...);
    const std::size_t hash = hashSignature(signature);
    if (const ASTNode* existant = chercherNoeud(signature, hash)) {
        return ExprPtr(const_cast<ASTNode*>(existant));
    }
    return enregistrerNoeud(new T(CleFabrique(), std::forward<Args>(args)...), hash);
}

// ============================================================================
// Classe de base
// ============================================================================

class ASTNode : public ObjetCompte {
public:
    ~ASTNode() override = default;

    /*
     * Nom : type
     * Description : Type concret du noeud.
     */
    TypeNoeud type() const { return m_type; }

    /*
     * Nom : hash
     * Description : Empreinte structurelle du noeud (type, valeurs, empreintes des enfants).
     */
    std::size_t hash() const { return m_hash; }

    /*
     * Nom : eval
     * Description : Évalue l'expression pour une valeur de la variable (lève std::logic_error
     *               si elle contient un paramètre ou un noeud non évalué).
     * Utilisation : double y = expr->eval(2.0);
     */
    virtual double eval(double x) const = 0;

    /*
     * Nom : derivee
     * Description : Dérivée symbolique (sous forme canonique). Les sous-expressions partagées
     *               ne sont dérivées qu'une fois par appel.
     * Utilisation : ExprPtr d = expr->derivee();
     */
    ExprPtr derivee() const;
    ExprPtr derivee(CacheDerivees& cache) const;

    /*
     * Nom : simplifier
     * Description : Simplification (la forme canonique est automatique ; simplifier()
     *               reconstruit l'expression à partir de ses parties simplifiées et applique
     *               les règles des fonctions). Le résultat est mémorisé dans le noeud.
     */
    ExprPtr simplifier() const;

    /*
     * Nom : integrer
     * Description : Primitive symbolique, ou IntegraleNonEvaluee si aucune règle ne s'applique.
     */
    ExprPtr integrer() const;

    /*
     * Nom : limite
     * Description : Limite quand x tend vers a, ou LimiteNonEvaluee si elle n'est pas déterminée.
     */
    ExprPtr limite(double a) const;

    /*
     * Nom : DL
     * Description : Développement limité de Taylor en a à l'ordre donné (arithmétique des
     *               séries tronquées). Lève std::domain_error si la fonction n'est pas
     *               développable en a.
     */
    ExprPtr DL(double a, int ordre) const;

    /*
     * Nom : afficher / texte
     * Description : Écriture lisible (« x^2 + 5*x + 6 », « sin(x)/x »).
     */
    void afficher(std::ostream& os) const;
    std::string texte() const;

    /*
     * Nom : clone
     * Description : Nouvelle référence vers ce noeud (l'arbre est immuable et partagé).
     */
    ExprPtr clone() const { return ExprPtr(const_cast<ASTNode*>(this)); }

    /*
     * Nom : estEgal
     * Description : Égalité mathématique structurelle. Grâce à la forme canonique et au
     *               hash-consing, c'est une simple comparaison d'adresses.
     */
    bool estEgal(const ASTNode& autre) const { return this == &autre; }

    /*
     * Nom : contientVariable
     * Description : Indique si l'expression dépend de la variable d'évaluation.
     */
    bool contientVariable() const { return m_contientVariable; }

    /*
     * Nom : estConstante / getValeurConstante
     * Description : Vrai pour un noeud Constante ; sa valeur approchée en double.
     */
    bool estConstante() const { return m_type == TypeNoeud::Constante; }
    double getValeurConstante() const;

protected:
    ASTNode(TypeNoeud type, bool contientVariable) : m_type(type), m_contientVariable(contientVariable) {}

    virtual ExprPtr calculerDerivee(CacheDerivees& cache) const = 0;
    virtual ExprPtr calculerSimplification() const = 0;
    virtual ExprPtr primitive() const = 0;
    virtual ExprPtr calculerLimite(double a) const = 0;

    ExprPtr integraleNonEvaluee() const;
    ExprPtr limiteNonEvaluee(double a) const;

    /*
     * Nom : detruire
     * Description : Retire le noeud de la table de hash-consing avant de le libérer.
     */
    void detruire() const override;

private:
    TypeNoeud m_type;
    bool m_contientVariable;
    std::size_t m_hash = 0;
    // Mémo de simplifier() : m_estSimplifie si le noeud est sa propre forme simplifiée,
    // sinon m_formeSimplifiee une fois calculée (elle ne peut pas contenir ce noeud)
    mutable bool m_estSimplifie = false;
    mutable ExprPtr m_formeSimplifiee;
    // Chaînage intrusif dans la table de hash-consing
    mutable const ASTNode* m_suivantTable = nullptr;

    friend ExprPtr enregistrerNoeud(ASTNode* nouveau, std::size_t hash);
    friend class TableNoeuds;
};

/*
 * Nom : comme
 * Description : Accès typé à un noeud : le noeud converti en T, ou nullptr s'il est d'un
 *               autre type.
 * Utilisation : if (const Sinus* s = comme<Sinus>(expr)) { ... s->m_argument ... }
 */
template <class T>
const T* comme(const ASTNode* n) {
    return n && T::correspond(n->type()) ? static_cast<const T*>(n) : nullptr;
}

template <class T>
const T* comme(const ExprPtr& e) { return comme<T>(e.get()); }

// Interdit sur un temporaire : le pointeur renvoyé survivrait à l'expression
// (ex. comme<Constante>(u->derivee()) pointerait vers un noeud déjà libéré)
template <class T>
const T* comme(ExprPtr&& e) = delete;

/*
 * Nom : comparer
 * Description : Ordre total et déterministe sur les expressions (tri des termes et des
 *               facteurs de la forme canonique). Renvoie <0, 0 ou >0.
 */
int comparer(const ASTNode& a, const ASTNode& b);

// ============================================================================
// Noeuds terminaux
// ============================================================================

/*
 * CLASSE CONSTANTE : nombre exact (entier ou rationnel de taille arbitraire) ou réel.
 */
class Constante : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Constante;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const Nombre& valeur) {
        Signature s{TYPE};
        s.nombre = &valeur;
        return s;
    }

    Constante(CleFabrique, Nombre valeur);

    double eval(double x) const override;
    const Nombre& getNombre() const { return m_valeur; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    Nombre m_valeur;
    double m_approx; // valeur en double, pour l'évaluation
};

/*
 * CLASSE VARIABLE : la variable d'évaluation (quel que soit son nom).
 */
class Variable : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Variable;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const std::string& nom) {
        Signature s{TYPE};
        s.nom = &nom;
        return s;
    }

    Variable(CleFabrique, const std::string& nom);

    double eval(double x) const override;
    const std::string& getNom() const { return m_nom; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    std::string m_nom;
};

/*
 * CLASSE PARAMETRE : constante symbolique sans valeur (ex. C1, C2 des solutions d'EDO) :
 * dérivée nulle, évaluation impossible (std::logic_error).
 */
class Parametre : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Parametre;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const std::string& nom) {
        Signature s{TYPE};
        s.nom = &nom;
        return s;
    }

    Parametre(CleFabrique, const std::string& nom);

    double eval(double x) const override;
    const std::string& getNom() const { return m_nom; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    std::string m_nom;
};

/*
 * CLASSE PI : la constante exacte pi. Les fonctions trigonométriques en donnent les valeurs
 * remarquables exactes (sin(pi/6) = 1/2, cos(3*pi/4) = -2^(1/2)/2...).
 */
class Pi : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Pi;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature() { return Signature{TYPE}; }

    explicit Pi(CleFabrique);

    double eval(double x) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

// ============================================================================
// Sommes, produits, puissances (forme canonique)
// ============================================================================

/*
 * CLASSE SOMME : constante + sum c_i * t_i, avec au moins deux éléments. Les termes t_i
 * sont distincts, triés, ni constants ni sommes, et leurs coefficients c_i non nuls.
 */
class Somme : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Somme;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const Nombre& constante, const std::vector<Terme>& termes) {
        Signature s{TYPE};
        s.nombre = &constante;
        s.termes = &termes;
        return s;
    }

    Somme(CleFabrique, Nombre constante, std::vector<Terme> termes);

    double eval(double x) const override;
    const Nombre& getConstante() const { return m_constante; }
    const std::vector<Terme>& getTermes() const { return m_termes; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    Nombre m_constante;
    std::vector<Terme> m_termes;
};

/*
 * CLASSE PRODUIT : coefficient * prod b_i ^ e_i. Les bases b_i sont distinctes, triées,
 * ni produits ni (sauf exposant non entier) constantes ; les exposants sont non nuls.
 */
class Produit : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Produit;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const Nombre& coefficient, const std::vector<Facteur>& facteurs) {
        Signature s{TYPE};
        s.nombre = &coefficient;
        s.facteurs = &facteurs;
        return s;
    }

    Produit(CleFabrique, Nombre coefficient, std::vector<Facteur> facteurs);

    double eval(double x) const override;
    const Nombre& getCoefficient() const { return m_coefficient; }
    const std::vector<Facteur>& getFacteurs() const { return m_facteurs; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    Nombre m_coefficient;
    std::vector<Facteur> m_facteurs;
};

/*
 * CLASSE PUISSANCE : base ^ exposant (un seul facteur ; exposant ni 0 ni 1).
 */
class Puissance : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Puissance;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const ExprPtr& base, const ExprPtr& exposant) {
        Signature s{TYPE};
        s.enfants[0] = base.get();
        s.enfants[1] = exposant.get();
        return s;
    }

    Puissance(CleFabrique, ExprPtr base, ExprPtr exposant);

    double eval(double x) const override;
    const ExprPtr& getBase() const { return m_base; }
    const ExprPtr& getExposant() const { return m_exposant; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    ExprPtr m_base;
    ExprPtr m_exposant;
    long long m_exposantEntier; // pour l'évaluation rapide, ou PAS_ENTIER
};

// ============================================================================
// Fonctions unaires
// ============================================================================

class FonctionUnaire : public ASTNode {
public:
    static constexpr bool correspond(TypeNoeud t) {
        return t >= TypeNoeud::Sinus && t <= TypeNoeud::ArcTangente;
    }

    ExprPtr m_argument;

protected:
    FonctionUnaire(TypeNoeud type, ExprPtr arg);
};

#define SYMALGO_FONCTION_UNAIRE(NOM)                                                           \
    class NOM : public FonctionUnaire {                                                        \
    public:                                                                                    \
        static constexpr TypeNoeud TYPE = TypeNoeud::NOM;                                      \
        static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }                    \
        static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); } \
        NOM(CleFabrique, ExprPtr arg) : FonctionUnaire(TYPE, std::move(arg)) {}                 \
        double eval(double x) const override;                                                  \
                                                                                               \
    protected:                                                                                 \
        ExprPtr calculerDerivee(CacheDerivees& cache) const override;                          \
        ExprPtr calculerSimplification() const override;                                       \
        ExprPtr primitive() const override;                                                    \
        ExprPtr calculerLimite(double a) const override;                                       \
    };

SYMALGO_FONCTION_UNAIRE(Sinus)
SYMALGO_FONCTION_UNAIRE(Cosinus)
SYMALGO_FONCTION_UNAIRE(Tangente)
SYMALGO_FONCTION_UNAIRE(Exponentielle)
SYMALGO_FONCTION_UNAIRE(Logarithme)
SYMALGO_FONCTION_UNAIRE(ArcSinus)
SYMALGO_FONCTION_UNAIRE(ArcCosinus)
SYMALGO_FONCTION_UNAIRE(ArcTangente)

#undef SYMALGO_FONCTION_UNAIRE

// ============================================================================
// Noeuds non évalués
// ============================================================================

/*
 * CLASSE INTEGRALENONEVALUEE : primitive qu'aucune règle ne sait calculer (au lieu d'un
 * résultat faux). Sa dérivée redonne l'intégrande ; l'évaluer lève std::logic_error.
 */
class IntegraleNonEvaluee : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::IntegraleNonEvaluee;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const ExprPtr& integrande) { return Signature::unaire(TYPE, integrande.get()); }

    IntegraleNonEvaluee(CleFabrique, ExprPtr integrande);

    double eval(double x) const override;
    const ExprPtr& getIntegrande() const { return m_integrande; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    ExprPtr m_integrande;
};

/*
 * CLASSE LIMITENONEVALUEE : limite non déterminée (forme indéterminée non résolue, ou
 * limite inexistante comme 1/x en 0). C'est un nombre (s'il existe) : dérivée nulle ;
 * l'évaluer lève std::logic_error.
 */
class LimiteNonEvaluee : public ASTNode {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::LimiteNonEvaluee;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }
    static Signature signature(const ExprPtr& expression, double point) {
        Signature s = Signature::unaire(TYPE, expression.get());
        s.reel = point;
        return s;
    }

    LimiteNonEvaluee(CleFabrique, ExprPtr expression, double point);

    double eval(double x) const override;
    const ExprPtr& getExpression() const { return m_expression; }
    double getPoint() const { return m_point; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;

private:
    ExprPtr m_expression;
    double m_point;
};

// ============================================================================
// Construction des expressions (toujours sous forme canonique)
// ============================================================================

/*
 * Nom : cst / frac / nombre
 * Description : Constantes. cst(2.0) est l'entier exact 2, cst(0.5) le réel 0.5 ;
 *               frac(1, 3) le rationnel exact 1/3 ; nombre(n) une Constante quelconque.
 */
ExprPtr cst(double valeur);
ExprPtr frac(std::int64_t num, std::int64_t den = 1);
ExprPtr nombre(const Nombre& valeur);

/*
 * Nom : var / param
 * Description : Variable d'évaluation et constante symbolique (paramètre).
 */
ExprPtr var(const std::string& nom = "x");
ExprPtr param(const std::string& nom);

/*
 * Nom : pi
 * Description : La constante exacte pi.
 */
ExprPtr pi();

/*
 * Nom : somme / produit
 * Description : Somme et produit canoniques d'une liste d'expressions.
 */
ExprPtr somme(const std::vector<ExprPtr>& termes);
ExprPtr produit(const std::vector<ExprPtr>& facteurs);

/*
 * Opérateurs : a - b est construit comme a + (-1)*b et a / b comme a * b^(-1).
 */
ExprPtr operator+(const ExprPtr& gauche, const ExprPtr& droite);
ExprPtr operator+(const ExprPtr& gauche, double droite);
ExprPtr operator+(double gauche, const ExprPtr& droite);
ExprPtr operator-(const ExprPtr& gauche, const ExprPtr& droite);
ExprPtr operator-(const ExprPtr& gauche, double droite);
ExprPtr operator-(double gauche, const ExprPtr& droite);
ExprPtr operator-(const ExprPtr& e);
ExprPtr operator*(const ExprPtr& gauche, const ExprPtr& droite);
ExprPtr operator*(const ExprPtr& gauche, double droite);
ExprPtr operator*(double gauche, const ExprPtr& droite);
ExprPtr operator/(const ExprPtr& gauche, const ExprPtr& droite);
ExprPtr operator/(const ExprPtr& gauche, double droite);
ExprPtr operator/(double gauche, const ExprPtr& droite);

/*
 * Nom : ast_pow
 * Description : Puissance canonique (évalue exactement 2^10 = 1024, 4^(1/2) = 2, (x^2)^3 = x^6).
 */
ExprPtr ast_pow(const ExprPtr& base, const ExprPtr& exposant);
ExprPtr ast_pow(const ExprPtr& base, double exposant);
ExprPtr ast_pow(double base, const ExprPtr& exposant);

/*
 * Nom : ast_sin / ast_cos / ast_tan / ast_exp / ast_ln
 * Description : Fonctions usuelles (valeurs exactes remarquables : sin(0) = 0, exp(0) = 1,
 *               ln(1) = 0, sin(pi/6) = 1/2, tan(pi/3) = 3^(1/2) ; exp(ln(u)) = u,
 *               ln(exp(u)) = u ; argument réel évalué).
 */
ExprPtr ast_sin(const ExprPtr& arg);
ExprPtr ast_cos(const ExprPtr& arg);
ExprPtr ast_tan(const ExprPtr& arg);
ExprPtr ast_exp(const ExprPtr& arg);
ExprPtr ast_ln(const ExprPtr& arg);

/*
 * Nom : ast_asin / ast_acos / ast_atan
 * Description : Fonctions réciproques (valeurs exactes remarquables : asin(1/2) = pi/6,
 *               acos(0) = pi/2, atan(1) = pi/4...).
 */
ExprPtr ast_asin(const ExprPtr& arg);
ExprPtr ast_acos(const ExprPtr& arg);
ExprPtr ast_atan(const ExprPtr& arg);

std::ostream& operator<<(std::ostream& os, const ExprPtr& e);

} // namespace symalgo
