#pragma once
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "Ref.hpp"

namespace symalgo {

class ASTNode;
class TableNoeuds;

/*
 * Nom : TypeNoeud
 * Description : Type concret d'un noeud, stocké sur un octet. Permet l'aiguillage par
 *               switch et les conversions typées sans RTTI (voir comme<T>). Les opérateurs
 *               binaires et les fonctions unaires occupent des plages contiguës.
 */
enum class TypeNoeud : std::uint8_t {
    Constante,
    Fraction,
    Variable,
    Parametre,
    // Opérateurs binaires
    Addition,
    Soustraction,
    Multiplication,
    Division,
    Puissance,
    // Fonctions unaires
    Sinus,
    Cosinus,
    Tangente,
    Exponentielle,
    Logarithme,
    // Noeuds non évalués
    IntegraleNonEvaluee,
    LimiteNonEvaluee,
};

/*
 * ExprPtr : pointeur partagé vers un noeud immuable de l'AST (compteur intrusif, voir Ref.hpp).
 */
using ExprPtr = Ref<ASTNode>;

/*
 * Nom : CacheDerivees
 * Description : Dérivées déjà calculées pendant un appel à derivee(), par noeud.
 */
using CacheDerivees = std::unordered_map<const ASTNode*, ExprPtr>;

/*
 * Nom : CleFabrique
 * Description : Clé de construction (idiome « passkey ») : seuls les helpers de fabrication
 *               peuvent créer des noeuds, qui sont donc toujours gérés par un ExprPtr.
 *               Un noeud ne peut pas être créé sur la pile.
 */
class CleFabrique {
    CleFabrique() = default;
    template <class T, class... Args>
    friend ExprPtr fabriquer(Args&&... args);
};

/*
 * Nom : Signature
 * Description : Description d'un noeud (type, valeurs, enfants) sans le construire. Sert au
 *               hash-consing : on cherche d'abord un noeud identique existant, et l'on
 *               n'alloue un nouveau noeud que s'il n'existe pas encore.
 */
struct Signature {
    TypeNoeud type;
    const ASTNode* enfants[2] = {nullptr, nullptr};
    double reel = 0.0;                // Constante, point d'une LimiteNonEvaluee
    std::int64_t num = 0, den = 1;    // Fraction (normalisée)
    const std::string* nom = nullptr; // Variable, Parametre

    static Signature unaire(TypeNoeud type, const ASTNode* a) {
        Signature s{type};
        s.enfants[0] = a;
        return s;
    }
    static Signature binaire(TypeNoeud type, const ASTNode* a, const ASTNode* b) {
        Signature s{type};
        s.enfants[0] = a;
        s.enfants[1] = b;
        return s;
    }
};

/*
 * Nom : hashSignature
 * Description : Empreinte d'une signature (les enfants contribuent par leur propre empreinte).
 */
std::size_t hashSignature(const Signature& s);

/*
 * Nom : chercherNoeud / enregistrerNoeud
 * Description : Table de hash-consing : recherche d'un noeud existant de même signature, et
 *               enregistrement d'un noeud nouvellement créé.
 */
const ASTNode* chercherNoeud(const Signature& s, std::size_t hash);
ExprPtr enregistrerNoeud(ASTNode* nouveau, std::size_t hash);

/*
 * Nom : nombreNoeudsVivants
 * Description : Nombre de noeuds distincts actuellement en mémoire (diagnostic, tests).
 */
std::size_t nombreNoeudsVivants();

/*
 * Nom : fabriquer
 * Description : Crée un noeud de type T, ou renvoie le noeud identique déjà existant
 *               (hash-consing) : chaque expression n'existe qu'une fois en mémoire et deux
 *               expressions identiques sont le même pointeur.
 * Utilisation : ExprPtr s = fabriquer<Sinus>(argument);
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

// Classe abstraite de base pour tous les noeuds de l'arbre
class ASTNode : public ObjetCompte {
public:
    /*
     * Nom : type
     * Description : Type concret du noeud.
     * Utilisation : if (noeud->type() == TypeNoeud::Sinus) { ... }
     */
    TypeNoeud type() const { return m_type; }

    /*
     * Nom : hash
     * Description : Empreinte structurelle du noeud (type, valeurs, empreintes des enfants).
     */
    std::size_t hash() const { return m_hash; }

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
     * Description : Calcule la dérivée symbolique de l'expression. Les sous-expressions
     *               partagées (fréquentes grâce au hash-consing) ne sont dérivées qu'une
     *               fois par appel.
     * Utilisation : ExprPtr d = noeud->derivee();
     */
    ExprPtr derivee() const;

    /*
     * Nom : derivee (avec cache)
     * Description : Variante utilisée par les règles de dérivation : réutilise les dérivées
     *               déjà calculées pendant l'appel en cours.
     */
    ExprPtr derivee(CacheDerivees& cache) const;

    /*
     * Nom : simplifier
     * Description : Simplifie mathématiquement l'expression. Le résultat est mémorisé dans le
     *               noeud : simplifier deux fois la même expression (ou une expression déjà
     *               simplifiée) est immédiat.
     * Utilisation : ExprPtr simp = noeud->simplifier();
     */
    ExprPtr simplifier() const;

    /*
     * Nom : afficher
     * Description : Affiche le contenu textuel de l'expression sur un flux.
     * Utilisation : noeud->afficher(std::cout);
     */
    virtual void afficher(std::ostream& os) const = 0;

    /*
     * Nom : clone
     * Description : Renvoie un pointeur partagé vers ce noeud : l'arbre est immuable, les
     *               sous-arbres sont donc partagés sans jamais être copiés.
     * Utilisation : ExprPtr copie = noeud->clone();
     */
    ExprPtr clone() const { return ExprPtr(const_cast<ASTNode*>(this)); }

    /*
     * Nom : integrer
     * Description : Calcule une primitive symbolique de l'expression. Une expression
     *               indépendante de x s'intègre en c*x ; sinon le calcul est délégué au
     *               noeud (primitive()). Si aucune règle ne s'applique, renvoie un noeud
     *               IntegraleNonEvaluee plutôt qu'un résultat faux.
     * Utilisation : ExprPtr p = noeud->integrer();
     */
    ExprPtr integrer() const;

    /*
     * Nom : limite
     * Description : Calcule la limite quand x tend vers a. Si elle ne peut pas être
     *               déterminée, renvoie un noeud LimiteNonEvaluee.
     * Utilisation : ExprPtr l = noeud->limite(a);
     */
    ExprPtr limite(double a) const;

    /*
     * Nom : contientVariable
     * Description : Indique si l'expression dépend de la variable d'évaluation x.
     *               Faux pour les constantes numériques et les paramètres symboliques.
     * Utilisation : bool dep = noeud->contientVariable();
     */
    virtual bool contientVariable() const { return false; }

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

protected:
    explicit ASTNode(TypeNoeud type) : m_type(type) {}

    /*
     * Nom : calculerDerivee
     * Description : Règle de dérivation propre au noeud ; les enfants se dérivent par
     *               enfant->derivee(cache).
     */
    virtual ExprPtr calculerDerivee(CacheDerivees& cache) const = 0;

    /*
     * Nom : calculerSimplification
     * Description : Règles de simplification propres au noeud, appelées par simplifier()
     *               lorsque le résultat n'est pas déjà connu.
     */
    virtual ExprPtr calculerSimplification() const = 0;

    /*
     * Nom : primitive
     * Description : Règles d'intégration propres au noeud, appelées par integrer()
     *               lorsque l'expression dépend de x.
     */
    virtual ExprPtr primitive() const = 0;

    /*
     * Nom : calculerLimite
     * Description : Règles de calcul de limite propres au noeud, appelées par limite().
     */
    virtual ExprPtr calculerLimite(double a) const = 0;

    /*
     * Nom : integraleNonEvaluee
     * Description : Renvoie le noeud IntegraleNonEvaluee portant sur ce noeud.
     */
    ExprPtr integraleNonEvaluee() const;

    /*
     * Nom : limiteNonEvaluee
     * Description : Renvoie le noeud LimiteNonEvaluee portant sur ce noeud au point a.
     */
    ExprPtr limiteNonEvaluee(double a) const;

    /*
     * Nom : detruire
     * Description : Retire le noeud de la table de hash-consing avant de le libérer.
     */
    void detruire() const override;

private:
    TypeNoeud m_type;
    std::size_t m_hash = 0;
    // Mémo de simplifier() : m_estSimplifie si le noeud est sa propre forme simplifiée,
    // sinon m_formeSimplifiee une fois calculée (elle ne peut pas contenir ce noeud)
    mutable bool m_estSimplifie = false;
    mutable ExprPtr m_formeSimplifiee;
    // Chaînage intrusif dans la table de hash-consing (aucune allocation par insertion)
    mutable const ASTNode* m_suivantTable = nullptr;

    friend ExprPtr enregistrerNoeud(ASTNode* nouveau, std::size_t hash);
    friend class TableNoeuds;
};

/*
 * Nom : comme
 * Description : Accès typé à un noeud : renvoie le noeud converti en T, ou nullptr s'il est
 *               d'un autre type (équivalent de std::dynamic_pointer_cast).
 * Utilisation : if (const Sinus* s = comme<Sinus>(expr)) { ... s->m_argument ... }
 */
template <class T>
const T* comme(const ASTNode* n) {
    return n && T::correspond(n->type()) ? static_cast<const T*>(n) : nullptr;
}

template <class T>
const T* comme(const ExprPtr& e) { return comme<T>(e.get()); }

// --- Noeuds Terminaux ---

class Constante : public ASTNode {
    double m_valeur;
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Constante;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Constante
     * Description : Constructeur initialisant la constante avec sa valeur numérique.
     * Utilisation : ExprPtr c = cst(5.0);
     */
    Constante(CleFabrique, double valeur);
    static Signature signature(double valeur) {
        Signature s{TYPE};
        s.reel = valeur;
        return s;
    }

    /*
     * Nom : eval
     * Description : Renvoie toujours la valeur de la constante, indépendamment de x.
     * Utilisation : double val = c.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche la valeur de la constante.
     * Utilisation : c.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est une constante de même valeur.
     * Utilisation : bool eq = c.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

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

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

/*
 * ============================================================================
 * CLASSE FRACTION
 * Représente un nombre rationnel exact A/B
 * ============================================================================
 */
class Fraction : public ASTNode {
private:
    int64_t m_num;
    int64_t m_den;
    double m_valeur_eval; // Cache pour eval()

public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Fraction;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    Fraction(CleFabrique, int64_t num, int64_t den);
    static Signature signature(int64_t num, int64_t den);

    double eval(double x) const override;
    void afficher(std::ostream& os) const override;
    bool estEgal(const ASTNode& autre) const override;

    bool estConstante() const override { return true; }
    double getValeurConstante() const override { return m_valeur_eval; }
    int64_t getNum() const { return m_num; }
    int64_t getDen() const { return m_den; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Variable : public ASTNode {
    std::string m_nom;
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Variable;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Variable
     * Description : Constructeur d'une variable mathématique avec un nom (par défaut "x").
     * Utilisation : ExprPtr v = var("y");
     */
    Variable(CleFabrique, const std::string& nom);
    static Signature signature(const std::string& nom) {
        Signature s{TYPE};
        s.nom = &nom;
        return s;
    }

    /*
     * Nom : eval
     * Description : Renvoie la valeur x passée en paramètre (valeur de la variable).
     * Utilisation : double val = v.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche le nom textuel de la variable.
     * Utilisation : v.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si l'autre noeud est une variable avec le même nom.
     * Utilisation : bool eq = v.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    bool contientVariable() const override { return true; }

    const std::string& getNom() const { return m_nom; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

/*
 * ============================================================================
 * CLASSE PARAMETRE
 * Constante symbolique sans valeur numérique (ex : les constantes C1, C2...
 * des solutions générales d'EDO). Contrairement à Variable, elle ne dépend
 * pas de x : sa dérivée est nulle. L'évaluer lève std::logic_error.
 * ============================================================================
 */
class Parametre : public ASTNode {
    std::string m_nom;
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Parametre;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Parametre
     * Description : Construit une constante symbolique nommée.
     * Utilisation : ExprPtr c = param("C1");
     */
    Parametre(CleFabrique, const std::string& nom);
    static Signature signature(const std::string& nom) {
        Signature s{TYPE};
        s.nom = &nom;
        return s;
    }

    /*
     * Nom : eval
     * Description : Un paramètre n'a pas de valeur numérique : lève std::logic_error.
     * Utilisation : Appelé indirectement par l'évaluation d'une expression.
     */
    double eval(double x) const override;

void afficher(std::ostream& os) const override;

    /*
     * Nom : estEgal
     * Description : Vérifie si l'autre noeud est un paramètre de même nom.
     * Utilisation : bool eq = p->estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

    const std::string& getNom() const { return m_nom; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

// --- Opérateurs Binaires ---

class OperateurBinaire : public ASTNode {
public:
    static constexpr bool correspond(TypeNoeud t) {
        return t >= TypeNoeud::Addition && t <= TypeNoeud::Puissance;
    }

    ExprPtr m_gauche;
    ExprPtr m_droite;

    /*
     * Nom : OperateurBinaire
     * Description : Constructeur de base pour tous les opérateurs prenant deux opérandes.
     * Utilisation : Appelé par les constructeurs des classes filles.
     */
    OperateurBinaire(TypeNoeud type, ExprPtr gauche, ExprPtr droite);

    bool contientVariable() const override {
        return m_gauche->contientVariable() || m_droite->contientVariable();
    }
};

class Addition : public OperateurBinaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Addition;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Addition
     * Description : Construit un noeud d'addition de deux expressions.
     * Utilisation : ExprPtr add = expr1 + expr2;
     */
    Addition(CleFabrique, ExprPtr gauche, ExprPtr droite);
    static Signature signature(const ExprPtr& g, const ExprPtr& d) { return Signature::binaire(TYPE, g.get(), d.get()); }

    /*
     * Nom : eval
     * Description : Evalue la somme des deux sous-arbres pour une valeur de x.
     * Utilisation : double val = add.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche l'addition au format (gauche + droite).
     * Utilisation : add.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie l'égalité commutative (a+b = a+b ou a+b = b+a).
     * Utilisation : bool eq = add.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Soustraction : public OperateurBinaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Soustraction;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Soustraction
     * Description : Construit un noeud de soustraction.
     * Utilisation : ExprPtr sub = expr1 - expr2;
     */
    Soustraction(CleFabrique, ExprPtr gauche, ExprPtr droite);
    static Signature signature(const ExprPtr& g, const ExprPtr& d) { return Signature::binaire(TYPE, g.get(), d.get()); }

    /*
     * Nom : eval
     * Description : Evalue la différence gauche - droite.
     * Utilisation : double val = sub.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche la soustraction au format (gauche - droite).
     * Utilisation : sub.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie que le noeud actuel est une soustraction de termes identiques.
     * Utilisation : bool eq = sub.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Multiplication : public OperateurBinaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Multiplication;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Multiplication
     * Description : Construit un noeud de multiplication de deux expressions.
     * Utilisation : ExprPtr mul = expr1 * expr2;
     */
    Multiplication(CleFabrique, ExprPtr gauche, ExprPtr droite);
    static Signature signature(const ExprPtr& g, const ExprPtr& d) { return Signature::binaire(TYPE, g.get(), d.get()); }

    /*
     * Nom : eval
     * Description : Evalue le produit des deux sous-arbres pour x.
     * Utilisation : double val = mul.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche le produit au format gauche * droite sans parenthèses.
     * Utilisation : mul.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie l'égalité commutative pour la multiplication.
     * Utilisation : bool eq = mul.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Division : public OperateurBinaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Division;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Division
     * Description : Construit un noeud de division de deux expressions.
     * Utilisation : ExprPtr div = expr1 / expr2;
     */
    Division(CleFabrique, ExprPtr gauche, ExprPtr droite);
    static Signature signature(const ExprPtr& g, const ExprPtr& d) { return Signature::binaire(TYPE, g.get(), d.get()); }

    /*
     * Nom : eval
     * Description : Evalue le quotient (gauche / droite).
     * Utilisation : double val = div.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche la division au format (gauche / droite).
     * Utilisation : div.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est la même division exacte.
     * Utilisation : bool eq = div.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Puissance : public OperateurBinaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Puissance;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Puissance
     * Description : Construit un noeud de puissance (base^exposant).
     * Utilisation : ExprPtr p = ast_pow(base, exposant);
     */
    Puissance(CleFabrique, ExprPtr base, ExprPtr exposant);
    static Signature signature(const ExprPtr& g, const ExprPtr& d) { return Signature::binaire(TYPE, g.get(), d.get()); }

    /*
     * Nom : eval
     * Description : Evalue (gauche) élevé à la puissance (droite).
     * Utilisation : double val = p.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche la puissance au format (base)^(exposant).
     * Utilisation : p.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est une puissance avec mêmes base et exposant.
     * Utilisation : bool eq = p.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

// --- Fonctions Unaires ---

class FonctionUnaire : public ASTNode {
public:
    static constexpr bool correspond(TypeNoeud t) {
        return t >= TypeNoeud::Sinus && t <= TypeNoeud::Logarithme;
    }

    ExprPtr m_argument;

    /*
     * Nom : FonctionUnaire
     * Description : Constructeur de base pour les fonctions mathématiques à un paramètre (ex: sin, cos).
     * Utilisation : Appelé par les constructeurs de Sinus, Cosinus, etc.
     */
    FonctionUnaire(TypeNoeud type, ExprPtr arg);

    bool contientVariable() const override { return m_argument->contientVariable(); }
};

class Sinus : public FonctionUnaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Sinus;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Sinus
     * Description : Construit un noeud pour la fonction sinus.
     * Utilisation : ExprPtr s = ast_sin(expr);
     */
    Sinus(CleFabrique, ExprPtr arg);
    static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); }

    /*
     * Nom : eval
     * Description : Evalue sin(argument).
     * Utilisation : double val = s.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche sous la forme sin(argument).
     * Utilisation : s.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est un sinus avec le même argument.
     * Utilisation : bool eq = s.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Cosinus : public FonctionUnaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Cosinus;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    /*
     * Nom : Cosinus
     * Description : Construit un noeud pour la fonction cosinus.
     * Utilisation : ExprPtr c = ast_cos(expr);
     */
    Cosinus(CleFabrique, ExprPtr arg);
    static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); }

    /*
     * Nom : eval
     * Description : Evalue cos(argument).
     * Utilisation : double val = c.eval(x);
     */
    double eval(double x) const override;

/*
     * Nom : afficher
     * Description : Affiche sous la forme cos(argument).
     * Utilisation : c.afficher(std::cout);
     */
    void afficher(std::ostream& os) const override;


    /*
     * Nom : estEgal
     * Description : Vérifie si un autre noeud est un cosinus avec le même argument.
     * Utilisation : bool eq = c.estEgal(autre);
     */
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Tangente : public FonctionUnaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Tangente;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    Tangente(CleFabrique, ExprPtr arg);
    static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); }
    double eval(double x) const override;
    void afficher(std::ostream& os) const override;
    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Exponentielle : public FonctionUnaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Exponentielle;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    Exponentielle(CleFabrique, ExprPtr arg);
    static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); }
    double eval(double x) const override;
    void afficher(std::ostream& os) const override;

    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

class Logarithme : public FonctionUnaire {
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::Logarithme;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    Logarithme(CleFabrique, ExprPtr arg);
    static Signature signature(const ExprPtr& a) { return Signature::unaire(TYPE, a.get()); }
    double eval(double x) const override;
    void afficher(std::ostream& os) const override;

    bool estEgal(const ASTNode& autre) const override;

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

// --- Noeuds non évalués ---

/*
 * ============================================================================
 * CLASSE INTEGRALENONEVALUEE
 * Représente la primitive d'une expression qu'aucune règle ne sait intégrer,
 * au lieu de renvoyer un résultat faux. Sa dérivée redonne l'intégrande
 * (théorème fondamental de l'analyse) ; l'évaluer lève std::logic_error.
 * ============================================================================
 */
class IntegraleNonEvaluee : public ASTNode {
    ExprPtr m_integrande;
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::IntegraleNonEvaluee;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    IntegraleNonEvaluee(CleFabrique, ExprPtr integrande);
    static Signature signature(const ExprPtr& integrande) { return Signature::unaire(TYPE, integrande.get()); }

    double eval(double x) const override;
    void afficher(std::ostream& os) const override;
    bool estEgal(const ASTNode& autre) const override;
    bool contientVariable() const override { return true; }

    const ExprPtr& getIntegrande() const { return m_integrande; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

/*
 * ============================================================================
 * CLASSE LIMITENONEVALUEE
 * Représente une limite qui n'a pas pu être déterminée (forme indéterminée
 * non résolue, ou limite inexistante comme 1/x en 0). C'est un nombre (s'il
 * existe) : sa dérivée est nulle ; l'évaluer lève std::logic_error.
 * ============================================================================
 */
class LimiteNonEvaluee : public ASTNode {
    ExprPtr m_expression;
    double m_point;
public:
    static constexpr TypeNoeud TYPE = TypeNoeud::LimiteNonEvaluee;
    static constexpr bool correspond(TypeNoeud t) { return t == TYPE; }

    LimiteNonEvaluee(CleFabrique, ExprPtr expression, double point);
    static Signature signature(const ExprPtr& expression, double point) {
        Signature s = Signature::unaire(TYPE, expression.get());
        s.reel = point;
        return s;
    }

    double eval(double x) const override;
    void afficher(std::ostream& os) const override;
    bool estEgal(const ASTNode& autre) const override;

    const ExprPtr& getExpression() const { return m_expression; }
    double getPoint() const { return m_point; }

protected:
    ExprPtr calculerDerivee(CacheDerivees& cache) const override;
    ExprPtr calculerSimplification() const override;
    ExprPtr primitive() const override;
    ExprPtr calculerLimite(double a) const override;
};

// --- Fonctions Helpers & Surcharge d'Opérateurs ---

/*
 * Nom : cst
 * Description : Helper pour créer rapidement un noeud constante.
 * Utilisation : ExprPtr c = cst(5.0);
 */
ExprPtr cst(double valeur);
ExprPtr frac(int64_t num, int64_t den = 1);

/*
 * Nom : var
 * Description : Helper pour créer rapidement un noeud variable.
 * Utilisation : ExprPtr v = var("x");
 */
ExprPtr var(const std::string& nom = "x");

/*
 * Nom : param
 * Description : Helper pour créer une constante symbolique (paramètre) nommée.
 * Utilisation : ExprPtr c1 = param("C1");
 */
ExprPtr param(const std::string& nom);

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
 * Description : Fonction pour créer le noeud tangente d'une expression.
 * Utilisation : ExprPtr t = ast_tan(expr);
 */
ExprPtr ast_tan(ExprPtr arg);

/*
 * Nom : ast_exp
 * Description : Fonction pour créer le noeud exponentielle d'une expression.
 * Utilisation : ExprPtr e = ast_exp(expr);
 */
ExprPtr ast_exp(ExprPtr arg);

/*
 * Nom : ast_ln
 * Description : Fonction pour créer le noeud logarithme d'une expression.
 * Utilisation : ExprPtr l = ast_ln(expr);
 */
ExprPtr ast_ln(ExprPtr arg);

} // namespace symalgo
