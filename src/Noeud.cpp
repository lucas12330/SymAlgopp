/**
 * @file Noeud.cpp
 * @brief Noyau de l'AST : table de hash-consing, ordre total, points d'entrée des
 *        opérations et construction canonique des expressions.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <sstream>
#include <unordered_map>

namespace symalgo {

// ============================================================================
// Hash-consing
// ============================================================================

namespace {

void combiner(std::size_t& graine, std::size_t valeur) {
    graine ^= valeur + 0x9e3779b97f4a7c15ULL + (graine << 6) + (graine >> 2);
}

/*
 * Nom : signatureDe
 * Description : Signature d'un noeud existant (même description que T::signature).
 */
Signature signatureDe(const ASTNode& n) {
    Signature s{n.type()};
    switch (n.type()) {
        case TypeNoeud::Constante: s.nombre = &static_cast<const Constante&>(n).getNombre(); break;
        case TypeNoeud::Variable: s.nom = &static_cast<const Variable&>(n).getNom(); break;
        case TypeNoeud::Parametre: s.nom = &static_cast<const Parametre&>(n).getNom(); break;
        case TypeNoeud::Pi: break;
        case TypeNoeud::Somme:
            s.nombre = &static_cast<const Somme&>(n).getConstante();
            s.termes = &static_cast<const Somme&>(n).getTermes();
            break;
        case TypeNoeud::Produit:
            s.nombre = &static_cast<const Produit&>(n).getCoefficient();
            s.facteurs = &static_cast<const Produit&>(n).getFacteurs();
            break;
        case TypeNoeud::Puissance:
            s.enfants[0] = static_cast<const Puissance&>(n).getBase().get();
            s.enfants[1] = static_cast<const Puissance&>(n).getExposant().get();
            break;
        case TypeNoeud::IntegraleNonEvaluee:
            s.enfants[0] = static_cast<const IntegraleNonEvaluee&>(n).getIntegrande().get();
            break;
        case TypeNoeud::LimiteNonEvaluee:
            s.enfants[0] = static_cast<const LimiteNonEvaluee&>(n).getExpression().get();
            s.reel = static_cast<const LimiteNonEvaluee&>(n).getPoint();
            break;
        default: s.enfants[0] = static_cast<const FonctionUnaire&>(n).m_argument.get(); break;
    }
    return s;
}

/*
 * Nom : memeSignature
 * Description : Égalité structurelle en un seul niveau : les enfants étant eux-mêmes uniques
 *               (hash-consing), il suffit de comparer leurs adresses.
 */
bool memeSignature(const Signature& a, const Signature& b) {
    if (a.type != b.type || a.enfants[0] != b.enfants[0] || a.enfants[1] != b.enfants[1]) return false;
    if ((a.nombre == nullptr) != (b.nombre == nullptr) || (a.nombre && !(*a.nombre == *b.nombre))) return false;
    if (a.reel != b.reel) return false;
    if ((a.nom == nullptr) != (b.nom == nullptr) || (a.nom && *a.nom != *b.nom)) return false;
    if (a.termes) {
        if (!b.termes || a.termes->size() != b.termes->size()) return false;
        for (std::size_t i = 0; i < a.termes->size(); ++i) {
            const Terme& ta = (*a.termes)[i];
            const Terme& tb = (*b.termes)[i];
            if (ta.expression != tb.expression || !(ta.coefficient == tb.coefficient)) return false;
        }
    }
    if (a.facteurs) {
        if (!b.facteurs || a.facteurs->size() != b.facteurs->size()) return false;
        for (std::size_t i = 0; i < a.facteurs->size(); ++i) {
            const Facteur& fa = (*a.facteurs)[i];
            const Facteur& fb = (*b.facteurs)[i];
            if (fa.base != fb.base || fa.exposant != fb.exposant) return false;
        }
    }
    return true;
}

} // namespace

std::size_t hashSignature(const Signature& s) {
    std::size_t h = static_cast<std::size_t>(s.type);
    for (const ASTNode* enfant : s.enfants) {
        if (enfant) combiner(h, enfant->hash());
    }
    if (s.nombre) combiner(h, s.nombre->hash());
    if (s.type == TypeNoeud::LimiteNonEvaluee) combiner(h, std::hash<double>()(s.reel == 0.0 ? 0.0 : s.reel));
    if (s.nom) combiner(h, std::hash<std::string>()(*s.nom));
    if (s.termes) {
        for (const Terme& t : *s.termes) {
            combiner(h, t.expression->hash());
            combiner(h, t.coefficient.hash());
        }
    }
    if (s.facteurs) {
        for (const Facteur& f : *s.facteurs) {
            combiner(h, f.base->hash());
            combiner(h, f.exposant->hash());
        }
    }
    return h;
}

/*
 * Nom : TableNoeuds
 * Description : Table de hachage des noeuds vivants, à chaînage intrusif (le lien vers le
 *               noeud suivant du seau est stocké dans le noeud) : aucune allocation à
 *               l'insertion ; la taille des seaux double quand la table se remplit.
 */
class TableNoeuds {
public:
    const ASTNode* trouver(const Signature& signature, std::size_t hash) const {
        if (m_seaux.empty()) return nullptr;
        for (const ASTNode* n = m_seaux[indice(hash)]; n; n = n->m_suivantTable) {
            if (n->m_hash == hash && memeSignature(signatureDe(*n), signature)) return n;
        }
        return nullptr;
    }

    void inserer(ASTNode* n) {
        if (m_taille + 1 > m_seaux.size()) redimensionner(m_seaux.empty() ? 1024 : 2 * m_seaux.size());
        const ASTNode*& tete = m_seaux[indice(n->m_hash)];
        n->m_suivantTable = tete;
        tete = n;
        ++m_taille;
    }

    void retirer(const ASTNode* n) {
        if (m_seaux.empty()) return;
        const ASTNode** lien = &m_seaux[indice(n->m_hash)];
        while (*lien && *lien != n) lien = &(*lien)->m_suivantTable;
        if (*lien) {
            *lien = n->m_suivantTable;
            --m_taille;
        }
    }

    std::size_t taille() const { return m_taille; }

private:
    std::size_t indice(std::size_t h) const { return h & (m_seaux.size() - 1); }

    void redimensionner(std::size_t nouvelleTaille) {
        std::vector<const ASTNode*> anciens(nouvelleTaille, nullptr);
        anciens.swap(m_seaux);
        for (const ASTNode* tete : anciens) {
            while (tete) {
                const ASTNode* suivant = tete->m_suivantTable;
                const ASTNode*& nouvelleTete = m_seaux[indice(tete->m_hash)];
                tete->m_suivantTable = nouvelleTete;
                nouvelleTete = tete;
                tete = suivant;
            }
        }
    }

    std::vector<const ASTNode*> m_seaux; // taille puissance de 2
    std::size_t m_taille = 0;
};

namespace {

// Table volontairement jamais détruite : des expressions globales peuvent être
// libérées après la fin de main(), elles doivent encore pouvoir s'en retirer
TableNoeuds& table() {
    static TableNoeuds* t = new TableNoeuds();
    return *t;
}

} // namespace

const ASTNode* chercherNoeud(const Signature& signature, std::size_t hash) {
    return table().trouver(signature, hash);
}

ExprPtr enregistrerNoeud(ASTNode* nouveau, std::size_t hash) {
    nouveau->m_hash = hash;
    table().inserer(nouveau);
    return ExprPtr(nouveau);
}

std::size_t nombreNoeudsVivants() { return table().taille(); }

void ASTNode::detruire() const {
    table().retirer(this);
    delete this;
}

// ============================================================================
// Points d'entrée des opérations
// ============================================================================

double ASTNode::getValeurConstante() const {
    return m_type == TypeNoeud::Constante ? static_cast<const Constante*>(this)->getNombre().versDouble() : 0.0;
}

ExprPtr ASTNode::derivee() const {
    CacheDerivees cache;
    return derivee(cache);
}

ExprPtr ASTNode::derivee(CacheDerivees& cache) const {
    // Un noeud référencé une seule fois ne peut être atteint qu'une fois : inutile de le
    // mémoriser (évite le coût de la table pour les expressions sans partage)
    if (nombreReferences() <= 1) return calculerDerivee(cache);
    const auto it = cache.find(this);
    if (it != cache.end()) return it->second;
    ExprPtr d = calculerDerivee(cache);
    cache.emplace(this, d);
    return d;
}

ExprPtr ASTNode::simplifier() const {
    if (m_estSimplifie) return clone();
    if (m_formeSimplifiee) return m_formeSimplifiee;
    ExprPtr resultat = calculerSimplification();
    resultat->m_estSimplifie = true;
    if (resultat.get() != this) m_formeSimplifiee = resultat;
    return resultat;
}

ExprPtr ASTNode::integrer() const {
    if (!contientVariable()) return clone() * var("x");
    return primitive();
}

ExprPtr ASTNode::limite(double a) const { return calculerLimite(a); }

ExprPtr ASTNode::integraleNonEvaluee() const { return fabriquer<IntegraleNonEvaluee>(clone()); }

ExprPtr ASTNode::limiteNonEvaluee(double a) const { return fabriquer<LimiteNonEvaluee>(clone(), a); }

std::string ASTNode::texte() const {
    std::ostringstream os;
    afficher(os);
    return os.str();
}

std::ostream& operator<<(std::ostream& os, const ExprPtr& e) {
    e->afficher(os);
    return os;
}

// ============================================================================
// Ordre total
// ============================================================================

namespace {

int rang(TypeNoeud t) {
    switch (t) {
        case TypeNoeud::Constante: return 0;
        case TypeNoeud::Pi: return 1;        // pi*x
        case TypeNoeud::Parametre: return 2; // C1*x plutôt que x*C1
        case TypeNoeud::Variable: return 3;
        case TypeNoeud::Puissance: return 4;
        case TypeNoeud::Produit: return 5;
        case TypeNoeud::Somme: return 6;
        default: return 7 + static_cast<int>(t) - static_cast<int>(TypeNoeud::Sinus);
    }
}

int signeDe(int c) { return (c > 0) - (c < 0); }

} // namespace

int comparer(const ASTNode& a, const ASTNode& b) {
    if (&a == &b) return 0;
    const int r = rang(a.type()) - rang(b.type());
    if (r != 0) return signeDe(r);
    switch (a.type()) {
        case TypeNoeud::Constante:
            return static_cast<const Constante&>(a).getNombre().comparer(static_cast<const Constante&>(b).getNombre());
        case TypeNoeud::Variable:
            return signeDe(static_cast<const Variable&>(a).getNom().compare(static_cast<const Variable&>(b).getNom()));
        case TypeNoeud::Parametre:
            return signeDe(static_cast<const Parametre&>(a).getNom().compare(static_cast<const Parametre&>(b).getNom()));
        case TypeNoeud::Pi: return 0; // noeud unique
        case TypeNoeud::Puissance: {
            const Puissance& pa = static_cast<const Puissance&>(a);
            const Puissance& pb = static_cast<const Puissance&>(b);
            if (const int c = comparer(*pa.getBase(), *pb.getBase())) return c;
            return comparer(*pa.getExposant(), *pb.getExposant());
        }
        case TypeNoeud::Produit: {
            const Produit& pa = static_cast<const Produit&>(a);
            const Produit& pb = static_cast<const Produit&>(b);
            const auto& fa = pa.getFacteurs();
            const auto& fb = pb.getFacteurs();
            for (std::size_t i = 0; i < std::min(fa.size(), fb.size()); ++i) {
                if (const int c = comparer(*fa[i].base, *fb[i].base)) return c;
                if (const int c = comparer(*fa[i].exposant, *fb[i].exposant)) return c;
            }
            if (fa.size() != fb.size()) return fa.size() < fb.size() ? -1 : 1;
            return pa.getCoefficient().comparer(pb.getCoefficient());
        }
        case TypeNoeud::Somme: {
            const Somme& sa = static_cast<const Somme&>(a);
            const Somme& sb = static_cast<const Somme&>(b);
            const auto& ta = sa.getTermes();
            const auto& tb = sb.getTermes();
            for (std::size_t i = 0; i < std::min(ta.size(), tb.size()); ++i) {
                if (const int c = comparer(*ta[i].expression, *tb[i].expression)) return c;
                if (const int c = ta[i].coefficient.comparer(tb[i].coefficient)) return c;
            }
            if (ta.size() != tb.size()) return ta.size() < tb.size() ? -1 : 1;
            return sa.getConstante().comparer(sb.getConstante());
        }
        case TypeNoeud::IntegraleNonEvaluee:
            return comparer(*static_cast<const IntegraleNonEvaluee&>(a).getIntegrande(),
                            *static_cast<const IntegraleNonEvaluee&>(b).getIntegrande());
        case TypeNoeud::LimiteNonEvaluee: {
            const LimiteNonEvaluee& la = static_cast<const LimiteNonEvaluee&>(a);
            const LimiteNonEvaluee& lb = static_cast<const LimiteNonEvaluee&>(b);
            if (const int c = comparer(*la.getExpression(), *lb.getExpression())) return c;
            return (la.getPoint() > lb.getPoint()) - (la.getPoint() < lb.getPoint());
        }
        default:
            return comparer(*static_cast<const FonctionUnaire&>(a).m_argument,
                            *static_cast<const FonctionUnaire&>(b).m_argument);
    }
}

// ============================================================================
// Constructeurs des noeuds
// ============================================================================

namespace {

// Exposant entier (pour l'évaluation rapide), ou PAS_ENTIER
long long exposantEntier(const ExprPtr& e) {
    long long k;
    if (const Constante* c = comme<Constante>(e); c && c->getNombre().versEntier(k)) return k;
    return PAS_ENTIER;
}

} // namespace

Constante::Constante(CleFabrique, Nombre valeur)
    : ASTNode(TYPE, false), m_valeur(std::move(valeur)), m_approx(m_valeur.versDouble()) {}

Variable::Variable(CleFabrique, const std::string& nom) : ASTNode(TYPE, true), m_nom(nom) {}

Parametre::Parametre(CleFabrique, const std::string& nom) : ASTNode(TYPE, false), m_nom(nom) {}

Pi::Pi(CleFabrique) : ASTNode(TYPE, false) {}

namespace {

bool termesContiennentVariable(const std::vector<Terme>& termes) {
    return std::any_of(termes.begin(), termes.end(), [](const Terme& t) { return t.expression->contientVariable(); });
}

bool facteursContiennentVariable(const std::vector<Facteur>& facteurs) {
    return std::any_of(facteurs.begin(), facteurs.end(), [](const Facteur& f) {
        return f.base->contientVariable() || f.exposant->contientVariable();
    });
}

} // namespace

Somme::Somme(CleFabrique, Nombre constante, std::vector<Terme> termes)
    : ASTNode(TYPE, termesContiennentVariable(termes)), m_constante(std::move(constante)), m_termes(std::move(termes)) {}

Produit::Produit(CleFabrique, Nombre coefficient, std::vector<Facteur> facteurs)
    : ASTNode(TYPE, facteursContiennentVariable(facteurs)),
      m_coefficient(std::move(coefficient)),
      m_facteurs(std::move(facteurs)) {}

Puissance::Puissance(CleFabrique, ExprPtr base, ExprPtr exposant)
    : ASTNode(TYPE, base->contientVariable() || exposant->contientVariable()),
      m_base(std::move(base)),
      m_exposant(std::move(exposant)),
      m_exposantEntier(exposantEntier(m_exposant)) {}

FonctionUnaire::FonctionUnaire(TypeNoeud type, ExprPtr arg)
    : ASTNode(type, arg->contientVariable()), m_argument(std::move(arg)) {}

IntegraleNonEvaluee::IntegraleNonEvaluee(CleFabrique, ExprPtr integrande)
    : ASTNode(TYPE, true), m_integrande(std::move(integrande)) {}

LimiteNonEvaluee::LimiteNonEvaluee(CleFabrique, ExprPtr expression, double point)
    : ASTNode(TYPE, false), m_expression(std::move(expression)), m_point(point) {}

// ============================================================================
// Construction canonique
// ============================================================================

void AccumulateurSomme::ajouter(const ExprPtr& e, const Nombre& facteur) {
    switch (e->type()) {
        case TypeNoeud::Constante:
            m_constante = m_constante + facteur * static_cast<const Constante&>(*e).getNombre();
            break;
        case TypeNoeud::Somme: {
            const Somme& s = static_cast<const Somme&>(*e);
            m_constante = m_constante + facteur * s.getConstante();
            for (const Terme& t : s.getTermes()) ajouterTerme(t.expression, facteur * t.coefficient);
            break;
        }
        case TypeNoeud::Produit: {
            // Le coefficient numérique d'un produit devient le coefficient du terme
            const Produit& p = static_cast<const Produit&>(*e);
            if (p.getCoefficient().estUn()) {
                ajouterTerme(e, facteur);
            } else {
                ajouterTerme(produitSansCoefficient(p), facteur * p.getCoefficient());
            }
            break;
        }
        default: ajouterTerme(e, facteur); break;
    }
}

namespace {

/*
 * Nom : trouverIndice
 * Description : Position de l'élément de clé donnée (adresse du noeud), ou taille si absent.
 *               Recherche linéaire pour les petites collections (cas courant : 2 ou 3
 *               éléments), index par adresse au-delà de SEUIL_INDEX.
 */
template <class Element, class Cle>
std::size_t trouverIndice(const std::vector<Element>& elements, std::unordered_map<const ASTNode*, std::size_t>& index,
                          const ASTNode* cle, Cle extraire) {
    if (elements.size() <= SEUIL_INDEX) {
        for (std::size_t i = 0; i < elements.size(); ++i) {
            if (extraire(elements[i]) == cle) return i;
        }
        return elements.size();
    }
    if (index.empty()) {
        for (std::size_t i = 0; i < elements.size(); ++i) index.emplace(extraire(elements[i]), i);
    }
    const auto it = index.find(cle);
    return it == index.end() ? elements.size() : it->second;
}

} // namespace

void AccumulateurSomme::ajouterTerme(const ExprPtr& terme, const Nombre& coefficient) {
    const std::size_t i =
        trouverIndice(m_termes, m_indices, terme.get(), [](const Terme& t) { return t.expression.get(); });
    if (i < m_termes.size()) {
        Nombre& c = m_termes[i].coefficient;
        c = c + coefficient;
        return;
    }
    if (!m_indices.empty()) m_indices.emplace(terme.get(), m_termes.size());
    m_termes.push_back({terme, coefficient});
}

ExprPtr AccumulateurSomme::construire() {
    std::vector<Terme> termes;
    termes.reserve(m_termes.size());
    for (Terme& t : m_termes) {
        if (!t.coefficient.estZero()) termes.push_back(std::move(t));
    }
    std::sort(termes.begin(), termes.end(),
              [](const Terme& a, const Terme& b) { return comparer(*a.expression, *b.expression) < 0; });
    if (termes.empty()) return nombre(m_constante);
    if (termes.size() == 1 && m_constante.estZero()) {
        if (termes[0].coefficient.estUn()) return termes[0].expression;
        AccumulateurProduit p;
        p.multiplier(nombre(termes[0].coefficient));
        p.multiplier(termes[0].expression);
        return p.construire();
    }
    if (m_constante.estZero()) m_constante = Nombre(0); // un 0.0 réel n'est pas conservé
    return fabriquer<Somme>(std::move(m_constante), std::move(termes));
}

void AccumulateurProduit::multiplier(const ExprPtr& e) {
    switch (e->type()) {
        case TypeNoeud::Constante:
            m_coefficient = m_coefficient * static_cast<const Constante&>(*e).getNombre();
            break;
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(*e);
            m_coefficient = m_coefficient * p.getCoefficient();
            for (const Facteur& f : p.getFacteurs()) ajouterFacteur(f.base, f.exposant);
            break;
        }
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(*e);
            ajouterFacteur(p.getBase(), p.getExposant());
            break;
        }
        default: ajouterFacteur(e, un()); break;
    }
}

void AccumulateurProduit::ajouterFacteur(const ExprPtr& base, const ExprPtr& exposant) {
    const std::size_t i =
        trouverIndice(m_facteurs, m_indices, base.get(), [](const Facteur& f) { return f.base.get(); });
    if (i < m_facteurs.size()) {
        ExprPtr& e = m_facteurs[i].exposant;
        e = e + exposant; // x^a * x^b = x^(a+b)
        return;
    }
    if (!m_indices.empty()) m_indices.emplace(base.get(), m_facteurs.size());
    m_facteurs.push_back({base, exposant});
}

ExprPtr AccumulateurProduit::construire() {
    std::vector<Facteur> facteurs;
    facteurs.reserve(m_facteurs.size());
    bool divisionParZero = false;
    for (Facteur& f : m_facteurs) {
        // Cas courant : base ordinaire et exposant numérique non nul, rien à réduire
        const Constante* e = comme<Constante>(f.exposant);
        const TypeNoeud tb = f.base->type();
        if (e && !e->getNombre().estZero() && tb != TypeNoeud::Constante && tb != TypeNoeud::Puissance &&
            tb != TypeNoeud::Produit) {
            facteurs.push_back(std::move(f));
            continue;
        }
        // Sinon le facteur est remis sous forme canonique : il peut se réduire à une
        // constante (2^3 = 8, 4^(1/2) = 2) ou disparaître (x^0 = 1)
        const ExprPtr p = ast_pow(f.base, f.exposant);
        switch (p->type()) {
            case TypeNoeud::Constante:
                m_coefficient = m_coefficient * static_cast<const Constante&>(*p).getNombre();
                break;
            case TypeNoeud::Puissance: {
                const Puissance& q = static_cast<const Puissance&>(*p);
                const Constante* b = comme<Constante>(q.getBase());
                divisionParZero = divisionParZero || (b && b->getNombre().estZero());
                facteurs.push_back({q.getBase(), q.getExposant()});
                break;
            }
            case TypeNoeud::Produit: {
                const Produit& q = static_cast<const Produit&>(*p);
                m_coefficient = m_coefficient * q.getCoefficient();
                for (const Facteur& g : q.getFacteurs()) facteurs.push_back(g);
                break;
            }
            default: facteurs.push_back({p, un()}); break;
        }
    }
    if (m_coefficient.estZero()) {
        // 0 * x^(-1) = 0, mais 0 * 0^(-1) (0/0) est indéfini
        return divisionParZero ? nombre(Nombre::reel(std::nan(""))) : nombre(m_coefficient);
    }
    std::sort(facteurs.begin(), facteurs.end(),
              [](const Facteur& a, const Facteur& b) { return comparer(*a.base, *b.base) < 0; });
    if (facteurs.empty()) return nombre(m_coefficient);
    if (facteurs.size() == 1) {
        const Facteur& f = facteurs[0];
        const bool exposantUn = f.exposant.get() == un().get();
        const ExprPtr seul = exposantUn ? f.base : fabriquer<Puissance>(f.base, f.exposant);
        if (m_coefficient.estUn()) return seul;
        // Un coefficient numérique se distribue sur une somme : 2*(x + 1) = 2*x + 2
        if (exposantUn && f.base->type() == TypeNoeud::Somme) {
            AccumulateurSomme s;
            s.ajouter(f.base, m_coefficient);
            return s.construire();
        }
    }
    return fabriquer<Produit>(std::move(m_coefficient), std::move(facteurs));
}

ExprPtr produitSansCoefficient(const Produit& p) {
    const auto& facteurs = p.getFacteurs();
    if (facteurs.size() == 1) return ast_pow(facteurs[0].base, facteurs[0].exposant);
    return fabriquer<Produit>(Nombre(1), facteurs);
}

const ExprPtr& un() {
    static const ExprPtr* valeur = new ExprPtr(fabriquer<Constante>(Nombre(1)));
    return *valeur;
}

ExprPtr nombre(const Nombre& valeur) { return fabriquer<Constante>(valeur); }

ExprPtr cst(double valeur) { return fabriquer<Constante>(Nombre::depuisDouble(valeur)); }

ExprPtr frac(std::int64_t num, std::int64_t den) {
    return fabriquer<Constante>(Nombre::rationnel(num, den));
}

ExprPtr var(const std::string& nom) { return fabriquer<Variable>(nom); }

ExprPtr param(const std::string& nom) { return fabriquer<Parametre>(nom); }

ExprPtr pi() { return fabriquer<Pi>(); }

ExprPtr somme(const std::vector<ExprPtr>& termes) {
    AccumulateurSomme acc;
    for (const ExprPtr& t : termes) acc.ajouter(t, Nombre(1));
    return acc.construire();
}

ExprPtr produit(const std::vector<ExprPtr>& facteurs) {
    AccumulateurProduit acc;
    for (const ExprPtr& f : facteurs) acc.multiplier(f);
    return acc.construire();
}

ExprPtr operator+(const ExprPtr& gauche, const ExprPtr& droite) {
    AccumulateurSomme acc;
    acc.ajouter(gauche, Nombre(1));
    acc.ajouter(droite, Nombre(1));
    return acc.construire();
}

ExprPtr operator-(const ExprPtr& gauche, const ExprPtr& droite) {
    AccumulateurSomme acc;
    acc.ajouter(gauche, Nombre(1));
    acc.ajouter(droite, Nombre(-1));
    return acc.construire();
}

ExprPtr operator-(const ExprPtr& e) {
    AccumulateurSomme acc;
    acc.ajouter(e, Nombre(-1));
    return acc.construire();
}

ExprPtr operator*(const ExprPtr& gauche, const ExprPtr& droite) {
    AccumulateurProduit acc;
    acc.multiplier(gauche);
    acc.multiplier(droite);
    return acc.construire();
}

ExprPtr operator/(const ExprPtr& gauche, const ExprPtr& droite) {
    AccumulateurProduit acc;
    acc.multiplier(gauche);
    acc.multiplier(ast_pow(droite, nombre(Nombre(-1))));
    return acc.construire();
}

ExprPtr operator+(const ExprPtr& g, double d) { return g + cst(d); }
ExprPtr operator+(double g, const ExprPtr& d) { return cst(g) + d; }
ExprPtr operator-(const ExprPtr& g, double d) { return g - cst(d); }
ExprPtr operator-(double g, const ExprPtr& d) { return cst(g) - d; }
ExprPtr operator*(const ExprPtr& g, double d) { return g * cst(d); }
ExprPtr operator*(double g, const ExprPtr& d) { return cst(g) * d; }
ExprPtr operator/(const ExprPtr& g, double d) { return g / cst(d); }
ExprPtr operator/(double g, const ExprPtr& d) { return cst(g) / d; }

namespace {

/*
 * Nom : extrairePuissances
 * Description : Écrit l'entier b > 0 (au plus 10^12) sous la forme exterieur^q * interieur
 *               avec exterieur maximal (division par les nombres premiers jusqu'à b^(1/q)).
 *               Renvoie faux s'il n'y a rien à extraire.
 */
bool extrairePuissances(const Nombre& b, long long q, Nombre& exterieur, Nombre& interieur) {
    long long v;
    if (!b.versEntier(v) || v <= 1 || v > 1000000000000LL || q < 2) return false;
    long long ext = 1, reste = v;
    for (long long d = 2; d <= 1000000 && d * d <= reste; ++d) {
        long long multiplicite = 0;
        while (reste % d == 0) {
            reste /= d;
            ++multiplicite;
        }
        for (long long i = 0; i < multiplicite / q; ++i) ext *= d;
    }
    // Le reste, s'il dépasse 1, est premier (multiplicité 1) : rien à extraire
    if (ext == 1) return false;
    exterieur = Nombre(ext);
    interieur = Nombre(v) / Nombre(ext).puissanceEntiere(q);
    return true;
}

} // namespace

ExprPtr ast_pow(const ExprPtr& base, const ExprPtr& exposant) {
    if (const Constante* ce = comme<Constante>(exposant)) {
        const Nombre& e = ce->getNombre();
        if (e.estZero()) return un();          // u^0 = 1
        if (e.estUn()) return base;            // u^1 = u
        if (const Constante* cb = comme<Constante>(base)) {
            const Nombre& b = cb->getNombre();
            if (b.estExact() && e.estExact()) {
                // 0^e (e < 0) est infini : une seule forme, 1/0 (0^(-3) et 0^(-1/2) aussi)
                if (b.estZero() && e.signe() < 0) {
                    return e.estMoinsUn() ? fabriquer<Puissance>(base, exposant) : ast_pow(base, nombre(Nombre(-1)));
                }
                long long k;
                if (e.versEntier(k)) {
                    return nombre(b.puissanceEntiere(k));
                } else if (b.signe() >= 0) {
                    // Exposant p/q : exact si b est une puissance q-ième parfaite (4^(1/2) = 2)
                    long long p, q;
                    Nombre racine;
                    if (e.versFraction(p, q) && b.racineExacte(q, racine)) {
                        return nombre(racine.puissanceEntiere(p));
                    }
                    // Forme unique des radicaux : b^(p/q) = b^n * b^(p/q - n), exposant dans ]0, 1[
                    // (3^(-1/2) s'écrit 3^(1/2)/3, 2^(3/2) s'écrit 2*2^(1/2))
                    if (e.versFraction(p, q) && !b.estZero()) {
                        const long long n = p >= 0 ? p / q : -((-p + q - 1) / q); // partie entière
                        if (n != 0) return nombre(b.puissanceEntiere(n)) * ast_pow(base, nombre(Nombre::rationnel(p - n * q, q)));
                        // Dénominateur rationalisé : (a/d)^(p/q) = (a d^(q-1))^(p/q) / d^p
                        const Nombre d = b.denominateurNombre();
                        if (!d.estUn()) {
                            return ast_pow(nombre(b * d.puissanceEntiere(q)), exposant) * nombre(d.puissanceEntiere(-p));
                        }
                        // Puissances q-ièmes parfaites sorties du radical : 8^(1/2) = 2*2^(1/2)
                        Nombre exterieur, interieur;
                        if (extrairePuissances(b, q, exterieur, interieur)) {
                            return nombre(exterieur.puissanceEntiere(p)) * ast_pow(nombre(interieur), exposant);
                        }
                    }
                }
            } else {
                const double v = std::pow(b.versDouble(), e.versDouble());
                if (std::isfinite(v)) return nombre(Nombre::reel(v));
            }
            return fabriquer<Puissance>(base, exposant);
        }
        long long k;
        if (e.versEntier(k)) {
            // (u^a)^k = u^(a k) pour k entier
            if (const Puissance* p = comme<Puissance>(base)) return ast_pow(p->getBase(), p->getExposant() * exposant);
            // (c * prod b_i^e_i)^k = c^k * prod b_i^(e_i k)
            if (const Produit* p = comme<Produit>(base)) {
                AccumulateurProduit acc;
                acc.multiplier(nombre(p->getCoefficient().puissanceEntiere(k)));
                for (const Facteur& f : p->getFacteurs()) acc.multiplier(ast_pow(f.base, f.exposant * exposant));
                return acc.construire();
            }
        }
    } else if (const Constante* cb = comme<Constante>(base); cb && cb->getNombre().estUn()) {
        return un(); // 1^u = 1
    }
    return fabriquer<Puissance>(base, exposant);
}

ExprPtr ast_pow(const ExprPtr& base, double exposant) { return ast_pow(base, cst(exposant)); }
ExprPtr ast_pow(double base, const ExprPtr& exposant) { return ast_pow(cst(base), exposant); }

namespace {

// Valeur numérique d'un argument réel (non exact), sinon faux
bool argumentReel(const ExprPtr& arg, double& v) {
    const Constante* c = comme<Constante>(arg);
    if (!c || c->getNombre().estExact()) return false;
    v = c->getNombre().versDouble();
    return true;
}

bool argumentExactNul(const ExprPtr& arg) {
    const Constante* c = comme<Constante>(arg);
    return c && c->getNombre().estExact() && c->getNombre().estZero();
}

} // namespace

namespace {

// Si arg = r*pi avec r rationnel exact (ou arg = 0), écrit r
bool multipleDePi(const ExprPtr& arg, Nombre& r) {
    if (argumentExactNul(arg)) {
        r = Nombre(0);
        return true;
    }
    if (arg->type() == TypeNoeud::Pi) {
        r = Nombre(1);
        return true;
    }
    const Produit* p = comme<Produit>(arg);
    if (p && p->getCoefficient().estExact() && p->getFacteurs().size() == 1 &&
        p->getFacteurs()[0].base->type() == TypeNoeud::Pi && p->getFacteurs()[0].exposant.get() == un().get()) {
        r = p->getCoefficient();
        return true;
    }
    return false;
}

/*
 * Nom : sinusExact
 * Description : Valeur exacte de sin(r*pi) quand r est un multiple de 1/6 ou de 1/4 :
 *               0, ±1/2, ±2^(1/2)/2, ±3^(1/2)/2, ±1.
 */
bool sinusExact(const Nombre& r, ExprPtr& valeur) {
    long long p, q;
    if (!r.versFraction(p, q) || 12 % q != 0) return false;
    long long k = (p % (2 * q)) * (12 / q); // r*12 modulo 24 (période 2*pi)
    k = ((k % 24) + 24) % 24;
    const ExprPtr demi = frac(1, 2);
    auto racine = [&](long long n) { return ast_pow(cst(static_cast<double>(n)), frac(1, 2)) * demi; };
    switch (k) {
        case 0: case 12: valeur = nombre(Nombre(0)); return true;
        case 2: case 10: valeur = demi; return true;
        case 14: case 22: valeur = -demi; return true;
        case 3: case 9: valeur = racine(2); return true;
        case 15: case 21: valeur = -racine(2); return true;
        case 4: case 8: valeur = racine(3); return true;
        case 16: case 20: valeur = -racine(3); return true;
        case 6: valeur = un(); return true;
        case 18: valeur = nombre(Nombre(-1)); return true;
        default: return false;
    }
}

bool cosinusExact(const Nombre& r, ExprPtr& valeur) { return sinusExact(r + Nombre::rationnel(1, 2), valeur); }

// Multiples r de pi dans [-1/2, 1/2] dont le sinus est tabulé (réciproques exactes)
const std::vector<Nombre>& anglesRemarquables() {
    static const std::vector<Nombre>* angles = new std::vector<Nombre>{
        Nombre::rationnel(-1, 2), Nombre::rationnel(-1, 3), Nombre::rationnel(-1, 4), Nombre::rationnel(-1, 6), Nombre(0),
        Nombre::rationnel(1, 6),  Nombre::rationnel(1, 4),  Nombre::rationnel(1, 3),  Nombre::rationnel(1, 2)};
    return *angles;
}

ExprPtr foisPi(const Nombre& r) { return nombre(r) * pi(); }

} // namespace

ExprPtr ast_sin(const ExprPtr& arg) {
    double v;
    Nombre r;
    ExprPtr exact;
    if (multipleDePi(arg, r) && sinusExact(r, exact)) return exact;
    if (argumentReel(arg, v)) return nombre(Nombre::reel(std::sin(v)));
    return fabriquer<Sinus>(arg);
}

ExprPtr ast_cos(const ExprPtr& arg) {
    double v;
    Nombre r;
    ExprPtr exact;
    if (multipleDePi(arg, r) && cosinusExact(r, exact)) return exact;
    if (argumentReel(arg, v)) return nombre(Nombre::reel(std::cos(v)));
    return fabriquer<Cosinus>(arg);
}

ExprPtr ast_tan(const ExprPtr& arg) {
    double v;
    Nombre r;
    ExprPtr s, c;
    if (multipleDePi(arg, r) && sinusExact(r, s) && cosinusExact(r, c) && !argumentExactNul(c)) return s / c;
    if (argumentReel(arg, v)) return nombre(Nombre::reel(std::tan(v)));
    return fabriquer<Tangente>(arg);
}

ExprPtr ast_asin(const ExprPtr& arg) {
    double v;
    for (const Nombre& r : anglesRemarquables()) {
        ExprPtr s;
        if (sinusExact(r, s) && s.get() == arg.get()) return foisPi(r); // asin(1/2) = pi/6
    }
    if (argumentReel(arg, v) && std::abs(v) <= 1.0) return nombre(Nombre::reel(std::asin(v)));
    return fabriquer<ArcSinus>(arg);
}

ExprPtr ast_acos(const ExprPtr& arg) {
    double v;
    for (const Nombre& r : anglesRemarquables()) {
        ExprPtr s;
        if (sinusExact(r, s) && s.get() == arg.get()) return foisPi(Nombre::rationnel(1, 2) - r); // pi/2 - asin
    }
    if (argumentReel(arg, v) && std::abs(v) <= 1.0) return nombre(Nombre::reel(std::acos(v)));
    return fabriquer<ArcCosinus>(arg);
}

ExprPtr ast_atan(const ExprPtr& arg) {
    double v;
    for (const Nombre& r : anglesRemarquables()) {
        ExprPtr s, c;
        if (r.valeurAbsolue() == Nombre::rationnel(1, 2)) continue; // tan non définie
        if (sinusExact(r, s) && cosinusExact(r, c) && (s / c).get() == arg.get()) return foisPi(r);
    }
    if (argumentReel(arg, v)) return nombre(Nombre::reel(std::atan(v)));
    return fabriquer<ArcTangente>(arg);
}

ExprPtr ast_exp(const ExprPtr& arg) {
    double v;
    if (argumentExactNul(arg)) return un();
    if (argumentReel(arg, v)) return nombre(Nombre::reel(std::exp(v)));
    if (const Logarithme* l = comme<Logarithme>(arg)) return l->m_argument; // exp(ln u) = u
    return fabriquer<Exponentielle>(arg);
}

ExprPtr ast_ln(const ExprPtr& arg) {
    double v;
    if (const Constante* c = comme<Constante>(arg); c && c->getNombre().estExact()) {
        const Nombre& b = c->getNombre();
        if (b.estUn()) return nombre(Nombre(0));
        // ln(p/q) = ln(p) - ln(q) et ln(a^n) = n ln(a) : ln(8) = 3*ln(2), ln(1/4) = -2*ln(2)
        if (b.signe() > 0 && !b.estEntier()) return ast_ln(nombre(b.numerateurNombre())) - ast_ln(nombre(b.denominateurNombre()));
        long long v64;
        if (b.versEntier(v64) && v64 > 1) {
            for (long long n = 62; n >= 2; --n) {
                Nombre racine;
                if (b.racineExacte(n, racine)) return nombre(Nombre(n)) * ast_ln(nombre(racine));
            }
        }
    }
    if (argumentReel(arg, v) && v > 0.0) return nombre(Nombre::reel(std::log(v)));
    if (const Exponentielle* e = comme<Exponentielle>(arg)) return e->m_argument; // ln(exp u) = u (u réel)
    return fabriquer<Logarithme>(arg);
}

// ============================================================================
// Substitution
// ============================================================================

ExprPtr appliquer(TypeNoeud fonction, const ExprPtr& arg) {
    switch (fonction) {
        case TypeNoeud::Sinus: return ast_sin(arg);
        case TypeNoeud::Cosinus: return ast_cos(arg);
        case TypeNoeud::Tangente: return ast_tan(arg);
        case TypeNoeud::Exponentielle: return ast_exp(arg);
        case TypeNoeud::Logarithme: return ast_ln(arg);
        case TypeNoeud::ArcSinus: return ast_asin(arg);
        case TypeNoeud::ArcCosinus: return ast_acos(arg);
        case TypeNoeud::ArcTangente: return ast_atan(arg);
        default: throw std::invalid_argument("appliquer : type de fonction inconnu");
    }
}

namespace {

class Substitution {
public:
    Substitution(const ExprPtr& cible, const ExprPtr& remplacement) : m_cible(cible), m_remplacement(remplacement) {}

    ExprPtr appliquerA(const ExprPtr& e) {
        if (e.get() == m_cible.get()) return m_remplacement;
        const auto it = m_memo.find(e.get());
        if (it != m_memo.end()) return it->second;
        ExprPtr r = e;
        switch (e->type()) {
            case TypeNoeud::Somme: {
                const Somme& s = static_cast<const Somme&>(*e);
                AccumulateurSomme acc;
                acc.ajouter(nombre(s.getConstante()), Nombre(1));
                for (const Terme& t : s.getTermes()) acc.ajouter(appliquerA(t.expression), t.coefficient);
                r = acc.construire();
                break;
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(*e);
                AccumulateurProduit acc;
                acc.multiplier(nombre(p.getCoefficient()));
                for (const Facteur& f : p.getFacteurs()) acc.multiplier(ast_pow(appliquerA(f.base), appliquerA(f.exposant)));
                r = acc.construire();
                break;
            }
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(*e);
                r = ast_pow(appliquerA(p.getBase()), appliquerA(p.getExposant()));
                break;
            }
            case TypeNoeud::IntegraleNonEvaluee:
                r = fabriquer<IntegraleNonEvaluee>(appliquerA(static_cast<const IntegraleNonEvaluee&>(*e).getIntegrande()));
                break;
            default:
                if (const FonctionUnaire* f = comme<FonctionUnaire>(e)) r = appliquer(e->type(), appliquerA(f->m_argument));
        }
        m_memo.emplace(e.get(), r);
        return r;
    }

private:
    ExprPtr m_cible, m_remplacement;
    std::unordered_map<const ASTNode*, ExprPtr> m_memo;
};

bool contientRecursif(const ASTNode& e, const ASTNode& cible, std::unordered_map<const ASTNode*, bool>& memo) {
    if (&e == &cible) return true;
    const auto it = memo.find(&e);
    if (it != memo.end()) return it->second;
    bool r = false;
    switch (e.type()) {
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(e).getTermes()) r = r || contientRecursif(*t.expression, cible, memo);
            break;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) {
                r = r || contientRecursif(*f.base, cible, memo) || contientRecursif(*f.exposant, cible, memo);
            }
            break;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            r = contientRecursif(*p.getBase(), cible, memo) || contientRecursif(*p.getExposant(), cible, memo);
            break;
        }
        case TypeNoeud::IntegraleNonEvaluee:
            r = contientRecursif(*static_cast<const IntegraleNonEvaluee&>(e).getIntegrande(), cible, memo);
            break;
        case TypeNoeud::LimiteNonEvaluee:
            r = contientRecursif(*static_cast<const LimiteNonEvaluee&>(e).getExpression(), cible, memo);
            break;
        default:
            if (const FonctionUnaire* f = comme<FonctionUnaire>(&e)) r = contientRecursif(*f->m_argument, cible, memo);
    }
    memo.emplace(&e, r);
    return r;
}

} // namespace

ExprPtr substituer(const ExprPtr& e, const ExprPtr& cible, const ExprPtr& remplacement) {
    return Substitution(cible, remplacement).appliquerA(e);
}

bool contient(const ExprPtr& e, const ExprPtr& cible) {
    std::unordered_map<const ASTNode*, bool> memo;
    return contientRecursif(*e, *cible, memo);
}

} // namespace symalgo
