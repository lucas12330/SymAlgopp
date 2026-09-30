/**
 * @file Affichage.cpp
 * @brief Écriture lisible des expressions canoniques (« x^2 + 5*x + 6 », « sin(x)/x »).
 *
 * Règles : parenthèses selon la précédence (somme < produit < puissance < atome), les
 * facteurs d'exposant négatif et les dénominateurs des coefficients rationnels forment un
 * quotient, et les termes d'une somme sont écrits par degré décroissant, la constante en
 * dernier.
 */

#include "ASTNode.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <sstream>

namespace symalgo {

namespace {

enum Precedence { SOMME = 1, PRODUIT = 2, PUISSANCE = 3, ATOME = 4 };

void ecrire(std::ostream& os, const ASTNode& n, int precedenceMin);

bool exposantNegatif(const ExprPtr& e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().signe() < 0;
}

int precedence(const ASTNode& n) {
    switch (n.type()) {
        case TypeNoeud::Somme: return SOMME;
        case TypeNoeud::Produit: return PRODUIT;
        case TypeNoeud::Constante: {
            const Nombre& v = static_cast<const Constante&>(n).getNombre();
            if (v.signe() < 0 || (v.estExact() && !v.estEntier())) return PRODUIT; // « -3 », « 1/3 »
            return ATOME;
        }
        case TypeNoeud::Puissance:
            return exposantNegatif(static_cast<const Puissance&>(n).getExposant()) ? PRODUIT : PUISSANCE;
        default: return ATOME;
    }
}

// Texte d'une sous-expression, avec le formatage (précision...) du flux de destination
std::string texte(const std::ostream& modele, const ASTNode& n, int precedenceMin) {
    std::ostringstream os;
    os.copyfmt(modele);
    ecrire(os, n, precedenceMin);
    return os.str();
}

std::string texteNombre(const std::ostream& modele, const Nombre& v) {
    std::ostringstream os;
    os.copyfmt(modele);
    v.afficher(os);
    return os.str();
}

// « b » ou « b^e » (e positif)
std::string textePuissance(const std::ostream& modele, const ExprPtr& base, const ExprPtr& exposant) {
    if (exposant.get() == un().get()) return texte(modele, *base, PUISSANCE);
    std::string t = texte(modele, *base, ATOME) + "^";
    // Exposant atomique (entier positif, variable, fonction) : pas de parenthèses (x^2, x^x)
    if (precedence(*exposant) == ATOME) return t + texte(modele, *exposant, ATOME);
    return t + "(" + texte(modele, *exposant, 0) + ")";
}

std::string joindre(const std::vector<std::string>& parties) {
    std::string r;
    for (std::size_t i = 0; i < parties.size(); ++i) r += (i ? "*" : "") + parties[i];
    return r;
}

/*
 * Nom : ecrireProduit
 * Description : Écrit coefficient * prod base^exposant (coefficient positif ou nul) en
 *               regroupant au dénominateur les exposants négatifs et le dénominateur du
 *               coefficient : « 2*x/3 », « sin(x)/x », « 1/(x + 1) ».
 */
void ecrireProduit(std::ostream& os, const Nombre& coefficient, const std::vector<Facteur>& facteurs) {
    std::vector<std::string> numerateur, denominateur;
    if (coefficient.estExact()) {
        if (coefficient.numerateur() != "1") numerateur.push_back(coefficient.numerateur());
        if (coefficient.denominateur() != "1") denominateur.push_back(coefficient.denominateur());
    } else if (!coefficient.estUn()) {
        numerateur.push_back(texteNombre(os, coefficient));
    }
    for (const Facteur& f : facteurs) {
        if (exposantNegatif(f.exposant)) {
            denominateur.push_back(textePuissance(os, f.base, -f.exposant));
        } else {
            numerateur.push_back(textePuissance(os, f.base, f.exposant));
        }
    }
    os << (numerateur.empty() ? std::string("1") : joindre(numerateur));
    if (!denominateur.empty()) {
        os << '/' << (denominateur.size() == 1 ? denominateur[0] : "(" + joindre(denominateur) + ")");
    }
}

// Facteurs d'un terme de somme (dont le coefficient est porté par la somme)
std::vector<Facteur> facteursDe(const ExprPtr& e) {
    if (const Produit* p = comme<Produit>(e)) return p->getFacteurs();
    if (const Puissance* p = comme<Puissance>(e)) return {{p->getBase(), p->getExposant()}};
    return {{e, un()}};
}

// Degré en la variable, pour ordonner l'affichage des sommes (x^2 + 5*x + 6)
double degre(const ASTNode& e) {
    auto degreFacteur = [](const ExprPtr& base, const ExprPtr& exposant) {
        const Constante* c = comme<Constante>(exposant);
        return base->type() == TypeNoeud::Variable && c ? c->getValeurConstante() : 0.0;
    };
    switch (e.type()) {
        case TypeNoeud::Variable: return 1.0;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(e);
            return degreFacteur(p.getBase(), p.getExposant());
        }
        case TypeNoeud::Produit: {
            double d = 0.0;
            for (const Facteur& f : static_cast<const Produit&>(e).getFacteurs()) d += degreFacteur(f.base, f.exposant);
            return d;
        }
        default: return 0.0;
    }
}

// Monôme : produit de variables et de paramètres à exposants numériques
bool estMonome(const ExprPtr& e) {
    auto facteurSimple = [](const ExprPtr& base, const ExprPtr& exposant) {
        return (base->type() == TypeNoeud::Variable || base->type() == TypeNoeud::Parametre) &&
               exposant->type() == TypeNoeud::Constante;
    };
    switch (e->type()) {
        case TypeNoeud::Variable:
        case TypeNoeud::Parametre: return true;
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(*e);
            return facteurSimple(p.getBase(), p.getExposant());
        }
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(*e).getFacteurs()) {
                if (!facteurSimple(f.base, f.exposant)) return false;
            }
            return true;
        default: return false;
    }
}

void ecrireSomme(std::ostream& os, const Somme& s) {
    std::vector<const Terme*> termes;
    for (const Terme& t : s.getTermes()) termes.push_back(&t);
    // Les monômes d'abord, par degré décroissant (x^2 + 5*x + 6), puis les autres termes
    // dans l'ordre canonique (C1*exp(-x) + C2*x*exp(-x))
    const auto finMonomes = std::stable_partition(termes.begin(), termes.end(),
                                                  [](const Terme* t) { return estMonome(t->expression); });
    std::stable_sort(termes.begin(), finMonomes,
                     [](const Terme* a, const Terme* b) { return degre(*a->expression) > degre(*b->expression); });
    bool premier = true;
    auto signe = [&](const Nombre& c) {
        const bool negatif = c.signe() < 0;
        if (premier) {
            if (negatif) os << '-';
        } else {
            os << (negatif ? " - " : " + ");
        }
        premier = false;
        return negatif ? -c : c;
    };
    for (const Terme* t : termes) {
        const Nombre valeurAbsolue = signe(t->coefficient);
        ecrireProduit(os, valeurAbsolue, facteursDe(t->expression));
    }
    if (!s.getConstante().estZero()) os << texteNombre(os, signe(s.getConstante()));
}

void ecrireBrut(std::ostream& os, const ASTNode& n) {
    switch (n.type()) {
        case TypeNoeud::Constante: static_cast<const Constante&>(n).getNombre().afficher(os); break;
        case TypeNoeud::Variable: os << static_cast<const Variable&>(n).getNom(); break;
        case TypeNoeud::Parametre: os << static_cast<const Parametre&>(n).getNom(); break;
        case TypeNoeud::Somme: ecrireSomme(os, static_cast<const Somme&>(n)); break;
        case TypeNoeud::Produit: {
            const Produit& p = static_cast<const Produit&>(n);
            if (p.getCoefficient().signe() < 0) {
                os << '-';
                ecrireProduit(os, -p.getCoefficient(), p.getFacteurs());
            } else {
                ecrireProduit(os, p.getCoefficient(), p.getFacteurs());
            }
            break;
        }
        case TypeNoeud::Puissance: {
            const Puissance& p = static_cast<const Puissance&>(n);
            ecrireProduit(os, Nombre(1), {{p.getBase(), p.getExposant()}});
            break;
        }
        case TypeNoeud::Sinus: os << "sin("; break;
        case TypeNoeud::Cosinus: os << "cos("; break;
        case TypeNoeud::Tangente: os << "tan("; break;
        case TypeNoeud::Exponentielle: os << "exp("; break;
        case TypeNoeud::Logarithme: os << "ln("; break;
        case TypeNoeud::IntegraleNonEvaluee:
            os << "integrale(";
            ecrire(os, *static_cast<const IntegraleNonEvaluee&>(n).getIntegrande(), 0);
            os << ')';
            break;
        case TypeNoeud::LimiteNonEvaluee: {
            const LimiteNonEvaluee& l = static_cast<const LimiteNonEvaluee&>(n);
            os << "lim(x->" << l.getPoint() << ", ";
            ecrire(os, *l.getExpression(), 0);
            os << ')';
            break;
        }
    }
    if (const FonctionUnaire* f = comme<FonctionUnaire>(&n)) {
        ecrire(os, *f->m_argument, 0);
        os << ')';
    }
}

void ecrire(std::ostream& os, const ASTNode& n, int precedenceMin) {
    if (precedence(n) < precedenceMin) {
        os << '(';
        ecrireBrut(os, n);
        os << ')';
    } else {
        ecrireBrut(os, n);
    }
}

} // namespace

void ASTNode::afficher(std::ostream& os) const { ecrire(os, *this, 0); }

} // namespace symalgo
