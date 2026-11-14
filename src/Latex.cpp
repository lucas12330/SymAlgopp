/**
 * @file Latex.cpp
 * @brief Rendu LaTeX des expressions, des solutions et du journal des étapes.
 *
 * Un seul parcours de l'expression, écrit directement à la fin d'une std::string : aucune
 * chaîne intermédiaire par sous-expression (les décisions de parenthèses ne regardent que
 * le noeud et ses enfants immédiats), donc un temps linéaire en la taille du texte produit.
 * La précédence est celle de Affichage.cpp (somme < produit < puissance < atome).
 */

#include "Latex.hpp"

#include "Canonique.hpp"
#include "Multivariable.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

namespace symalgo {

namespace {

enum Precedence { SOMME = 1, PRODUIT = 2, PUISSANCE = 3, ATOME = 4 };

// ----------------------------------------------------------------------------
// Nombres
// ----------------------------------------------------------------------------

// Écriture la plus courte d'un double positif fini, qui se relit au même double
struct TexteReel {
    char tampon[40];
    int taille = 0;
    int e = -1; // position du 'e' (notation scientifique), -1 sinon

    explicit TexteReel(double v) {
        const auto r = std::to_chars(tampon, tampon + sizeof(tampon) - 1, v);
        taille = static_cast<int>(r.ptr - tampon);
        tampon[taille] = '\0';
        for (int i = 0; i < taille; ++i) {
            if (tampon[i] == 'e') e = i;
        }
    }
    bool scientifique() const { return e >= 0; }
};

void ajouterEntier(std::string& s, long long v) {
    char t[24];
    const auto r = std::to_chars(t, t + sizeof(t), v);
    s.append(t, static_cast<std::size_t>(r.ptr - t));
}

// Entier exact positif (de n'importe quelle taille)
void ajouterEntierExact(std::string& s, const Nombre& n) {
    long long v;
    if (n.versEntier(v)) {
        ajouterEntier(s, v);
    } else {
        s += n.numerateur();
    }
}

// Réel positif ou nul, fini ou non : 1.5, 10^{-10}, 2.5 \times 10^{20}, \infty
void ajouterReel(std::string& s, double v) {
    if (std::isnan(v)) {
        s += "\\mathrm{NaN}";
        return;
    }
    if (std::isinf(v)) {
        s += "\\infty";
        return;
    }
    const TexteReel t(v);
    if (!t.scientifique()) {
        s.append(t.tampon, static_cast<std::size_t>(t.taille));
        return;
    }
    if (!(t.e == 1 && t.tampon[0] == '1')) {
        s.append(t.tampon, static_cast<std::size_t>(t.e));
        s += " \\times ";
    }
    s += "10^{";
    ajouterEntier(s, std::strtol(t.tampon + t.e + 1, nullptr, 10));
    s += '}';
}

// Valeur absolue d'un nombre : entier, fraction ou réel
void ajouterNombreAbsolu(std::string& s, const Nombre& n, bool fractions) {
    if (!n.estExact()) {
        ajouterReel(s, std::abs(n.versDouble()));
        return;
    }
    const Nombre a = n.valeurAbsolue();
    if (a.estEntier()) {
        ajouterEntierExact(s, a);
        return;
    }
    if (fractions) s += "\\frac{";
    ajouterEntierExact(s, a.numerateurNombre());
    s += fractions ? "}{" : "/";
    ajouterEntierExact(s, a.denominateurNombre());
    if (fractions) s += '}';
}

// ----------------------------------------------------------------------------
// Noms (variables, paramètres)
// ----------------------------------------------------------------------------

const char* const GRECS[] = {"Delta",   "Gamma",   "Lambda",     "Omega",   "Phi",     "Pi",      "Psi",    "Sigma",
                             "Theta",   "Upsilon", "Xi",         "alpha",   "beta",    "chi",     "delta",  "epsilon",
                             "eta",     "gamma",   "iota",       "kappa",   "lambda",  "mu",      "nu",     "omega",
                             "phi",     "pi",      "psi",        "rho",     "sigma",   "tau",     "theta",  "upsilon",
                             "varepsilon", "varphi", "vartheta",  "xi",      "zeta"};

bool estGrec(const char* debut, std::size_t n) {
    for (const char* g : GRECS) {
        if (std::strlen(g) == n && std::strncmp(g, debut, n) == 0) return true;
    }
    return false;
}

// Symboles unicode qui n'existent pas en mode texte LaTeX : écrits en mode mathématique
struct SymboleUnicode {
    const char* utf8;
    const char* latex;
};
const SymboleUnicode SYMBOLES[] = {
    {"\xCE\x94", "$\\Delta$"}, {"\xCF\x80", "$\\pi$"},     {"\xE2\x88\x92", "$-$"},          {"\xE2\x88\x9A", "$\\surd$"},
    {"\xE2\x89\xA5", "$\\geq$"}, {"\xE2\x89\xA4", "$\\leq$"}, {"\xC2\xB1", "$\\pm$"},        {"\xC2\xB2", "$^{2}$"},
    {"\xC2\xB3", "$^{3}$"},       {"\xC3\x97", "$\\times$"},  {"\xE2\x86\x92", "$\\rightarrow$"}, {"\xE2\x89\x88", "$\\approx$"},
    {"\xE2\x89\xA0", "$\\neq$"}, {"\xE2\x88\x88", "$\\in$"},  {"\xE2\x88\x9E", "$\\infty$"},   {"\xE2\x84\x9D", "$\\mathbb{R}$"},
    {"\xE2\x84\xA4", "$\\mathbb{Z}$"},
};

// Texte (mode texte) ou nom (mode mathématique) : caractères spéciaux de LaTeX échappés
void ajouterEchappe(std::string& s, const char* debut, std::size_t n, bool mathematique) {
    for (std::size_t i = 0; i < n; ++i) {
        const char c = debut[i];
        switch (c) {
            case '\\': s += mathematique ? "\\backslash " : "\\textbackslash{}"; break;
            case '&': case '%': case '$': case '#': case '_': case '{': case '}': s += '\\'; s += c; break;
            case '~': s += "\\textasciitilde{}"; break;
            case '^': s += "\\textasciicircum{}"; break;
            case '<': s += mathematique ? "<" : "$<$"; break;
            case '>': s += mathematique ? ">" : "$>$"; break;
            default:
                if (!mathematique && (static_cast<unsigned char>(c) & 0x80)) {
                    bool trouve = false;
                    for (const SymboleUnicode& u : SYMBOLES) {
                        const std::size_t t = std::strlen(u.utf8);
                        if (i + t <= n && std::memcmp(debut + i, u.utf8, t) == 0) {
                            s += u.latex;
                            i += t - 1;
                            trouve = true;
                            break;
                        }
                    }
                    if (trouve) break;
                }
                s += c;
        }
    }
}

// x -> x, alpha -> \alpha, theta1 -> \theta_{1}, C1 -> C_{1}, vitesse -> \mathrm{vitesse} ;
// variable interne (variable dont le nom commence par « _ », ex. « _t » du solveur) : _t -> t.
// Un paramètre « _a » garde son « _ » (\mathrm{\_a}) : il ne doit pas s'écrire comme « a ».
void ajouterNom(std::string& s, const std::string& nom, bool interne = false) {
    std::size_t i = 0;
    while (interne && i < nom.size() && nom[i] == '_') ++i;
    std::size_t j = i;
    while (j < nom.size() && std::isalpha(static_cast<unsigned char>(nom[j]))) ++j;
    if (j == i) { // pas de lettre ASCII après les « _ » (« _1 », « 2a », symbole unicode...) : nom entier, droit
        s += "\\mathrm{";
        ajouterEchappe(s, nom.c_str(), nom.size(), true);
        s += '}';
        return;
    }
    if (j - i == 1) {
        s += nom[i];
    } else if (estGrec(nom.c_str() + i, j - i)) {
        s += '\\';
        s.append(nom, i, j - i);
    } else {
        s += "\\mathrm{";
        s.append(nom, i, j - i);
        s += '}';
    }
    if (j < nom.size()) {
        std::size_t d = j;
        if (nom[d] == '_') ++d;
        if (d < nom.size()) {
            s += "_{";
            ajouterEchappe(s, nom.c_str() + d, nom.size() - d, true);
            s += '}';
        }
    }
}

// ----------------------------------------------------------------------------
// Expressions
// ----------------------------------------------------------------------------

bool exposantNegatif(const ASTNode* e) {
    const Constante* c = comme<Constante>(e);
    return c && c->getNombre().signe() < 0;
}

// Fonctions écrites avec l'exposant avant l'argument : \sin^{2}\left(x\right)
bool fonctionAExposant(TypeNoeud t) {
    switch (t) {
        case TypeNoeud::Sinus:
        case TypeNoeud::Cosinus:
        case TypeNoeud::Tangente:
        case TypeNoeud::Logarithme:
        case TypeNoeud::SinusHyperbolique:
        case TypeNoeud::CosinusHyperbolique:
        case TypeNoeud::TangenteHyperbolique: return true;
        default: return false;
    }
}

const char* nomFonction(TypeNoeud t) {
    switch (t) {
        case TypeNoeud::Sinus: return "\\sin";
        case TypeNoeud::Cosinus: return "\\cos";
        case TypeNoeud::Tangente: return "\\tan";
        case TypeNoeud::Logarithme: return "\\ln";
        case TypeNoeud::ArcSinus: return "\\arcsin";
        case TypeNoeud::ArcCosinus: return "\\arccos";
        case TypeNoeud::ArcTangente: return "\\arctan";
        case TypeNoeud::Signe: return "\\operatorname{sgn}";
        case TypeNoeud::SinusHyperbolique: return "\\sinh";
        case TypeNoeud::CosinusHyperbolique: return "\\cosh";
        case TypeNoeud::TangenteHyperbolique: return "\\tanh";
        default: return "";
    }
}

class Ecrivain {
public:
    Ecrivain(std::string& s, const OptionsLatex& o) : m_s(s), m_o(o) {}

    void ecrire(const ASTNode& n) { brut(n); }

private:
    std::string& m_s;
    const OptionsLatex& m_o;

    // Pour chaque facteur base^exposant de e (exposant nul : 1)
    template <class Fonction>
    static void chaqueFacteur(const ASTNode& e, Fonction f) {
        if (const Produit* p = comme<Produit>(&e)) {
            for (const Facteur& g : p->getFacteurs()) f(*g.base, g.exposant.get());
        } else if (const Puissance* p = comme<Puissance>(&e)) {
            f(*p->getBase(), p->getExposant().get());
        } else {
            f(e, static_cast<const ASTNode*>(nullptr));
        }
    }

    int precedence(const ASTNode& n) const {
        switch (n.type()) {
            case TypeNoeud::Somme: return SOMME;
            case TypeNoeud::Produit: return PRODUIT;
            case TypeNoeud::Constante: {
                const Nombre& v = static_cast<const Constante&>(n).getNombre();
                if (v.signe() < 0) return PRODUIT;
                if (v.estExact()) return v.estEntier() ? ATOME : PRODUIT;
                const double d = v.versDouble();
                return std::isfinite(d) && TexteReel(d).scientifique() ? PRODUIT : ATOME;
            }
            case TypeNoeud::Puissance:
                return exposantNegatif(static_cast<const Puissance&>(n).getExposant().get()) ? PRODUIT : PUISSANCE;
            case TypeNoeud::Exponentielle: return PUISSANCE; // e^{u}
            case TypeNoeud::IntegraleNonEvaluee:
            case TypeNoeud::LimiteNonEvaluee: return PRODUIT;
            default: return ATOME;
        }
    }

    void noeud(const ASTNode& n, int precedenceMin) {
        if (precedence(n) < precedenceMin) {
            m_s += "\\left(";
            brut(n);
            m_s += "\\right)";
        } else {
            brut(n);
        }
    }

    // Exposant n (positif) : \sqrt{b}, \sqrt[n]{b}
    bool racine(const ASTNode& base, const Nombre& n) {
        long long num, den;
        if (!m_o.racines || !n.versFraction(num, den) || num != 1 || den < 2) return false;
        if (den == 2) {
            m_s += "\\sqrt{";
        } else {
            m_s += "\\sqrt[";
            ajouterEntier(m_s, den);
            m_s += "]{";
        }
        brut(base);
        m_s += '}';
        return true;
    }

    // Le facteur base^exposant (inverse : exposant changé de signe, pour un dénominateur)
    // seul : unique élément d'un numérateur ou d'un dénominateur de \frac, dont les accolades
    // délimitent déjà une somme
    void puissance(const ASTNode& base, const ASTNode* exposant, bool inverse, bool seul) {
        if (!exposant) {
            noeud(base, seul ? 0 : PUISSANCE);
            return;
        }
        if (const Constante* c = comme<Constante>(exposant)) {
            const Nombre n = inverse ? -c->getNombre() : c->getNombre();
            if (n.estExact() && n.estUn()) {
                noeud(base, seul ? 0 : PUISSANCE);
                return;
            }
            if (racine(base, n)) return;
            long long k;
            if (fonctionAExposant(base.type()) && n.versEntier(k) && k >= 2) {
                m_s += nomFonction(base.type());
                m_s += "^{";
                ajouterEntier(m_s, k);
                m_s += "}\\left(";
                brut(*static_cast<const FonctionUnaire&>(base).m_argument);
                m_s += "\\right)";
                return;
            }
            noeud(base, ATOME);
            m_s += "^{";
            ajouterNombreAbsolu(m_s, n, m_o.fractions);
            m_s += '}';
            return;
        }
        noeud(base, ATOME);
        m_s += "^{";
        brut(*exposant);
        m_s += '}';
    }

    // Vrai si le facteur s'écrit en commençant par un chiffre (« 2 \cdot 3^{x} » et non « 2 3^{x} »)
    bool commenceParChiffre(const ASTNode& base, const ASTNode* exposant, bool inverse) const {
        if (base.type() != TypeNoeud::Constante || precedence(base) != ATOME) return false;
        const Constante* c = exposant ? comme<Constante>(exposant) : nullptr;
        if (!c) return true;
        const Nombre n = inverse ? -c->getNombre() : c->getNombre();
        long long num, den;
        return !(m_o.racines && n.versFraction(num, den) && num == 1 && den >= 2);
    }

    // coefficient (positif ou nul) * produit des facteurs de corps, regroupé en une fraction
    void produit(const Nombre& coefficient, const ASTNode& corps) {
        bool numerateurCoef, denominateurCoef;
        if (coefficient.estExact()) {
            denominateurCoef = !coefficient.estEntier();
            numerateurCoef = denominateurCoef ? !coefficient.numerateurNombre().estUn() : !coefficient.estUn();
        } else {
            denominateurCoef = false;
            numerateurCoef = !coefficient.estUn();
        }
        int nombreNumerateur = numerateurCoef ? 1 : 0, nombreDenominateur = denominateurCoef ? 1 : 0;
        chaqueFacteur(corps, [&](const ASTNode&, const ASTNode* e) {
            if (exposantNegatif(e)) ++nombreDenominateur; else ++nombreNumerateur;
        });
        if (nombreDenominateur == 0) {
            liste(coefficient, numerateurCoef, false, corps, false, false);
            if (nombreNumerateur == 0) m_s += '1';
            return;
        }
        if (m_o.fractions) m_s += "\\frac{";
        liste(coefficient, numerateurCoef, false, corps, false, m_o.fractions && nombreNumerateur == 1);
        if (nombreNumerateur == 0) m_s += '1';
        m_s += m_o.fractions ? "}{" : "/";
        const bool parentheses = !m_o.fractions && nombreDenominateur > 1;
        if (parentheses) m_s += "\\left(";
        liste(coefficient, denominateurCoef, true, corps, true, m_o.fractions && nombreDenominateur == 1);
        if (parentheses) m_s += "\\right)";
        if (m_o.fractions) m_s += '}';
    }

    // Numérateur ou dénominateur : nombre éventuel puis facteurs de même signe d'exposant
    void liste(const Nombre& coefficient, bool avecNombre, bool nombreEstDenominateur, const ASTNode& corps,
               bool denominateur, bool seul) {
        bool premier = true;
        if (avecNombre) {
            if (nombreEstDenominateur) {
                ajouterEntierExact(m_s, coefficient.denominateurNombre());
            } else if (coefficient.estExact()) {
                ajouterEntierExact(m_s, coefficient.numerateurNombre());
            } else {
                ajouterReel(m_s, coefficient.versDouble());
            }
            premier = false;
        }
        // pi en dernier (2 k \pi) : deux passes, les autres facteurs puis pi
        for (int passe = 0; passe < 2; ++passe) {
            chaqueFacteur(corps, [&](const ASTNode& base, const ASTNode* e) {
                if (exposantNegatif(e) != denominateur || (base.type() == TypeNoeud::Pi) != (passe == 1)) return;
                if (!premier) m_s += m_o.multiplicationExplicite || commenceParChiffre(base, e, denominateur) ? " \\cdot " : " ";
                premier = false;
                puissance(base, e, denominateur, seul && !avecNombre);
            });
        }
    }

    static double degre(const ASTNode& e) {
        auto degreFacteur = [](const ASTNode& base, const ASTNode* exposant) {
            const Constante* c = comme<Constante>(exposant);
            return base.type() == TypeNoeud::Variable && c ? c->getNombre().versDouble() : 0.0;
        };
        double d = 0.0;
        if (e.type() == TypeNoeud::Variable) return 1.0;
        if (e.type() == TypeNoeud::Puissance || e.type() == TypeNoeud::Produit) {
            chaqueFacteur(e, [&](const ASTNode& base, const ASTNode* ex) { d += degreFacteur(base, ex); });
        }
        return d;
    }

    // Produit de variables et de paramètres à exposants numériques
    static bool estMonome(const ASTNode& e) {
        auto simple = [](const ASTNode& base, const ASTNode* exposant) {
            return (base.type() == TypeNoeud::Variable || base.type() == TypeNoeud::Parametre) &&
                   (!exposant || exposant->type() == TypeNoeud::Constante);
        };
        switch (e.type()) {
            case TypeNoeud::Variable:
            case TypeNoeud::Parametre: return true;
            case TypeNoeud::Puissance:
            case TypeNoeud::Produit: {
                bool oui = true;
                chaqueFacteur(e, [&](const ASTNode& base, const ASTNode* ex) { oui = oui && simple(base, ex); });
                return oui;
            }
            default: return false;
        }
    }

    struct TermeOrdonne {
        const Terme* terme;
        double degre;
        bool monome;
    };

    // x^2 + 5 x + 6 : monômes par degré décroissant, puis les autres termes dans l'ordre canonique
    void somme(const Somme& s) {
        const std::vector<Terme>& termes = s.getTermes();
        TermeOrdonne local[24];
        std::vector<TermeOrdonne> grand;
        TermeOrdonne* t = local;
        if (termes.size() > 24) {
            grand.resize(termes.size());
            t = grand.data();
        }
        std::size_t n = 0;
        for (const Terme& terme : termes) {
            const bool monome = estMonome(*terme.expression);
            t[n++] = {&terme, monome ? degre(*terme.expression) : 0.0, monome};
        }
        TermeOrdonne* finMonomes = std::stable_partition(t, t + n, [](const TermeOrdonne& x) { return x.monome; });
        std::stable_sort(t, finMonomes, [](const TermeOrdonne& a, const TermeOrdonne& b) { return a.degre > b.degre; });
        bool premier = true;
        auto signe = [&](bool negatif) {
            if (premier) {
                if (negatif) m_s += '-';
            } else {
                m_s += negatif ? " - " : " + ";
            }
            premier = false;
        };
        for (std::size_t i = 0; i < n; ++i) {
            const Nombre& c = t[i].terme->coefficient;
            signe(c.signe() < 0);
            produit(c.signe() < 0 ? -c : c, *t[i].terme->expression);
        }
        if (!s.getConstante().estZero()) {
            signe(s.getConstante().signe() < 0);
            ajouterNombreAbsolu(m_s, s.getConstante(), m_o.fractions);
        }
    }

    void brut(const ASTNode& n) {
        switch (n.type()) {
            case TypeNoeud::Constante: {
                const Nombre& v = static_cast<const Constante&>(n).getNombre();
                if (v.signe() < 0) m_s += '-';
                ajouterNombreAbsolu(m_s, v, m_o.fractions);
                break;
            }
            case TypeNoeud::Variable: ajouterNom(m_s, static_cast<const Variable&>(n).getNom(), true); break;
            case TypeNoeud::Parametre: ajouterNom(m_s, static_cast<const Parametre&>(n).getNom()); break;
            case TypeNoeud::Pi: m_s += "\\pi"; break;
            case TypeNoeud::Somme: somme(static_cast<const Somme&>(n)); break;
            case TypeNoeud::Produit: {
                const Nombre& c = static_cast<const Produit&>(n).getCoefficient();
                if (c.signe() < 0) m_s += '-';
                produit(c.signe() < 0 ? -c : c, n);
                break;
            }
            case TypeNoeud::Puissance: produit(Nombre(1), n); break;
            case TypeNoeud::Exponentielle:
                m_s += "e^{";
                brut(*static_cast<const FonctionUnaire&>(n).m_argument);
                m_s += '}';
                break;
            case TypeNoeud::ValeurAbsolue:
                m_s += "\\left|";
                brut(*static_cast<const FonctionUnaire&>(n).m_argument);
                m_s += "\\right|";
                break;
            case TypeNoeud::IntegraleNonEvaluee: {
                const ExprPtr& f = static_cast<const IntegraleNonEvaluee&>(n).getIntegrande();
                m_s += "\\int ";
                noeud(*f, PRODUIT);
                m_s += "\\,d";
                const std::vector<ExprPtr> v = variables(f);
                if (v.empty()) {
                    m_s += 'x';
                } else {
                    brut(*v[0]);
                }
                break;
            }
            case TypeNoeud::LimiteNonEvaluee: {
                const LimiteNonEvaluee& l = static_cast<const LimiteNonEvaluee&>(n);
                const std::vector<ExprPtr> v = variables(l.getExpression());
                m_s += "\\lim_{";
                if (v.empty()) {
                    m_s += 'x';
                } else {
                    brut(*v[0]);
                }
                m_s += " \\to ";
                if (l.getPoint() < 0.0) m_s += '-';
                ajouterReel(m_s, std::abs(l.getPoint()));
                m_s += "} ";
                noeud(*l.getExpression(), PRODUIT);
                break;
            }
            default: { // fonctions usuelles
                m_s += nomFonction(n.type());
                m_s += "\\left(";
                brut(*static_cast<const FonctionUnaire&>(n).m_argument);
                m_s += "\\right)";
            }
        }
    }
};

// ----------------------------------------------------------------------------
// Solutions
// ----------------------------------------------------------------------------

void ajouterApproximation(std::string& s, double v) {
    if (v < 0.0) s += '-';
    const double a = std::abs(v);
    if (!std::isfinite(a)) {
        ajouterReel(s, a);
        return;
    }
    char t[40];
    std::snprintf(t, sizeof(t), "%.5g", a); // 5 chiffres significatifs, zéros finaux omis
    const double arrondi = std::strtod(t, nullptr);
    ajouterReel(s, arrondi);
}

} // namespace

std::string latex(const ExprPtr& e, const OptionsLatex& options) {
    if (!e) throw std::invalid_argument("latex : expression nulle");
    std::string s;
    s.reserve(128);
    Ecrivain(s, options).ecrire(*e);
    return s;
}

std::string latex(const Solutions& solutions, const std::string& inconnue, const OptionsLatex& options) {
    std::string x;
    ajouterNom(x, inconnue, true);
    std::string s;
    s.reserve(64 + 48 * solutions.liste.size());
    if (solutions.toutReel) return x + " \\in \\mathbb{R}";
    if (solutions.liste.empty()) {
        s = solutions.complet ? x + " \\in \\emptyset" : "\\text{aucune solution trouv\xC3\xA9" "e (r\xC3\xA9sultat incomplet)}";
        return s;
    }
    bool premier = true;
    for (const Solution& sol : solutions.liste) {
        if (!premier) s += " \\text{ ou } ";
        premier = false;
        s += x;
        if (!sol.exacte && !sol.estFamille()) {
            s += " \\approx ";
            ajouterApproximation(s, sol.approximation);
        } else {
            s += " = ";
            Ecrivain(s, options).ecrire(*sol.valeur);
        }
        if (sol.multiplicite == 2) {
            s += " \\; (\\text{double})";
        } else if (sol.multiplicite == 3) {
            s += " \\; (\\text{triple})";
        } else if (sol.multiplicite > 3) {
            s += " \\; (\\text{multiplicit\xC3\xA9 } ";
            ajouterEntier(s, sol.multiplicite);
            s += ')';
        }
        if (sol.estFamille()) {
            s += ",\\ ";
            for (std::size_t i = 0; i < sol.entiers.size(); ++i) {
                if (i) s += ", ";
                ajouterNom(s, sol.entiers[i]);
            }
            s += " \\in \\mathbb{Z}";
        }
    }
    if (!solutions.complet) s += " \\text{ ou d'autres solutions non d\xC3\xA9termin\xC3\xA9" "es}";
    return s;
}

std::string latex(const Journal& journal, const OptionsLatex& options) {
    std::string s;
    s.reserve(256 + 192 * journal.etapes().size());
    constexpr int NIVEAUX_MAX = 4; // profondeur maximale des listes imbriquées de LaTeX
    int ouvertes = 0;
    for (const Etape& etape : journal.etapes()) {
        const int voulues = std::min(etape.niveau, NIVEAUX_MAX - 1) + 1;
        while (ouvertes < voulues) {
            s.append(static_cast<std::size_t>(2 * ouvertes), ' ');
            s += "\\begin{itemize}\n";
            ++ouvertes;
        }
        while (ouvertes > voulues) {
            --ouvertes;
            s.append(static_cast<std::size_t>(2 * ouvertes), ' ');
            s += "\\end{itemize}\n";
        }
        s.append(static_cast<std::size_t>(2 * ouvertes), ' ');
        s += "\\item ";
        ajouterEchappe(s, etape.description.data(), etape.description.size(), false);
        if (etape.avant || etape.apres) {
            const char fin = etape.description.empty() ? ' ' : etape.description.back();
            s += fin == ':' ? " " : " : ";
            s += "$\\displaystyle ";
            Ecrivain ecrivain(s, options);
            if (etape.avant) ecrivain.ecrire(*etape.avant);
            if (etape.avant && etape.apres) s += etape.egalite ? (etape.approchee ? " \\approx " : " = ") : " \\;\\longrightarrow\\; ";
            if (etape.apres) ecrivain.ecrire(*etape.apres);
            s += '$';
        }
        s += '\n';
    }
    while (ouvertes > 0) {
        --ouvertes;
        s.append(static_cast<std::size_t>(2 * ouvertes), ' ');
        s += "\\end{itemize}\n";
    }
    return s;
}

std::string documentLatex(const std::string& corps) {
    std::string s;
    s.reserve(corps.size() + 256);
    s += "\\documentclass{article}\n"
         "\\usepackage[utf8]{inputenc}\n"
         "\\usepackage[T1]{fontenc}\n"
         "\\usepackage{amsmath, amssymb}\n"
         "\\begin{document}\n";
    s += corps;
    if (corps.empty() || corps.back() != '\n') s += '\n';
    s += "\\end{document}\n";
    return s;
}

} // namespace symalgo
