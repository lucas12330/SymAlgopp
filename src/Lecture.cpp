/**
 * @file Lecture.cpp
 * @brief Analyse lexicale puis descente récursive (voir la grammaire dans Lecture.hpp).
 */

#include "Lecture.hpp"

#include <cstdint>
#include <cstring>
#include <limits>
#include <string_view>
#include <vector>

namespace symalgo {

namespace {

// Indice (à partir de 1) du caractère commençant à l'octet donné, en UTF-8
std::size_t indiceCaractere(const std::string& texte, std::size_t octet) {
    std::size_t indice = 1;
    for (std::size_t i = 0; i < octet && i < texte.size(); ++i) {
        if ((static_cast<unsigned char>(texte[i]) & 0xC0) != 0x80) ++indice;
    }
    return indice;
}

std::string message(const std::string& texte, std::size_t position, const std::string& raison) {
    return "lecture, position " + std::to_string(position) + " : " + raison + "\n  " + texte + "\n  " +
           std::string(position - 1, ' ') + "^";
}

enum class Jeton { Nombre, Identifiant, Plus, Moins, Fois, Divise, Puissance, Exposant, Ouvrante, Fermante, Egal, Fin };

// Un lexème désigne un morceau du texte (12 octets) : les nombres et les noms ne sont
// convertis qu'à l'analyse syntaxique
struct Lexeme {
    Jeton type;
    std::uint32_t octet;     // position dans le texte
    std::uint32_t longueur;
};

// Symboles Unicode reconnus (UTF-8) et leur jeton
struct SymboleUnicode {
    const char* octets;
    Jeton type;
};

constexpr const char* PI_UNICODE = "\xCF\x80";   // π
constexpr const char* CARRE = "\xC2\xB2";        // ²

constexpr SymboleUnicode SYMBOLES[] = {
    {"\xC3\x97", Jeton::Fois},          // ×
    {"\xC2\xB7", Jeton::Fois},          // ·
    {"\xC3\xB7", Jeton::Divise},        // ÷
    {"\xE2\x88\x92", Jeton::Moins},     // −
    {CARRE, Jeton::Exposant},
    {"\xC2\xB3", Jeton::Exposant},      // ³
    {PI_UNICODE, Jeton::Identifiant},
};

// Au-delà, 10^exposant n'a plus de sens pour une saisie et coûterait cher en mémoire
constexpr long long EXPOSANT_DECIMAL_MAX = 4096;

// Imbrication maximale (parenthèses, fonctions, exposants), pour ne jamais épuiser la pile
constexpr int PROFONDEUR_MAX = 500;

bool estChiffre(char c) { return c >= '0' && c <= '9'; }
bool debutIdentifiant(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_'; }
bool suiteIdentifiant(char c) { return debutIdentifiant(c) || estChiffre(c); }

class Lecteur {
public:
    Lecteur(const std::string& texte, const OptionsLecture& options) : m_texte(texte), m_options(options) {
        if (texte.size() >= std::numeric_limits<std::uint32_t>::max()) erreurA(0, "texte trop long");
        decouper();
    }

    EgaliteLue equation() {
        if (courant().type == Jeton::Fin) erreur("expression vide");
        EgaliteLue resultat{expression(), frac(0)};
        if (accepter(Jeton::Egal)) {
            if (courant().type == Jeton::Fin) erreur("expression attendue apres '='");
            resultat.droite = expression();
        }
        if (courant().type == Jeton::Egal) erreur("un seul '=' est permis");
        terminer();
        return resultat;
    }

    ExprPtr expressionSeule() {
        if (courant().type == Jeton::Fin) erreur("expression vide");
        ExprPtr e = expression();
        if (courant().type == Jeton::Egal) erreur("'=' inattendu (utiliser lireEquation)");
        terminer();
        return e;
    }

private:
    const std::string& m_texte;
    const OptionsLecture& m_options;
    std::vector<Lexeme> m_lexemes;
    std::size_t m_indice = 0;
    int m_profondeur = 0;

    [[noreturn]] void erreurA(std::size_t octet, const std::string& raison) const {
        throw ErreurLecture(m_texte, indiceCaractere(m_texte, octet), raison);
    }
    [[noreturn]] void erreur(const std::string& raison) const { erreurA(courant().octet, raison); }

    // ----- Analyse lexicale -----

    void ajouter(Jeton type, std::size_t debut, std::size_t fin) {
        m_lexemes.push_back({type, static_cast<std::uint32_t>(debut), static_cast<std::uint32_t>(fin - debut)});
    }

    void decouper() {
        std::size_t i = 0;
        const std::size_t n = m_texte.size();
        m_lexemes.reserve(n / 2 + 2);
        while (i < n) {
            const char c = m_texte[i];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                ++i;
                continue;
            }
            const std::size_t debut = i;
            if (estChiffre(c) || (c == '.' && i + 1 < n && estChiffre(m_texte[i + 1]))) {
                i = finNombre(i);
                ajouter(Jeton::Nombre, debut, i);
                continue;
            }
            if (debutIdentifiant(c)) {
                while (i < n && suiteIdentifiant(m_texte[i])) ++i;
                ajouter(Jeton::Identifiant, debut, i);
                continue;
            }
            Jeton type = Jeton::Fin;
            switch (c) {
                case '+': type = Jeton::Plus; break;
                case '-': type = Jeton::Moins; break;
                case '/': type = Jeton::Divise; break;
                case '^': type = Jeton::Puissance; break;
                case '(': type = Jeton::Ouvrante; break;
                case ')': type = Jeton::Fermante; break;
                case '=': type = Jeton::Egal; break;
                case '*':
                    type = Jeton::Fois;
                    if (i + 1 < n && m_texte[i + 1] == '*') { // ** (notation Python)
                        type = Jeton::Puissance;
                        ++i;
                    }
                    break;
                default: break;
            }
            if (type != Jeton::Fin) {
                ajouter(type, debut, ++i);
                continue;
            }
            if (!lireSymbole(i)) erreurA(debut, "caractere inattendu");
        }
        ajouter(Jeton::Fin, n, n);
    }

    bool lireSymbole(std::size_t& i) {
        for (const SymboleUnicode& s : SYMBOLES) {
            const std::size_t longueur = std::strlen(s.octets);
            if (m_texte.compare(i, longueur, s.octets) != 0) continue;
            ajouter(s.type, i, i + longueur);
            i += longueur;
            return true;
        }
        return false;
    }

    // Fin d'un nombre : chiffres [. chiffres] [e|E [+|-] chiffres]. « e » n'introduit un
    // exposant que s'il est suivi d'un chiffre (« 2e » vaut 2*e)
    std::size_t finNombre(std::size_t i) const {
        const std::size_t n = m_texte.size();
        while (i < n && estChiffre(m_texte[i])) ++i;
        if (i < n && m_texte[i] == '.') {
            ++i;
            while (i < n && estChiffre(m_texte[i])) ++i;
        }
        if (i < n && (m_texte[i] == 'e' || m_texte[i] == 'E')) {
            std::size_t j = i + 1;
            if (j < n && (m_texte[j] == '+' || m_texte[j] == '-')) ++j;
            if (j < n && estChiffre(m_texte[j])) {
                while (j < n && estChiffre(m_texte[j])) ++j;
                i = j;
            }
        }
        return i;
    }

    // Valeur exacte d'un lexème nombre : mantisse entière * 10^(exposant - décimales)
    Nombre valeurNombre(const Lexeme& l) const {
        std::size_t i = l.octet;
        const std::size_t fin = l.octet + l.longueur;
        std::string chiffres;
        long long decimales = 0;
        bool apresPoint = false;
        for (; i < fin && m_texte[i] != 'e' && m_texte[i] != 'E'; ++i) {
            if (m_texte[i] == '.') {
                apresPoint = true;
                continue;
            }
            chiffres += m_texte[i];
            if (apresPoint) ++decimales;
        }
        long long exposant = 0;
        if (i < fin) {
            const bool negatif = m_texte[++i] == '-';
            if (m_texte[i] == '+' || m_texte[i] == '-') ++i;
            for (; i < fin; ++i) {
                exposant = exposant * 10 + (m_texte[i] - '0');
                if (exposant > EXPOSANT_DECIMAL_MAX) erreurA(l.octet, "exposant decimal trop grand");
            }
            if (negatif) exposant = -exposant;
        }
        const long long puissance = exposant - decimales;
        if (puissance > EXPOSANT_DECIMAL_MAX || puissance < -EXPOSANT_DECIMAL_MAX) {
            erreurA(l.octet, "nombre trop long");
        }
        return Nombre::depuisTexte(chiffres) * Nombre(10).puissanceEntiere(puissance);
    }

    // ----- Analyse syntaxique -----

    const Lexeme& courant() const { return m_lexemes[m_indice]; }

    bool accepter(Jeton type) {
        if (courant().type != type) return false;
        ++m_indice;
        return true;
    }

    void attendre(Jeton type, const char* description) {
        if (!accepter(type)) erreur(std::string(description) + " attendue");
    }

    void terminer() const {
        if (courant().type == Jeton::Fermante) erreur("parenthese fermante sans ouvrante");
        if (courant().type != Jeton::Fin) erreur("operateur attendu");
    }

    // Garde de profondeur, pour les constructions récursives
    struct Profondeur {
        Lecteur& lecteur;
        explicit Profondeur(Lecteur& l) : lecteur(l) {
            if (++lecteur.m_profondeur > PROFONDEUR_MAX) lecteur.erreur("expression trop imbriquee");
        }
        ~Profondeur() { --lecteur.m_profondeur; }
        Profondeur(const Profondeur&) = delete;
        Profondeur& operator=(const Profondeur&) = delete;
    };

    ExprPtr expression() {
        ExprPtr resultat = terme();
        for (;;) {
            if (accepter(Jeton::Plus)) {
                resultat = resultat + terme();
            } else if (accepter(Jeton::Moins)) {
                resultat = resultat - terme();
            } else {
                return resultat;
            }
        }
    }

    // Les facteurs d'un terme sont rassemblés puis multipliés en une fois : le signe et les
    // coefficients portent sur tout le produit, comme dans la forme canonique. Construire
    // de gauche à droite développerait « -(a + b)*c » ou « 5*(x - 1)/x » dès le premier
    // facteur, et le texte affiché ne se relirait pas à l'identique.
    ExprPtr terme() {
        std::vector<ExprPtr> facteurs;
        if (signes()) facteurs.push_back(frac(-1));
        facteurs.push_back(puissance());
        for (;;) {
            if (accepter(Jeton::Fois)) {
                facteurs.push_back(unaire());
            } else if (accepter(Jeton::Divise)) {
                facteurs.push_back(ast_pow(unaire(), frac(-1)));
            } else if (courant().type == Jeton::Identifiant || courant().type == Jeton::Ouvrante) {
                facteurs.push_back(puissance()); // multiplication implicite
            } else {
                return facteurs.size() == 1 ? facteurs[0] : produit(facteurs);
            }
        }
    }

    // Consomme une suite de signes « + » et « - » (sans récursion) ; vrai si elle est négative
    bool signes() {
        bool negatif = false;
        for (;;) {
            if (accepter(Jeton::Moins)) {
                negatif = !negatif;
            } else if (!accepter(Jeton::Plus)) {
                return negatif;
            }
        }
    }

    ExprPtr unaire() {
        Profondeur garde(*this);
        const bool negatif = signes();
        ExprPtr e = puissance();
        return negatif ? -e : e;
    }

    ExprPtr puissance() {
        ExprPtr base = primaire();
        if (accepter(Jeton::Puissance)) return ast_pow(base, unaire()); // associative à droite
        if (courant().type == Jeton::Exposant) {
            const bool carre = m_texte.compare(courant().octet, courant().longueur, CARRE) == 0;
            ++m_indice;
            return ast_pow(base, frac(carre ? 2 : 3));
        }
        return base;
    }

    ExprPtr primaire() {
        const Lexeme& l = courant();
        switch (l.type) {
            case Jeton::Nombre: ++m_indice; return nombre(valeurNombre(l));
            case Jeton::Identifiant: return identifiant();
            case Jeton::Ouvrante: {
                Profondeur garde(*this);
                ++m_indice;
                if (courant().type == Jeton::Fermante) erreur("expression attendue entre les parentheses");
                ExprPtr e = expression();
                attendre(Jeton::Fermante, "parenthese fermante");
                return e;
            }
            case Jeton::Fin: erreur("expression incomplete");
            case Jeton::Fermante: erreur("expression attendue avant ')'");
            case Jeton::Egal: erreur("expression attendue avant '='");
            default: erreur("expression attendue");
        }
    }

    ExprPtr identifiant() {
        const std::size_t octet = courant().octet;
        const std::string_view nom(m_texte.data() + octet, courant().longueur);
        ++m_indice;
        if (nom == m_options.variable) return var(m_options.variable);
        for (const std::string& autre : m_options.autresVariables) {
            if (nom == autre) return var(autre);
        }

        struct Fonction {
            const char* nom;
            ExprPtr (*appliquer)(const ExprPtr&);
        };
        static const Fonction FONCTIONS[] = {
            {"sin", ast_sin},   {"cos", ast_cos},     {"tan", ast_tan},     {"exp", ast_exp},
            {"ln", ast_ln},     {"log", ast_ln},      {"asin", ast_asin},   {"acos", ast_acos},
            {"atan", ast_atan}, {"arcsin", ast_asin}, {"arccos", ast_acos}, {"arctan", ast_atan},
            {"sqrt", [](const ExprPtr& u) { return ast_pow(u, frac(1, 2)); }},
        };
        for (const Fonction& f : FONCTIONS) {
            if (nom != f.nom) continue;
            if (courant().type != Jeton::Ouvrante) erreur("parenthese ouvrante attendue apres " + std::string(nom));
            Profondeur garde(*this);
            ++m_indice;
            if (courant().type == Jeton::Fermante) erreur("argument attendu");
            ExprPtr argument = expression();
            attendre(Jeton::Fermante, "parenthese fermante");
            return f.appliquer(argument);
        }
        if (nom == "pi" || nom == PI_UNICODE) return pi();
        if (nom == "e") return ast_exp(frac(1));
        // Un paramètre suivi d'une parenthèse est presque toujours une fonction inconnue
        // (sinh(x), f(x)) : la lire comme un produit donnerait un résultat faux
        if (courant().type == Jeton::Ouvrante && courant().octet == octet + nom.size()) {
            erreurA(octet, "fonction inconnue : " + std::string(nom));
        }
        return param(std::string(nom));
    }
};

} // namespace

ErreurLecture::ErreurLecture(const std::string& texte, std::size_t position, const std::string& raison)
    : std::invalid_argument(message(texte, position, raison)), m_position(position), m_raison(raison) {}

ExprPtr lire(const std::string& texte, const OptionsLecture& options) {
    return Lecteur(texte, options).expressionSeule();
}

EgaliteLue lireEquation(const std::string& texte, const OptionsLecture& options) {
    return Lecteur(texte, options).equation();
}

} // namespace symalgo
