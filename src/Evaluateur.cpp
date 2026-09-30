/**
 * @file Evaluateur.cpp
 * @brief Compilation des expressions en programme linéaire et interprétation (point par
 *        point ou par blocs vectorisés).
 */

#include "Evaluateur.hpp"
#include "Canonique.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace symalgo {

/*
 * Nom : CompilateurExpression
 * Description : Traduit le graphe de l'expression en instructions (une valeur par noeud
 *               distinct, les noeuds partagés n'étant compilés qu'une fois), puis attribue
 *               les registres en réutilisant ceux des valeurs qui ne servent plus.
 */
class CompilateurExpression {
public:
    explicit CompilateurExpression(ProgrammeEvaluation& programme) : m_programme(programme) {}

    // Une seule variable, quel que soit son nom (les expressions à plusieurs variables sont refusées)
    void compiler(const ExprPtr& racine) {
        m_uneVariable = true;
        compiler(std::vector<ExprPtr>{racine});
    }

    // Les variables sont les entrées, numérotées dans l'ordre de la liste
    void compiler(const std::vector<ExprPtr>& racines, const std::vector<ExprPtr>& entrees) {
        for (std::size_t i = 0; i < entrees.size(); ++i) {
            const Variable* v = comme<Variable>(entrees[i]);
            if (!v) throw std::invalid_argument("ProgrammeEvaluation : une entree n'est pas une variable");
            if (!m_entrees.emplace(v->identifiantVariable(), i).second) {
                throw std::invalid_argument("ProgrammeEvaluation : variable '" + v->getNom() + "' en double");
            }
        }
        m_programme.m_nombreEntrees = static_cast<std::uint32_t>(entrees.size());
        compiler(racines);
    }

private:
    void compiler(const std::vector<ExprPtr>& racines) {
        try {
            std::vector<std::uint32_t> resultats;
            for (const ExprPtr& racine : racines) {
                if (m_uneVariable && racine->plusieursVariables()) {
                    throw std::logic_error("Impossible d'evaluer en un seul reel l'expression '" + racine->texte() +
                                           "' : elle a plusieurs variables");
                }
                resultats.push_back(valeur(*racine));
            }
            allouerRegistres(resultats);
            std::unordered_map<const ASTNode*, double> visites;
            double total = 0.0;
            for (const ExprPtr& racine : racines) total += nombreVisites(*racine, visites);
            m_programme.m_gainPartage = total / static_cast<double>(m_ssa.size());
        } catch (const std::logic_error& e) {
            m_programme.m_instructions.clear();
            m_programme.m_erreur = e.what();
        }
    }

    using Instruction = ProgrammeEvaluation::Instruction;
    using Code = ProgrammeEvaluation::Code;

    // Émet une instruction ; sa destination est son numéro de valeur (forme SSA)
    std::uint32_t emettre(Instruction ins) {
        ins.destination = static_cast<std::uint32_t>(m_ssa.size());
        m_ssa.push_back(ins);
        return ins.destination;
    }

    std::uint32_t constante(double v) {
        Instruction ins{Code::Constante};
        ins.valeur = v;
        return emettre(ins);
    }

    std::uint32_t facteur(const ExprPtr& base, const ExprPtr& exposant) {
        const std::uint32_t b = valeur(*base);
        long long k;
        if (const Constante* e = comme<Constante>(exposant); e && e->getNombre().versEntier(k)) {
            if (k == 1) return b;
            Instruction ins{Code::PuissanceEntiere};
            ins.a = b;
            ins.entier = k;
            return emettre(ins);
        }
        Instruction ins{Code::Puissance};
        ins.a = b;
        ins.b = valeur(*exposant);
        return emettre(ins);
    }

    std::uint32_t valeur(const ASTNode& n) {
        const auto it = m_memo.find(&n);
        if (it != m_memo.end()) return it->second;
        std::uint32_t v = 0;
        switch (n.type()) {
            case TypeNoeud::Constante: v = constante(n.getValeurConstante()); break;
            case TypeNoeud::Variable: {
                Instruction ins{Code::Variable};
                if (!m_uneVariable) {
                    const auto it = m_entrees.find(n.identifiantVariable());
                    if (it == m_entrees.end()) {
                        throw std::logic_error("Impossible d'evaluer la variable '" + static_cast<const Variable&>(n).getNom() +
                                               "' : elle n'a pas de valeur");
                    }
                    ins.entier = static_cast<long long>(it->second);
                }
                v = emettre(ins);
                break;
            }
            case TypeNoeud::Somme: {
                const Somme& s = static_cast<const Somme&>(n);
                bool premier = s.getConstante().estZero();
                if (!premier) v = constante(s.getConstante().versDouble());
                for (const Terme& t : s.getTermes()) {
                    const std::uint32_t tv = valeur(*t.expression);
                    Instruction ins{premier ? Code::Echelle : Code::Axpy};
                    ins.valeur = t.coefficient.versDouble();
                    if (premier) {
                        ins.a = tv; // v = c * t
                        v = t.coefficient.estUn() ? tv : emettre(ins);
                        premier = false;
                    } else {
                        ins.a = v; // v = v + c * t
                        ins.b = tv;
                        v = emettre(ins);
                    }
                }
                break;
            }
            case TypeNoeud::Produit: {
                const Produit& p = static_cast<const Produit&>(n);
                bool premier = true;
                for (const Facteur& f : p.getFacteurs()) {
                    const std::uint32_t fv = facteur(f.base, f.exposant);
                    if (premier) {
                        v = fv;
                        premier = false;
                    } else {
                        Instruction ins{Code::Produit};
                        ins.a = v;
                        ins.b = fv;
                        v = emettre(ins);
                    }
                }
                if (!p.getCoefficient().estUn()) {
                    Instruction ins{Code::Echelle};
                    ins.a = v;
                    ins.valeur = p.getCoefficient().versDouble();
                    v = emettre(ins);
                }
                break;
            }
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(n);
                v = facteur(p.getBase(), p.getExposant());
                break;
            }
            case TypeNoeud::Sinus: v = unaire(Code::Sinus, n); break;
            case TypeNoeud::Cosinus: v = unaire(Code::Cosinus, n); break;
            case TypeNoeud::Tangente: v = unaire(Code::Tangente, n); break;
            case TypeNoeud::Exponentielle: v = unaire(Code::Exponentielle, n); break;
            case TypeNoeud::Logarithme: v = unaire(Code::Logarithme, n); break;
            case TypeNoeud::ArcSinus: v = unaire(Code::ArcSinus, n); break;
            case TypeNoeud::ArcCosinus: v = unaire(Code::ArcCosinus, n); break;
            case TypeNoeud::ArcTangente: v = unaire(Code::ArcTangente, n); break;
            case TypeNoeud::Pi: v = constante(n.eval(0.0)); break;
            default:
                // Paramètre ou noeud non évalué : même erreur que l'évaluation directe
                n.eval(0.0);
                throw std::logic_error("Expression non evaluable");
        }
        m_memo.emplace(&n, v);
        return v;
    }

    std::uint32_t unaire(Code code, const ASTNode& n) {
        Instruction ins{code};
        ins.a = valeur(*static_cast<const FonctionUnaire&>(n).m_argument);
        return emettre(ins);
    }

    // Noeuds visités par l'évaluation récursive de l'arbre (partage compté à chaque usage)
    static double nombreVisites(const ASTNode& n, std::unordered_map<const ASTNode*, double>& memo) {
        const auto it = memo.find(&n);
        if (it != memo.end()) return it->second;
        double v = 1.0;
        switch (n.type()) {
            case TypeNoeud::Somme:
                for (const Terme& t : static_cast<const Somme&>(n).getTermes()) v += nombreVisites(*t.expression, memo);
                break;
            case TypeNoeud::Produit:
                for (const Facteur& f : static_cast<const Produit&>(n).getFacteurs()) {
                    v += nombreVisites(*f.base, memo) + visitesExposant(*f.exposant, memo);
                }
                break;
            case TypeNoeud::Puissance: {
                const Puissance& p = static_cast<const Puissance&>(n);
                v += nombreVisites(*p.getBase(), memo) + visitesExposant(*p.getExposant(), memo);
                break;
            }
            default:
                if (const FonctionUnaire* f = comme<FonctionUnaire>(&n)) v += nombreVisites(*f->m_argument, memo);
        }
        memo.emplace(&n, v);
        return v;
    }

    // L'arbre lit directement un exposant constant, sans le visiter
    static double visitesExposant(const ASTNode& e, std::unordered_map<const ASTNode*, double>& memo) {
        return e.type() == TypeNoeud::Constante ? 0.0 : nombreVisites(e, memo);
    }

    static int nombreOperandes(Code c) {
        switch (c) {
            case Code::Constante:
            case Code::Variable: return 0;
            case Code::Axpy:
            case Code::Produit:
            case Code::Puissance: return 2;
            default: return 1;
        }
    }

    /*
     * Nom : allouerRegistres
     * Description : Remplace les numéros de valeurs par des registres : un registre est
     *               libéré après la dernière utilisation de sa valeur et réutilisé. Écrire
     *               dans le registre d'un opérande qui meurt est sûr, car chaque instruction
     *               lit ses opérandes avant d'écrire (y compris élément par élément).
     */
    void allouerRegistres(const std::vector<std::uint32_t>& resultats) {
        const std::size_t n = m_ssa.size();
        std::vector<std::size_t> derniereUtilisation(n, 0);
        for (std::size_t i = 0; i < n; ++i) {
            const int k = nombreOperandes(m_ssa[i].code);
            if (k >= 1) derniereUtilisation[m_ssa[i].a] = i;
            if (k >= 2) derniereUtilisation[m_ssa[i].b] = i;
        }
        for (const std::uint32_t resultat : resultats) derniereUtilisation[resultat] = std::numeric_limits<std::size_t>::max();

        std::vector<std::uint32_t> registre(n, 0), libres;
        std::uint32_t nombreRegistres = 0;
        for (std::size_t i = 0; i < n; ++i) {
            Instruction ins = m_ssa[i];
            const int k = nombreOperandes(ins.code);
            const std::uint32_t va = ins.a, vb = ins.b;
            if (k >= 1) ins.a = registre[va];
            if (k >= 2) ins.b = registre[vb];
            if (k >= 1 && derniereUtilisation[va] == i) libres.push_back(registre[va]);
            if (k >= 2 && vb != va && derniereUtilisation[vb] == i) libres.push_back(registre[vb]);
            if (libres.empty()) {
                ins.destination = nombreRegistres++;
            } else {
                ins.destination = libres.back();
                libres.pop_back();
            }
            registre[i] = ins.destination;
            m_programme.m_instructions.push_back(ins);
        }
        for (const std::uint32_t resultat : resultats) m_programme.m_resultats.push_back(registre[resultat]);
        m_programme.m_nombreRegistres = nombreRegistres;
    }

    ProgrammeEvaluation& m_programme;
    bool m_uneVariable = false; // toute variable est l'entrée 0
    std::unordered_map<std::uint16_t, std::size_t> m_entrees; // identifiant de variable -> entrée
    std::vector<Instruction> m_ssa;
    std::unordered_map<const ASTNode*, std::uint32_t> m_memo;
};

ProgrammeEvaluation::ProgrammeEvaluation(const ExprPtr& expression) {
    CompilateurExpression(*this).compiler(expression);
}

ProgrammeEvaluation::ProgrammeEvaluation(const std::vector<ExprPtr>& expressions, const std::vector<ExprPtr>& entrees) {
    CompilateurExpression(*this).compiler(expressions, entrees);
}

void ProgrammeEvaluation::verifierValidite() const {
    if (!m_erreur.empty()) throw std::logic_error(m_erreur);
}

void ProgrammeEvaluation::executer(const double* entrees, double* r) const {
    for (const Instruction& ins : m_instructions) {
        double& d = r[ins.destination];
        switch (ins.code) {
            case Code::Constante: d = ins.valeur; break;
            case Code::Variable: d = entrees[ins.entier]; break;
            case Code::Axpy: d = r[ins.a] + ins.valeur * r[ins.b]; break;
            case Code::Echelle: d = ins.valeur * r[ins.a]; break;
            case Code::Produit: d = r[ins.a] * r[ins.b]; break;
            case Code::PuissanceEntiere: d = puissanceEntiereReelle(r[ins.a], ins.entier); break;
            case Code::Puissance: d = std::pow(r[ins.a], r[ins.b]); break;
            case Code::Sinus: d = std::sin(r[ins.a]); break;
            case Code::Cosinus: d = std::cos(r[ins.a]); break;
            case Code::Tangente: d = std::tan(r[ins.a]); break;
            case Code::Exponentielle: d = std::exp(r[ins.a]); break;
            case Code::Logarithme: d = std::log(r[ins.a]); break;
            case Code::ArcSinus: d = std::asin(r[ins.a]); break;
            case Code::ArcCosinus: d = std::acos(r[ins.a]); break;
            case Code::ArcTangente: d = std::atan(r[ins.a]); break;
        }
    }
}

double ProgrammeEvaluation::evaluer(double x) const {
    verifierValidite();
    if (m_nombreEntrees != 1 || m_resultats.size() != 1) {
        throw std::logic_error("evaluer(x) : le programme a plusieurs entrees ou plusieurs sorties, utiliser evaluerEn");
    }
    constexpr std::size_t REGISTRES_SUR_PILE = 64;
    double pile[REGISTRES_SUR_PILE] = {}; // toujours écrit avant lecture (ordre des instructions)
    std::vector<double> tas;
    double* r = pile;
    if (m_nombreRegistres > REGISTRES_SUR_PILE) {
        tas.resize(m_nombreRegistres);
        r = tas.data();
    }
    executer(&x, r);
    return r[m_resultats[0]];
}

void ProgrammeEvaluation::evaluerEn(const double* valeurs, double* sorties) const {
    verifierValidite();
    constexpr std::size_t REGISTRES_SUR_PILE = 64;
    double pile[REGISTRES_SUR_PILE] = {};
    std::vector<double> tas;
    double* r = pile;
    if (m_nombreRegistres > REGISTRES_SUR_PILE) {
        tas.resize(m_nombreRegistres);
        r = tas.data();
    }
    executer(valeurs, r);
    for (std::size_t i = 0; i < m_resultats.size(); ++i) sorties[i] = r[m_resultats[i]];
}

std::vector<double> ProgrammeEvaluation::evaluerEn(const std::vector<double>& valeurs) const {
    if (valeurs.size() != m_nombreEntrees) {
        throw std::invalid_argument("evaluerEn : " + std::to_string(m_nombreEntrees) + " valeur(s) attendue(s), " +
                                    std::to_string(valeurs.size()) + " recue(s)");
    }
    std::vector<double> sorties(m_resultats.size());
    evaluerEn(valeurs.data(), sorties.data());
    return sorties;
}

void ProgrammeEvaluation::evaluer(const double* xs, double* ys, std::size_t n) const {
    verifierValidite();
    if (m_nombreEntrees != 1 || m_resultats.size() != 1) {
        throw std::logic_error("evaluer(xs) : le programme a plusieurs entrees ou plusieurs sorties, utiliser evaluerEn");
    }
    constexpr std::size_t BLOC = 256;
    std::vector<double> registres(static_cast<std::size_t>(m_nombreRegistres) * BLOC);
    for (std::size_t debut = 0; debut < n; debut += BLOC) {
        const std::size_t m = std::min(BLOC, n - debut);
        const double* x = xs + debut;
        for (const Instruction& ins : m_instructions) {
            double* d = &registres[ins.destination * BLOC];
            const double* a = &registres[ins.a * BLOC];
            const double* b = &registres[ins.b * BLOC];
            const double c = ins.valeur;
            switch (ins.code) {
                case Code::Constante: std::fill(d, d + m, c); break;
                case Code::Variable: std::copy(x, x + m, d); break;
                case Code::Axpy:
                    for (std::size_t j = 0; j < m; ++j) d[j] = a[j] + c * b[j];
                    break;
                case Code::Echelle:
                    for (std::size_t j = 0; j < m; ++j) d[j] = c * a[j];
                    break;
                case Code::Produit:
                    for (std::size_t j = 0; j < m; ++j) d[j] = a[j] * b[j];
                    break;
                case Code::PuissanceEntiere:
                    switch (ins.entier) {
                        case 2:
                            for (std::size_t j = 0; j < m; ++j) d[j] = a[j] * a[j];
                            break;
                        case 3:
                            for (std::size_t j = 0; j < m; ++j) d[j] = a[j] * a[j] * a[j];
                            break;
                        case -1:
                            for (std::size_t j = 0; j < m; ++j) d[j] = 1.0 / a[j];
                            break;
                        default:
                            for (std::size_t j = 0; j < m; ++j) d[j] = puissanceEntiereReelle(a[j], ins.entier);
                    }
                    break;
                case Code::Puissance:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::pow(a[j], b[j]);
                    break;
                case Code::Sinus:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::sin(a[j]);
                    break;
                case Code::Cosinus:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::cos(a[j]);
                    break;
                case Code::Tangente:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::tan(a[j]);
                    break;
                case Code::Exponentielle:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::exp(a[j]);
                    break;
                case Code::Logarithme:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::log(a[j]);
                    break;
                case Code::ArcSinus:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::asin(a[j]);
                    break;
                case Code::ArcCosinus:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::acos(a[j]);
                    break;
                case Code::ArcTangente:
                    for (std::size_t j = 0; j < m; ++j) d[j] = std::atan(a[j]);
                    break;
            }
        }
        const double* resultat = &registres[m_resultats[0] * BLOC];
        std::copy(resultat, resultat + m, ys + debut);
    }
}

std::vector<double> ProgrammeEvaluation::evaluer(const std::vector<double>& xs) const {
    std::vector<double> ys(xs.size());
    evaluer(xs.data(), ys.data(), xs.size());
    return ys;
}

} // namespace symalgo
