/**
 * @file Multivariable.cpp
 * @brief Variables d'une expression, gradient, jacobienne, hessienne et évaluation
 *        numérique à plusieurs variables.
 */

#include "Multivariable.hpp"
#include "Evaluateur.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

namespace symalgo {

namespace {

// Parcourt le graphe une fois par noeud partagé ; les sous-arbres sans variable sont ignorés
void collecter(const ExprPtr& e, std::unordered_set<const ASTNode*>& vus, std::vector<ExprPtr>& trouvees) {
    if (!e->contientVariable() || !vus.insert(e.get()).second) return;
    switch (e->type()) {
        case TypeNoeud::Variable: trouvees.push_back(e); break;
        case TypeNoeud::Somme:
            for (const Terme& t : static_cast<const Somme&>(*e).getTermes()) collecter(t.expression, vus, trouvees);
            break;
        case TypeNoeud::Produit:
            for (const Facteur& f : static_cast<const Produit&>(*e).getFacteurs()) {
                collecter(f.base, vus, trouvees);
                collecter(f.exposant, vus, trouvees);
            }
            break;
        case TypeNoeud::Puissance:
            collecter(static_cast<const Puissance&>(*e).getBase(), vus, trouvees);
            collecter(static_cast<const Puissance&>(*e).getExposant(), vus, trouvees);
            break;
        case TypeNoeud::IntegraleNonEvaluee:
            collecter(static_cast<const IntegraleNonEvaluee&>(*e).getIntegrande(), vus, trouvees);
            break;
        default:
            if (const FonctionUnaire* f = comme<FonctionUnaire>(e)) collecter(f->m_argument, vus, trouvees);
    }
}

const std::vector<ExprPtr>& variablesOuDefaut(const ExprPtr& f, const std::vector<ExprPtr>& vars,
                                              std::vector<ExprPtr>& defaut) {
    if (!vars.empty()) return vars;
    defaut = variables(f);
    return defaut;
}

} // namespace

std::vector<ExprPtr> variables(const std::vector<ExprPtr>& expressions) {
    std::unordered_set<const ASTNode*> vus;
    std::vector<ExprPtr> trouvees;
    for (const ExprPtr& e : expressions) collecter(e, vus, trouvees);
    std::sort(trouvees.begin(), trouvees.end(), [](const ExprPtr& a, const ExprPtr& b) {
        return static_cast<const Variable&>(*a).getNom() < static_cast<const Variable&>(*b).getNom();
    });
    return trouvees;
}

std::vector<ExprPtr> variables(const ExprPtr& e) { return variables(std::vector<ExprPtr>{e}); }

std::vector<ExprPtr> gradient(const ExprPtr& f, const std::vector<ExprPtr>& vars) {
    std::vector<ExprPtr> defaut;
    const std::vector<ExprPtr>& v = variablesOuDefaut(f, vars, defaut);
    std::vector<ExprPtr> g;
    g.reserve(v.size());
    for (const ExprPtr& x : v) g.push_back(f->derivee(x)->simplifier());
    return g;
}

std::vector<std::vector<ExprPtr>> jacobienne(const std::vector<ExprPtr>& fs, const std::vector<ExprPtr>& vars) {
    const std::vector<ExprPtr> v = vars.empty() ? variables(fs) : vars;
    std::vector<std::vector<ExprPtr>> J;
    J.reserve(fs.size());
    for (const ExprPtr& f : fs) J.push_back(gradient(f, v));
    return J;
}

std::vector<std::vector<ExprPtr>> hessienne(const ExprPtr& f, const std::vector<ExprPtr>& vars) {
    std::vector<ExprPtr> defaut;
    const std::vector<ExprPtr>& v = variablesOuDefaut(f, vars, defaut);
    const std::vector<ExprPtr> g = gradient(f, v);
    const std::size_t n = v.size();
    std::vector<std::vector<ExprPtr>> H(n, std::vector<ExprPtr>(n));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            H[i][j] = g[i]->derivee(v[j])->simplifier();
            H[j][i] = H[i][j];
        }
    }
    return H;
}

ExprPtr laplacien(const ExprPtr& f, const std::vector<ExprPtr>& vars) {
    std::vector<ExprPtr> defaut;
    const std::vector<ExprPtr>& v = variablesOuDefaut(f, vars, defaut);
    ExprPtr somme = nombre(Nombre(0));
    for (const ExprPtr& x : v) somme = somme + f->derivee(x)->derivee(x);
    return somme->simplifier();
}

ExprPtr divergence(const std::vector<ExprPtr>& champ, const std::vector<ExprPtr>& vars) {
    if (champ.size() != vars.size()) {
        throw std::invalid_argument("divergence : " + std::to_string(champ.size()) + " composante(s) pour " +
                                    std::to_string(vars.size()) + " variable(s)");
    }
    ExprPtr somme = nombre(Nombre(0));
    for (std::size_t i = 0; i < champ.size(); ++i) somme = somme + champ[i]->derivee(vars[i]);
    return somme->simplifier();
}

ExprPtr deriveeMixte(const ExprPtr& f, const std::vector<ExprPtr>& ordre) {
    ExprPtr d = f;
    for (const ExprPtr& x : ordre) d = d->derivee(x);
    return d->simplifier();
}

std::vector<double> evaluer(const std::vector<ExprPtr>& expressions, const Valeurs& valeurs) {
    const std::vector<ExprPtr> vars = variables(expressions);
    std::vector<double> point;
    point.reserve(vars.size());
    for (const ExprPtr& v : vars) {
        const std::string& nom = static_cast<const Variable&>(*v).getNom();
        const auto it = valeurs.find(nom);
        if (it == valeurs.end()) throw std::invalid_argument("evaluer : pas de valeur pour la variable '" + nom + "'");
        point.push_back(it->second);
    }
    return ProgrammeEvaluation(expressions, vars).evaluerEn(point);
}

double evaluer(const ExprPtr& e, const Valeurs& valeurs) { return evaluer(std::vector<ExprPtr>{e}, valeurs)[0]; }

} // namespace symalgo
