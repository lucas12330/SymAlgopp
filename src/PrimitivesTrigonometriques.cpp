/**
 * @file PrimitivesTrigonometriques.cpp
 * @brief Intégration des fractions rationnelles de sinus, cosinus et tangentes.
 *
 * Avec c = cos u et s = sin u (u affine), toute fraction rationnelle R(sin k u, cos k u,
 * tan k u) a une forme unique
 *
 *        f = A(c) + s B(c),    A et B fractions rationnelles d'une variable,
 *
 * dans l'extension quadratique Q(c)[s]/(s^2 + c^2 - 1) (ElementQuadratique). Les multiples
 * de l'angle se ramènent à u par cos((k+1)u) = cos(ku) c - sin(ku) s et
 * sin((k+1)u) = sin(ku) c + cos(ku) s. Les règles de Bioche se lisent sur A et B :
 *
 *   A = 0                    f du invariant par u -> -u     t = cos u    int -B(t) dt
 *   A et B impaires          invariant par u -> pi - u       t = sin u    int A1(1 - t^2) + t B1(1 - t^2) dt,
 *                                                                         A = c A1(c^2), B = c B1(c^2)
 *   A paire, B impaire       invariant par u -> pi + u       t = tan u    c^2 = 1/(1 + t^2), dt/(1 + t^2)
 *   sinon                                                    t = tan(u/2) c = (1 - t^2)/(1 + t^2),
 *                                                                         s = 2t/(1 + t^2), du = 2 dt/(1 + t^2)
 *
 * La fraction rationnelle en t est intégrée par primitiveFraction, avec l'expression de t à
 * la place de la variable. Pour t = tan(u/2), la primitive n'est continue que sur les
 * intervalles où tan(u/2) l'est.
 */

#include "Canonique.hpp"
#include "Journal.hpp"
#include "Primitives.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace symalgo {

namespace {

constexpr int DEGRE_MAX = 96;
constexpr long long MULTIPLE_MAX = 24;

bool estTrigonometrique(TypeNoeud t) { return t == TypeNoeud::Sinus || t == TypeNoeud::Cosinus || t == TypeNoeud::Tangente; }

struct ReglesTrigonometriques {
    const Multiples& m;
    const std::vector<ElementQuadratique>& cosinus; // cos(k u), k >= 0
    const std::vector<ElementQuadratique>& sinus;   // sin(k u)

    Feuille feuille(const ASTNode& e, ElementQuadratique& r) const {
        const TypeNoeud type = e.type();
        if (!estTrigonometrique(type)) return Feuille::Descendre;
        const long long k = m.multipleDe(static_cast<const FonctionUnaire&>(e).m_argument.get());
        const std::size_t n = static_cast<std::size_t>(k < 0 ? -k : k);
        const ElementQuadratique sin = k < 0 ? sinus[n] * Nombre(-1) : sinus[n];
        if (type == TypeNoeud::Cosinus) {
            r = cosinus[n];
        } else if (type == TypeNoeud::Sinus) {
            r = sin;
        } else {
            r = sin * cosinus[n].puissance(-1);
        }
        return Feuille::Convertie;
    }

    bool puissance(const ASTNode&, const ASTNode&, ElementQuadratique&) const { return false; }
};

ExprPtr integrer(const ExprPtr& f, bool impaireSeulement) {
    const TypeNoeud type = f->type();
    if (type != TypeNoeud::Produit && type != TypeNoeud::Puissance) return nullptr;
    Multiples m;
    if (!collecter(*f, m, estTrigonometrique) || !m.resoudre(impaireSeulement ? 1 : MULTIPLE_MAX)) return nullptr;

    // Forme A(c) + s B(c)
    const Fraction c = Fraction::variable(), un(Nombre(1));
    const Fraction carre = un - c * c; // s^2
    const ElementQuadratique unite{un, Fraction(), &carre};
    long long multipleMax = 1;
    for (long long k : m.multiples) multipleMax = std::max(multipleMax, k < 0 ? -k : k);
    std::vector<ElementQuadratique> cosinus{unite}, sinus{{Fraction(), Fraction(), &carre}};
    const ElementQuadratique cosU{c, Fraction(), &carre}, sinU{Fraction(), un, &carre};
    for (long long k = 0; k < multipleMax; ++k) {
        const ElementQuadratique& ck = cosinus.back();
        const ElementQuadratique& sk = sinus.back();
        ElementQuadratique suivantC = ck * cosU + sk * sinU * Nombre(-1), suivantS = sk * cosU + ck * sinU;
        cosinus.push_back(std::move(suivantC));
        sinus.push_back(std::move(suivantS));
    }
    ElementQuadratique e = unite;
    ReglesTrigonometriques regles{m, cosinus, sinus};
    try {
        if (!convertirRationnel(*f, e, regles, unite, DEGRE_MAX)) return nullptr;
    } catch (const std::domain_error&) { // division par une expression identiquement nulle
        return nullptr;
    }

    // Règles de Bioche
    const Fraction t = Fraction::variable();
    const int pariteA = e.a.parite(), pariteB = e.b.estNulle() ? -1 : e.b.parite();
    Fraction h;
    ExprPtr variable, angle; // angle : atan(variable) quand la variable est une tangente
    const char* regle;
    const auto parSinus = [&]() { // f = c (A1(c^2) + s B1(c^2)), c^2 = 1 - t^2
        const Fraction y = un - t * t;
        return (e.a / c).enPuissance(2).composer(y) + t * (e.b / c).enPuissance(2).composer(y);
    };
    const auto taille = [](const Fraction& g) {
        int n = g.denominateur().degre();
        for (int i = 0; i <= g.numerateur().degre(); ++i) n += g.numerateur().coefficient(i).estZero() ? 0 : 1;
        return n;
    };
    if (e.a.estNulle()) {
        h = -e.b;
        variable = ast_cos(m.u);
        regle = "f(-u) = -f(u) : t = cos(u)";
        // Les deux règles s'appliquent (sin^3 cos) : la fraction la plus courte
        if (pariteB == -1) {
            const Fraction hSinus = parSinus();
            if (taille(hSinus) < taille(h)) {
                h = hSinus;
                variable = ast_sin(m.u);
                regle = "f(pi - u) = -f(u) : t = sin(u)";
            }
        }
    } else if (pariteA == -1 && pariteB == -1) {
        h = parSinus();
        variable = ast_sin(m.u);
        regle = "f(pi - u) = -f(u) : t = sin(u)";
    } else if (impaireSeulement) {
        return nullptr;
    } else if (pariteA == 1 && pariteB == -1) {
        const Fraction y = un / (un + t * t); // c^2
        h = (e.a.enPuissance(2).composer(y) + t * (e.b * c).enPuissance(2).composer(y)) * y;
        variable = ast_tan(m.u);
        angle = m.u;
        regle = "f(pi + u) = f(u) : t = tan(u)";
    } else {
        const Fraction d = un + t * t;
        const Fraction cosT = (un - t * t) / d, sinT = t * Nombre(2) / d;
        h = (e.a.composer(cosT) + sinT * e.b.composer(cosT)) * Nombre(2) / d;
        angle = m.u / 2.0;
        variable = ast_tan(angle);
        regle = "t = tan(u/2)";
    }
    if (impaireSeulement && !h.estPolynome()) return nullptr;
    if (h.depasse(2 * DEGRE_MAX)) return nullptr;

    Journal* const j = journal::actif();
    const std::size_t etapes = j ? j->etapes().size() : 0;
    journal::noter([&](Journal& jr) {
        jr.ajouter("integration.bioche", std::string("Règle de Bioche avec u = ") + m.u->texte() + ", " + regle, f,
                   h.versExpression(var("t")) / m.pente);
    });
    bool complete;
    ExprPtr F;
    {
        const journal::SousNiveau sous;
        F = primitiveFraction(h, variable, &complete);
    }
    if (!complete) {
        if (j) j->tronquer(etapes);
        return nullptr;
    }
    // atan(tan(w)) = w à une constante près sur chaque intervalle : w est la primitive continue
    if (angle) F = substituer(F, ast_atan(variable), angle);
    return F / m.pente;
}

} // namespace

ExprPtr parTrigonometrieImpaire(const ExprPtr& f) { return integrer(f, true); }

ExprPtr parTrigonometrie(const ExprPtr& f) { return integrer(f, false); }

} // namespace symalgo
