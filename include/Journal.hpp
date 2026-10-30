/**
 * @file Journal.hpp
 * @brief Journal des étapes d'un calcul (résolution, intégration...), pour les expliquer
 *        ou les exporter en LaTeX.
 *
 * Le journal est facultatif et ne coûte rien quand on ne le demande pas : les algorithmes
 * notent leurs étapes par journal::noter(), qui ne fait qu'un test de pointeur tant
 * qu'aucun journal n'est actif (le texte de l'étape n'est alors jamais construit).
 *
 *   Journal j;
 *   Solutions s = resoudre(gauche, droite, j);   // surcharge publique avec journal
 *   for (const Etape& e : j.etapes()) std::cout << e.description << '\n';
 *
 * Côté algorithme :
 *
 *   journal::noter([&](Journal& j) {
 *       j.ajouter("resolution.produit_nul", "Un produit est nul si l'un de ses facteurs l'est",
 *                 equation, facteur);
 *   });
 *   journal::SousNiveau sous; // étapes suivantes en retrait, jusqu'à la fin du bloc
 *
 * Le journal actif est propre au thread (comme les expressions, qui ne se partagent pas
 * entre threads). Les caches (mémo de simplification, dérivées partagées) peuvent éviter
 * un calcul déjà fait : un algorithme qui veut des étapes complètes doit les contourner
 * quand journal::actif() n'est pas nul.
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "ASTNode.hpp"

namespace symalgo {

/*
 * Nom : Etape
 * Description : Une étape d'un calcul. regle est un identifiant stable (« module.regle »,
 *               ex. « integration.substitution ») pour les tests et la traduction ;
 *               description est une phrase en français ; avant / apres sont les
 *               expressions avant et après l'étape (nulles si sans objet) ; niveau est la
 *               profondeur d'imbrication (0 : étape principale).
 */
struct Etape {
    std::string regle;
    std::string description;
    ExprPtr avant;
    ExprPtr apres;
    int niveau = 0;
};

/*
 * Nom : Journal
 * Description : Suite ordonnée des étapes d'un ou plusieurs calculs.
 */
class Journal {
public:
    void ajouter(std::string regle, std::string description, ExprPtr avant = ExprPtr(),
                 ExprPtr apres = ExprPtr()) {
        m_etapes.push_back({std::move(regle), std::move(description), std::move(avant), std::move(apres), m_niveau});
    }

    const std::vector<Etape>& etapes() const noexcept { return m_etapes; }
    bool vide() const noexcept { return m_etapes.empty(); }
    void effacer() noexcept {
        m_etapes.clear();
        m_niveau = 0;
    }

    int niveau() const noexcept { return m_niveau; }
    void entrer() noexcept { ++m_niveau; }
    void sortir() noexcept {
        if (m_niveau > 0) --m_niveau;
    }

private:
    std::vector<Etape> m_etapes;
    int m_niveau = 0;
};

namespace journal {

/*
 * Nom : actif
 * Description : Journal du calcul en cours sur ce thread, nul si aucune étape n'est demandée.
 */
inline Journal*& actif() noexcept {
    thread_local Journal* courant = nullptr;
    return courant;
}

/*
 * Nom : noter
 * Description : Appelle fabrique(journal) seulement si un journal est actif : le texte de
 *               l'étape n'est construit que s'il sert.
 */
template <class Fabrique>
inline void noter(Fabrique&& fabrique) {
    if (Journal* j = actif()) fabrique(*j);
}

/*
 * Nom : Activation
 * Description : Rend un journal actif jusqu'à la fin du bloc (puis restaure le précédent).
 *               À utiliser dans les surcharges publiques qui reçoivent un Journal&.
 */
class Activation {
public:
    explicit Activation(Journal& j) noexcept : m_precedent(actif()) { actif() = &j; }
    ~Activation() { actif() = m_precedent; }
    Activation(const Activation&) = delete;
    Activation& operator=(const Activation&) = delete;

private:
    Journal* m_precedent;
};

/*
 * Nom : SousNiveau
 * Description : Met en retrait les étapes notées jusqu'à la fin du bloc (sans effet si
 *               aucun journal n'est actif).
 */
class SousNiveau {
public:
    SousNiveau() noexcept : m_journal(actif()) {
        if (m_journal) m_journal->entrer();
    }
    ~SousNiveau() {
        if (m_journal) m_journal->sortir();
    }
    SousNiveau(const SousNiveau&) = delete;
    SousNiveau& operator=(const SousNiveau&) = delete;

private:
    Journal* m_journal;
};

} // namespace journal

} // namespace symalgo
