/**
 * @file Ref.hpp
 * @brief Pointeur intelligent à compteur de références intrusif.
 *
 * Remplace std::shared_ptr pour les noeuds de l'AST :
 *  - le compteur est stocké dans l'objet : une seule allocation, pas de bloc de
 *    contrôle séparé ni de weak_ptr (8 octets par pointeur au lieu de 16) ;
 *  - l'objet est détruit dès que le dernier Ref disparaît, ce qui permet à la
 *    table de hash-consing de le retirer immédiatement (voir ObjetCompte::detruire).
 *
 * Le compteur n'est pas atomique : comme GiNaC, SymAlgo++ ne permet pas de partager
 * une même expression entre plusieurs threads.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <utility>

namespace symalgo {

template <class T>
class Ref;

/*
 * Nom : ObjetCompte
 * Description : Classe de base des objets gérés par Ref : porte le compteur de références.
 */
class ObjetCompte {
public:
    ObjetCompte(const ObjetCompte&) = delete;
    ObjetCompte& operator=(const ObjetCompte&) = delete;

    /*
     * Nom : nombreReferences
     * Description : Nombre de Ref pointant actuellement vers l'objet (diagnostic, tests).
     */
    std::uint32_t nombreReferences() const { return m_references; }

protected:
    ObjetCompte() = default;
    virtual ~ObjetCompte() = default;

    /*
     * Nom : detruire
     * Description : Appelée quand le dernier Ref disparaît. Les classes dérivées peuvent la
     *               redéfinir pour se désinscrire d'une table avant la destruction.
     */
    virtual void detruire() const { delete this; }

private:
    mutable std::uint32_t m_references = 0;

    template <class>
    friend class Ref;
};

/*
 * Nom : Ref
 * Description : Pointeur partagé vers un ObjetCompte (sémantique proche de std::shared_ptr).
 * Utilisation : Ref<Noeud> p = ...; p->methode(); Ref<Base> b = p;
 */
template <class T>
class Ref {
public:
    Ref() noexcept = default;
    Ref(std::nullptr_t) noexcept {}

    /*
     * Nom : Ref
     * Description : Prend une référence sur un objet existant (le compteur est incrémenté).
     */
    explicit Ref(T* objet) noexcept : m_objet(objet) { acquerir(); }

    Ref(const Ref& autre) noexcept : m_objet(autre.m_objet) { acquerir(); }
    Ref(Ref&& autre) noexcept : m_objet(autre.m_objet) { autre.m_objet = nullptr; }

    // Conversion Ref<Derivee> -> Ref<Base>
    template <class U, class = std::enable_if_t<std::is_convertible_v<U*, T*>>>
    Ref(const Ref<U>& autre) noexcept : m_objet(autre.get()) { acquerir(); }

    template <class U, class = std::enable_if_t<std::is_convertible_v<U*, T*>>>
    Ref(Ref<U>&& autre) noexcept : m_objet(autre.relacher()) {}

    ~Ref() { liberer(); }

    Ref& operator=(Ref autre) noexcept {
        std::swap(m_objet, autre.m_objet);
        return *this;
    }

    T* get() const noexcept { return m_objet; }
    T& operator*() const noexcept { return *m_objet; }
    T* operator->() const noexcept { return m_objet; }
    explicit operator bool() const noexcept { return m_objet != nullptr; }

    /*
     * Nom : relacher
     * Description : Abandonne l'objet sans décrémenter le compteur (transfert de propriété).
     */
    T* relacher() noexcept {
        T* objet = m_objet;
        m_objet = nullptr;
        return objet;
    }

    void reset() noexcept { Ref().swap(*this); }
    void swap(Ref& autre) noexcept { std::swap(m_objet, autre.m_objet); }

private:
    void acquerir() const noexcept {
        if (m_objet) ++static_cast<const ObjetCompte*>(m_objet)->m_references;
    }

    void liberer() noexcept {
        if (m_objet) {
            const ObjetCompte* objet = m_objet;
            m_objet = nullptr;
            if (--objet->m_references == 0) objet->detruire();
        }
    }

    T* m_objet = nullptr;
};

template <class T, class U>
bool operator==(const Ref<T>& a, const Ref<U>& b) noexcept { return a.get() == b.get(); }
template <class T, class U>
bool operator!=(const Ref<T>& a, const Ref<U>& b) noexcept { return a.get() != b.get(); }
template <class T>
bool operator==(const Ref<T>& a, std::nullptr_t) noexcept { return a.get() == nullptr; }
template <class T>
bool operator!=(const Ref<T>& a, std::nullptr_t) noexcept { return a.get() != nullptr; }
template <class T>
bool operator==(std::nullptr_t, const Ref<T>& a) noexcept { return a.get() == nullptr; }
template <class T>
bool operator!=(std::nullptr_t, const Ref<T>& a) noexcept { return a.get() != nullptr; }

} // namespace symalgo

// Permet d'utiliser Ref comme clé de std::unordered_map / std::unordered_set
template <class T>
struct std::hash<symalgo::Ref<T>> {
    std::size_t operator()(const symalgo::Ref<T>& r) const noexcept { return std::hash<T*>()(r.get()); }
};
