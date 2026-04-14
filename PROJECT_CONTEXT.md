# Contexte du Projet et Lignes Directrices : SymAlgo++
**Document de Référence pour Agents IA**

Ce document sert de pilier et de contexte fondamental pour tout assistant IA travaillant sur le projet **SymAlgo++**. Son but est d'informer l'IA des choix technologiques, de l'architecture et de la philosophie de cette librairie afin de garantir des propositions de code cohérentes, performantes et structurées.

---

## 1. Présentation Générale
* **Nom du Projet** : SymAlgo++
* **Type de Projet** : Librairie Gratuite et Open Source.
* **Langage** : C++ (choisi spécifiquement pour ses performances et ses puissantes capacités de Programmation Orientée Objet).
* **Mission** : Fournir une interface efficace pour représenter, manipuler et résoudre des équations mathématiques de toutes sortes.

## 2. Architecture Orientée Objet (Le Cœur du Projet)
L'architecture de base repose sur l'héritage et le polymorphisme C++. L'IA doit toujours respecter cette structure en arbre lors de l'ajout de fonctionnalités.

* **Classe Abstraite/Super Classe `Equation`** : 
  Sert de fondation. Elle définit l'interface virtuelle que toutes les équations devront implémenter.
* **Classe `EquationClassique` (Hérite de `Equation`)** : 
  Gère les expressions algébriques et les fonctions classiques (ex: $y = 3x^2 + 2$).
* **Classe `EquationDifferentielle` (Hérite de `Equation`)** : 
  Dédiée aux équations faisant intervenir des dérivées (ex: $y' + 2y = 0$).

## 3. Fonctionnalités Clés à Implémenter
Les classes du projet doivent être capables de s'interfacer avec les méthodes suivantes, que ce soit via des méthodes intrinsèques ou des fonctions utilitaires :
1. **Résolution / Évaluation** : Capacité à évaluer l'équation en un point précis (ex: pour $x = A$), ou sur une plage finie de points.
2. **Dérivation** : Méthodes permettant de dériver l'équation.
3. **Intégration** : Méthodes permettant d'intégrer l'équation sur un domaine.
4. **Calcul de Limites** : Fonctionnalités prévues pour évaluer la limite d'une équation en un point.
5. **Tracé et Affichage** : Outils facilitant la récolte de données de la fonction en vue d'être tracées (génération de tableaux de points).

## 4. Écosystème et Dépendances
L'intégration de bibliothèques tierces doit être justifiée par un besoin d'efficacité :
* **STL (Standard Template Library)** : L'utilisation de types comme `std::vector` est fortement recommandée. Elle servira d'alternative locale aux "arrays" que l'on aurait pu trouver dans des outils Python comme *Numpy*, permettant une manipulation souple des listes de points.
* **Eigen** : Pour tout besoin futur de calcul matriciel intense ou d'algèbre linéaire avancée, **Eigen** peut être utilisé.
* Le but reste de garder la librairie la plus autonome possible sans la surcharger de trop de dépendances.

## 5. Directives de Qualité de Code pour l'IA
* **Performance** : Le temps d'exécution est un critère critique et sera évalué. Veillez à utiliser les structures de données adéquates et à optimiser les passages coûteux en ressource (ex: éviter les copies inutiles).
* **Paradigmes C++** : Utilisez du C++ moderne. Pensez à sécuriser la mémoire (utilisation pertinente de `std::unique_ptr` ou `std::shared_ptr` dans la hiérarchie polymorphique), utilisez `const` par défaut lorsque c'est pertinent, et le passage par référence constante pour les objets lourds.
* **Extensibilité** : La conception doit permettre l'ajout futur de nouvelles familles d'équations avec un minimum de modifications du code existant.

---
*Note aux Agents : Lisez ce document avant de modifier ou lire l'architecture du projet.*
