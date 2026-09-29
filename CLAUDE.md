# Directives et Contexte pour Claude Code : SymAlgo++

> 🛑 **RÈGLE SYSTÈME ABSOLUE : PRIORITÉ AU CONTEXTE**
> Le fichier [`PROJECT_CONTEXT.md`](PROJECT_CONTEXT.md) constitue la source de vérité fondamentale du projet. Toutes les démarches, architectures et propositions techniques doivent obligatoirement s'y aligner sans exception.

---

## 🚀 Présentation Rapide
**SymAlgo++** est une bibliothèque C++17 open-source moderne pour la représentation, la dérivation, l'intégration, la simplification et la résolution (analytique et numérique) d'expressions algébriques et d'équations différentielles ordinaires (EDO).

---

## 🛠️ Commandes de Build, Test et Benchmark

* **Compiler la bibliothèque et les démos/tests** :
  ```bash
  make
  ```
* **Exécuter la suite de tests automatisés** (échoue au premier test en échec) :
  ```bash
  make run_tests
  ```
* **Exécuter les benchmarks comparatifs (-O3 / Google Benchmark / GiNaC)** :
  ```bash
  make bench
  ```
* **Nettoyer les artéfacts de compilation** :
  ```bash
  make clean
  ```
* **Exécution directe des binaires individuels** :
  * `./bin/test_ast` : Tests unitaires de l'Arbre de Syntaxe Abstraite (AST).
  * `./bin/test_differentielle` : Tests des équations différentielles (matrice compagnon, RK4, solutions littérale et de Cauchy).
  * `./bin/demo` : Programme de démonstration.
  * `./bin/bench_suite` : Exécutable du benchmark de performance.

---

## 🏛️ Architecture et Hiérarchie des Classes

Le projet repose sur la Programmation Orientée Objet et le polymorphisme C++17. Tout le code est dans `namespace symalgo`.

1. **`Equation`** (`include/Equation.hpp`) : Classe abstraite fondamentale.
   * `virtual double eval(double x) const = 0;`
   * `virtual std::unique_ptr<Equation> deriveeGenerique() const = 0;` (dérivée polymorphe ; chaque classe dérivée expose aussi `derivee()` qui renvoie son propre type **par valeur**)

2. **`EquationClassique`** (`include/EquationClassique.hpp`, `src/EquationClassique.cpp`) :
   * Wrapper orienté objet autour de la racine d'un AST (`ExprPtr`).
   * Méthodes clés (résultats renvoyés **par valeur**, jamais de `new`) : `eval()`, `derivee()`, `simplifier()`, `integrer()`, `limite(x0)` (L'Hôpital, formes indéterminées), `DL(x0, ordre)` (Taylor par arithmétique des séries), `genererPointsTrace(xMin, xMax, tolerance)` (échantillonnage adaptatif), `getExpression()`.

3. **`EquationDifferentielle`** (`include/EquationDifferentielle.hpp`, `src/EquationDifferentielle.cpp`) :
   * Représentation linéaire d'EDO : $\sum a_i y^{(i)} = 0$.
   * Méthodes clés : `ajouterTerme(rang, coeff)`, `setConditionsInitiales(ci)`, `getMatriceCompagnon()` (Eigen state-space matrix), `resoudreLitteral()` (solution générale avec paramètres C1..Cn, racines multiples gérées), `resoudreProblemeCauchy()` (solution exacte avec conditions initiales), `derivee()` (EDO dont la solution est y'), `eval(x)` (solveur numérique RK4).

4. **`ASTNode`** (`include/ASTNode.hpp`, `src/ASTNode.cpp`) :
   * Heritage de `std::enable_shared_from_this<ASTNode>` pour l'optimisation mémoire du `clone()`.
   * Alias de pointeur intelligent : `using ExprPtr = std::shared_ptr<ASTNode>;`.
   * **Nœuds terminaux** : `Constante`, `Fraction` (rationnel exact $A/B$, arithmétique exacte lors de la simplification), `Variable` (variable d'évaluation, quel que soit son nom), `Parametre` (constante symbolique : dérivée nulle, évaluation impossible).
   * **Nœuds non évalués** : `IntegraleNonEvaluee`, `LimiteNonEvaluee` (renvoyés quand aucune règle ne s'applique, au lieu d'un résultat faux).
   * **Points d'entrée** : `integrer()` et `limite(a)` sont non virtuels ; les règles propres à chaque nœud sont dans les méthodes protégées `primitive()` et `calculerLimite(a)`.
   * **Opérateurs binaires** : `Addition`, `Soustraction`, `Multiplication`, `Division`, `Puissance` (support $u(x)^{v(x)}$ et exposants constants).
   * **Fonctions unaires** : `Sinus`, `Cosinus`, `Tangente`, `Exponentielle`, `Logarithme`.
   * **Helpers & Surcharges** : `cst()`, `frac()`, `var()`, `param()`, `ast_pow()`, `ast_sin()`, `ast_cos()`, `ast_tan()`, `ast_exp()`, `ast_ln()`, surcharges d'opérateurs `+`, `-`, `*`, `/`.

---

## 🎨 Normes de Code & Bonnes Pratiques C++

* **Langage & Standard** : C++17 moderne.
* **Sécurité Mémoire** : Utiliser `std::shared_ptr` / `std::unique_ptr`. Éviter les `raw pointers` pour la gestion de propriété.
* **Const-Correctness** : Marquer `const` toutes les méthodes de consultation et passer les types non primitifs par référence constante (`const std::string&`, `const std::vector<double>&`).
* **Tests** : chaque correction ou fonctionnalité s'accompagne de cas dans `tests/test_*.cpp` (mini-framework `tests/test_framework.hpp` : `TEST_CASE`, `CHECK`, `CHECK_EQ`, `CHECK_NEAR`, `CHECK_THROWS`). Les nœuds de l'AST se créent toujours via les helpers (jamais sur la pile).
* **Dépendances** :
  * **STL** (`std::vector`, `std::map`, `std::shared_ptr`) comme alternative native et performante.
  * **Eigen** (`vendor/eigen`) pour l'algèbre linéaire / matrices d'espace d'état.
  * **Google Benchmark** & **GiNaC** pour les mesures comparatives de performance dans `benchmarks/`.

---

## 📌 Gestion du Dépôt Git (Règles Importantes)

* **Branches** :
  * Travailler **toujours** sur une branche de fonctionnalité (`feature/nom-de-la-tache`).
  * **Interdiction stricte** de commit ou push directement sur `main`.
* **Convention de Commits (en français)** :
  * `feat: ...` (ajout de fonctionnel)
  * `fix: ...` (correction de bug)
  * `refactor: ...` (restructuration ou optimisation sans changement d'API)
  * `docs: ...` (documentation)
  * `chore: ...` (Makefile, .gitignore, configurations)
* **Propreté du dépôt** :
  * Ne jamais commiter de fichiers générés (`.o`, `build/`, `bin/`, fichiers temporaires).
