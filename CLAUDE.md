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
* **Exécuter les tests sous AddressSanitizer + UBSan** (à lancer avant chaque commit touchant la mémoire) :
  ```bash
  make check
  ```
* **Exécuter les benchmarks comparatifs (-O3 / Google Benchmark / GiNaC)** — résultats consignés dans `benchmarks/RESULTATS.md` :
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
  * `./bin/test_nombre` : Tests des nombres exacts (`Nombre`, repli GMP).
  * `./bin/test_evaluateur` : Tests de l'évaluation compilée.
  * `./bin/test_polynome` : Tests de `developper`, `Polynome` (Sturm, racines certifiées) et `factoriser`.
  * `./bin/test_solveur` : Tests du solveur d'équations (exact, familles, intervalles, numérique).
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
   * Méthodes clés (résultats renvoyés **par valeur**, jamais de `new`) : `eval(x)` (compilé automatiquement si le partage le rend rentable), `eval(xs)` (tableau, évaluation compilée par blocs), `derivee()`, `simplifier()`, `integrer()`, `limite(x0)`, `DL(x0, ordre)`, `genererPointsTrace(xMin, xMax, tolerance)`, `developper()`, `factoriser()`, `resoudre()`, `resoudre(a, b)`, `getExpression()`.

3. **`EquationDifferentielle`** (`include/EquationDifferentielle.hpp`, `src/EquationDifferentielle.cpp`) :
   * Représentation linéaire d'EDO : $\sum a_i y^{(i)} = 0$.
   * Méthodes clés : `ajouterTerme(rang, coeff)`, `setConditionsInitiales(ci)`, `getMatriceCompagnon()` (Eigen state-space matrix), `resoudreLitteral()` (solution générale avec paramètres C1..Cn, racines multiples gérées), `resoudreProblemeCauchy()` (solution exacte avec conditions initiales), `derivee()` (EDO dont la solution est y'), `eval(x)` (solveur numérique RK4).

4. **Expressions (AST)** — schéma complet dans `docs/ARCHITECTURE.md` :
   * `ExprPtr = Ref<ASTNode>` (`include/Ref.hpp`) : compteur de références **intrusif**, non atomique (pas de partage d'une expression entre threads). Les nœuds ne se créent que via `fabriquer<T>()` et les helpers (clé `CleFabrique`) ; `comme<T>(e)` remplace `dynamic_cast` (interdit sur un temporaire).
   * **Hash-consing** : chaque expression est unique en mémoire (`Signature` + table intrusive dans `src/Noeud.cpp`) ; `estEgal` = comparaison d'adresses.
   * **Forme canonique automatique** (`src/Noeud.cpp`, `include/Canonique.hpp`) : `Somme` (constante + termes triés à coefficients `Nombre`), `Produit` (coefficient + facteurs `base^exposant` triés), `Puissance`. Pas de nœuds Soustraction/Division : `a - b = a + (-1)*b`, `a / b = a * b^(-1)`. Ordre total : `comparer()`.
   * **`Nombre`** (`include/Nombre.hpp`, `src/Nombre.cpp`) : rationnel exact (int64 avec repli GMP) ou réel ; `cst(2.0)` est exact, `cst(0.5)` réel.
   * **Nœuds** : `Constante`, `Variable`, `Parametre`, `Pi`, `Somme`, `Produit`, `Puissance`, `Sinus`, `Cosinus`, `Tangente`, `Exponentielle`, `Logarithme`, `ArcSinus`, `ArcCosinus`, `ArcTangente`, `IntegraleNonEvaluee`, `LimiteNonEvaluee`. Valeurs remarquables exactes (tables de `src/Noeud.cpp`) et radicaux canoniques.
   * **Points d'entrée non virtuels** (`derivee`, `simplifier`, `integrer`, `limite`) ; règles par nœud dans les méthodes protégées `calculerDerivee`, `calculerSimplification`, `primitive`, `calculerLimite`. Modules : `Regles.cpp` (éval, dérivées, primitives), `Limites.cpp`, `Series.cpp` (DL), `Affichage.cpp`, `Evaluateur.cpp` (programme compilé).
   * **Helpers** : `cst()`, `frac()`, `nombre()`, `var()`, `param()`, `pi()`, `somme()`, `produit()`, `ast_pow()`, `ast_sin()`, `ast_cos()`, `ast_tan()`, `ast_exp()`, `ast_ln()`, `ast_asin()`, `ast_acos()`, `ast_atan()`, `appliquer()`, `substituer()`, `contient()`, opérateurs `+ - * /` et moins unaire.

5. **Algèbre** :
   * `include/Polynome.hpp`, `src/Polynome.cpp` : `developper()`, `Polynome` (coefficients `Nombre` exacts, `estExact()` faux si un coefficient était réel ; PGCD, `sansCarre()` de Yun, `racinesReelles()` certifiées par Sturm), `factoriser()` sur Q.
   * `include/Solveur.hpp`, `src/Solveur.cpp` : `resoudre(gauche, droite)`, `resoudreSurIntervalle(e, a, b)`, `resoudreNumerique(f, a, b)` ; `Solutions { liste, complet, toutReel }`, `Solution { valeur, approximation, multiplicite, exacte, entiers }` (familles trigonométriques paramétrées par `k`, `k2`...).

---

## 🎨 Normes de Code & Bonnes Pratiques C++

* **Langage & Standard** : C++17 moderne.
* **Sécurité Mémoire** : Utiliser `std::unique_ptr` / `std::shared_ptr`, et `Ref` (`ExprPtr`) pour les nœuds de l'AST. Éviter les `raw pointers` pour la gestion de propriété ; ne jamais conserver le pointeur renvoyé par `comme<T>()` au-delà de la vie de l'`ExprPtr` source.
* **Const-Correctness** : Marquer `const` toutes les méthodes de consultation et passer les types non primitifs par référence constante (`const std::string&`, `const std::vector<double>&`).
* **Tests** : chaque correction ou fonctionnalité s'accompagne de cas dans `tests/test_*.cpp` (mini-framework `tests/test_framework.hpp` : `TEST_CASE`, `CHECK`, `CHECK_EQ`, `CHECK_NEAR`, `CHECK_THROWS`). Les nœuds de l'AST se créent toujours via les helpers (jamais sur la pile).
* **Dépendances** :
  * **STL** (`std::vector`, `std::map`, `std::shared_ptr`) comme alternative native et performante.
  * **Eigen** (`vendor/eigen`) pour l'algèbre linéaire / matrices d'espace d'état.
  * **GMP** (`-lgmpxx -lgmp`) pour les rationnels exacts de taille arbitraire, confiné à `src/Nombre.cpp`.
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
