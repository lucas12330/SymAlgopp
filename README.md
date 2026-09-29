# SymAlgo++

SymAlgo++ est une bibliothèque C++17 moderne et performante conçue pour représenter, manipuler, évaluer, dériver, intégrer et simplifier des expressions algébriques classiques ainsi que des équations différentielles ordinaires (EDO). Elle s'appuie sur une modélisation orientée objet robuste utilisant un Arbre de Syntaxe Abstraite (AST) pour le calcul symbolique et des solveurs algébriques/numériques intégrés.

---

## 🏛️ Architecture de la bibliothèque

Le projet est structuré autour d'une hiérarchie de classes exploitant le polymorphisme et le C++17 moderne :

* **`Equation`** : Classe abstraite de base définissant l'interface commune de toutes les équations du système.
  * `eval(double x)` : Évalue l'équation pour une valeur rélle donnée.
  * `deriveeGenerique()` : Dérivée polymorphe, renvoyée sous forme de `std::unique_ptr<Equation>`. Chaque classe dérivée propose aussi `derivee()`, qui renvoie son propre type par valeur.
* **`EquationClassique`** (hérite d'`Equation`) : Encapsule un arbre de syntaxe abstraite (`ASTNode`) et fournit des opérations algébriques avancées (dérivation formelle, intégration symbolique, calcul de limites par L'Hôpital, développements limités et tracé adaptatif).
* **`EquationDifferentielle`** (hérite d'`Equation`) : Représente des équations différentielles linéaires de la forme $\sum a_i y^{(i)} = 0$. Offre à la fois une résolution analytique/littérale exacte et une résolution numérique via Runge-Kutta 4 (RK4) basée sur la matrice compagnon d'état (Eigen).

---

## 🌲 Arbre de Syntaxe Abstraite (AST)

La manipulation symbolique repose sur la classe abstraite `ASTNode` bénéficiant de l'optimisation mémoire `std::enable_shared_from_this<ASTNode>` et de l'alias de pointeur intelligent :

```cpp
using ExprPtr = std::shared_ptr<ASTNode>;
```

Toute la bibliothèque est déclarée dans l'espace de noms `symalgo` (les opérateurs `+ - * /` sur `ExprPtr` sont trouvés automatiquement par ADL).

### Nœuds disponibles dans l'AST

1. **Nœuds terminaux** :
   * `Constante` : Contient une valeur réelle (`double`).
   * `Fraction` : Représente un nombre rationnel exact $A/B$ (`int64_t`) avec un cache d'évaluation float haute performance.
   * `Variable` : Représente l'inconnue symbolique (par défaut `"x"`).
2. **Opérateurs binaires** :
   * `Addition` ($+$)
   * `Soustraction` ($-$)
   * `Multiplication` ($*$)
   * `Division` ($/$)
   * `Puissance` (base ^ exposant, optimisée pour exposants constants et formes générales $u(x)^{v(x)}$).
3. **Fonctions unaires** :
   * `Sinus` ($\sin$)
   * `Cosinus` ($\cos$)
   * `Tangente` ($\tan$)
   * `Exponentielle` ($\exp$)
   * `Logarithme` ($\ln$)

### Simplification et Fonctions Helpers

Des constructeurs d'aide (`cst`, `frac`, `var`, `ast_pow`, `ast_sin`, `ast_cos`, `ast_tan`, `ast_exp`, `ast_ln`) ainsi que des surcharges d'opérateurs arithmétiques permettent une écriture proche des mathématiques :

```cpp
auto X = var("x");
// Représente (1/3)*x^2 + sin(x) + ln(x)
auto expr = frac(1, 3) * ast_pow(X, 2) + ast_sin(X) + ast_ln(X);
```

---

## ⚡ Fonctionnalités Avancées

* **Intégration Symbolique (`integrer()`)** : Moteur formel de recherche de primitives par reconnaissance de motifs (polynômes, séries trigonométriques, linéarité).
* **Limites Symboliques (`limite(x0)`)** : Calcul des limites formelles intégrant l'application récursive de la règle de L'Hôpital pour lever les formes indéterminées $0/0$.
* **Développements Limités (`DL(x0, ordre)`)** : Approximation symbolique de Taylor / Maclaurin à n'importe quel ordre.
* **Tracé Adaptatif (`genererPointsTrace(xMin, xMax, tolerance)`)** : Génération de courbes par échantillonnage adaptatif récursif, réduisant le nombre de points requis en zone linéaire tout en affinant les zones de forte courbure.
* **Résolution Littérale d'EDO (`resoudreLitteral()`)** : Résolution formelle exacte des équations différentielles linéaires homogènes à coefficients constants, fournissant la combinaison linéaire des solutions (réelles ou complexes conjuguées).
* **Solveur Numérique EDO (RK4)** : Simulation numérique temporelle d'ordre 4 basée sur le schéma d'espace d'état et la matrice compagnon (Eigen).

---

## 💡 Exemples d'utilisation

### 1. Expressions classiques, dérivation, intégration et limite
```cpp
#include <iostream>
#include "EquationClassique.hpp"
#include "ASTNode.hpp"

using namespace symalgo;

int main() {
    auto X = var("x");
    // (x^2 - 1) / (x - 1)
    auto expr = (ast_pow(X, 2) - 1.0) / (X - 1.0);
    EquationClassique eq(expr);

    std::cout << "Equation : ";
    eq.afficher();

    // Calcul de la limite en x = 1 (L'Hôpital) -> 2
    EquationClassique eq_lim = eq.limite(1.0);
    std::cout << "Limite en x->1 : ";
    eq_lim.afficher();

    // Dérivation et Intégration d'un polynôme (résultats renvoyés par valeur)
    EquationClassique poly(ast_pow(X, 2) + X * 5 + 6);
    EquationClassique poly_der = poly.derivee();
    EquationClassique poly_int = poly.integrer();

    std::cout << "Derivee   : "; poly_der.afficher(); // (2*x + 5)
    std::cout << "Integrale : "; poly_int.afficher(); // (1/3*x^3 + 2.5*x^2 + 6*x)
    return 0;
}
```

### 2. Équations différentielles (Littérale et Numérique RK4)
```cpp
#include <iostream>
#include "EquationDifferentielle.hpp"
#include "EquationClassique.hpp"

using namespace symalgo;

int main() {
    // Oscillateur Harmonique : y'' + 4y = 0
    EquationDifferentielle eq;
    eq.ajouterTerme(2, 1.0); // y''
    eq.ajouterTerme(0, 4.0); // 4y

    std::cout << "EDO : ";
    eq.afficher(); // y'' + 4*y = 0

    // Solution générale : C1, C2 sont des paramètres symboliques
    EquationClassique sol_generale = eq.resoudreLitteral();
    std::cout << "Solution analytique : ";
    sol_generale.afficher(); // C1*cos(2*x) + C2*sin(2*x)

    // Avec les conditions initiales y(0)=1, y'(0)=0
    eq.setConditionsInitiales({1.0, 0.0});
    EquationClassique sol_exacte = eq.resoudreProblemeCauchy();
    sol_exacte.afficher(); // cos(2*x)

    // Évaluation numérique via RK4
    std::cout << "Evaluation RK4 en x=pi/4 : " << eq.eval(3.141592653589793 / 4.0) << std::endl; // ~ 0
    return 0;
}
```

---

## 🛠️ Compilation et Tests

Le projet est doté d'un `Makefile` complet et portable :

* **Compiler l'ensemble du projet** :
  ```bash
  make
  ```
* **Lancer la suite de tests automatisés** :
  ```bash
  make run_tests
  ```
* **Nettoyer les fichiers de build** :
  ```bash
  make clean
  ```

---

## 📊 Benchmarks de Performances

Les performances de SymAlgo++ sont mesurées à l'aide de **Google Benchmark** et comparées à **GiNaC** (bibliothèque C++ formelle de référence).

### 1. Exécution des Benchmarks
```bash
sudo pacman -Syu benchmark ginac # Arch / EndeavourOS
make bench
```

### 2. Résultats des Mesures

Les tests ont été réalisés sur un processeur Intel Core i5 @ 2.60 GHz (EndeavourOS).

#### A. Évaluation et Dérivation (SymAlgo++ vs GiNaC)

| Opération / Scénario | SymAlgo++ | GiNaC | Comparaison |
| :--- | :--- | :--- | :--- |
| **Évaluation numérique** ($x=5$) | **~ 274 ns** | ~ 33 074 ns | 🚀 **SymAlgo++ est ~120x plus rapide** |
| **Dérivation symbolique** ($f'(x)$) | ~ 41 696 ns | **~ 28 757 ns** | ⚠️ GiNaC est ~1.45x plus rapide |

![Comparaison SymAlgo++ vs GiNaC](docs/images/bench_comparison.svg)

* **Analyse** :
  * **Évaluation** : Grâce au cache de valeurs réelles dans le nœud `Fraction` et au polymorphisme direct, SymAlgo++ offre un débit d'évaluation exceptionnel (~120x supérieur à GiNaC).
  * **Dérivation** : Les récentes optimisations (`std::enable_shared_from_this`, réutilisation de pointeurs immuables) ont réduit l'écart avec GiNaC de 1.7x à 1.45x.

#### B. Scalabilité EDO (Génération Matrice Compagnon)

![Scalabilité de l'Équation Différentielle](docs/images/bench_ode_scaling.svg)

* **Analyse** : La génération de la matrice compagnon conserve une complexité $O(N)$ strictement linéaire par rapport à l'ordre $N$ de l'équation (mesurée jusqu'à l'ordre 15+).

---

## 📂 Conventions et Contribution

* **Arborescence** :
  * `include/` : En-têtes de la bibliothèque (`.hpp`).
  * `src/` : Fichiers sources (`.cpp`).
  * `tests/` : Suite de tests automatisés.
  * `benchmarks/` : Suite de benchmarks Google Benchmark et scripts d'analyse.
  * `docs/` : Documentation et visuels SVG.
* **Stratégie de branche** :
  Travail exclusivement sur branches dédiées (`feature/nom-de-la-tache`). Direct commit sur `main` proscrit.
* **Conventions de commit (en français)** :
  * `feat:` Nouvelle fonctionnalité.
  * `fix:` Correction de bogue.
  * `refactor:` Amélioration de structure/code.
  * `docs:` Documentation.
  * `chore:` Configuration, Makefile, .gitignore.
