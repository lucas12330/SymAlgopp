# SymAlgo++

SymAlgo++ est une bibliothèque C++17 moderne et performante conçue pour représenter, manipuler, évaluer, dériver, intégrer et simplifier des expressions algébriques classiques ainsi que des équations différentielles ordinaires (EDO). Elle s'appuie sur une modélisation orientée objet robuste utilisant un Arbre de Syntaxe Abstraite (AST) pour le calcul symbolique et des solveurs algébriques/numériques intégrés.

---

## 🏛️ Architecture de la bibliothèque

Le projet est structuré autour d'une hiérarchie de classes exploitant le polymorphisme et le C++17 moderne :

* **`Equation`** : Classe abstraite de base définissant l'interface commune de toutes les équations du système.
  * `eval(double x)` : Évalue l'équation pour une valeur réelle donnée.
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
   * `Fraction` : Représente un nombre rationnel exact $A/B$ (`int64_t`) avec un cache d'évaluation float haute performance. La simplification calcule exactement sur les fractions ($1/3 + 1/6 = 1/2$), avec détection des débordements.
   * `Variable` : Représente la variable d'évaluation (par défaut `"x"` ; `var("v")` fonctionne de la même façon).
   * `Parametre` : Constante symbolique sans valeur (ex. les constantes $C_1, C_2$ des solutions d'EDO) : dérivée nulle, évaluation impossible (`std::logic_error`).
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
4. **Nœuds non évalués** : `IntegraleNonEvaluee` et `LimiteNonEvaluee` représentent une primitive ou une limite que la bibliothèque ne sait pas déterminer. Ils remplacent tout résultat faux (la dérivée d'une intégrale non évaluée redonne l'intégrande ; les évaluer lève `std::logic_error`).

### Simplification et Fonctions Helpers

Des constructeurs d'aide (`cst`, `frac`, `var`, `param`, `ast_pow`, `ast_sin`, `ast_cos`, `ast_tan`, `ast_exp`, `ast_ln`) ainsi que des surcharges d'opérateurs arithmétiques permettent une écriture proche des mathématiques :

```cpp
auto X = var("x");
// Représente (1/3)*x^2 + sin(x) + ln(x)
auto expr = frac(1, 3) * ast_pow(X, 2) + ast_sin(X) + ast_ln(X);
```

---

## ⚡ Fonctionnalités Avancées

* **Intégration Symbolique (`integrer()`)** : Moteur formel de recherche de primitives par reconnaissance de motifs (polynômes, fonctions usuelles, linéarité, facteurs constants, substitution linéaire $\int f(ax+b)\,dx = F(ax+b)/a$). Hors de ces règles, le résultat est une `IntegraleNonEvaluee`.
* **Limites (`limite(x0)`)** : Règle de L'Hôpital pour $0/0$ et $\infty/\infty$ (profondeur bornée), formes $0 \cdot \infty$, $1^\infty$, $0^0$, $\infty^0$, signe de l'infini déterminé par le développement du dénominateur. Une limite inexistante (ex. $1/x$ en 0, où les limites à gauche et à droite diffèrent) ou non déterminée est une `LimiteNonEvaluee`.
* **Développements Limités (`DL(x0, ordre)`)** : Taylor / Maclaurin à n'importe quel ordre, calculé par arithmétique des séries tronquées (technique de la différentiation automatique) : l'ordre 20 de $e^{\sin x}$ s'obtient en quelques dizaines de microsecondes.
* **Tracé Adaptatif (`genererPointsTrace(xMin, xMax, tolerance)`)** : Génération de courbes par échantillonnage adaptatif récursif, réduisant le nombre de points requis en zone linéaire tout en affinant les zones de forte courbure.
* **Résolution Littérale d'EDO (`resoudreLitteral()`)** : Résolution formelle exacte des équations différentielles linéaires homogènes à coefficients constants, fournissant la combinaison linéaire des solutions (racines réelles, complexes conjuguées et multiples : $x^k e^{rx}$).
* **Problème de Cauchy (`resoudreProblemeCauchy()`)** : Solution exacte satisfaisant les conditions initiales, directement évaluable.
* **Solveur Numérique EDO (RK4)** : Simulation numérique temporelle d'ordre 4 basée sur le schéma d'espace d'état et la matrice compagnon (Eigen) ; chaque pas se réduit à un produit matrice-vecteur précalculé.

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
* **Lancer la suite de tests automatisés** (mini-framework sans dépendance, `tests/test_framework.hpp` ; code de retour non nul en cas d'échec) :
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

Les mesures ci-dessous (médiane de 3 répétitions) ont été réalisées sur un processeur Intel Celeron N4120 @ 1.10 GHz (Arch Linux, GCC 16, `-O3`).

#### A. Évaluation et Dérivation (SymAlgo++ vs GiNaC)

| Opération / Scénario | SymAlgo++ | GiNaC | Comparaison |
| :--- | :--- | :--- | :--- |
| **Évaluation numérique** ($x=5$) | **~ 143 ns** | ~ 30 641 ns | 🚀 **SymAlgo++ est ~210x plus rapide** |
| **Dérivation symbolique** ($f'(x)$, simplifiée) | **~ 6 816 ns** | ~ 28 021 ns | 🚀 **SymAlgo++ est ~4x plus rapide** |

![Comparaison SymAlgo++ vs GiNaC](docs/images/bench_comparison.svg)

* **Analyse** :
  * **Évaluation** : L'évaluation parcourt directement l'arbre en `double` (polymorphisme, aucune allocation), là où GiNaC substitue puis évalue symboliquement.
  * **Dérivation** : Le partage des sous-arbres immuables (`std::enable_shared_from_this`) évite toute copie. Les anciennes mesures (GiNaC 1.45x plus rapide) étaient faussées : le `Makefile` liait le benchmark à une bibliothèque compilée sans optimisation.
  * *Comparaison à nuancer* : GiNaC est un système de calcul formel complet (forme canonique, arithmétique exacte généralisée), dont les opérations font davantage de travail.

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
