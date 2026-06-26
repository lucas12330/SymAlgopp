# SymAlgo++

SymAlgo++ est une bibliothèque C++ conçue pour représenter, manipuler, évaluer, dériver et simplifier des expressions algébriques classiques ainsi que des équations différentielles ordinaires. Elle s'appuie sur une modélisation orientée objet moderne utilisant un arbre de syntaxe abstraite (AST) pour la partie symbolique.

---

## Architecture de la bibliothèque

Le projet est structuré autour d'une hiérarchie de classes exploitant le polymorphisme en C++ :

* **`Equation`** : Classe abstraite de base définissant l'interface commune de toutes les équations du système.
  * `eval(double x)` : Évalue l'équation pour une valeur de variable donnée.
  * `derivee()` : Retourne un pointeur vers une nouvelle instance d'équation représentant sa dérivée.
* **`EquationClassique`** (hérite d'`Equation`) : Représente des expressions algébriques simples. Elle encapsule un arbre de syntaxe abstraite (`ASTNode`) et implémente la résolution formelle (dérivation symbolique, simplification mathématique).
* **`EquationDifferentielle`** (hérite d'`Equation`) : Dédiée à la représentation linéaire d'équations différentielles de la forme $\sum a_i y^{(i)} = 0$. Elle stocke en interne une table associative associant l'ordre de dérivation à son coefficient multiplicateur.

---

## Arbre de Syntaxe Abstraite (AST)

La manipulation symbolique des expressions classiques repose sur la classe abstraite **`ASTNode`** et l'alias de pointeur intelligent associé :
```cpp
using ExprPtr = std::shared_ptr<ASTNode>;
```

### Nœuds disponibles
1. **Nœuds terminaux** :
   * `Constante` : Contient une valeur réelle (`double`).
   * `Variable` : Représente la variable (par défaut `"x"`).
2. **Opérateurs binaires** :
   * `Addition` ($+$)
   * `Soustraction` ($-$)
   * `Multiplication` ($*$)
   * `Division` ($/$)
   * `Puissance` (base ^ exposant)
3. **Fonctions unaires** :
   * `Sinus` ($\sin$)
   * `Cosinus` ($\cos$)

### Simplification et Surcharges
Des fonctions d'aide (`cst`, `var`) et des surcharges d'opérateurs arithmétiques permettent d'écrire des expressions de manière proche de l'écriture mathématique usuelle.

Par exemple, au lieu d'instancier manuellement chaque nœud, l'expression $x^2 + 5x + 6$ s'écrit :
```cpp
auto X = var("x");
auto expr = ast_pow(X, 2) + (X * 5) + 6;
```

---

## Exemples d'utilisation

### 1. Expressions classiques et dérivation symbolique
Cet exemple montre comment construire l'expression $x^2 + 5x + 6$, l'évaluer pour $x = 2$, calculer sa dérivée symbolique ($2x + 5$) et évaluer cette dernière.

```cpp
#include <iostream>
#include "EquationClassique.hpp"
#include "ASTNode.hpp"

int main() {
    auto X = var("x");
    // Construction de x^2 + 5x + 6
    auto expr = ast_pow(X, 2) + (X * 5) + 6;
    
    EquationClassique eq(expr);
    std::cout << "Equation originale : ";
    eq.afficher(); // Affiche : (((x)^(2) + x * 5) + 6) = 0
    
    std::cout << "Evaluation (x = 2) : " << eq.eval(2.0) << std::endl; // Affiche : 20

    // Dérivation formelle
    EquationClassique* eq_derivee = eq.derivee();
    std::cout << "Derivee simplifiee : ";
    eq_derivee->afficher(); // Affiche : (2 * x + 5) = 0
    
    std::cout << "Evaluation derivee (x = 2) : " << eq_derivee->eval(2.0) << std::endl; // Affiche : 9

    delete eq_derivee;
    return 0;
}
```

### 2. Équations différentielles
Cet exemple montre comment définir et manipuler des équations différentielles sous forme linéaire en y ajoutant des termes.

```cpp
#include <iostream>
#include "EquationDifferentielle.hpp"

int main() {
    // Exemple 1 : Oscillateur Harmonique (y'' + 4y = 0)
    EquationDifferentielle eqHarmonique;
    eqHarmonique.ajouterTerme(2, 1.0); // Ajoute 1.0 * y''
    eqHarmonique.ajouterTerme(0, 4.0); // Ajoute 4.0 * y
    
    std::cout << "Oscillateur Harmonique : ";
    eqHarmonique.afficher(); // Affiche : 4*y^0 + 1*y^2 + 0 = 0

    // Exemple 2 : Équation d'ordre supérieur (3y''' - 2y' + 5y = 0)
    EquationDifferentielle eqComplexe;
    eqComplexe.ajouterTerme(3, 3.0);  // Ajoute 3.0 * y'''
    eqComplexe.ajouterTerme(1, -2.0); // Ajoute -2.0 * y'
    eqComplexe.ajouterTerme(0, 5.0);  // Ajoute 5.0 * y
    
    std::cout << "Equation complexe : ";
    eqComplexe.afficher(); // Affiche : 5*y^0 + -2*y^1 + 3*y^3 + 0 = 0

    return 0;
}
```

---

## Compilation et Tests

Le projet fournit un `Makefile` permettant de gérer la compilation de la bibliothèque et des exécutables de test.

* **Compiler l'ensemble du projet** :
  ```bash
  make
  ```
* **Lancer la suite de tests automatisés** (recommandé pour valider les changements) :
  ```bash
  make run_tests
  ```
* **Nettoyer les fichiers de build (fichiers objets et binaires)** :
  ```bash
  make clean
  ```

Les fichiers intermédiaires sont compilés dans le dossier `build/` et les exécutables finaux sont générés dans le dossier `bin/`.

---

## Organisation du dépôt et conventions

* **Arborescence** :
  * `include/` : Contient les fichiers d'en-tête de la bibliothèque (`.hpp`).
  * `src/` : Contient les fichiers sources d'implémentation (`.cpp`).
  * `tests/` : Contient les programmes de validation et tests unitaires.
* **Stratégie de branche** :
  Toute modification doit être effectuée sur une branche de fonctionnalité dédiée (`feature/nom-de-la-tache`) et soumise via Pull Request. Les commits directs sur la branche `main` sont interdits.
* **Conventions de messages de commits** :
  Les messages doivent respecter les préfixes suivants :
  * `feat:` (nouvelle fonctionnalité)
  * `fix:` (correction de bug)
  * `refactor:` (restructuration du code sans changement de comportement)
  * `docs:` (documentation)
  * `chore:` (fichiers de configuration, Makefile, .gitignore)
