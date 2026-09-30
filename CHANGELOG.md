# Changelog

Ce document répertorie tous les ajouts et correctifs majeurs de la bibliothèque SymAlgo++.

## [Unreleased] - Résolution d'équations (branches feature/solveur-equations → feature/resolution-equations)

### Ajouts (Additions)
* **`developper()`** : distribution des produits et puissances entières de sommes (exponentiation rapide), mémorisée sur le graphe partagé.
* **`Polynome`** : polynômes à une variable à coefficients `Nombre` exacts ; opérations, division euclidienne, PGCD, dérivée, décomposition sans carré (Yun), évaluation exacte.
* **Racines réelles certifiées** : isolement par suites de Sturm et bissection rationnelle exacte ; racines rationnelles reconnues exactement, degré 2 par radicaux, les autres arrondies au `double` le plus proche (polynôme de Wilkinson de degré 20 : 11,6 ms).
* **`factoriser()`** : factorisation sur $\mathbb{Q}$ (facteurs linéaires rationnels avec multiplicité, parties primitives restantes).
* **Constante `pi` exacte** et fonctions **`asin`, `acos`, `atan`** (évaluation, dérivées, primitives, limites, DL, affichage, évaluation compilée).
* **Valeurs remarquables exactes** : sin/cos/tan des multiples de $\pi/6$ et $\pi/4$ et leurs réciproques ; radicaux canoniques (`8^(1/2) = 2*2^(1/2)`, dénominateurs rationalisés) ; `ln(8) = 3*ln(2)`.
* **`resoudre()`** : résolution exacte par isolement de l'inconnue, règle du produit nul, changement de variable, polynômes ; solutions avec approximation, multiplicité et indicateur de complétude ; vérification du domaine (`x*ln(x) = 0` donne seulement 1).
* **Familles trigonométriques** : `sin(x) = 1/2` donne `pi/6 + 2*pi*k` et `5*pi/6 + 2*pi*k`.
* **`resoudreSurIntervalle()`** et **`resoudreNumerique()`** : familles dépliées sur un intervalle ; recherche numérique (échantillonnage, méthode de Brent, racines doubles aux points critiques).
* `EquationClassique` : `resoudre()`, `resoudre(a, b)`, `developper()`, `factoriser()`.
* Helpers `substituer(e, cible, remplacement)` et `contient(e, cible)`.

## [Unreleased] - Optimisation de fond (branches feature/mesure-performances → feature/evaluation-compilee)

### Performances (mesurées contre GiNaC, détail dans `benchmarks/RESULTATS.md`)
* **Dérivée 10e** : de 4,4 s et 901 Mo alloués à 1,4 ms et 0,56 Mo (GiNaC : 7,5 ms).
* **Mémoire** : le résultat d'une dérivée 6e occupe 4,2 Ko au lieu de 375 Ko (GiNaC : 4,9 Ko).
* **Évaluation** : grande expression 1 450 fois plus rapide que GiNaC, tableau de 10 000 points 494 fois plus rapide.
* SymAlgo++ est désormais devant GiNaC sur tous les scénarios mesurés.

### Ajouts (Additions)
* **Forme canonique automatique** : sommes et produits n-aires triés, regroupement des termes semblables (`x + x = 2*x`, `x*x^2 = x^3`, `2*(x+1) = 2*x + 2`), simplification des quotients.
* **Hash-consing** : chaque expression est unique en mémoire ; égalité en temps constant.
* **`Nombre`** : rationnels exacts de taille arbitraire (repli GMP) ou réels ; `2^100`, `4^(1/2) = 2` exacts.
* **Évaluation compilée** (`ProgrammeEvaluation`, `eq.eval(std::vector<double>)`) : sous-expressions partagées calculées une fois, blocs vectorisés, compilation automatique quand elle est rentable.
* Limites désormais déterminées : $x \ln x \to 0$ et $x^x \to 1$ en 0.
* Cible `make check` (AddressSanitizer + UBSan) ; benchmark élargi avec mesure de la mémoire.

### Correctifs (Patches)
* Use-after-free dans la détection des arguments linéaires de l'intégration (détecté par AddressSanitizer) ; `comme<T>()` est désormais interdit sur un temporaire.
* `0/0` reste indéfini (NaN) sous forme canonique.

### Changements d'API (Breaking Changes)
* `ExprPtr` devient `Ref<ASTNode>` (compteur intrusif) ; `std::dynamic_pointer_cast` est remplacé par `comme<T>()`.
* Les classes `Fraction`, `Addition`, `Soustraction`, `Multiplication`, `Division` et `OperateurBinaire` sont remplacées par `Constante` (nombre exact), `Somme`, `Produit` et `Puissance`.
* Affichage : `2*x`, `x^2`, `sin(x)/x`, `x^3/3` (au lieu de `2 * x`, `(x)^(2)`, `(sin(x) / x)`).
* `cst(2.0)` est l'entier exact 2 ; `sin(2)` reste symbolique (évalué par `eval`), alors que `sin(2.5)` (argument réel) est évalué.
* Les nœuds ne se construisent plus sur la pile (erreur de compilation) ; le partage d'une expression entre threads n'est pas pris en charge.

## Branche feature/corrections-robustesse

### Correctifs majeurs (Major Patches)
* **Constantes d'EDO** : les constantes $C_1, C_2$ de `resoudreLitteral()` étaient des variables (évaluées comme $x$, dérivée 1). Nouveau nœud `Parametre` (helper `param()`).
* **Intégration** : les primitives non calculables renvoyaient 0 ($\int x \cdot x = 0$, $\int \sin 2x = 0$). Nouveau nœud `IntegraleNonEvaluee`, et règle de substitution linéaire $\int f(ax+b) = F(ax+b)/a$.
* **Limites** : $x \cdot (1/x)$ en 0 donnait NaN, $-1/x$ en 0 donnait $+\infty$, $\ln$ d'un argument négatif donnait $-\infty$. Formes indéterminées, signe de l'infini et limites inexistantes (`LimiteNonEvaluee`) sont désormais gérés ; récursion de L'Hôpital bornée.
* **Racines multiples** : $y'' + 2y' + y = 0$ donnait une base dégénérée ($C_1 e^{-x} + C_2 e^{-x}$) au lieu de $(C_1 + C_2 x) e^{-x}$.
* **Coefficient dominant nul** d'une EDO : division par zéro dans la matrice compagnon.
* **Simplification** : toute constante inférieure à $10^{-9}$ était considérée nulle ($6.674 \cdot 10^{-11} x \to 0$) ; les fractions perdaient leur exactitude ($1/3 + 1/3 = 0.666667$) ; $0/0$ devenait 0.
* **Développements limités** : termes d'ordre $\geq 13$ perdus (seuil absolu sur $1/k!$), factorielle débordant à $21!$, NaN silencieux aux points non développables.
* **Affichage** des équations différentielles : « 1y » pour un coefficient $\pm 1$.
* **`Makefile`** : bibliothèque compilée sans optimisation (benchmarks faussés), headers non suivis.
* Nœud créé sur la pile : `std::logic_error` explicite au lieu de `std::bad_weak_ptr`.

### Ajouts (Additions)
* `EquationDifferentielle::resoudreProblemeCauchy()` : solution exacte avec conditions initiales.
* `EquationDifferentielle::derivee()` : EDO dont la solution est $y'$.
* Suite de tests à assertions (`tests/test_framework.hpp`, 63 cas) ; `make run_tests` échoue en cas d'erreur.
* Cible `make demo` (`bin/demo`).

### Performances
* `DL()` par arithmétique des séries de Taylor : ordre 10 de $e^{\sin x}$ de 2,9 s à 0,015 ms.
* RK4 : un pas = un produit matrice-vecteur précalculé (x6 à x8).
* Dérivation : ~4x plus rapide que GiNaC (mesures refaites avec la bibliothèque optimisée).

### Changements d'API (Breaking Changes)
* Tout le code est dans `namespace symalgo`.
* `derivee()`, `integrer()`, `limite()`, `DL()`, `resoudreLitteral()` renvoient un `EquationClassique` **par valeur** (plus de `delete`).
* `Equation::derivee()` devient `deriveeGenerique()` (renvoie `std::unique_ptr<Equation>`).

## Branche feature/equation-classique

### Ajouts majeurs (Major Additions)
* **Intégration Symbolique** : Implémentation d'un moteur de calcul formel d'intégrales par reconnaissance de motifs pour les expressions classiques (`ASTNode::integrer`). Prise en charge des polynômes, fonctions trigonométriques (sinus, cosinus), et de la linéarité.
* **Limites Symboliques (L'Hôpital)** : Ajout d'un système de calcul de limites formelles (`ASTNode::limite`). Pour le cas des fractions avec indétermination de type $\frac{0}{0}$, l'algorithme dérive symboliquement le numérateur et le dénominateur selon la célèbre règle de L'Hôpital.
* **Développements Limités (DL)** : Implémentation de la génération symbolique de développements de Taylor et Maclaurin (`ASTNode::DL`) via le calcul récursif de dérivées formelles (ex: approximation polynomiale).
* **Tracé et Tableaux de Points Adaptatifs** : Ajout de la méthode `genererPointsTrace` utilisant un algorithme récursif d'échantillonnage adaptatif. Cela permet de réduire radicalement la quantité de points calculés sur les sections linéaires des courbes tout en gardant une excellente précision dans les variations.
* **Dérivation Généralisée** : Le système dérive désormais les fonctions puissances de la forme $u(x)^{v(x)}$ en utilisant rigoureusement le logarithme népérien, tout en optimisant explicitement le comportement pour les exposants constants.
* **Nouvelle primitive mathématique** : Ajout du nœud `Logarithme` (`ast_ln`).

### Correctifs majeurs (Major Patches)
* *Aucun correctif majeur pour le moment.*
