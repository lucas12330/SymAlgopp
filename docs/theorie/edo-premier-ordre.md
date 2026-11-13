# Équations différentielles du premier ordre (non linéaires)

Ce document rassemble la théorie nécessaire à un solveur symbolique d'équations différentielles ordinaires (EDO) du premier ordre dans SymAlgo++ : pour chaque classe classique, comment la **reconnaître** sur une expression canonique, la **formule** de la solution, les solutions **perdues** par la méthode, le **domaine** de validité, et un exemple **vérifié par réinjection**. Le solveur existant (`EquationDifferentielle`) ne traite que les équations linéaires homogènes à coefficients constants ; tout ce qui suit est à construire (voir la section 15, « État dans le code »).

Le second ordre réductible au premier ordre fait l'objet de [`edo-second-ordre-reductible.md`](edo-second-ordre-reductible.md).

## Sommaire

1. [Notations](#1-notations)
2. [Trois écritures d'une EDO du premier ordre](#2-trois-écritures-dune-edo-du-premier-ordre)
3. [Existence, unicité et solutions singulières](#3-existence-unicité-et-solutions-singulières)
4. [Outils communs : test de nullité, gel de variable, réinjection](#4-outils-communs--test-de-nullité-gel-de-variable-réinjection)
5. [Variables séparables](#5-variables-séparables)
6. [Équation linéaire](#6-équation-linéaire)
7. [Équation de Bernoulli](#7-équation-de-bernoulli)
8. [Équation y' = F(ax + by + c)](#8-équation-y--fax--by--c)
9. [Équations homogènes et réductibles à homogènes](#9-équations-homogènes-et-réductibles-à-homogènes)
10. [Équations exactes et facteurs intégrants](#10-équations-exactes-et-facteurs-intégrants)
11. [Équation de Riccati](#11-équation-de-riccati)
12. [Équations de Clairaut et de Lagrange](#12-équations-de-clairaut-et-de-lagrange)
13. [Stratégie d'un solveur : ordre des tests et échec](#13-stratégie-dun-solveur--ordre-des-tests-et-échec)
14. [Pièges et cas de contrôle](#14-pièges-et-cas-de-contrôle)
15. [État dans le code](#15-état-dans-le-code)
16. [Annexe : script de vérification](#16-annexe--script-de-vérification)

## 1. Notations

| Notation | Sens |
| :--- | :--- |
| $x$ | variable indépendante (réelle) |
| $y = y(x)$ | fonction inconnue ; $p = y'$ |
| $f(x, y)$ | second membre de la forme résolue $y' = f(x, y)$ |
| $f_x, f_y, f_{xy}$ | dérivées partielles $\partial f/\partial x$, $\partial f/\partial y$, $\partial^2 f/\partial x \partial y$ |
| $M\,dx + N\,dy = 0$ | forme différentielle, $M = M(x, y)$, $N = N(x, y)$ |
| $F(x, y, p) = 0$ | forme implicite (non résolue en $p = y'$) |
| $C, K$ | constantes arbitraires réelles ; $K$ désigne une constante « absorbée » (par exemple $K = \pm e^{C}$) |
| $I$ | intervalle ouvert de définition d'une solution |
| $e \equiv 0$ | $e$ est identiquement nulle (comme fonction), décidé par le test de la section 4.1 |

Dans le code, on note `x = var("x")`, `y = var("y")`, `p = var("p")` (trois variables distinctes, voir `include/Multivariable.hpp`) et les constantes `param("C1")`, `param("C2")`… comme le fait déjà `EquationDifferentielle::resoudreLitteral` (`src/EquationDifferentielle.cpp:239`).

## 2. Trois écritures d'une EDO du premier ordre

**Forme résolue.** $y' = f(x, y)$, avec $f$ définie sur un ouvert $U \subset \mathbb{R}^2$. Une solution est un couple $(I, y)$ : $y$ dérivable sur l'intervalle $I$, $(x, y(x)) \in U$ et $y'(x) = f(x, y(x))$ pour tout $x \in I$.

**Forme différentielle.** $M(x, y)\,dx + N(x, y)\,dy = 0$. Une solution est une courbe le long de laquelle la forme s'annule. Là où $N \neq 0$, c'est $y' = -M/N$ ; là où $M \neq 0$, c'est $dx/dy = -N/M$, ce qui autorise des solutions verticales $x = x_0$ (si $N(x_0, y) \equiv 0$). Le passage de la forme résolue à la forme différentielle n'est pas unique : $y' = P/Q$ s'écrit $-P\,dx + Q\,dy = 0$ mais aussi $-\frac{P}{Q}\,dx + dy = 0$ ; les facteurs intégrants (section 10) dépendent de ce choix. **Recommandation** : si $f$ est un quotient, prendre $M = -\text{numérateur}$, $N = \text{dénominateur}$.

**Forme implicite.** $F(x, y, y') = 0$. Là où $F_p \neq 0$, le théorème des fonctions implicites donne localement une ou plusieurs formes résolues $y' = f_k(x, y)$ (une par branche, par exemple $y'^2 = x$ donne $y' = \pm\sqrt{x}$). Les points où $F = 0$ et $F_p = 0$ simultanément sont le lieu des **solutions singulières possibles** (section 3.3).

## 3. Existence, unicité et solutions singulières

### 3.1 Énoncés

**Théorème de Cauchy–Lipschitz.** Soit $U \subset \mathbb{R} \times \mathbb{R}^n$ ouvert et $f : U \to \mathbb{R}^n$ continue et localement lipschitzienne par rapport à $y$ (c'est le cas si $\partial f / \partial y$ existe et est continue sur $U$). Pour tout $(x_0, y_0) \in U$, le problème $y' = f(x, y)$, $y(x_0) = y_0$ admet une unique solution maximale ; toute solution du problème en est une restriction. Une solution maximale sort de tout compact de $U$ aux bornes de son intervalle (théorème des bouts).

**Théorème de Cauchy–Peano.** Si $f$ est seulement continue sur $U$, le problème admet au moins une solution locale ; l'unicité peut être perdue.

Conséquence pour $n = 1$ et $f$ de classe $C^1$ : les graphes de deux solutions maximales distinctes ne se coupent pas. En particulier, si $y \equiv y_0$ est une solution constante (un équilibre), aucune autre solution n'atteint la valeur $y_0$ en un point de $U$.

### 3.2 Solution générale, solutions perdues, solutions singulières

Le vocabulaire varie selon les manuels ; on fixe ici :

* **solution générale** : famille $y = \varphi(x, C)$ (ou $\Phi(x, y) = C$, ou une forme paramétrée) obtenue par une méthode, $C$ parcourant un ensemble de valeurs admissibles ;
* **solution perdue** : solution de l'équation qui n'est obtenue pour **aucune** valeur admissible de $C$ — typiquement une solution annulant un diviseur utilisé par la méthode ($h(y) = 0$ pour une séparable, $y = 0$ pour une Bernoulli, la solution particulière d'une Riccati qui correspond à « $C = \infty$ ») ;
* **solution singulière** (sens fort) : solution en chaque point de laquelle l'unicité est en défaut (par exemple l'enveloppe d'une famille de Clairaut). Une solution singulière est toujours perdue par la méthode qui produit la famille, mais une solution perdue n'est pas forcément singulière (exemple : $y = 0$ pour $y' = x y^2$ est perdue par la séparation des variables, mais l'unicité y est vraie).

Pour un solveur, la distinction utile est : **toute division par une expression qui peut s'annuler doit produire une liste de candidats** (les zéros du diviseur), chacun étant **réinjecté** dans l'équation de départ ; ceux qui la vérifient et n'appartiennent pas à la famille sont ajoutés au résultat.

### 3.3 Perte d'unicité et raccordements

Si $f$ n'est pas lipschitzienne en $y$ au voisinage d'un équilibre, des solutions peuvent atteindre l'équilibre en un temps fini et s'y raccorder.

**Exemple (vérifié).** $y' = 3 y^{2/3}$ avec la racine cubique réelle. La famille $y = (x + C)^3$ convient : $y' = 3(x + C)^2 = 3\big((x+C)^3\big)^{2/3}$. L'équilibre $y = 0$ convient. Pour $a \le b$, la fonction égale à $(x - a)^3$ sur $]-\infty, a]$, à $0$ sur $[a, b]$ et à $(x - b)^3$ sur $[b, +\infty[$ est encore une solution $C^1$ : par chaque point $(x_0, 0)$ passent une infinité de solutions. $y = 0$ est ici **singulière** au sens fort. Un solveur doit renvoyer la famille, l'équilibre, et signaler la non-unicité (il ne peut pas lister les raccordements).

> Note d'implémentation : `Puissance::eval` calcule `std::pow(b, e)` pour un exposant non entier (`src/Regles.cpp:168`), qui vaut NaN pour $b < 0$ et $e = 2/3$. La formule $(x + C)^3$ est correcte partout, mais la réinjection numérique dans $3 y^{2/3}$ doit se faire en des points où $y > 0$.

### 3.4 Forme implicite : le $p$-discriminant

Pour $F(x, y, p) = 0$ avec $F$ de classe $C^1$, si une solution $y$ est singulière en $x_0$, alors $F_p(x_0, y(x_0), y'(x_0)) = 0$ (sinon le théorème des fonctions implicites ramène localement à une forme résolue, régulière). Les candidats sont donc obtenus en **éliminant $p$** entre $F = 0$ et $F_p = 0$. Cette élimination donne une condition **nécessaire** seulement : la courbe obtenue peut être un lieu de points de rebroussement ou ne pas être solution. **Chaque candidat doit être réinjecté.**

## 4. Outils communs : test de nullité, gel de variable, réinjection

### 4.1 Test de nullité

Tous les tests de reconnaissance ci-dessous se ramènent à décider si une expression $e(x, y)$ est identiquement nulle. Le problème est indécidable en général dès que les expressions mêlent exponentielles, fonctions trigonométriques et valeur absolue (théorème de Richardson) ; on procède donc en trois étages :

```text
Nul(e) :
  1. si e est la constante 0 après forme canonique         -> NUL (certain)
  2. e2 <- simplifier(developper(e)) ; si e2 est 0           -> NUL (certain)
  3. évaluer e en k points rationnels tirés au hasard dans le domaine
     (k = 6 à 10, points où l'évaluation est finie) :
       si une valeur dépasse tol * échelle                  -> NON NUL (certain, aux arrondis près)
       sinon                                                -> PROBABLEMENT NUL
```

L'étage 3 est asymétrique : une valeur clairement non nulle **prouve** la non-nullité, des valeurs nulles ne font que la suggérer. Un test de reconnaissance qui conclut sur « probablement nul » doit être suivi d'une **réinjection** de la solution finale (section 4.4) ; si elle échoue, on traite le cas comme un échec de la méthode.

Coût : les dérivées partielles et la forme canonique sont linéaires en la taille du graphe partagé (une seule dérivation par nœud et par appel, voir `docs/ARCHITECTURE.md`) ; l'étage 3 coûte $k$ évaluations.

### 4.2 Dérivées partielles

`e->derivee("y")` (`include/ASTNode.hpp:202-203`, `src/Noeud.cpp:233`) donne $\partial e / \partial y$, les autres variables étant des constantes. Les paramètres (`Parametre`) ont une dérivée nulle (`src/Regles.cpp:71`).

### 4.3 Intégrer par rapport à une variable en « gelant » l'autre

`ASTNode::integrer()` refuse les expressions à plusieurs variables et renvoie une `IntegraleNonEvaluee` (`src/Noeud.cpp:265-268`). Pour calculer $\int M(x, y)\,dx$ à $y$ fixé :

```text
IntegrerEnX(M) :
  Y <- param("_y")                  # nom réservé, absent de M
  A <- substituer(M, y, Y)->integrer()
  renvoyer substituer(A, Y, y)
```

Pour intégrer **par rapport à $y$**, il faut aujourd'hui **renommer $y$ en $x$** : plusieurs règles de primitive écrivent en dur la variable `x` — `Constante::primitive` (`src/Regles.cpp:58`), `Parametre::primitive` (`src/Regles.cpp:73`), le terme constant d'une somme (`src/Regles.cpp:106`) et l'intégrale d'une expression sans variable (`src/Noeud.cpp:266`). Ainsi, intégrer `1 + y^2` « en $y$ » donnerait aujourd'hui `x + y^3/3`. D'où :

```text
IntegrerEnY(H) :                    # H peut contenir x (gelé) et y
  X <- param("_x")
  H1 <- substituer(substituer(H, x, X), y, x)     # y devient x
  A  <- H1->integrer()
  renvoyer substituer(substituer(A, x, y), X, x)
```

L'ordre des substitutions compte (geler $x$ avant de renommer $y$ en $x$, et l'inverse au retour).

### 4.4 Réinjection

Une solution explicite $y = \varphi(x, C)$ se vérifie en calculant $r = \varphi' - f(x, \varphi)$ et en appliquant `Nul(r)`. Une solution implicite $\Phi(x, y) = C$ se vérifie par $\Phi_x + \Phi_y\,f \equiv 0$ (forme résolue) ou $\Phi_x N - \Phi_y M \equiv 0$ (forme différentielle : les gradients sont colinéaires). Une solution paramétrée $(X(p), Y(p))$ de $F(x, y, p) = 0$ se vérifie par $Y'(p) - p\,X'(p) \equiv 0$ et $F(X(p), Y(p), p) \equiv 0$. Les constantes $C$ restent des paramètres symboliques ; pour l'étage numérique, on leur donne aussi des valeurs aléatoires.

**Toute solution renvoyée par le solveur doit avoir passé la réinjection** ; c'est ce qui garantit l'absence de résultat faux, quelle que soit la sophistication des tests de reconnaissance.

### 4.5 Constantes absorbées et domaine

Les étapes « $\ln|y| = G(x) + C \Rightarrow y = K e^{G(x)}$ » font apparaître une constante $K = \pm e^{C} \neq 0$. Si $K = 0$ redonne une solution perdue (ici $y \equiv 0$), on le dit et on autorise $K \in \mathbb{R}$. De même, une formule n'est valable que sur les intervalles où toutes les étapes ont un sens (argument d'un logarithme non nul, racine d'un nombre positif, dénominateur non nul) ; le résultat doit porter ces conditions, au moins sous forme de commentaire (voir l'exemple de Bernoulli $n = 1/2$, section 7.4, où la formule hors domaine **n'est pas** solution).

## 5. Variables séparables

### 5.1 Définition et théorème

L'équation est **à variables séparées** si $f(x, y) = g(x)\,h(y)$ sur un pavé $I \times J$, avec $g$ continue sur $I$ et $h$ continue sur $J$.

**Théorème.**
1. Pour tout $y_0 \in J$ tel que $h(y_0) = 0$, la constante $y \equiv y_0$ est solution.
2. Sur tout intervalle où $h(y(x)) \neq 0$, $y$ est solution si et seulement si $H(y(x)) = G(x) + C$, où $H' = 1/h$ et $G' = g$. Comme $H$ est strictement monotone sur chaque intervalle de $J$ où $h$ ne s'annule pas, $y = H^{-1}(G(x) + C)$.
3. Si $h$ est localement lipschitzienne (par exemple $C^1$), une solution non constante n'atteint jamais un zéro de $h$ (section 3.1) : les points 1 et 2 donnent alors **toutes** les solutions. Sinon, des raccordements sont possibles (section 3.3).

*Preuve du point 2.* $(H \circ y)' = y'/h(y) = g$, puis on intègre. ∎

### 5.2 Test de reconnaissance

**Proposition.** Soit $f$ de classe $C^2$ et sans zéro sur un pavé. Alors $f(x, y) = g(x) h(y)$ sur ce pavé si et seulement si
$$R = f\,f_{xy} - f_x\,f_y \equiv 0 .$$

*Preuve.* $(\ln|f|)_{xy} = R / f^2$. Si $R \equiv 0$, $\ln|f| = A(x) + B(y)$ sur le pavé (connexe), donc $f = \pm e^{A(x)} e^{B(y)}$ avec un signe constant. Réciproquement, si $f = g h$, $R = g h\,g' h' - g' h\,g h' = 0$. ∎

Le test reste une condition nécessaire là où $f$ s'annule ; la décision finale est donc prise par l'**extraction suivie d'une vérification** :

```text
Separer(f) :                                  # renvoie (g, h) ou ECHEC
  si Nul(f_y) : renvoyer (f, 1)                # quadrature y = ∫ f dx + C
  si Nul(f_x) : renvoyer (1, f)                # équation autonome
  si non Nul(f*f_xy - f_x*f_y) : ECHEC
  choisir (x0, y0) parmi des petits rationnels (1, 2, 1/2, 3, -1, ...)
      tels que f(x0, y0) soit fini et non nul
  h <- f(x0, y)                                # substituer(f, x, x0)
  g <- f(x, y0) / f(x0, y0)
  si Nul(f - g*h) : renvoyer (g, h) sinon ECHEC
```

Raccourci structurel (moins coûteux, à tenter d'abord) : si $f$ est un `Produit` canonique, répartir ses facteurs entre ceux qui ne contiennent que $x$, ceux qui ne contiennent que $y$ et les constantes (`contient(facteur, y)`, `include/ASTNode.hpp:762`). Si un facteur contient les deux variables, revenir au test.

**Pièges.**
* **Facteurs constants** : $f = 2xy$ donne $h(y) = f(1, y) = 2y$ et $g(x) = f(x, 1)/f(1, 1) = x$ : la constante est rangée dans $h$. Sans importance mathématique, mais le choix doit être stable (pour des tests reproductibles).
* **Somme qui est un produit** : $f = xy + x + y + 1 = (x+1)(y+1)$ arrive comme une `Somme` canonique ; le raccourci structurel échoue, le test réussit ($R = 0$).
* **Fonction d'une somme** : $e^{x+y} = e^x e^y$ est séparable sans être un produit syntaxique ; le test réussit.
* **Partie non séparable qui se simplifie** : $\sin(x+y) + \sin(x-y) = 2 \sin x \cos y$. Le test est vrai, mais l'extraction par points donne $h(y) = \sin(1+y) + \sin(1-y)$, correcte mais non réduite : l'intégration de $1/h$ échouera si $h$ n'est pas simplifiée en $2 \sin(1) \cos y$. À défaut de développement trigonométrique, signaler l'échec plutôt que renvoyer une intégrale non évaluée inutile.
* **Domaine** : $\sqrt{xy} = \sqrt{x}\sqrt{y}$ seulement pour $x, y > 0$ ; l'extraction en $(x_0, y_0) = (1, 1)$ est juste sur ce quadrant et fausse pour $x, y < 0$. Le domaine du résultat doit être celui du point d'extraction.
* **Faux négatif de nullité** : si `simplifier` ne prouve pas $R \equiv 0$, l'étage numérique décide ; la réinjection finale couvre le risque.

### 5.3 Algorithme de résolution et solutions perdues

```text
ResoudreSeparable(g, h) :
  equilibres <- racines réelles de h(y) = 0          # Solveur, en y seul
  H <- IntegrerEnY(1/h) ; G <- IntegrerEnX(g)
  si H ou G non évaluée : renvoyer forme implicite non évaluée (ou ECHEC)
  essayer d'isoler y dans H(y) = G(x) + C            # sinon : implicite H(y) - G(x) = C
  ajouter chaque équilibre y0 absent de la famille (vérifier : la famille
     le contient-elle pour une valeur de C ou de K ?)
  si h n'est pas lipschitzienne en un équilibre : signaler la non-unicité
```

Complexité : dominée par les deux intégrations et l'inversion de $H$.

### 5.4 Exemples travaillés (vérifiés)

**S1.** $y' = x y^2$. Ici $g = x$, $h = y^2$, équilibre $y = 0$. Pour $y \neq 0$ : $-1/y = x^2/2 + C$, d'où
$$y = -\frac{2}{x^2 + C}.$$
Réinjection : $y' = \dfrac{4x}{(x^2 + C)^2}$ et $x y^2 = x \cdot \dfrac{4}{(x^2+C)^2}$. ✓ La solution $y \equiv 0$ est **perdue** (elle correspond à « $C = \infty$ »).

**S4.** $y' = e^{x+y}$. Test : $f f_{xy} - f_x f_y = e^{2(x+y)} - e^{2(x+y)} = 0$. Puis $-e^{-y} = e^x + C$, soit $y = -\ln(K - e^x)$ sur $\{x : e^x < K\}$, $K > 0$. Réinjection : $y' = \dfrac{e^x}{K - e^x}$ et $e^{x+y} = e^x \cdot \dfrac{1}{K - e^x}$. ✓ Pas d'équilibre ; pour $K \le 0$, aucune solution (domaine vide) : les solutions explosent en temps fini $x = \ln K$.

## 6. Équation linéaire

### 6.1 Théorème

$y' + a(x)\,y = b(x)$ avec $a, b$ continues sur $I$. Soit $A$ une primitive de $a$. Les solutions sur $I$ sont exactement
$$y = e^{-A(x)}\left(\int b(x)\,e^{A(x)}\,dx + C\right), \qquad C \in \mathbb{R}.$$
Elles sont définies sur **tout** $I$ et il n'y a **aucune** solution perdue.

*Preuve.* $(y\,e^{A})' = (y' + a y)\,e^{A} = b\,e^{A}$. ∎ La **variation de la constante** donne la même formule : on cherche $y = c(x)\,e^{-A}$, d'où $c' = b\,e^{A}$.

### 6.2 Test de reconnaissance

$f$ est affine en $y$ si et seulement si $f_y$ ne dépend pas de $y$. Alors $a = -f_y$ et $b = f - f_y\,y$.

```text
Lineaire(f) :
  a <- -simplifier(f_y)
  si contient(a, y) : ECHEC                     # structurel ; NE PAS remplacer par Nul(a_y)
                                                # (faux positif sur |y|, voir le piège)
  b <- simplifier(f + a*y)
  si contient(b, y) : ECHEC                     # garde-fou
  renvoyer (a, b)
```

**Piège** : avec la valeur absolue (fonction `abs`, présente sur la branche `feature/fonctions-usuelles` et sur celle de Fable, pas encore sur `main`), $f = |y|$ a pour dérivée $\mathrm{sign}(y)$ dont la dérivée est nulle presque partout ; le test « $f_{yy} \equiv 0$ » seul conclurait à tort à la linéarité. Exiger que $a$ **ne contienne pas** $y$ (test structurel) évite ce faux positif.

### 6.3 Algorithme et coût

Deux intégrations : $A = \int a$, puis $\int b\,e^{A}$. La seconde peut être non élémentaire (par exemple $\int e^{x^2}\,dx$) : renvoyer alors la formule avec une intégrale non évaluée est **correct** (c'est une vraie solution), à condition de le signaler.

### 6.4 Exemples (vérifiés)

**L1.** $y' + y/x = x^2$ sur $x > 0$ ou $x < 0$. $A = \ln|x|$, $e^{A} = |x|$ ; sur chaque demi-droite $(x y)' = x^3$, d'où $y = \dfrac{x^3}{4} + \dfrac{C}{x}$. Réinjection : $y' + y/x = \tfrac{3x^2}{4} - \tfrac{C}{x^2} + \tfrac{x^2}{4} + \tfrac{C}{x^2} = x^2$. ✓

**L3.** $y' + y \tan x = \cos x$ sur $]-\pi/2, \pi/2[$. $e^{A} = 1/\cos x$, $(y/\cos x)' = 1$, d'où $y = (x + C)\cos x$. Réinjection : $\cos x - (x + C)\sin x + (x + C)\sin x = \cos x$. ✓

## 7. Équation de Bernoulli

### 7.1 Théorème

$y' + a(x)\,y = b(x)\,y^n$, $n \in \mathbb{R} \setminus \{0, 1\}$. Sur un intervalle où $y > 0$ (ou $y \neq 0$ si $n$ est entier), $z = y^{1-n}$ vérifie l'équation **linéaire**
$$z' + (1 - n)\,a(x)\,z = (1 - n)\,b(x).$$
*Preuve.* $z' = (1 - n)\,y^{-n} y' = (1-n)\,y^{-n}(b y^n - a y) = (1-n)(b - a z)$. ∎

### 7.2 Solutions perdues et domaine

* Si $n > 0$, $y \equiv 0$ est solution ; elle n'est pas dans la famille $y = z^{1/(1-n)}$ (sauf cas particulier) : **perdue**.
* Si $0 < n < 1$, $y \mapsto y^n$ n'est pas lipschitzienne en 0 : non-unicité, raccordements avec $y = 0$ (comme en 3.3).
* Si $n < 0$, $y = 0$ n'est pas dans le domaine.
* Retour à $y$ : $y = z^{1/(1-n)}$ exige $z > 0$ si $1/(1-n)$ n'est pas un entier ; si $1 - n$ est pair (par exemple $n = 3$, $z = y^{-2}$), il y a **deux branches** $y = \pm z^{-1/2}$, le signe de $y$ étant conservé le long d'une solution.

### 7.3 Test de reconnaissance

Pour $f = \alpha(x)\,y + \beta(x)\,y^n$, on a $Q = y\,f_y - f = (n - 1)\,\beta(x)\,y^n$, donc
$$n = \frac{y\,Q_y}{Q}.$$

```text
Bernoulli(f) :
  Q <- simplifier(y*f_y - f)
  si Nul(Q) : c'est f = alpha(x)*y (linéaire, déjà traité)
  n <- simplifier(y*Q_y / Q)
  si n contient x ou y, ou n ∈ {0, 1} : ECHEC
  beta  <- Q / ((n - 1) * y^n)
  alpha <- (f - beta*y^n) / y
  si contient(beta, y) ou contient(alpha, y) : ECHEC
  renvoyer (n, alpha, beta)                     # y' = alpha y + beta y^n
```

Pièges : $n$ doit rester **exact** (`Nombre` rationnel : $1/2$, pas $0{,}5$), sinon $y^{1-n}$ devient approché ; un troisième terme ($f = \alpha y + \beta y^n + \gamma$) rend $n$ non constant, et le test échoue comme il faut.

### 7.4 Exemples (vérifiés)

**B1.** $y' + y/x = x y^2$ ($n = 2$). $z = 1/y$ : $z' - z/x = -x$, $(z/x)' = -1$, $z = Cx - x^2$, d'où
$$y = \frac{1}{Cx - x^2}.$$
Réinjection : $y' = \dfrac{2x - C}{(Cx - x^2)^2}$ et
$$x y^2 - \frac{y}{x} = \frac{x}{(Cx - x^2)^2} - \frac{1}{x(Cx - x^2)} = \frac{x^2 - (Cx - x^2)}{x\,(Cx - x^2)^2} = \frac{2x - C}{(Cx - x^2)^2}. \;✓$$
Perdue : $y \equiv 0$ (limite $C \to \infty$).

**B2 (domaine).** $y' = y + \sqrt{y}$ ($n = 1/2$, $y \ge 0$). $z = \sqrt{y}$ : $z' = z/2 + 1/2$, $z = C e^{x/2} - 1$, d'où $y = (C e^{x/2} - 1)^2$ **seulement là où** $C e^{x/2} \ge 1$ (car $z = \sqrt{y} \ge 0$). Hors de ce domaine, la formule n'est pas solution : en $x = 0$, $C = 1/2$, le résidu $y' - y - \sqrt{y}$ vaut $-1$ (calculé). Perdue et singulière : $y \equiv 0$, avec raccordements $y = 0$ pour $x \le x_1$ puis $y = (e^{(x - x_1)/2} - 1)^2$.

**B3 (deux branches).** $y' = y - y^3$ ($n = 3$). $z = y^{-2}$ : $z' = -2z + 2$, $z = 1 + C e^{-2x}$, d'où $y = \pm \dfrac{1}{\sqrt{1 + C e^{-2x}}}$. Perdue : $y \equiv 0$. Les équilibres $y = \pm 1$ sont dans la famille ($C = 0$).

## 8. Équation $y' = F(ax + by + c)$

**Théorème.** Si $f(x, y) = F(\lambda x + y)$ ($\lambda$ constant), alors $z = \lambda x + y$ vérifie l'équation **autonome** $z' = \lambda + F(z)$, donc séparable. Les zéros $z_0$ de $\lambda + F$ donnent les droites $y = z_0 - \lambda x$, solutions **perdues** par la séparation.

**Test.** $f = F(\lambda x + y)$ avec $F$ dérivable et $f_y \neq 0$ si et seulement si $f_x / f_y$ est une constante $\lambda$ (alors $f$ est constante le long de la direction $(1, -\lambda)$, dont l'invariant est $\lambda x + y$).

```text
FonctionAffine(f) :
  si Nul(f_y) : ECHEC (quadrature, déjà traité)
  lambda <- simplifier(f_x / f_y)
  si lambda contient x ou y : ECHEC
  F(z) <- substituer(f, y, z - lambda*x) ; vérifier que F ne contient plus x
```

**Exemple S6 (vérifié).** $y' = \sin(x + y)$. Le test de séparabilité échoue ($f f_{xy} - f_x f_y = -\sin^2 - \cos^2 = -1$) ; ici $\lambda = 1$, $z = x + y$, $z' = 1 + \sin z$. Avec $t = \tan(z/2)$ : $\int \frac{dz}{1 + \sin z} = -\frac{2}{1 + t}$ (à une constante près), d'où, après regroupement des constantes,
$$y = 2 \arctan\frac{1 + x + C}{1 - x - C} - x \quad (\text{modulo } 2\pi).$$
Perdues : $1 + \sin z = 0$, soit $y = -x - \pi/2 + 2k\pi$ ; réinjection : $y' = -1 = \sin(-\pi/2)$. ✓

## 9. Équations homogènes et réductibles à homogènes

### 9.1 Homogène $y' = F(y/x)$

$f$ est **homogène de degré 0** si $f(tx, ty) = f(x, y)$ pour tout $t > 0$. Alors, pour $x > 0$, $f(x, y) = F(y/x)$ avec $F(u) = f(1, u)$ ; pour $x < 0$, $F(u) = f(-1, -u)$, qui peut différer si l'homogénéité n'a lieu que pour $t > 0$.

**Théorème.** Avec $y = u\,x$ : $x\,u' = F(u) - u$, équation séparable. Les racines $u_0$ de $F(u) = u$ donnent les droites $y = u_0 x$ (**perdues**). Si $F(u) \equiv u$, l'équation est $y' = y/x$ et $y = Cx$.

**Test.**
```text
Homogene(f) :
  g <- simplifier(substituer(f, y, u*x))        # u : nouvelle variable
  si non Nul(g_x) : ECHEC
  tester aussi numériquement en x < 0 et en x > 0 (piège |x|) ;
  F(u) <- substituer(g, x, 1)  (et substituer(g, x, -1) pour x < 0 si différent)
```
Variante : identité d'Euler $x f_x + y f_y \equiv 0$ (pour $f$ différentiable, équivalente à l'invariance pour $t > 0$ sur un domaine conique).

**Piège.** $f = \sqrt{x^2 + y^2}/x$ : $g = \sqrt{x^2(1 + u^2)}/x = \frac{|x|}{x}\sqrt{1 + u^2}$. La dérivée formelle de $|x|/x$ est nulle là où elle existe ; le test passe, mais $F$ dépend du signe de $x$. Il faut deux résolutions (une par demi-plan).

### 9.2 Exemples (vérifiés)

**H1.** $y' = (y^2 + xy)/x^2$, $F(u) = u^2 + u$, $x u' = u^2$, $-1/u = \ln|x| + C$ :
$$y = -\frac{x}{\ln|x| + C}.$$
Réinjection (sur $x > 0$ et sur $x < 0$) : ✓. Perdue : $u_0 = 0$, $y \equiv 0$.

**H2.** $y' = \dfrac{x + y}{x - y}$, $F(u) = \frac{1 + u}{1 - u}$, $x u' = \frac{1 + u^2}{1 - u}$ ; solution **implicite**
$$\arctan\frac{y}{x} - \frac{1}{2}\ln(x^2 + y^2) = C$$
(spirales logarithmiques). Réinjection : $\Phi_x + \Phi_y f \equiv 0$. ✓ Pas de racine réelle de $F(u) = u$ ($u^2 + 1 = 0$).

**H3 (branches).** $y' = y/x + \tan(y/x)$ : $\cot u\,du = dx/x$, $\sin(y/x) = Cx$, d'où $y = x \arcsin(Cx)$ **et** $y = x(\pi - \arcsin(Cx))$ (et leurs translatées de $2k\pi x$). Les deux sont vérifiées par réinjection. Les droites $y = k\pi x$ ($\sin u = 0$, divisées lors de la séparation) correspondent à $C = 0$ dans ces familles translatées : aucune n'est perdue **si** le résultat garde toutes les branches ; si l'on ne renvoie que $y = x\arcsin(Cx)$, toutes les droites $k \neq 0$ sont perdues.

### 9.3 Réductibles à homogènes : $y' = F\!\left(\frac{a x + b y + c}{d x + e y + f}\right)$

* Si $ae - bd \neq 0$ : soit $(x_0, y_0)$ l'intersection des droites $ax + by + c = 0$ et $dx + ey + f = 0$. Avec $X = x - x_0$, $Y = y - y_0$, l'équation devient $\frac{dY}{dX} = F\!\left(\frac{aX + bY}{dX + eY}\right)$, homogène.
* Si $ae - bd = 0$ : les deux formes affines sont proportionnelles à une même $\lambda x + y$ (si $b$ ou $e$ est non nul), et l'on est dans le cas de la section 8.

**Reconnaissance.** On ne cherche pas $F$ : on repère dans l'AST les formes affines (sous-expressions `Somme` de degré 1 en $x, y$), on calcule l'intersection de deux d'entre elles, on translate (`substituer(f, x, X + x0)`, idem pour $y$) et l'on applique le test `Homogene`. Si le déterminant est nul, le test de la section 8 s'applique directement à $f$.

**R1 (vérifié).** $y' = \dfrac{x + y - 2}{x - y}$ : intersection $(1, 1)$, puis H2 translaté :
$$\arctan\frac{y - 1}{x - 1} - \frac{1}{2}\ln\big((x-1)^2 + (y-1)^2\big) = C. \;✓$$

**R2 (vérifié, cas parallèle).** $y' = \dfrac{x + y + 1}{x + y - 1}$ : $z = x + y$, $z' = \frac{2z}{z - 1}$, $\frac{z - 1}{2z}\,dz = dx$, d'où la solution implicite $y - x - \ln|x + y| = C$. Perdue : $z = 0$, soit $y = -x$ ; réinjection : $y' = -1$ et $\frac{0 + 1}{0 - 1} = -1$. ✓

## 10. Équations exactes et facteurs intégrants

### 10.1 Équation exacte

**Théorème (lemme de Poincaré pour les 1-formes).** Si $M, N$ sont de classe $C^1$ sur un ouvert **étoilé** (ou simplement connexe) $U$ et $M_y = N_x$ sur $U$, il existe $\Phi$ de classe $C^2$ avec $\Phi_x = M$, $\Phi_y = N$. Les solutions de $M\,dx + N\,dy = 0$ sont alors les courbes de niveau $\Phi(x, y) = C$ ; là où $N \neq 0$, elles définissent $y$ comme fonction de $x$ (théorème des fonctions implicites).

Sur un ouvert non simplement connexe, $M_y = N_x$ ne donne qu'un potentiel **local** : $\frac{-y\,dx + x\,dy}{x^2 + y^2}$ est fermée sur $\mathbb{R}^2 \setminus \{0\}$ sans y être exacte, mais l'algorithme ci-dessous renvoie $\Phi = -\arctan(x/y)$ (valable pour $y \neq 0$), ce qui suffit pour une EDO (les solutions sont les demi-droites issues de l'origine).

```text
Potentiel(M, N) :
  si non Nul(M_y - N_x) : ECHEC
  A  <- IntegrerEnX(M)                      # y gelé
  hp <- simplifier(N - A_y)
  si contient(hp, x) : ECHEC                # (incohérence : ne devrait pas arriver)
  h  <- IntegrerEnY(hp)
  renvoyer Phi = A + h
```
Coût : une dérivée croisée, deux intégrations. Échec si l'une des intégrales est non évaluée.

**E1 (vérifié).** $(2xy + 3)\,dx + (x^2 + 4y)\,dy = 0$ : $M_y = 2x = N_x$, $A = x^2 y + 3x$, $h' = 4y$, $\Phi = x^2 y + 3x + 2y^2 = C$. Isoler $y$ donne deux branches $y = \frac{-x^2 \pm \sqrt{x^4 - 8(3x - C)}}{4}$ (trinôme en $y$).

### 10.2 Facteurs intégrants

Un **facteur intégrant** est une fonction $\mu \neq 0$ telle que $\mu M\,dx + \mu N\,dy$ soit exacte : $(\mu M)_y = (\mu N)_x$. Dans chaque cas ci-dessous, on suppose $\mu = \mu(s)$ pour une variable $s$ fixée, on écrit la condition et l'on obtient $\mu'/\mu$ ; le critère est que ce rapport ne dépende que de $s$ (tous vérifiés symboliquement, voir l'annexe).

| Forme de $\mu$ | Condition (doit dépendre de $s$ seul) | Facteur |
| :--- | :--- | :--- |
| $\mu(x)$, $s = x$ | $\phi = \dfrac{M_y - N_x}{N}$ | $\mu = \exp\int \phi\,dx$ |
| $\mu(y)$, $s = y$ | $\psi = \dfrac{N_x - M_y}{M}$ | $\mu = \exp\int \psi\,dy$ |
| $\mu(x + y)$, $s = x + y$ | $\phi = \dfrac{M_y - N_x}{N - M}$ | $\mu = \exp\int \phi(s)\,ds$ |
| $\mu(xy)$, $s = xy$ | $\phi = \dfrac{M_y - N_x}{yN - xM}$ | $\mu = \exp\int \phi(s)\,ds$ |
| $M, N$ homogènes de même degré | $xM + yN \not\equiv 0$ | $\mu = \dfrac{1}{xM + yN}$ |

*Dérivation pour $\mu(x + y)$* : $\mu_x = \mu_y = \mu'$, donc $\mu' M + \mu M_y = \mu' N + \mu N_x$, soit $\mu'/\mu = (M_y - N_x)/(N - M)$. *Pour $\mu(xy)$* : $\mu_y = x\mu'$, $\mu_x = y\mu'$, donc $\mu'/\mu = (M_y - N_x)/(yN - xM)$. Le cas homogène a été vérifié pour $M, N$ homogènes de degré 2 à coefficients symboliques ; c'est un résultat classique pour tout degré commun.

**Test « dépend de $s$ seul ».** Substituer $y = s - x$ (resp. $y = s/x$, avec $x \neq 0$) dans le rapport, simplifier et tester `Nul` de la dérivée en $x$. Pour $s = x$ : le rapport ne doit pas contenir $y$.

**Solutions perdues ou ajoutées.** Multiplier par $\mu$ peut **ajouter** des courbes ($\mu = 0$) et **perdre** des solutions (là où $\mu$ est infini ou non défini). Les candidats $\{\mu = 0\} \cup \{1/\mu = 0\}$ sont réinjectés dans l'équation **de départ**.

### 10.3 Exemples (vérifiés)

**E2, $\mu(x)$.** $(3xy + y^2)\,dx + (x^2 + xy)\,dy = 0$ : $\frac{M_y - N_x}{N} = \frac{x + y}{x(x + y)} = \frac{1}{x}$, $\mu = x$, $\Phi = x^3 y + \frac{x^2 y^2}{2} = C$. $\mu = 0$ : la droite $x = 0$ est solution de la forme de départ ($N(0, y) = 0$) et correspond à $C = 0$.

**E3, $\mu(y)$.** $y\,dx + (2x - y e^y)\,dy = 0$ : $\frac{N_x - M_y}{M} = \frac{1}{y}$, $\mu = y$, $\Phi = x y^2 - (y^2 - 2y + 2)\,e^y = C$. $y = 0$ est solution de départ et correspond à $C = -2$.

**E4, $\mu(xy)$.** $y(1 + xy)\,dx + x(1 - xy)\,dy = 0$ : $\frac{M_y - N_x}{yN - xM} = \frac{4xy}{-2x^2y^2} = -\frac{2}{xy}$, donc $\mu = (xy)^{-2}$ et
$$\Phi = -\frac{1}{xy} + \ln\left|\frac{x}{y}\right| = C.$$
**Perdues** : $y = 0$ (et la verticale $x = 0$), où $\mu$ est infini ; $y \equiv 0$ vérifie bien la forme de départ.

**E5, homogène.** $(x^2 + y^2)\,dx - 2xy\,dy = 0$ : $\mu = \frac{1}{x(x^2 - y^2)}$ ; on trouve $\frac{x^2 - y^2}{x} = C$. Les droites $y = \pm x$ (où $\mu$ est infini) sont solutions et correspondent à $C = 0$.

**Contrôles du critère $\mu(x + y)$** (exemples construits pour le test) : $M = \frac{2x}{x + y}$, $N = \frac{1}{x + y}$ donne $\phi = \frac{1}{x + y}$, $\mu = x + y$ ; $M = y e^{-x-y}$, $N = x e^{-x-y}$ donne $\phi = 1$, $\mu = e^{x + y}$.

## 11. Équation de Riccati

### 11.1 Réduction par une solution particulière

$y' = q_0(x) + q_1(x)\,y + q_2(x)\,y^2$, $q_2 \not\equiv 0$. Si $y_1$ est une solution particulière, alors $y = y_1 + 1/v$ est solution si et seulement si
$$v' = -\big(q_1 + 2 q_2 y_1\big)\,v - q_2,$$
équation **linéaire** en $v$ (vérifié symboliquement). La solution $y_1$ elle-même correspond à $v = \infty$ : elle est **perdue** par la formule et doit être ajoutée.

**Deux solutions particulières** $y_1 \neq y_2$ : $w = \frac{y - y_1}{y - y_2}$ vérifie $\frac{w'}{w} = q_2\,(y_1 - y_2)$, donc la solution générale s'obtient par **une quadrature**. **Trois solutions** : le birapport de quatre solutions est constant (aucune quadrature). Ces deux énoncés sont vérifiés symboliquement (annexe).

### 11.2 Transformation en EDO linéaire du second ordre

Avec $y = -\dfrac{u'}{q_2\,u}$ :
$$y' - (q_0 + q_1 y + q_2 y^2) = -\frac{1}{q_2 u}\left(u'' - \Big(q_1 + \frac{q_2'}{q_2}\Big) u' + q_0 q_2\,u\right),$$
donc $y$ est solution si et seulement si $u'' - (q_1 + q_2'/q_2)\,u' + q_0 q_2\,u = 0$ (vérifié). Si $q_1 + q_2'/q_2$ et $q_0 q_2$ sont **constants**, cette équation est exactement du ressort d'`EquationDifferentielle` (coefficients constants). La solution générale $u = C_1 u_1 + C_2 u_2$ donne $y$ avec **une seule** constante effective ($C_1 : C_2$), les cas $C_1 = 0$ ou $C_2 = 0$ redonnant les solutions particulières.

### 11.3 Heuristiques pour une solution particulière

Aucune n'est complète ; chaque candidat est **vérifié par réinjection**.
1. **Constante** $y_1 = c$ : résoudre $q_0 + q_1 c + q_2 c^2 = 0$ en un point $x_1$ (trinôme en $c$), puis vérifier l'identité en $x$.
2. **Monôme** $y_1 = c\,x^k$ (en particulier $c/x$) : substituer, égaler les puissances de $x$.
3. **Polynôme** de degré $d$ si les $q_i$ sont des polynômes : l'équilibre des degrés dominants parmi $\deg q_2 + 2d$, $\deg q_1 + d$, $\deg q_0$ et $d - 1$ fournit les $d$ candidats (typiquement $d = \frac{\deg q_0 - \deg q_2}{2}$ ou $d = \deg q_1 - \deg q_2$) ; puis coefficients indéterminés (système polynomial, résolu du coefficient dominant vers le bas).
4. Si $q_0, q_1, q_2$ sont constants, l'équation est **autonome**, donc séparable : ne pas passer par Riccati.

### 11.4 Exemples (vérifiés)

**Ri1.** $y' = y^2 - 2/x^2$. Monôme $y_1 = a/x$ : $-a = a^2 - 2$, donc $a = 1$ ou $a = -2$. Avec $y_1 = 1/x$ : $v' = -\frac{2}{x}v - 1$, $(x^2 v)' = -x^2$, $v = \frac{C}{x^2} - \frac{x}{3}$, d'où
$$y = \frac{1}{x} + \frac{3x^2}{K - x^3}.$$
$K = 0$ redonne $y = -2/x$ ; $y = 1/x$ est **perdue** ($K = \infty$). Par la section 11.2 : $u'' - \frac{2}{x^2}u = 0$ (Euler), $u = C_1 x^2 + C_2/x$, et $y = -u'/u$ redonne la même famille. Toutes les formes vérifiées par réinjection.

**Ri2 (intégrale non élémentaire).** $y' = y^2 - xy + 1$, $y_1 = x$ ($1 = x^2 - x^2 + 1$). $v' = -xv - 1$, $\big(v e^{x^2/2}\big)' = -e^{x^2/2}$ :
$$y = x + \frac{e^{x^2/2}}{C - \int_0^x e^{t^2/2}\,dt}, \qquad \int_0^x e^{t^2/2}\,dt = \sqrt{\pi/2}\;\mathrm{erfi}(x/\sqrt{2}).$$
Résultat **correct** avec une intégrale non évaluée (vérifié avec la forme $\mathrm{erfi}$). Perdue : $y = x$.

**Ri3 (plusieurs classes).** $y' = (y - x)^2 + 1$ est une Riccati ($y_1 = x$, $v' = -1$) et aussi une $F(\lambda x + y)$ avec $\lambda = -1$ : $y = x + \dfrac{1}{C - x}$. Perdue : $y = x$.

**Échec attendu.** $y' = x^2 + y^2$ n'a pas de solution élémentaire (exemple classique de Liouville ; la solution s'exprime par des fonctions de Bessel d'ordre $1/4$). Aucune heuristique ne trouve de solution particulière : le solveur doit renvoyer « non évalué ».

## 12. Équations de Clairaut et de Lagrange

### 12.1 Clairaut $y = x\,p + g(p)$

En dérivant : $p = p + (x + g'(p))\,p'$, donc $p' = 0$ ou $x + g'(p) = 0$.
* $p' = 0$ : **famille de droites** $y = Cx + g(C)$, pour $C$ dans le domaine de $g$.
* $x = -g'(p)$ : **solution singulière** paramétrée $x = -g'(p)$, $y = g(p) - p\,g'(p)$ (si $g'' \neq 0$), enveloppe de la famille. Elle coïncide avec le $p$-discriminant ($F_p = x + g'(p) = 0$).
* Des solutions $C^1$ **composées** existent (un arc d'enveloppe prolongé par ses tangentes) ; le solveur ne les liste pas mais peut les mentionner.

**Cl1 (vérifié).** $y = xp + p^2$ : $y = Cx + C^2$ ; singulière $x = -2p$, $y = -p^2$, soit $y = -x^2/4$. Réinjection : $p = -x/2$, $xp + p^2 = -x^2/2 + x^2/4 = -x^2/4$. ✓

**Cl2 (vérifié).** $y = xp + \sqrt{1 + p^2}$ : $y = Cx + \sqrt{1 + C^2}$ ; singulière $y = \sqrt{1 - x^2}$, $\lvert x \rvert < 1$ (demi-cercle).

**Cl3 (vérifié).** $y = xp - \ln p$ : $y = Cx - \ln C$ ($C > 0$) ; singulière $y = 1 + \ln x$ ($x > 0$).

### 12.2 Lagrange (d'Alembert) $y = x\,f(p) + g(p)$, $f(p) \not\equiv p$

En dérivant : $p = f(p) + \big(x f'(p) + g'(p)\big)\,p'$. Là où $p' \neq 0$, on prend $p$ comme paramètre et $x$ comme fonction de $p$ :
$$\frac{dx}{dp} - \frac{f'(p)}{p - f(p)}\,x = \frac{g'(p)}{p - f(p)},$$
équation **linéaire** en $x(p)$. Sa solution $x = X(p, C)$ donne la solution générale **paramétrée** $x = X(p, C)$, $y = X(p, C)\,f(p) + g(p)$ ; on élimine $p$ quand c'est possible.

**Solutions perdues.** Les racines $p_0$ de $p = f(p)$ donnent les droites $y = p_0 x + g(p_0)$ (cas $p' = 0$), à réinjecter.

**La1 (vérifié).** $y = 2xp + p^2$ : $\frac{dx}{dp} = -\frac{2x}{p} - 2$, $(p^2 x)' = -2p^2$, d'où
$$x = \frac{C}{p^2} - \frac{2p}{3}, \qquad y = \frac{2C}{p} - \frac{p^2}{3}.$$
Vérifié par $Y'(p) = p\,X'(p)$ et $Y - 2Xp - p^2 \equiv 0$. Perdue : $p_0 = 0$, $y \equiv 0$.

**La2 (vérifié, élimination possible).** $y = x p^2$ ($f = p^2$, $g = 0$) : $(1 - p) = 2x\,p'$ pour $p \neq 0$, d'où $p = 1 - C/\sqrt{x}$ et $y = (\sqrt{x} - C)^2$ ($x > 0$). Racines de $p = p^2$ : $p_0 = 0$ ($y \equiv 0$, **perdue**) et $p_0 = 1$ ($y = x$, dans la famille avec $C = 0$).

### 12.3 Reconnaissance

```text
ClairautLagrange(F) :                        # F(x, y, p) = 0
  si Nul(F_y) ou contient(F_y, y) : ECHEC    # F doit être affine en y
  G <- simplifier(y - F / F_y)               # y = G(x, p) (diviser par F_y peut perdre
                                             # des solutions là où F_y = 0 : les réinjecter)
  si contient(G, y) : ECHEC
  fp <- simplifier(G_x)
  si contient(fp, x) : ECHEC                 # G n'est pas affine en x
  g <- simplifier(G - x*fp)
  si Nul(fp - p) : Clairaut(g) sinon Lagrange(fp, g)
```

## 13. Stratégie d'un solveur : ordre des tests et échec

Principe : essayer d'abord les classes **les moins coûteuses** et celles qui donnent une solution **explicite et complète**, puis les plus générales.

```text
ResoudreEDO1(equation) :
  0. Normaliser. Forme implicite F(x, y, p) = 0 :
       - si F est affine en p : f <- -F(x, y, 0) / F_p
       - sinon : essayer ClairautLagrange ; sinon, résoudre en p (branches) si possible ;
         sinon ECHEC.
  1. Nul(f_y)                         -> quadrature y = ∫ f dx + C
  2. Lineaire(f)                      -> formule 6.1 (aucune solution perdue)
  3. Bernoulli(f)                     -> linéaire en z = y^(1-n) ; + y = 0 si n > 0
  4. Separer(f)                       -> 5.3 ; + équilibres
  5. FonctionAffine(f) (y' = F(λx+y)) -> séparable en z ; + droites perdues
  6. Homogene(f), puis translation    -> séparable en u ; + droites y = u0 x
  7. Forme M dx + N dy (M = -num, N = den) :
       exacte ; μ(x) ; μ(y) ; μ(x+y) ; μ(xy) ; homogène 1/(xM + yN)
  8. Riccati (Nul(f_yyy)) :
       heuristiques 11.3 ; sinon transformation 11.2 si coefficients constants
  9. ECHEC -> résultat « non évalué »
  Pour chaque succès :
    - réinjecter la solution (4.4) ; en cas d'échec : traiter comme ECHEC de la méthode
    - réinjecter chaque candidat perdu (zéros des diviseurs, équilibres, p0)
    - joindre domaine et indication de non-unicité
```

Justifications :
* la **linéaire** passe avant la séparable : $y' = xy$ est des deux classes, mais la formule linéaire $y = C e^{x^2/2}$ contient déjà $y = 0$ ($C = 0$), alors que la séparation donne $\ln|y| = \ldots$ et une solution perdue à recoller ;
* **Bernoulli** passe avant la séparable : pour $y' = y - y^3$, Bernoulli donne une forme explicite, la séparation une forme implicite à inverser ;
* les tests 1 à 6 ne coûtent que quelques dérivées et substitutions ; les étapes 7 et 8 demandent des intégrations et des résolutions.

**Quand rien ne marche** : renvoyer un résultat « non évalué » qui garde l'équation (sur le modèle d'`IntegraleNonEvaluee`), jamais une solution non vérifiée, jamais d'exception pour une équation simplement hors de portée. Une résolution numérique (RK4 sur $y' = f(x, y)$) peut être proposée en complément, marquée approchée.

**Forme du résultat (proposition, à discuter)** : `{ forme : explicite | implicite | paramétrée | non évaluée ; expression(s) ; constantes ; solutionsPerdues ; domaine ; unicité ; vérifiée }`.

## 14. Pièges et cas de contrôle

Batterie destinée aux tests du solveur. Pour chaque entrée : la ou les classes attendues, la solution générale exacte, les solutions perdues ou singulières. **Toutes** les solutions ci-dessous ont été vérifiées par réinjection avec SymPy 1.14 (script en annexe, 70 contrôles, 0 échec) : simplification symbolique du résidu, ou, quand `simplify` ne conclut pas, évaluation du résidu à 30 chiffres en 6 points rationnels aléatoires du domaine (indiqué « num. »).

**Séparables**

1. $y' = x y^2$ — séparable, Bernoulli ($n = 2$). $y = -\dfrac{2}{x^2 + C}$. Perdue : $y = 0$.
2. $y' = (1 + y^2)\,e^x$ — séparable. $y = \tan(e^x + C)$. Aucune perdue ($1 + y^2 > 0$). Domaine : un intervalle où $\cos(e^x + C) \neq 0$ (explosion en temps fini).
3. $y' = xy + x + y + 1$ — **piège** : `Somme` canonique qui est $(x + 1)(y + 1)$ ; linéaire et séparable. $y = C e^{x^2/2 + x} - 1$ (contient $y = -1$, $C = 0$).
4. $y' = e^{x + y}$ — **piège** : séparable sans produit syntaxique. $y = -\ln(K - e^x)$, $K > 0$, $e^x < K$.
5. $y' = 3 y^{2/3}$ — séparable, **non-unicité**. $y = (x + C)^3$ (num., points où $y > 0$). Singulière : $y = 0$ et raccordements.
6. $y' = \sin(x + y)$ — **non séparable** (le test vaut $-1$) ; classe $F(\lambda x + y)$, $\lambda = 1$. $y = 2\arctan\frac{1 + x + C}{1 - x - C} - x$. Perdues : $y = -x - \pi/2 + 2k\pi$.
7. $y' = \sin(x + y) + \sin(x - y)$ — **piège** : séparable après simplification ($2 \sin x \cos y$ ; test $R = 0$). $\int \sec y\,dy = -2\cos x + C$, d'où $y = \arcsin\big(\tanh(C - 2\cos x)\big)$ et la symétrique $y = \pi - \arcsin\big(\tanh(C - 2\cos x)\big)$ (modulo $2\pi$). Perdues : $y = \pi/2 + k\pi$.

**Linéaires**

8. $y' + y/x = x^2$ — $y = \dfrac{x^3}{4} + \dfrac{C}{x}$, sur $x > 0$ ou $x < 0$.
9. $y' - 2xy = x$ — $y = C e^{x^2} - \dfrac{1}{2}$.
10. $y' + y \tan x = \cos x$ — $y = (x + C)\cos x$ sur $]-\pi/2, \pi/2[$.
11. $y' = y + x$ — $y = C e^{x} - x - 1$ (aussi : coefficients constants avec second membre, hors de portée de l'actuel `EquationDifferentielle`).
12. $y' = \lvert y \rvert$ — **piège** (avec `abs`) : **non linéaire** ; $f_{yy} = 0$ presque partout mais $f_y = \mathrm{sign}(y)$ contient $y$. Solutions : $y = K e^{x}$ pour $K \ge 0$, $y = K e^{-x}$ pour $K < 0$ (num. ; autonome, donc séparable : le signe de $y$ est conservé).

**Bernoulli**

13. $y' + y/x = x y^2$ — $y = \dfrac{1}{Cx - x^2}$. Perdue : $y = 0$.
14. $y' = y + \sqrt{y}$ — $n = 1/2$. $y = (C e^{x/2} - 1)^2$ **pour $C e^{x/2} \ge 1$** (num.) ; hors domaine, résidu $-1$ en $x = 0$, $C = 1/2$. Singulière : $y = 0$.
15. $y' = y - y^3$ — $n = 3$, deux branches $y = \pm\dfrac{1}{\sqrt{1 + C e^{-2x}}}$. Perdue : $y = 0$.
16. $x y' + y = y^2 \ln x$ — $n = 2$ ; $y = \dfrac{1}{\ln x + 1 + Cx}$, $x > 0$. Perdue : $y = 0$.

**Homogènes et réductibles**

17. $y' = (y^2 + xy)/x^2$ — $y = -\dfrac{x}{\ln\lvert x \rvert + C}$ (vérifié sur $x > 0$ et $x < 0$). Perdue : $y = 0$.
18. $y' = (x + y)/(x - y)$ — implicite $\arctan(y/x) - \tfrac{1}{2}\ln(x^2 + y^2) = C$.
19. $y' = y/x + \tan(y/x)$ — $y = x\arcsin(Cx)$ et $y = x(\pi - \arcsin(Cx))$, plus leurs translatées de $2k\pi x$. Les droites $y = k\pi x$ sont obtenues pour $C = 0$ si toutes les branches sont gardées ; sinon elles sont perdues.
20. $y' = y/x + x/y$ — $y = x\sqrt{2\ln x + C}$ ($x > 0$ ; et $y = -x\sqrt{\ldots}$ par symétrie $y \mapsto -y$).
21. $y' = \dfrac{x + y - 2}{x - y}$ — translation en $(1, 1)$ ; implicite $\arctan\frac{y-1}{x-1} - \tfrac{1}{2}\ln\big((x-1)^2 + (y-1)^2\big) = C$.
22. $y' = \dfrac{x + y + 1}{x + y - 1}$ — déterminant nul ; implicite $y - x - \ln\lvert x + y \rvert = C$. Perdue : $y = -x$.
23. $y' = \sqrt{x^2 + y^2}/x$ — **piège** : homogène pour $t > 0$ seulement ; $F$ dépend du signe de $x$ (deux résolutions). Non vérifié en tant que solution (contrôle de reconnaissance uniquement).

**Exactes et facteurs intégrants**

24. $(2xy + 3)\,dx + (x^2 + 4y)\,dy = 0$ — exacte ; $x^2 y + 3x + 2y^2 = C$.
25. $(3xy + y^2)\,dx + (x^2 + xy)\,dy = 0$ — $\mu = x$ ; $x^3 y + \tfrac{1}{2}x^2 y^2 = C$.
26. $y\,dx + (2x - y e^y)\,dy = 0$ — $\mu = y$ ; $x y^2 - (y^2 - 2y + 2)e^y = C$.
27. $y(1 + xy)\,dx + x(1 - xy)\,dy = 0$ — $\mu = (xy)^{-2}$ (critère $-2/(xy)$) ; $-\frac{1}{xy} + \ln\lvert x/y \rvert = C$. Perdue : $y = 0$.
28. $(x^2 + y^2)\,dx - 2xy\,dy = 0$ — homogène, $\mu = 1/(xM + yN)$ ; $(x^2 - y^2)/x = C$ ; $y = \pm x$ pour $C = 0$.
29. $\frac{2x}{x+y}\,dx + \frac{1}{x+y}\,dy = 0$ — critère $\mu(x + y)$ : $\phi = 1/(x+y)$, $\mu = x + y$ (exacte après multiplication).

**Riccati**

30. $y' = y^2 - 2/x^2$ — $y = \dfrac{1}{x} + \dfrac{3x^2}{K - x^3}$ ; $K = 0$ : $y = -2/x$. Perdue : $y = 1/x$.
31. $y' = y^2 - xy + 1$ — $y = x + \dfrac{e^{x^2/2}}{C - \sqrt{\pi/2}\,\mathrm{erfi}(x/\sqrt{2})}$. Perdue : $y = x$.
32. $y' = (y - x)^2 + 1$ — $y = x + \dfrac{1}{C - x}$. Perdue : $y = x$.
33. $y' = x^2 + y^2$ — **échec attendu** : « non évalué » (pas de solution élémentaire ; résultat classique, non revérifié ici).

**Clairaut et Lagrange** (forme $F(x, y, p) = 0$)

34. $y = xp + p^2$ — $y = Cx + C^2$ ; singulière $y = -x^2/4$.
35. $y = xp + \sqrt{1 + p^2}$ — $y = Cx + \sqrt{1 + C^2}$ ; singulière $y = \sqrt{1 - x^2}$ (num.).
36. $y = xp - \ln p$ — $y = Cx - \ln C$ ; singulière $y = 1 + \ln x$ (num.).
37. $y = 2xp + p^2$ — paramétrée $x = \frac{C}{p^2} - \frac{2p}{3}$, $y = \frac{2C}{p} - \frac{p^2}{3}$. Perdue : $y = 0$.
38. $y = x p^2$ — $y = (\sqrt{x} - C)^2$. Perdue : $y = 0$ ; $y = x$ pour $C = 0$.

**Contrôles des tests de reconnaissance** (valeurs exactes calculées)

* $f f_{xy} - f_x f_y$ : $0$ pour $(x+1)(y+1)$ développé, $e^{x+y}$, $\sin(x+y) + \sin(x-y)$ ; $-1$ pour $\sin(x+y)$ et pour $x + y$.
* Bernoulli : $y\,Q_y / Q = n$ pour $f = a(x)\,y + b(x)\,y^n$ général.
* Critères de facteur intégrant : $-2/(xy)$ (cas 27), $1/(x+y)$ (cas 29), $1$ pour $M = y e^{-x-y}$, $N = x e^{-x-y}$.

## 15. État dans le code

Références vérifiées dans le worktree `theorie` (commit `569e686`, identique à `origin/main`).

| Élément | État | Où |
| :--- | :--- | :--- |
| Classe `EquationDifferentielle` : $\sum a_i y^{(i)} = 0$, coefficients `double` | implémenté (linéaire, homogène, coefficients constants uniquement) | `include/EquationDifferentielle.hpp:26`, stockage `m_terme` `:118` |
| Solution générale littérale (constantes `C1..Cn`) | implémenté pour le cas ci-dessus ; racines du polynôme caractéristique numériques (Eigen) | `src/EquationDifferentielle.cpp:207` (`baseDeSolutions`), `:235` (`resoudreLitteral`) |
| Problème de Cauchy en $x = 0$ | implémenté pour le cas ci-dessus (système linéaire numérique, QR) | `src/EquationDifferentielle.cpp:244` |
| Résolution numérique RK4 | implémentée pour le cas ci-dessus seulement (pas pour $y' = f(x, y)$ général) | `src/EquationDifferentielle.cpp:276` |
| Représentation d'une EDO non linéaire $y' = f(x, y)$, $M\,dx + N\,dy$, $F(x, y, p)$ | **absent** | — |
| Classes 5 à 12 (séparable, linéaire à coefficients variables, Bernoulli, homogène, exacte, facteurs intégrants, Riccati, Clairaut, Lagrange) | **absent** | — |
| Dérivées partielles (tests de reconnaissance) | implémenté | `include/ASTNode.hpp:202-203`, `src/Noeud.cpp:233` |
| Plusieurs variables, `variables(e)`, `evaluer(e, valeurs)` (étage numérique du test de nullité) | implémenté | `include/Multivariable.hpp:40`, `:94` |
| `substituer`, `contient` | implémenté | `include/ASTNode.hpp:756`, `:762` ; `src/Noeud.cpp:1079`, `:1083` |
| Constantes symboliques `param("C")`, de dérivée nulle | implémenté | `include/ASTNode.hpp:683`, `src/Regles.cpp:71` |
| Intégration d'une expression à plusieurs variables | **refusée** (`IntegraleNonEvaluee`) : utiliser le gel de variable (4.3) | `src/Noeud.cpp:265-268` |
| Intégration par rapport à une variable autre que `x` | **partiel / piège** : `x` écrit en dur pour les constantes et le terme constant d'une somme | `src/Regles.cpp:58`, `:73`, `:106` ; `src/Noeud.cpp:266` |
| Isoler $y$ dans $H(y) = G(x) + C$ | **absent** : `resoudre` exige une seule variable et des solutions numériques | `src/Solveur.cpp:530` (`exigerUneVariable`), `:112` (solution non numérique => incomplet), moteur d'isolement `:92` (`Resolveur`), `:318` (`isoler`) |
| Racines de $h(y) = 0$ (équilibres) quand $h$ ne contient que $y$ | implémenté : `resoudre` accepte une variable unique quel que soit son nom (polynôme exact : Sturm ; sinon isolement, familles trigonométriques) | `src/Solveur.cpp:539`, `include/Polynome.hpp:58` |
| Polynôme en $y$ à coefficients fonctions de $x$ (Riccati, Clairaut) | **absent** : `Polynome::depuisExpression` exige des coefficients numériques | `include/Polynome.hpp:58` |
| Test de nullité : égalité structurelle (hash-consing) | implémenté | `include/ASTNode.hpp:251` (`estEgal`) |
| Simplification rationnelle à plusieurs variables | **absent** (ni sur `main`, ni sur la branche de Fable, qui se limite à une variable) | voir ci-dessous |

### Branche en cours (feature/primitives-substitution, commit 96de779)

Lu en lecture seule le 13/11/2026 ; une partie du travail est indexée mais non commitée (`src/Primitives.cpp`, `include/Primitives.hpp` non suivis).
* **Aucune EDO non linéaire** pour l'instant (recherche de « Bernoulli », « Riccati », « séparable », « Clairaut » : aucun résultat dans `src/`, `include/`, `tests/`).
* Intégration par méthodes successives : `ASTNode::integrer` déplacé dans `src/Primitives.cpp:559` (non commité).
* Le terme constant d'une somme est intégré par rapport à la variable de l'expression (`variableDe`, `src/Regles.cpp:99` de l'arbre de travail, non commité) : corrige en partie le piège de la section 4.3 ; `Constante::primitive` et l'intégrale d'une expression sans variable gardent `var("x")`, le renommage $y \to x$ reste donc nécessaire.
* Simplification rationnelle (`include/Simplification.hpp`, commitée) limitée aux expressions **à une variable** : elle n'aide pas les tests à deux variables.

## 16. Annexe : script de vérification

Lancé avec Python 3 et SymPy 1.14 ; sortie : « 70 contrôles, 0 échec(s) ». Les fonctions `explicite`, `implicite`, `forme_diff`, `implicite_F` et `parametrique` réalisent les réinjections de la section 4.4.

<details>
<summary>verifier_premier_ordre.py</summary>

```python
# Vérification par réinjection des exemples de docs/theorie/edo-premier-ordre.md
# Lancement : python3 -I verifier_premier_ordre.py (SymPy 1.14)
import random
from sympy import (tanh, symbols, Function, sqrt, exp, log, sin, cos, tan, atan, asin, erfi, pi,
                   simplify, diff, Rational, Abs, S, N)

x, y, C, K, p = symbols('x y C K p')
resultats = []


def nul(expr, substitutions_aleatoires, nom):
    """Vrai si expr est nulle : simplify, sinon contrôle numérique en plusieurs points."""
    e = simplify(expr)
    if e == 0:
        return "symbolique"
    for _ in range(6):
        vals = {s: f() for s, f in substitutions_aleatoires.items()}
        v = complex(N(expr.subs(vals), 30))
        if abs(v) > 1e-12:
            return f"ECHEC ({vals} -> {v})"
    return "numérique (6 points)"


def explicite(nom, f, sol, dom):
    """y' = f(x, y) ; sol = y(x) ; dom : {symbole: tirage}."""
    r = diff(sol, x) - f.subs(y, sol)
    resultats.append((nom, "y=" + str(sol), nul(r, dom, nom)))


def implicite(nom, f, phi, dom):
    """Solution Phi(x, y) = C de y' = f : Phi_x + Phi_y f = 0."""
    r = diff(phi, x) + diff(phi, y) * f
    resultats.append((nom, "Phi=" + str(phi), nul(r, dom, nom)))


def forme_diff(nom, M, N_, phi, dom):
    """M dx + N dy = 0, solution Phi = C : Phi_x N - Phi_y M = 0 (gradients colinéaires)."""
    r = diff(phi, x) * N_ - diff(phi, y) * M
    resultats.append((nom, "Phi=" + str(phi), nul(r, dom, nom)))


def implicite_F(nom, F, sol, dom):
    """F(x, y, p) = 0 avec y = sol(x)."""
    r = F.subs({y: sol, p: diff(sol, x)})
    resultats.append((nom, "y=" + str(sol), nul(r, dom, nom)))


def parametrique(nom, F, X, Y, dom):
    """Solution paramétrée par p : dY/dp = p dX/dp et F(X, Y, p) = 0."""
    r1 = diff(Y, p) - p * diff(X, p)
    r2 = F.subs({x: X, y: Y}, simultaneous=True)
    resultats.append((nom, "param", nul(r1, dom, nom) + " / " + nul(r2, dom, nom)))


random.seed(1)
u = lambda a, b: (lambda: Rational(random.randint(int(a * 100), int(b * 100)), 100))
dx_pos = {x: u(0.3, 2.0), C: u(3, 5), K: u(3, 5)}

# --- Séparables
explicite("S1 y'=x y^2", x * y**2, -2 / (x**2 + C), dx_pos)
explicite("S1 perdue y=0", x * y**2, S(0), dx_pos)
explicite("S2 y'=(1+y^2)e^x", (1 + y**2) * exp(x), tan(exp(x) + C), {x: u(-2, -1), C: u(-0.3, 0.3)})
explicite("S3 y'=xy+x+y+1", x * y + x + y + 1, C * exp(x**2 / 2 + x) - 1, dx_pos)
explicite("S4 y'=exp(x+y)", exp(x + y), -log(K - exp(x)), {x: u(-1, 1), K: u(3, 5)})
explicite("S5 y'=3y^(2/3) (y>=0)", 3 * y**Rational(2, 3), (x + C)**3, {x: u(0.1, 2), C: u(0.1, 1)})
explicite("S6 y'=sin(x+y)", sin(x + y), 2 * atan((1 + x + C) / (1 - x - C)) - x, {x: u(-0.4, 0.4), C: u(-0.3, 0.3)})
explicite("S6 perdue y=-x-pi/2", sin(x + y), -x - pi / 2, dx_pos)
explicite("S7 y'=sin(x+y)+sin(x-y)", sin(x + y) + sin(x - y), asin(tanh(C - 2 * cos(x))), {x: u(-2, 2), C: u(-1, 1)})
explicite("S7 symétrique pi-y", sin(x + y) + sin(x - y), pi - asin(tanh(C - 2 * cos(x))), {x: u(-2, 2), C: u(-1, 1)})
explicite("S7 perdue y=pi/2", sin(x + y) + sin(x - y), pi / 2, dx_pos)
# --- Linéaires
explicite("L1 y'+y/x=x^2", x**2 - y / x, x**3 / 4 + C / x, dx_pos)
explicite("L2 y'-2xy=x", x + 2 * x * y, C * exp(x**2) - Rational(1, 2), dx_pos)
explicite("L3 y'+y tan x=cos x", cos(x) - y * tan(x), (x + C) * cos(x), {x: u(-1.4, 1.4), C: u(-2, 2)})
explicite("L4 y'=y+x", y + x, C * exp(x) - x - 1, dx_pos)
explicite("L5 y'=|y| (K>0)", Abs(y), K * exp(x), dx_pos)
explicite("L5 y'=|y| (K<0)", Abs(y), -K * exp(-x), dx_pos)
# --- Bernoulli
explicite("B1 y'+y/x=xy^2", x * y**2 - y / x, 1 / (C * x - x**2), dx_pos)
explicite("B2 y'=y+sqrt(y) (Ce^(x/2)>=1)", y + sqrt(y), (C * exp(x / 2) - 1)**2, {x: u(0, 2), C: u(1.1, 3)})
# Hors du domaine : la formule n'est PAS solution (C e^(x/2) < 1)
r_hors = N((diff((C * exp(x / 2) - 1)**2, x) - ((C * exp(x / 2) - 1)**2 + sqrt((C * exp(x / 2) - 1)**2))).subs({x: 0, C: Rational(1, 2)}))
resultats.append(("B2 hors domaine (x=0, C=1/2)", "résidu attendu non nul", f"résidu = {r_hors}"))
explicite("B3 y'=y-y^3 (+)", y - y**3, 1 / sqrt(1 + C * exp(-2 * x)), dx_pos)
explicite("B3 y'=y-y^3 (-)", y - y**3, -1 / sqrt(1 + C * exp(-2 * x)), dx_pos)
explicite("B4 xy'+y=y^2 ln x", (y**2 * log(x) - y) / x, 1 / (log(x) + 1 + C * x), dx_pos)
# --- Homogènes
explicite("H1 y'=(y^2+xy)/x^2", (y**2 + x * y) / x**2, -x / (log(x) + C), dx_pos)
explicite("H1 x<0", (y**2 + x * y) / x**2, -x / (log(-x) + C), {x: u(-2, -0.3), C: u(3, 5)})
implicite("H2 y'=(x+y)/(x-y)", (x + y) / (x - y), atan(y / x) - log(x**2 + y**2) / 2, {x: u(0.3, 2), y: u(-2, 2)})
explicite("H3 y'=y/x+tan(y/x)", y / x + tan(y / x), x * asin(C * x), {x: u(0.1, 0.9), C: u(-1, 1)})
explicite("H3 autre branche", y / x + tan(y / x), x * (pi - asin(C * x)), {x: u(0.1, 0.9), C: u(-1, 1)})
explicite("H4 y'=y/x+x/y", y / x + x / y, x * sqrt(2 * log(x) + C), dx_pos)
# --- Réductibles
implicite("R1 y'=(x+y-2)/(x-y)", (x + y - 2) / (x - y),
          atan((y - 1) / (x - 1)) - log((x - 1)**2 + (y - 1)**2) / 2, {x: u(1.3, 3), y: u(-2, 2)})
implicite("R2 y'=(x+y+1)/(x+y-1)", (x + y + 1) / (x + y - 1), y - x - log(x + y), {x: u(0.5, 2), y: u(0.5, 2)})
explicite("R2 perdue y=-x", (x + y + 1) / (x + y - 1), -x, dx_pos)
# --- Exactes et facteurs intégrants
forme_diff("E1 exacte", 2 * x * y + 3, x**2 + 4 * y, x**2 * y + 3 * x + 2 * y**2, {x: u(-2, 2), y: u(-2, 2)})
forme_diff("E2 mu=x", 3 * x * y + y**2, x**2 + x * y, x**3 * y + x**2 * y**2 / 2, {x: u(-2, 2), y: u(-2, 2)})
forme_diff("E3 mu=y", y, 2 * x - y * exp(y), x * y**2 - (y**2 - 2 * y + 2) * exp(y), {x: u(-2, 2), y: u(-2, 2)})
forme_diff("E4 mu=1/(xy)^2", y * (1 + x * y), x * (1 - x * y), -1 / (x * y) + log(x) - log(y), {x: u(0.3, 2), y: u(0.3, 2)})
forme_diff("E5 homogène, mu=1/(xM+yN)", x**2 + y**2, -2 * x * y, (x**2 - y**2) / x, {x: u(0.3, 2), y: u(-2, 2)})
# Les facteurs intégrants rendent la forme exacte
def exacte(M, N_):
    return simplify(diff(M, y) - diff(N_, x)) == 0
resultats.append(("E2 x*forme exacte", "", str(exacte(x * (3 * x * y + y**2), x * (x**2 + x * y)))))
resultats.append(("E3 y*forme exacte", "", str(exacte(y * y, y * (2 * x - y * exp(y))))))
m4 = 1 / (x * y)**2
resultats.append(("E4 mu*forme exacte", "", str(exacte(m4 * y * (1 + x * y), m4 * x * (1 - x * y)))))
m5 = 1 / (x * (x**2 + y**2) + y * (-2 * x * y))
resultats.append(("E5 mu*forme exacte", "", str(exacte(m5 * (x**2 + y**2), m5 * (-2 * x * y)))))
m6 = x + y  # contrôle mu(x+y) : M = 2x/(x+y), N = 1/(x+y)
resultats.append(("E6 mu=(x+y) exacte", "", str(exacte(m6 * 2 * x / (x + y), m6 / (x + y)))))
M6, N6 = 2 * x / (x + y), 1 / (x + y)
resultats.append(("E6 critère (My-Nx)/(N-M)", "", str(simplify((diff(M6, y) - diff(N6, x)) / (N6 - M6)))))
M7, N7 = y * exp(-x - y), x * exp(-x - y)
resultats.append(("E7 critère (My-Nx)/(N-M)", "", str(simplify((diff(M7, y) - diff(N7, x)) / (N7 - M7)))))
M4, N4 = y * (1 + x * y), x * (1 - x * y)
resultats.append(("E4 critère (My-Nx)/(yN-xM)", "", str(simplify((diff(M4, y) - diff(N4, x)) / (y * N4 - x * M4)))))
# --- Riccati
explicite("Ri1 y'=y^2-2/x^2", y**2 - 2 / x**2, 1 / x + 3 * x**2 / (K - x**3), {x: u(0.3, 1), K: u(3, 5)})
explicite("Ri1 perdue y=1/x", y**2 - 2 / x**2, 1 / x, dx_pos)
explicite("Ri1 K=0 : y=-2/x", y**2 - 2 / x**2, -2 / x, dx_pos)
explicite("Ri2 y'=y^2-xy+1", y**2 - x * y + 1,
          x + exp(x**2 / 2) / (C - sqrt(pi / 2) * erfi(x / sqrt(2))), {x: u(-0.5, 0.5), C: u(3, 5)})
explicite("Ri3 y'=(y-x)^2+1", (y - x)**2 + 1, x + 1 / (C - x), dx_pos)
# Transformation en EDO linéaire d'ordre 2 : y = -u'/(q2 u)
q0, q1, q2 = Function('q0')(x), Function('q1')(x), Function('q2')(x)
U = Function('u')(x)
Y = -diff(U, x) / (q2 * U)
eq2 = diff(U, x, 2) - (q1 + diff(q2, x) / q2) * diff(U, x) + q2 * q0 * U
res = simplify((diff(Y, x) - (q0 + q1 * Y + q2 * Y**2)) * q2 * U + eq2)
resultats.append(("Riccati -> ordre 2", "y=-u'/(q2 u)", "symbolique" if res == 0 else f"ECHEC {res}"))
# Réduction y = y1 + 1/v
y1, v = Function('y1')(x), Function('v')(x)
ric = lambda z: q0 + q1 * z + q2 * z**2
res = simplify((diff(y1 + 1 / v, x) - ric(y1 + 1 / v)).subs(diff(y1, x), ric(y1))
               - (-(diff(v, x) + (q1 + 2 * q2 * y1) * v + q2) / v**2))
resultats.append(("Riccati y=y1+1/v", "v'=-(q1+2q2y1)v-q2", "symbolique" if res == 0 else f"ECHEC {res}"))
# Bernoulli z = y^(1-n) : z' = (1-n)(-a z + b) pour y' = a y + b y^n  (forme y' + P y = Q y^n)
n, a, b = symbols('n'), Function('a')(x), Function('b')(x)
yy = Function('y')(x)
z = yy**(1 - n)
res = simplify((diff(z, x).subs(diff(yy, x), a * yy + b * yy**n)) - (1 - n) * (a * z + b))
resultats.append(("Bernoulli z=y^(1-n)", "z'=(1-n)(a z+b)", "symbolique" if res == 0 else f"ECHEC {res}"))
# Test de Bernoulli : Q = y f_y - f = (n-1) b y^n, n = y Q_y / Q
fB = a * y + b * y**n
Q = y * diff(fB, y) - fB
resultats.append(("Test Bernoulli", "y Q_y/Q", str(simplify(y * diff(Q, y) / Q))))
# Test de séparabilité
for nom, f in [("sep sin(x+y)+sin(x-y)", sin(x + y) + sin(x - y)), ("sep (x+1)(y+1)", x * y + x + y + 1), ("sep exp(x+y)", exp(x + y)), ("non sep sin(x+y)", sin(x + y)),
               ("non sep x+y", x + y)]:
    R = f * diff(f, x, y) - diff(f, x) * diff(f, y)
    resultats.append(("Test séparable " + nom, "f f_xy - f_x f_y", str(simplify(R))))
# --- Clairaut et Lagrange (F(x, y, p) = 0)
F1 = y - x * p - p**2
implicite_F("Cl1 famille", F1, C * x + C**2, dx_pos)
implicite_F("Cl1 singulière", F1, -x**2 / 4, dx_pos)
F2 = y - x * p - sqrt(1 + p**2)
implicite_F("Cl2 famille", F2, C * x + sqrt(1 + C**2), dx_pos)
implicite_F("Cl2 singulière", F2, sqrt(1 - x**2), {x: u(-0.9, 0.9)})
F3 = y - x * p + log(p)
implicite_F("Cl3 famille", F3, C * x - log(C), dx_pos)
implicite_F("Cl3 singulière", F3, 1 + log(x), dx_pos)
F4 = y - 2 * x * p - p**2
parametrique("La1 y=2xp+p^2", F4, C / p**2 - 2 * p / 3, 2 * C / p - p**2 / 3, {p: u(0.3, 2), C: u(-2, 2)})
implicite_F("La1 singulière y=0", F4, S(0), dx_pos)
F5 = y - x * p**2
implicite_F("La2 y=x p^2", F5, (sqrt(x) - C)**2, {x: u(0.3, 2), C: u(-2, 2)})
implicite_F("La2 singulière y=0", F5, S(0), dx_pos)
implicite_F("La2 y=x (C=0)", F5, x, dx_pos)

random.seed(1)
largeur = max(len(r[0]) for r in resultats)
echecs = 0
for nom, sol, verdict in resultats:
    if verdict.startswith("ECHEC") or "ECHEC" in verdict:
        echecs += 1
    print(f"{nom:<{largeur}}  {verdict}")
print(f"\n{len(resultats)} contrôles, {echecs} échec(s)")
```

</details>

<details>
<summary>Réductions générales de Riccati et facteur intégrant homogène</summary>

```python
# Riccati : deux solutions particulières y1, y2 => w = (y-y1)/(y-y2) vérifie w'/w = q2 (y1 - y2)
from sympy import symbols, Function, simplify, diff
x = symbols('x')
q0, q1, q2 = [Function(n)(x) for n in ('q0', 'q1', 'q2')]
y, y1, y2, y3 = [Function(n)(x) for n in ('y', 'y1', 'y2', 'y3')]
R = lambda z: q0 + q1 * z + q2 * z**2
w = (y - y1) / (y - y2)
dw = diff(w, x).subs({diff(y, x): R(y), diff(y1, x): R(y1), diff(y2, x): R(y2)})
print("w'/w - q2 (y1-y2) =", simplify(dw / w - q2 * (y1 - y2)))
# Birapport de quatre solutions constant
br = (y - y1) * (y3 - y2) / ((y - y2) * (y3 - y1))
dbr = diff(br, x).subs({diff(y, x): R(y), diff(y1, x): R(y1), diff(y2, x): R(y2), diff(y3, x): R(y3)})
print("d(birapport)/dx =", simplify(dbr))
# Facteur intégrant 1/(xM+yN) pour M, N homogènes de même degré (cas général degré 2 symbolique)
X, Y = symbols('X Y')
a, b, c, d, e, f = symbols('a b c d e f')
M = a*X**2 + b*X*Y + c*Y**2
N = d*X**2 + e*X*Y + f*Y**2
mu = 1 / (X*M + Y*N)
print("homogène deg 2 : (mu M)_y - (mu N)_x =", simplify(diff(mu*M, Y) - diff(mu*N, X)))
# Euler : x f_x + y f_y = 0 pour f homogène de degré 0 ; contre-exemple t<0
f0 = (X**2 + Y**2) / (X*Y)
print("Euler deg 0 :", simplify(X*diff(f0, X) + Y*diff(f0, Y)))
```

</details>
