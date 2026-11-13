# Équations différentielles du second ordre réductibles au premier ordre

Certaines équations du second ordre non linéaires se ramènent, par un changement d'inconnue, à une équation du premier ordre suivie d'une quadrature. Ce document décrit ces réductions (variable $y$ absente, variable $x$ absente, $y'' = f(y)$, homogénéité en $y$, invariance d'échelle en $x$, dérivée exacte), les solutions qu'elles perdent, et la détermination des constantes dans un problème de Cauchy. Les méthodes du premier ordre auxquelles on se ramène sont celles de [`edo-premier-ordre.md`](edo-premier-ordre.md), dont on reprend les outils (test de nullité, gel de variable, réinjection : section 4 de ce document).

## Sommaire

1. [Notations](#1-notations)
2. [Cadre, existence et unicité](#2-cadre-existence-et-unicité)
3. [Variable y absente : F(x, y', y'') = 0](#3-variable-y-absente--fx-y-y--0)
4. [Variable x absente (autonome) : F(y, y', y'') = 0](#4-variable-x-absente-autonome--fy-y-y--0)
5. [Cas y'' = f(y) : intégrale première de l'énergie](#5-cas-y--fy--intégrale-première-de-lénergie)
6. [Équation homogène en y](#6-équation-homogène-en-y)
7. [Équation équidimensionnelle en x](#7-équation-équidimensionnelle-en-x)
8. [Dérivée exacte](#8-dérivée-exacte)
9. [Problème de Cauchy : déterminer les constantes](#9-problème-de-cauchy--déterminer-les-constantes)
10. [Stratégie d'un solveur](#10-stratégie-dun-solveur)
11. [Pièges et cas de contrôle](#11-pièges-et-cas-de-contrôle)
12. [État dans le code](#12-état-dans-le-code)
13. [Annexe : script de vérification](#13-annexe--script-de-vérification)

## 1. Notations

| Notation | Sens |
| :--- | :--- |
| $x$, $y = y(x)$ | variable indépendante, inconnue |
| $p = y'$, $q = y''$ | dérivées, traitées comme des variables indépendantes dans $F$ |
| $F(x, y, p, q) = 0$ | forme implicite ; forme résolue $y'' = f(x, y, y')$ |
| $F_y, F_p, F_q$ | dérivées partielles de $F$ |
| $C_1, C_2$, $K$, $a$ | constantes arbitraires |
| $(x_0, y_0, y_1)$ | données de Cauchy : $y(x_0) = y_0$, $y'(x_0) = y_1$ |
| $e \equiv 0$ | nullité identique, décidée par le test `Nul` (premier ordre, section 4.1) |

Dans le code : `x`, `y`, `p = var("p")`, `q = var("q")` (quatre variables distinctes) et `param("C1")`, `param("C2")`.

## 2. Cadre, existence et unicité

**Théorème (Cauchy–Lipschitz pour le second ordre).** Si $f(x, y, p)$ est continue sur un ouvert $U \subset \mathbb{R}^3$ et localement lipschitzienne par rapport à $(y, p)$, alors pour tout $(x_0, y_0, y_1) \in U$ il existe une unique solution maximale de $y'' = f(x, y, y')$, $y(x_0) = y_0$, $y'(x_0) = y_1$. *Preuve* : on applique le théorème du premier ordre au système $Y = (y, y')$, $Y' = (y', f(x, y, y'))$. ∎

Conséquences :
* la solution générale dépend de **deux** constantes ; une famille à une seule constante obtenue par une réduction n'est jamais la solution générale : c'est une famille de solutions perdues (ou singulières) qu'il faut lister à part ;
* pour une forme implicite $F(x, y, y', y'') = 0$, le théorème s'applique là où $F_q \neq 0$ (fonctions implicites). Si $F_q = 0$ aux données de Cauchy, il peut n'y avoir **aucune** solution (cas D15) ou plusieurs.

Pour un solveur, chaque réduction ci-dessous transforme $F$ en une équation du premier ordre, résolue par `ResoudreEDO1` (premier ordre, section 13), puis effectue une quadrature. Les solutions perdues au premier ordre (par exemple $p = p_0(x)$) deviennent des familles **à une constante** ($y = \int p_0 + C_2$), qu'il faut réinjecter dans l'équation de départ.

## 3. Variable $y$ absente : $F(x, y', y'') = 0$

**Théorème.** Si $F$ ne dépend pas de $y$, $y$ est solution si et seulement si $p = y'$ est solution de l'équation du premier ordre $F(x, p, p') = 0$. Si $p = P(x, C_1)$ est sa solution générale, $y = \int P(x, C_1)\,dx + C_2$.

**Test** : `contient(F, y)` est faux (structurel, coût linéaire).

**Solutions perdues** : chaque solution perdue $p_0(x)$ de l'équation en $p$ donne la famille $y = \int p_0\,dx + C_2$.

**Piège : discussion selon le signe d'une constante.** La quadrature $\int P(x, C_1)\,dx$ peut changer de forme selon la valeur de $C_1$ (cas D4 : $\int \frac{dx}{x^2 + C_1}$ est un arc tangente, un logarithme ou $-1/x$ selon le signe de $C_1$). Un intégrateur symbolique qui traite $C_1$ comme un symbole sans signe ne peut pas trancher ; il faut soit renvoyer les trois cas avec leurs conditions, soit poser $C_1 = \pm a^2$ explicitement.

### Exemples (vérifiés)

**D1.** $y'' = y'/x$ : $p' = p/x$, $p = 2C_1 x$, $y = C_1 x^2 + C_2$.

**D2.** $y'' = 1 + y'^2$ : $p' = 1 + p^2$, $p = \tan(x + C_1)$, $y = -\ln\lvert\cos(x + C_1)\rvert + C_2$.

**D3 (familles perdues).** $x y'' = y' + y'^2$ : $x p' = p(1 + p)$ (séparable). Pour $p \notin \{0, -1\}$ : $\frac{p}{1 + p} = Kx$, $p = \frac{Kx}{1 - Kx}$, et
$$y = -x - \frac{\ln\lvert 1 - Kx \rvert}{K} + C_2 \qquad (K \neq 0).$$
Perdues au premier ordre : $p = 0$ et $p = -1$, d'où les familles $y = C_2$ et $y = -x + C_2$ (vérifiées). La première est la limite $K \to 0$ de la formule (calculée : $-x - \ln(1 - Kx)/K \to 0$) mais la formule elle-même n'est pas définie en $K = 0$ ; la seconde correspond à $K = \infty$.

**D4 (discussion).** $y'' = 2x\,y'^2$ : $p' = 2x p^2$, $-1/p = x^2 + C_1$, $p = -\frac{1}{x^2 + C_1}$.
* $C_1 = a^2 > 0$ : $y = -\frac{1}{a}\arctan\frac{x}{a} + C_2$ ;
* $C_1 = -a^2 < 0$ : $y = \frac{1}{2a}\ln\left\lvert\frac{x + a}{x - a}\right\rvert + C_2$ ;
* $C_1 = 0$ : $y = \frac{1}{x} + C_2$ ;
* perdue ($p = 0$) : $y = C_2$.

## 4. Variable $x$ absente (autonome) : $F(y, y', y'') = 0$

**Théorème.** Soit $y$ une solution sur un intervalle où $y' \neq 0$. Alors $y$ est strictement monotone, donc inversible, et $P(y) := y'(x(y))$ vérifie
$$y'' = \frac{d}{dx}P(y(x)) = P'(y)\,y' = P\,\frac{dP}{dy},$$
donc $F(y, P, P\,P') = 0$, équation du **premier ordre** en $P(y)$. Si $P = \Pi(y, C_1)$, alors $y$ s'obtient par la séparation $\int \frac{dy}{\Pi(y, C_1)} = x + C_2$.

**Test** : `contient(F, x)` est faux.

**Solutions perdues.**
* Les **constantes** $y = c$ (où $y' \equiv 0$, exclues par l'hypothèse $y' \neq 0$) : elles sont solutions si et seulement si $F(c, 0, 0) = 0$. Il faut résoudre cette équation en $c$ (`resoudre`, une variable) et ajouter les solutions qui ne sont pas dans la famille.
* Les **points de rebroussement** : une solution non constante peut avoir $y'(x_1) = 0$ en des points isolés (par exemple un minimum de $y$). La réduction la découpe en deux branches $P = \pm\ldots$ qu'il faut **recoller** ; une formule close de la solution générale (quand on peut l'écrire) évite ce problème, voir D14.

### Exemples (vérifiés)

**D5.** $y\,y'' = y'^2$ : $y P P' = P^2$, donc (pour $P \neq 0$) $P = C_1 y$ et $y = C_2 e^{C_1 x}$. Les constantes ($F(c, 0, 0) = 0$ pour tout $c$) sont dans la famille ($C_1 = 0$).

**D6.** $2y\,y'' = 1 + y'^2$ : $2y P P' = 1 + P^2$, $\frac{2P\,dP}{1 + P^2} = \frac{dy}{y}$, $1 + P^2 = C_1 y$, puis
$$y = \frac{1}{C_1} + \frac{C_1 (x + C_2)^2}{4}$$
(paraboles). Réinjection : $y' = \frac{C_1(x + C_2)}{2}$, $y'' = \frac{C_1}{2}$, $2y\,y'' = C_1 y = 1 + \frac{C_1^2 (x + C_2)^2}{4} = 1 + y'^2$. ✓ Constantes : $F(c, 0, 0) = -1 \neq 0$, aucune.

**D7.** $y\,y'' + y'^2 = 0$ : c'est $(y\,y')' = 0$ (section 8), $y^2 = C_1 x + C_2$, soit $y = \pm\sqrt{C_1 x + C_2}$. Les constantes sont dans la famille ($C_1 = 0$).

## 5. Cas $y'' = f(y)$ : intégrale première de l'énergie

**Théorème.** Si $f$ est continue et $\mathcal{F}' = f$, toute solution de $y'' = f(y)$ vérifie
$$\frac{1}{2}\,y'^2 - \mathcal{F}(y) = E \quad (\text{constante}).$$
*Preuve* : la dérivée vaut $y'\,y'' - f(y)\,y' = 0$ (vérifié symboliquement). ∎

Sur un intervalle où $y' \neq 0$ : $y' = \sigma\sqrt{2(\mathcal{F}(y) + E)}$ avec $\sigma = \pm 1$ constant, et
$$x = \sigma\int \frac{dy}{\sqrt{2(\mathcal{F}(y) + E)}} + C_2.$$

**Solutions perdues** : les équilibres $y = y_0$ avec $f(y_0) = 0$ (cas $y' \equiv 0$). Si $f$ est lipschitzienne, ce sont les seules solutions atteignant $y' = 0$ avec $f(y) = 0$.

**Limites.** L'intégrale $\int \frac{dy}{\sqrt{2(\mathcal{F} + E)}}$ est en général non élémentaire : pour $\mathcal{F}$ polynôme de degré 3 ou 4, c'est une intégrale elliptique, et la solution générale s'exprime par des fonctions elliptiques. Le solveur doit alors renvoyer la **quadrature implicite** (intégrale non évaluée) ou les seules sous-familles élémentaires (valeurs particulières de $E$), en le disant.

### Exemples (vérifiés)

**D8.** $y'' = 2y^3$ : $\frac{1}{2}y'^2 = \frac{y^4}{2} + E$. Pour $E = 0$ : $y' = \pm y^2$, d'où $y = \frac{1}{C_2 - x}$ et $y = -\frac{1}{C_2 - x}$ (vérifiées). Pour $E \neq 0$ : intégrale elliptique, **non évaluée**. Équilibre : $y = 0$.

**D9 (entièrement élémentaire).** $y'' = -y^{-3}$ : $\mathcal{F} = \frac{1}{2y^2}$, $y'^2 = \frac{1}{y^2} + a$ ($a = 2E$), soit $y\,y' = \pm\sqrt{1 + a y^2}$ ; pour $a \neq 0$, $\frac{d}{dx}\sqrt{1 + a y^2} = \pm a$, d'où
$$y^2 = a(x + C_2)^2 - \frac{1}{a} \quad (a \neq 0), \qquad y^2 = 2(x + C_2) \;\text{ ou }\; y^2 = -2(x + C_2) \quad (a = 0).$$
Vérifié pour $a > 0$, $a < 0$ et $a = 0$. Pas d'équilibre ($-y^{-3}$ ne s'annule pas).

## 6. Équation homogène en $y$

**Définition.** $F$ est homogène de degré $m$ en $(y, p, q)$ si $F(x, ky, kp, kq) = k^m F(x, y, p, q)$ pour tout $k > 0$. Test (identité d'Euler) :
$$\frac{y F_y + p F_p + q F_q}{F} \text{ est une constante } m.$$

**Réduction.** Avec $y' = u\,y$ (soit $u = (\ln\lvert y\rvert)'$), on a $y'' = y\,(u' + u^2)$ et $F(x, y, uy, y(u' + u^2)) = y^m F(x, 1, u, u' + u^2)$ pour $y > 0$ : l'équation devient du premier ordre en $u$. Puis $y = C_2 \exp\int u\,dx$.

**Solution perdue** : $y \equiv 0$ (division par $y^m$), à réinjecter ; elle correspond souvent à $C_2 = 0$.

**D12 (vérifié).** $y\,y'' - y'^2 = x\,y^2$ (degré $m = 2$, calculé) : $u' + u^2 - u^2 = x$, $u = \frac{x^2}{2} + C_1$,
$$y = C_2 \exp\Big(\frac{x^3}{6} + C_1 x\Big).$$
$y \equiv 0$ est solution et correspond à $C_2 = 0$. D5 est aussi homogène (degré 2) : $u' = 0$.

## 7. Équation équidimensionnelle en $x$

**Définition.** $F$ est équidimensionnelle en $x$ si elle est invariante, à un facteur près, par le changement d'échelle $x \mapsto kx$ ($y$ fixe, donc $p \mapsto p/k$, $q \mapsto q/k^2$). Test infinitésimal :
$$\frac{x F_x - p F_p - 2q F_q}{F} \text{ est une constante}.$$
C'est le cas de toute équation de la forme $G(y, x y', x^2 y'') = 0$.

**Réduction.** Avec $x = e^t$ ($x > 0$ ; $x = -e^t$ pour $x < 0$) et $Y(t) = y(e^t)$ :
$$x\,y' = \dot Y, \qquad x^2 y'' = \ddot Y - \dot Y$$
(vérifié), et l'équation en $Y(t)$ est **autonome** (section 4).

**D11 (vérifié).** $x^2 y'' = (x y')^2$ (poids 0, calculé) : $\ddot Y - \dot Y = \dot Y^2$, autonome et sans $Y$ ; en revenant à $x$, ou directement par la section 3 ($p' = p^2$) :
$$y = -\ln\lvert x + C_1 \rvert + C_2.$$
Perdue : $y = C_2$ ($p = 0$). Avec le paramétrage équivalent $y = -\ln\lvert 1 - Kx \rvert + C_2$, la même famille $y = C_2$ correspond à $K = 0$ et n'est plus perdue : **le choix du paramétrage décide de ce qui est « perdu »** ; un test doit donc vérifier l'appartenance à la famille par réinjection, pas par comparaison de formules.

D1 ($y'' = y'/x$) est aussi équidimensionnelle (poids $-2$, calculé) ; D4 ne l'est pas.

## 8. Dérivée exacte

Si $F(x, y, y', y'')$ est la dérivée totale d'une expression $G(x, y, y')$, l'équation $F = 0$ équivaut à $G = C_1$, du premier ordre.

**Critère (opérateur d'Euler–Lagrange).** Pour $F$ régulière, $F$ est une dérivée totale si et seulement si
$$\mathcal{E}(F) = F_y - \frac{d}{dx}F_p + \frac{d^2}{dx^2}F_q \equiv 0$$
(énoncé classique du calcul des variations sur les « lagrangiens nuls », sur un domaine convenable ; je ne détaille pas ici la reconstruction de $G$, qui se fait par intégrations successives comme pour le potentiel du premier ordre). Calculé : $\mathcal{E}(y y'' + y'^2) = 0$ (c'est $(y y')'$) ; $\mathcal{E}(y y'' - y'^2) = 4y''$ (non exacte) mais $\mathcal{E}\big((y y'' - y'^2)/y^2\big) = 0$ (c'est $(y'/y)'$, le facteur $1/y^2$ jouant le rôle d'un facteur intégrant) ; $\mathcal{E}(y'' - 1 - y'^2) = 2y''$.

Ce test coûte trois dérivations totales ; il est utile surtout comme **heuristique** (essayer aussi $F/y^2$, $F/y'^2$…). Il n'est pas indispensable pour un premier solveur.

## 9. Problème de Cauchy : déterminer les constantes

### 9.1 Méthode recommandée : imposer les conditions à chaque étage

Plutôt que de résoudre le système non linéaire $\{\varphi(x_0, C_1, C_2) = y_0,\ \varphi'(x_0, C_1, C_2) = y_1\}$ :
1. **à l'étage du premier ordre**, déterminer $C_1$ à partir de $y_1$ : $P(x_0, C_1) = y_1$ (section 3), ou $\Pi(y_0, C_1) = y_1$ (section 4), ou $E = \frac{1}{2}y_1^2 - \mathcal{F}(y_0)$ (section 5, sans résolution du tout) ;
2. **à l'étage de la quadrature**, déterminer $C_2$ à partir de $y_0$, équation en général **affine** en $C_2$.

Chaque étape est une équation à **une** inconnue. Dans SymAlgo++, `resoudre` exige une variable (`Variable`), pas un paramètre : il faut substituer `var` à `param("C1")` le temps de la résolution.

### 9.2 Choix de branche

* Énergie : $\sigma = \mathrm{signe}(y_1)$. Si $y_1 = 0$ et $f(y_0) \neq 0$, $x_0$ est un point de rebroussement de $y'$ : les deux branches se recollent en $x_0$ ($y''(x_0) = f(y_0)$ fixe le sens). Si $y_1 = 0$ et $f(y_0) = 0$ : équilibre $y \equiv y_0$ (unique si $f$ est lipschitzienne).
* Constantes périodiques ($p = \tan(x + C_1)$) : $C_1 = \arctan(y_1) - x_0$ ; les autres valeurs ($+k\pi$) donnent la même solution. Vérifier que $x_0$ est dans l'intervalle maximal de la formule choisie.

### 9.3 Cas d'échec à reconnaître

1. **Donnée sur une solution perdue** : la condition $P(x_0, C_1) = y_1$ n'a pas de solution, ou seulement une valeur exclue ($K = 0$ hors domaine de la formule, D17). Il faut alors tester les solutions perdues listées.
2. **Données incompatibles** : $F_q = 0$ aux données et $F \neq 0$ (D15) : aucune solution.
3. **Équation transcendante en $C_1$** non résoluble exactement : solution approchée (`resoudreNumerique`), marquée comme telle.
4. **Sous-famille non élémentaire** : la valeur de $C_1$ imposée tombe dans le cas non évalué (D18c, elliptique).
5. Dans tous les cas, la solution renvoyée est **réinjectée** et les conditions initiales **vérifiées** ; sinon, résultat non évalué (ou intégration numérique, marquée approchée).

### 9.4 Exemples (vérifiés)

**D14.** $2y\,y'' = 1 + y'^2$, $y(0) = 1$, $y'(0) = 0$. La branche $P = \pm\sqrt{C_1 y - 1}$ est singulière en $x_0$ ($P = 0$) : mauvaise voie. Avec la formule close D6 : $y'(0) = \frac{C_1 C_2}{2} = 0$ et $C_1 \neq 0$ donnent $C_2 = 0$, puis $y(0) = 1/C_1 = 1$ :
$$y = 1 + \frac{x^2}{4}. \;✓$$

**D15 (aucune solution).** $y\,y'' = y'^2$, $y(0) = 0$, $y'(0) = 1$ : en $x = 0$, $y\,y'' - y'^2 = 0 \cdot y''(0) - 1 = -1 \neq 0$ quelle que soit $y''(0)$. Aucune solution deux fois dérivable. La famille $C_2 e^{C_1 x}$ le confirme ($y(0) = 0 \Rightarrow C_2 = 0 \Rightarrow y'(0) = 0$).

**D16.** $y'' = 1 + y'^2$, $y(0) = 0$, $y'(0) = 0$ : $\tan C_1 = 0$, $C_1 = 0$ ; $C_2 = 0$ ; $y = -\ln(\cos x)$ sur $]-\pi/2, \pi/2[$. ✓

**D17 (donnée sur une solution perdue).** $x y'' = y' + y'^2$, $y(1) = 0$, $y'(1) = 0$ : $p(1) = \frac{K}{1 - K} = 0$ impose $K = 0$, exclu de la formule D3. La famille perdue $y = C_2$ convient : $y \equiv 0$, unique car $x_0 = 1 \neq 0$ rend l'équation résolue lipschitzienne.

**D18 (énergie).** $y'' = 2y^3$, $y(0) = 1$ :
* (a) $y'(0) = 1$ : $E = \frac{1}{2} - \frac{1}{2} = 0$, $\sigma = +1$, $y' = y^2$, $y = \frac{1}{1 - x}$ ($x < 1$). ✓
* (b) $y'(0) = -1$ : $E = 0$, $\sigma = -1$, $y = \frac{1}{1 + x}$ ($x > -1$). ✓
* (c) $y'(0) = 0$ : $E = -\frac{1}{2} \neq 0$ : cas elliptique, **non évalué** (ou numérique).

**D19.** $y\,y'' = y'^2$, $y(0) = 1$, $y'(0) = 2$ : $C_2 = 1$, $C_1 = 2$, $y = e^{2x}$. ✓

## 10. Stratégie d'un solveur

```text
ResoudreEDO2(F) :                              # F(x, y, p, q) = 0
  0. Linéaire à coefficients constants, homogène  -> EquationDifferentielle existante
  1. non contient(F, y)      -> section 3 : ResoudreEDO1 en p, puis quadrature
  2. non contient(F, x)      -> si F = q - f(y) : énergie (section 5)
                                sinon section 4 : ResoudreEDO1 en P(y), puis séparation
                                + constantes c avec F(c, 0, 0) = 0
  3. homogène en (y, p, q)   -> section 6 ; + y = 0
  4. équidimensionnelle      -> x = e^t, puis étape 2
  5. dérivée exacte          -> G = C1 (section 8), puis ResoudreEDO1
  6. ECHEC                   -> résultat « non évalué »
  Après chaque succès : réinjecter (deux constantes symboliques, puis valeurs aléatoires),
  convertir chaque solution perdue du premier ordre en famille à une constante et la
  réinjecter ; joindre les conditions de domaine et la discussion de signe éventuelle.
```

Les tests 1 et 2 sont structurels (coût linéaire) ; 3 et 4 demandent trois dérivées partielles et un test « constante » ; 5 trois dérivées totales. Les réductions 3 et 4 se combinent (D1, D11) : il suffit de prendre la première qui réussit.

## 11. Pièges et cas de contrôle

Toutes les solutions ont été vérifiées par réinjection avec SymPy 1.14 (script en annexe : 47 contrôles, 0 échec, tous symboliques), sauf mention contraire. $F$ est écrite $= 0$.

**y absente**
1. $y'' = y'/x$ — $y = C_1 x^2 + C_2$. Aussi équidimensionnelle (poids $-2$).
2. $y'' = 1 + y'^2$ — $y = -\ln\lvert\cos(x + C_1)\rvert + C_2$.
3. $x y'' = y' + y'^2$ — $y = -x - \frac{1}{K}\ln\lvert 1 - Kx\rvert + C_2$ ; perdues $y = C_2$ (limite $K \to 0$) et $y = -x + C_2$.
4. $y'' = 2x\,y'^2$ — **piège** du signe de $C_1$ : trois formes (arc tangente, logarithme, $1/x$) ; perdue $y = C_2$.

**x absente**
5. $y\,y'' = y'^2$ — $y = C_2 e^{C_1 x}$ ; constantes incluses ($C_1 = 0$). Aussi homogène en $y$.
6. $2y\,y'' = 1 + y'^2$ — $y = \frac{1}{C_1} + \frac{C_1(x + C_2)^2}{4}$ ; **aucune** constante solution ($F(c, 0, 0) = -1$).
7. $y\,y'' + y'^2 = 0$ — $y^2 = C_1 x + C_2$ ; dérivée exacte de $y\,y'$.

**$y'' = f(y)$**
8. $y'' = 2y^3$ — sous-famille $E = 0$ : $y = \pm\frac{1}{C_2 - x}$ ; $E \neq 0$ : **non évalué** (elliptique) ; équilibre $y = 0$.
9. $y'' = -1/y^3$ — $y^2 = a(x + C_2)^2 - \frac{1}{a}$ ($a \neq 0$, vérifié pour $a > 0$ et $a < 0$), $y^2 = \pm 2(x + C_2)$ ($a = 0$, vérifié pour le signe $+$).
10. $y'' = -y$ — linéaire à coefficients constants : doit être routée vers `EquationDifferentielle` (étape 0), $y = C_1\cos x + C_2\sin x$ (immédiat ; le cas voisin $y'' + 4y = 0$ est testé dans `tests/test_differentielle.cpp:34`). Piège : la reconnaissance « $x$ absente » (étape 2) s'applique aussi ; elle doit passer **après** l'étape 0, sinon on obtient une quadrature $\int dy/\sqrt{a - y^2}$ au lieu de la forme linéaire.

**Homogène en y, équidimensionnelle**
11. $x^2 y'' = (x y')^2$ — $y = -\ln\lvert x + C_1\rvert + C_2$ (ou $-\ln\lvert 1 - Kx\rvert + C_2$) ; perdue $y = C_2$ selon le paramétrage.
12. $y\,y'' - y'^2 = x\,y^2$ — $y = C_2\exp\big(\frac{x^3}{6} + C_1 x\big)$ ; $y = 0$ pour $C_2 = 0$.

**Tests de reconnaissance (valeurs calculées)**
* Degré d'homogénéité en $(y, p, q)$ : 2 pour les cas 5 et 12 ; non constant pour le cas 6.
* Poids d'échelle $(x F_x - p F_p - 2q F_q)/F$ : 0 pour le cas 11, $-2$ pour le cas 1, non constant pour le cas 4.
* $\mathcal{E}(F)$ : 0 pour le cas 7 et pour $(y y'' - y'^2)/y^2$ ; $4y''$ pour le cas 5 tel quel ; $2y''$ pour le cas 2.

**Problèmes de Cauchy**
13. $2y\,y'' = 1 + y'^2$, $y(0) = 1$, $y'(0) = 0$ — $y = 1 + x^2/4$ (piège : branche $\pm\sqrt{\cdot}$ singulière en $x_0$).
14. $y\,y'' = y'^2$, $y(0) = 0$, $y'(0) = 1$ — **aucune solution** (vérifié à la main : résidu $-1$ en $x_0$).
15. $y'' = 1 + y'^2$, $y(0) = 0$, $y'(0) = 0$ — $y = -\ln(\cos x)$, $\lvert x\rvert < \pi/2$.
16. $x y'' = y' + y'^2$, $y(1) = 0$, $y'(1) = 0$ — $y \equiv 0$ (solution perdue ; vérification immédiate).
17. $y'' = 2y^3$, $y(0) = 1$, $y'(0) = 1$ — $y = \frac{1}{1 - x}$.
18. $y'' = 2y^3$, $y(0) = 1$, $y'(0) = -1$ — $y = \frac{1}{1 + x}$.
19. $y'' = 2y^3$, $y(0) = 1$, $y'(0) = 0$ — **non évalué** (elliptique) ou numérique marqué approché.
20. $y\,y'' = y'^2$, $y(0) = 1$, $y'(0) = 2$ — $y = e^{2x}$.

## 12. État dans le code

Références vérifiées dans le worktree `theorie` (commit `569e686`, identique à `origin/main`).

| Élément | État | Où |
| :--- | :--- | :--- |
| EDO linéaire homogène à coefficients constants (tout ordre) : solution littérale, Cauchy en $x = 0$, RK4 | implémenté | `include/EquationDifferentielle.hpp:26` ; `src/EquationDifferentielle.cpp:235` (`resoudreLitteral`), `:244` (`resoudreProblemeCauchy`), `:276` (`eval`, RK4) |
| Conditions initiales hors de $x_0 = 0$ | **absent** (Cauchy et RK4 partent de $x = 0$) | `src/EquationDifferentielle.cpp:254` (dérivées évaluées en 0), `:284` |
| Représentation d'une EDO non linéaire du second ordre $F(x, y, p, q)$ | **absent** | — |
| Réductions des sections 3 à 8 | **absent** | — |
| Détermination non linéaire des constantes | **absent** ; le cas linéaire passe par un système linéaire numérique (QR d'Eigen) | `src/EquationDifferentielle.cpp:262` |
| Résolution d'une équation à une inconnue (pour $C_1$, $C_2$, ou $F(c, 0, 0) = 0$) | implémenté pour une **variable** (pas un paramètre) et des solutions numériques | `src/Solveur.cpp:539` (`resoudre`), `:530` (une seule variable), `:112` |
| Recherche numérique de racines (constante transcendante) | implémenté | `src/Solveur.cpp:561` (`resoudreNumerique`) |
| Intégration numérique d'un système non linéaire $Y' = G(x, Y)$ | **absent** (le RK4 existant est spécialisé au cas linéaire : un pas = un produit matrice-vecteur) | `src/EquationDifferentielle.cpp:289-303` |
| Outils communs (dérivées partielles, `substituer`, `contient`, `param`, intégration) | voir premier ordre, section 15 | `edo-premier-ordre.md` |

### Branche en cours (feature/primitives-substitution, commit 96de779)

Lu en lecture seule le 13/11/2026 : aucune EDO non linéaire (ni du premier ni du second ordre) n'y est encore écrite ; `src/EquationDifferentielle.cpp` y est inchangé quant aux fonctionnalités décrites ci-dessus. Les progrès de l'intégration symbolique (`src/Primitives.cpp`, non commité) profiteront directement aux quadratures des sections 3 à 5.

## 13. Annexe : script de vérification

Lancé avec Python 3 et SymPy 1.14 ; sortie : « 47 contrôles, 0 échec(s) ». `verifier(nom, F, sol, dom)` substitue $y$, $y'$, $y''$ dans $F$ et teste la nullité (symbolique, sinon numérique en 6 points) ; `cauchy` vérifie les conditions initiales.

<details>
<summary>verifier_second_ordre.py</summary>

```python
# Vérification par réinjection des exemples de docs/theorie/edo-second-ordre-reductible.md
# Lancement : python3 -I verifier_second_ordre.py (SymPy 1.14)
import random
from sympy import (symbols, sqrt, exp, log, cos, atan, simplify, diff, Rational, S, N, limit, series,
                   Function, Abs)

random.seed(1)
x, C1, C2, K, a = symbols('x C1 C2 K a')
Y0, Y1, Y2 = symbols('Y0 Y1 Y2')  # y, y', y''
resultats = []


def u(lo, hi):
    return lambda: Rational(random.randint(int(lo * 100), int(hi * 100)), 100)


def verifier(nom, F, sol, dom):
    """F(x, y, y', y'') = 0 avec y = sol(x)."""
    r = F.subs({Y2: diff(sol, x, 2), Y1: diff(sol, x), Y0: sol})
    e = simplify(r)
    if e == 0:
        verdict = "symbolique"
    else:
        verdict = "numérique (6 points)"
        for _ in range(6):
            vals = {s: f() for s, f in dom.items()}
            v = complex(N(r.subs(vals), 30))
            if abs(v) > 1e-12:
                verdict = f"ECHEC ({vals} -> {v})"
                break
    resultats.append((nom, verdict))


def cauchy(nom, sol, x0, y0, y1):
    ok = simplify(sol.subs(x, x0) - y0) == 0 and simplify(diff(sol, x).subs(x, x0) - y1) == 0
    resultats.append((nom + " (conditions)", "ok" if ok else "ECHEC"))


pos = {x: u(0.2, 0.9), C1: u(0.1, 0.9), C2: u(-2, 2), K: u(0.1, 0.9)}

# y absent
F = Y2 - Y1 / x
verifier("D1 y''=y'/x", F, C1 * x**2 + C2, pos)
F = Y2 - 1 - Y1**2
verifier("D2 y''=1+y'^2", F, -log(cos(x + C1)) + C2, {x: u(-0.5, 0.5), C1: u(-0.5, 0.5), C2: u(-2, 2)})
F = x * Y2 - Y1 - Y1**2
verifier("D3 xy''=y'+y'^2", F, -x - log(1 - K * x) / K + C2, pos)
verifier("D3 perdue y=C", F, C2, pos)
verifier("D3 perdue y=-x+C", F, -x + C2, pos)
resultats.append(("D3 limite K->0", str(limit(-x - log(1 - K * x) / K, K, 0))))
F = Y2 - 2 * x * Y1**2
verifier("D4 C1=a^2>0", F, -atan(x / a) / a + C2, {x: u(-2, 2), a: u(0.5, 2), C2: u(-2, 2)})
verifier("D4 C1=-a^2<0", F, log((x + a) / (x - a)) / (2 * a) + C2, {x: u(2.5, 4), a: u(0.5, 2), C2: u(-2, 2)})
verifier("D4 C1=0", F, 1 / x + C2, pos)
verifier("D4 perdue y=C", F, C2, pos)
# x absent
F = Y0 * Y2 - Y1**2
verifier("D5 yy''=y'^2", F, C2 * exp(C1 * x), pos)
F = 2 * Y0 * Y2 - 1 - Y1**2
verifier("D6 2yy''=1+y'^2", F, 1 / C1 + C1 * (x + C2)**2 / 4, pos)
F = Y0 * Y2 + Y1**2
verifier("D7 yy''+y'^2=0", F, sqrt(C1 * x + C2), {x: u(0.2, 0.9), C1: u(0.1, 0.9), C2: u(1, 2)})
F = Y2 - 2 * Y0**3
verifier("D8 y''=2y^3 (+)", F, 1 / (C2 - x), {x: u(-2, 0), C2: u(1, 2)})
verifier("D8 y''=2y^3 (-)", F, -1 / (C2 - x), {x: u(-2, 0), C2: u(1, 2)})
F = Y2 + 1 / Y0**3
verifier("D9 y''=-1/y^3", F, sqrt(a * (x + C2)**2 - 1 / a), {x: u(1, 2), a: u(1, 2), C2: u(0, 1)})
verifier("D9 a<0", F, sqrt(a * (x + C2)**2 - 1 / a), {x: u(-0.3, 0.3), a: u(-1, -0.5), C2: u(-0.1, 0.1)})
verifier("D9 a=0", F, sqrt(2 * (x + C2)), {x: u(0.2, 2), C2: u(0, 1)})
# équidimensionnelle et homogène en y
F = x**2 * Y2 - (x * Y1)**2
verifier("D11 x^2y''=(xy')^2", F, -log(1 - K * x) + C2, pos)
verifier("D11 autre paramétrage", F, -log(x + C1) + C2, pos)
verifier("D11 perdue y=C", F, C2, pos)
F = Y0 * Y2 - Y1**2 - x * Y0**2
verifier("D12 yy''-y'^2=xy^2", F, C2 * exp(x**3 / 6 + C1 * x), pos)
verifier("D12 perdue y=0", F, S(0), pos)
# Euler : x = e^t
t = symbols('t')
g = Function('g')
yx = g(log(x))
lhs1 = simplify(x * diff(yx, x))
lhs2 = simplify(x**2 * diff(yx, x, 2))
resultats.append(("x y' = dy/dt", str(lhs1)))
resultats.append(("x^2 y'' = y_tt - y_t", str(lhs2)))
# Problèmes de Cauchy
F = 2 * Y0 * Y2 - 1 - Y1**2
verifier("D14 Cauchy 2yy''=1+y'^2", F, 1 + x**2 / 4, pos)
cauchy("D14", 1 + x**2 / 4, 0, 1, 0)
F = Y2 - 1 - Y1**2
verifier("D16 Cauchy y''=1+y'^2", F, -log(cos(x)), {x: u(-1.5, 1.5)})
cauchy("D16", -log(cos(x)), 0, 0, 0)
F = Y2 - 2 * Y0**3
verifier("D18a y''=2y^3", F, 1 / (1 - x), {x: u(-2, 0.9)})
cauchy("D18a", 1 / (1 - x), 0, 1, 1)
verifier("D18b y''=2y^3", F, 1 / (1 + x), {x: u(-0.9, 2)})
cauchy("D18b", 1 / (1 + x), 0, 1, -1)
F = Y0 * Y2 - Y1**2
verifier("D19 Cauchy yy''=y'^2", F, exp(2 * x), pos)
cauchy("D19", exp(2 * x), 0, 1, 2)
# Énergie : y'' = f(y) => (y'^2/2 - F(y))' = 0
yf = Function('y')(x)
f = Function('f')
Fp = Function('Fprim')
E = diff(yf, x)**2 / 2 - Fp(yf)
dE = diff(E, x).subs(diff(yf, x, 2), f(yf))
resultats.append(("dE/dx avec F'=f", str(simplify(dE.subs(diff(Fp(yf), yf), f(yf)).doit()))))
# Autonome : y'' = p dp/dy
P = Function('P')
yy = Function('y')(x)
d2 = diff(P(yy), x)
resultats.append(("(P(y))' = P'(y) y'", str(d2)))

# Tests de reconnaissance par opérateurs d'Euler (Y0 = y, Y1 = p, Y2 = q)
def degre_homogene_y(F):
    return simplify((Y0 * diff(F, Y0) + Y1 * diff(F, Y1) + Y2 * diff(F, Y2)) / F)
def equidimensionnel(F):
    return simplify((x * diff(F, x) - Y1 * diff(F, Y1) - 2 * Y2 * diff(F, Y2)) / F)
resultats.append(("homogène en y : D12", "degré " + str(degre_homogene_y(Y0 * Y2 - Y1**2 - x * Y0**2))))
resultats.append(("homogène en y : D5", "degré " + str(degre_homogene_y(Y0 * Y2 - Y1**2))))
resultats.append(("homogène en y : D6 (non)", str(degre_homogene_y(2 * Y0 * Y2 - 1 - Y1**2))))
resultats.append(("équidimensionnel : D11", "poids " + str(equidimensionnel(x**2 * Y2 - (x * Y1)**2))))
resultats.append(("équidimensionnel : D1", "poids " + str(equidimensionnel(Y2 - Y1 / x))))
resultats.append(("équidimensionnel : D4 (non)", "poids " + str(equidimensionnel(Y2 - 2 * x * Y1**2))))
# Opérateur d'Euler-Lagrange : E(F) = F_y - D(F_p) + D^2(F_q), nul ssi F est une dérivée totale
yfun = Function('y')(x)
def total(F):
    return F.subs({Y2: diff(yfun, x, 2), Y1: diff(yfun, x), Y0: yfun}, simultaneous=True)
def euler(F):
    Fy, Fp, Fq = (total(diff(F, v)) for v in (Y0, Y1, Y2))
    return simplify(Fy - diff(Fp, x) + diff(Fq, x, 2))
resultats.append(("Euler yy''+y'^2 (dérivée de yy')", str(euler(Y0 * Y2 + Y1**2))))
resultats.append(("Euler yy''-y'^2 (non exacte)", str(euler(Y0 * Y2 - Y1**2))))
resultats.append(("Euler (yy''-y'^2)/y^2 = (y'/y)'", str(euler((Y0 * Y2 - Y1**2) / Y0**2))))
resultats.append(("Euler y''-1-y'^2 (non exacte)", str(euler(Y2 - 1 - Y1**2))))

largeur = max(len(r[0]) for r in resultats)
echecs = sum(1 for r in resultats if "ECHEC" in r[1])
for nom, v in resultats:
    print(f"{nom:<{largeur}}  {v}")
print(f"\n{len(resultats)} contrôles, {echecs} échec(s)")
```

</details>
