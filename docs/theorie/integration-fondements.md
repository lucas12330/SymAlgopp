# Intégration symbolique — fondements : Liouville et Risch (cas transcendant)

Ce document pose le cadre algébrique de l'intégration en termes finis : corps différentiels,
extensions élémentaires, théorème de Liouville et structure de l'algorithme de Risch dans le
cas transcendant. Il explique **pourquoi** certaines primitives n'existent pas sous forme
élémentaire (et comment le démontrer), et **comment** un intégrateur complet décide. Les
méthodes pratiques sont dans [`integration-heuristiques.md`](integration-heuristiques.md), le cas
des fractions rationnelles (qui est le cas de base de Risch) dans
[`integration-fractions-rationnelles.md`](integration-fractions-rationnelles.md).

## Table des matières

1. [Notations](#notations)
2. [Corps différentiels](#1-corps-différentiels)
3. [Extensions élémentaires](#2-extensions-élémentaires)
4. [Le théorème de Liouville](#3-le-théorème-de-liouville)
5. [Conséquences : primitives non élémentaires](#4-conséquences--primitives-non-élémentaires)
6. [L'algorithme de Risch, cas transcendant](#5-lalgorithme-de-risch-cas-transcendant)
7. [Exemples travaillés dans une tour](#6-exemples-travaillés-dans-une-tour)
8. [Pièges et cas de contrôle](#7-pièges-et-cas-de-contrôle)
9. [État dans le code](#8-état-dans-le-code)
10. [Références](#9-références)

## Notations

| Symbole | Sens |
|---|---|
| $K$, $E$, $F$ | corps de caractéristique 0 munis d'une dérivation |
| $D$ ou $'$ | la dérivation ; sur $\mathbb{C}(x)$, $D = d/dx$ |
| $\mathrm{Const}(K)$ | sous-corps des constantes $\{c \in K : c' = 0\}$ |
| $C$ | corps des constantes (en pratique $\overline{\mathbb{Q}}$ ou $\mathbb{C}$, supposé algébriquement clos sauf mention) |
| $t$ | monôme : élément transcendant sur $K$ ajouté à $K$ (logarithme, exponentielle, tangente) |
| $K[t]$, $K(t)$ | polynômes et fractions en $t$ à coefficients dans $K$ |
| $\nu_a(y)$ | ordre de $y \in C(x)$ en $a$ : $y = (x-a)^{\nu_a(y)} \cdot u$ avec $u(a) \ne 0, \infty$ ; $\nu_a < 0$ : pôle d'ordre $-\nu_a$ |
| $\mathrm{Res}_t(P, Q)$ | résultant de $P$ et $Q$ par rapport à $t$ |
| « élémentaire » | appartenant à une extension élémentaire de $C(x)$ (section 2) |

Toutes les primitives sont définies à une constante additive près ; « $F$ primitive de $f$ »
signifie $F' = f$ dans le corps différentiel considéré (donc sur un intervalle où tout est
défini, pour les fonctions).

---

## 1. Corps différentiels

### 1.1 Définitions

**Définition.** Une *dérivation* sur un corps $K$ est une application $D : K \to K$ additive
qui vérifie la règle de Leibniz $D(uv) = u\,D(v) + v\,D(u)$. Le couple $(K, D)$ est un *corps
différentiel*. Les éléments de dérivée nulle forment un sous-corps $\mathrm{Const}(K)$.

Conséquences immédiates (récurrence et Leibniz) :

- $D(1) = 0$, $D(u/v) = (v\,Du - u\,Dv)/v^2$, $D(u^n) = n\,u^{n-1} Du$ pour $n \in \mathbb{Z}$ ;
- **dérivée logarithmique** : pour $u, v \ne 0$ et $m, n \in \mathbb{Z}$,
  $$\frac{D(u^m v^n)}{u^m v^n} = m\,\frac{Du}{u} + n\,\frac{Dv}{v}.$$
  C'est l'identité qui fait que $\log$ « transforme les produits en sommes » sans jamais
  écrire de logarithme : elle est au cœur de la preuve de Liouville.

**Exemple fondamental.** $(\mathbb{C}(x), d/dx)$, de constantes $\mathbb{C}$. Plus généralement
$(C(x), d/dx)$ pour tout corps $C$ de caractéristique 0.

### 1.2 Monômes ; polynômes normaux et spéciaux

Soit $(K, D)$ un corps différentiel et $E = K(t) \supset K$ une extension différentielle avec
$t$ transcendant sur $K$ et $D(t) \in K[t]$. On dit que $t$ est un *monôme* sur $K$
(terminologie de Bronstein). Pour $p \in K[t]$, $D(p) \in K[t]$.

**Définition.** Un polynôme $p \in K[t]$ est *normal* si $\gcd(p, Dp) = 1$, *spécial* si
$p \mid Dp$.

**Lemme 1 (cas logarithmique).** Si $Dt \in K$ et $\mathrm{Const}(K(t)) = \mathrm{Const}(K)$,
alors tout $p \in K[t]$ irréductible est normal.

*Preuve.* Écrivons $p$ unitaire de degré $n \ge 1$. Comme $D(t^k) = k\,t^{k-1}\,Dt$ est de
degré $< k$, $D(p)$ est de degré $< n$. Si $p \mid Dp$, alors $Dp = 0$, donc $p$ est une
constante de $K(t)$, donc un élément de $K$ : contradiction avec $n \ge 1$. Comme $p$ est
irréductible, $\gcd(p, Dp) \in \{1, p\}$, d'où $\gcd(p, Dp) = 1$. $\square$

**Lemme 2 (cas exponentiel).** Si $Dt = \eta'\,t$ avec $\eta \in K$ et
$\mathrm{Const}(K(t)) = \mathrm{Const}(K)$, alors $t$ est spécial, et tout $p \in K[t]$
irréductible non associé à $t$ est normal.

*Preuve.* $D(t) = \eta' t$ est divisible par $t$. Soit $p$ unitaire irréductible de degré
$n \ge 1$, $p \ne t$, donc $p(0) \ne 0$. Comme $D(a\,t^k) = (a' + k\eta' a)\,t^k$, $D$ conserve
le degré et $Dp = n\eta'\,t^n + \dots$. Si $p \mid Dp$, alors $Dp = n\eta'\,p$ (comparaison des
termes dominants), donc
$D(p/t^n) = (Dp\cdot t^n - p\cdot n\eta' t^n)/t^{2n} = 0$ : $p/t^n$ est une constante, donc
$p = c\,t^n$, contradiction avec $p(0) \ne 0$. $\square$

Les monômes **tangents** $Dt = \eta'(1 + t^2)$ (cas de $\tan \eta$) ont pour polynôme spécial
$1 + t^2$ ; Bronstein les traite directement, ce qui évite de passer par $e^{i\eta}$.

**Pourquoi ces notions comptent.** Si $p$ est normal et $y \in K(t)$ a un pôle d'ordre
$k \ge 1$ en $p$ (c'est-à-dire $y = a/(p^k b)$ avec $p \nmid ab$), alors $Dy$ a un pôle
d'ordre exactement $k + 1$ en $p$ : la dérivation **augmente** strictement l'ordre des pôles
normaux. Une dérivée n'a donc jamais de pôle simple en un polynôme normal, alors qu'une
dérivée logarithmique $Du/u$ n'a **que** des pôles simples. Toute la théorie (Liouville,
Hermite, Rothstein–Trager) repose sur cette dichotomie.

*Preuve de l'augmentation d'ordre.* $D(a/(p^k b)) = \dfrac{Da}{p^k b} - k\dfrac{a\,Dp}{p^{k+1} b}
- \dfrac{a\,Db}{p^k b^2}$ ; le terme du milieu a un pôle d'ordre $k+1$ car $p \nmid a\,Dp$
(normalité) et $p \nmid b$, les deux autres un pôle d'ordre $\le k$. $\square$

---

## 2. Extensions élémentaires

### 2.1 Définitions

Soit $(K, D)$ un corps différentiel et $E \supset K$ une extension différentielle. Un élément
$t \in E$ est :

- **logarithmique** sur $K$ s'il existe $u \in K^*$ tel que $Dt = Du/u$ (« $t = \log u$ ») ;
- **exponentiel** sur $K$ s'il existe $u \in K$ tel que $Dt = t\,Du$ (« $t = \exp u$ ») ;
- **algébrique** sur $K$ s'il annule un polynôme non nul de $K[T]$.

$E$ est une **extension élémentaire** de $K$ s'il existe une tour
$K = K_0 \subset K_1 \subset \dots \subset K_n = E$ avec $K_i = K_{i-1}(t_i)$ et chaque $t_i$
logarithmique, exponentiel ou algébrique sur $K_{i-1}$. Une **fonction élémentaire** est un
élément d'une extension élémentaire de $C(x)$ ; $f$ a une **primitive élémentaire** s'il
existe $F$ élémentaire (dans une extension élémentaire de $C(x)$ contenant $f$) avec $F' = f$.

On exige en outre, pour les énoncés qui suivent, que l'extension **n'ajoute pas de
constantes** : $\mathrm{Const}(E) = \mathrm{Const}(K)$. Sans cette hypothèse, on pourrait
« intégrer » n'importe quoi en adjoignant formellement une primitive qui serait une nouvelle
constante d'un autre corps ; et des extensions naturelles en ajoutent (adjoindre $e^{x}$ et
$e^{x+1}$ séparément fait apparaître la constante $e^{x+1}/e^x = e$ si elle n'était pas dans
$C$).

### 2.2 Les fonctions usuelles sont élémentaires (sur $\mathbb{C}$)

Avec des constantes complexes :

$$
\sin u = \frac{e^{iu} - e^{-iu}}{2i},\quad \cos u = \frac{e^{iu} + e^{-iu}}{2},\quad
\tan u = -i\,\frac{e^{2iu} - 1}{e^{2iu} + 1},\quad
\arctan u = \frac{i}{2}\log\frac{1 - iu}{1 + iu},
$$
$$
\arcsin u = -i\log\!\left(iu + \sqrt{1-u^2}\right),\qquad u^{r} = \exp(r\log u),\qquad
\sqrt{u}\ \text{algébrique}\ (T^2 - u).
$$

(Les identités en $\log$ s'entendent comme égalités de dérivées, c'est-à-dire à une constante
additive près sur chaque intervalle, ce qui suffit pour l'intégration.) Les fonctions
hyperboliques sont des fractions rationnelles en $e^u$. Une primitive « élémentaire réelle »
est donc en particulier élémentaire sur $\mathbb{C}$ : **prouver la non-existence sur
$\mathbb{C}$ suffit**.

### 2.3 Indépendance : le théorème de structure

Une même fonction peut être décrite par des tours différentes, et une tour naïve peut être
illégitime : $e^{2x} = (e^x)^2$ n'est pas transcendant sur $C(x, e^x)$ ;
$\log(x^2) = 2\log x$ ; $e^{\log x} = x$ ; $\exp(\tfrac12 \log x) = \sqrt{x}$ est algébrique.
L'algorithme de Risch exige des monômes **transcendants sans nouvelles constantes** ; avant de
l'appliquer, il faut donc réécrire l'intégrande dans une tour réduite. Le *théorème de
structure* de Risch (traité chez Bronstein) donne un critère effectif : un nouveau logarithme
$\log u$ (resp. une nouvelle exponentielle $\exp u$) est algébrique sur la tour courante si et
seulement si $u$ (resp. $u'$) s'exprime par une relation linéaire à coefficients rationnels
avec les arguments des exponentielles et logarithmes déjà présents. En pratique : décider si
$\exp(u)$ est une puissance rationnelle d'un produit d'exponentielles déjà présentes, si
$\log(u)$ est une combinaison rationnelle de logarithmes présents, etc. La forme canonique de
SymAlgo++ fait une petite partie de ce travail (exposants rationnels de $e^{v}$ reconnus par
la branche en cours, voir section 8).

---

## 3. Le théorème de Liouville

### 3.1 Énoncé

**Théorème (Liouville, forme forte, telle que l'énonce Bronstein).** Soit $(K, D)$ un corps
différentiel de caractéristique 0, de corps des constantes $C$, et $f \in K$. S'il existe une
extension élémentaire $E$ de $K$ avec $\mathrm{Const}(E) = C$ et $g \in E$ tel que $g' = f$,
alors il existe $v \in K$, des constantes $c_1, \dots, c_m \in \overline{C}$ et
$u_1, \dots, u_m \in K(c_1, \dots, c_m)^*$ tels que
$$
f = v' + \sum_{i=1}^{m} c_i\,\frac{u_i'}{u_i}.
$$
Autrement dit, si $f$ a une primitive élémentaire, elle en a une de la forme
$$
\int f = v + \sum_i c_i \log u_i,
$$
où **tout est déjà dans $K$**, aux constantes algébriques $c_i$ près : la seule chose qu'une
primitive élémentaire peut ajouter à $K$, ce sont des **logarithmes à coefficients
constants**. En particulier aucune nouvelle exponentielle, aucun nouveau radical.

Si $C$ est algébriquement clos (cas $\overline{\mathbb{Q}}$, $\mathbb{C}$), on peut prendre
$c_i \in C$ et $u_i \in K^*$.

### 3.2 Idée de la preuve

On raisonne par récurrence sur la longueur de la tour $K = K_0 \subset \dots \subset K_n = E$.
L'hypothèse de récurrence, appliquée au corps $K_1 = K(t)$ (qui contient $f$) et à la tour
plus courte $K_1 \subset \dots \subset K_n$, donne
$f = v' + \sum c_i u_i'/u_i$ avec $v, u_i \in K(t)$. Il reste à « redescendre » dans $K$. On
peut supposer les $c_i$ linéairement indépendants sur $\mathbb{Q}$ (regrouper les termes :
si $c_3 = 2c_1 - c_2$, remplacer $u_1$ par $u_1 u_3^2$, $u_2$ par $u_2/u_3$, grâce à
l'identité de la dérivée logarithmique), et $C$ algébriquement clos pour simplifier.

**Cas $t$ logarithmique ou exponentiel (transcendant).** Écrivons chaque $u_i = w_i \prod_p
p^{e_{i,p}}$ avec $w_i \in K^*$, $p$ parcourant les irréductibles unitaires de $K[t]$ et
$e_{i,p} \in \mathbb{Z}$. Alors
$$
\sum_i c_i \frac{u_i'}{u_i} = \sum_i c_i\frac{w_i'}{w_i} + \sum_p \Big(\sum_i c_i e_{i,p}\Big)\frac{p'}{p}.
$$
Soit $p$ normal. Le membre de gauche $f \in K$ n'a pas de pôle en $p$. $v'$ a en $p$ soit aucun
pôle, soit un pôle d'ordre $\ge 2$ (section 1.2) ; le terme en $p'/p$ a un pôle simple, de
coefficient $\sum_i c_i e_{i,p}$, et aucun autre terme n'a de pôle en $p$. Donc $v$ n'a pas de
pôle en $p$, puis $\sum_i c_i e_{i,p} = 0$, ce qui force $e_{i,p} = 0$ pour tout $i$
(indépendance des $c_i$ sur $\mathbb{Q}$). Conclusion :

- cas logarithmique (tout irréductible est normal) : $v \in K[t]$ et $u_i \in K^*$ ;
- cas exponentiel ($t$ seul spécial) : $v \in K[t, t^{-1}]$ et $u_i = w_i\,t^{e_i}$.

*Logarithmique,* $t' = a'/a$, $a \in K$. Si $v = \sum_{j \le d} v_j t^j$ avec $d \ge 2$ et
$v_d \ne 0$ : le coefficient de $t^d$ dans $v'$ est $v_d'$, il doit être nul ($f$ et les
$w_i'/w_i$ sont dans $K$), donc $v_d \in C$ ; celui de $t^{d-1}$ est
$v_{d-1}' + d\,v_d\,t'$, nul lui aussi, donc $t' = (-v_{d-1}/(d v_d))'$ et
$t + v_{d-1}/(d v_d)$ est une constante : $t \in K$, contradiction. Donc $v = v_1 t + v_0$ avec
$v_1 \in C$ (même argument pour $d = 1$), et
$f = v_0' + v_1\,a'/a + \sum_i c_i w_i'/w_i$ : c'est la forme voulue dans $K$, avec un
logarithme de plus.

*Exponentiel,* $t' = \eta' t$, $\eta \in K$. Comme $(v_j t^j)' = (v_j' + j\eta' v_j)\,t^j$, les
coefficients de $t^j$, $j \ne 0$, donnent $(v_j t^j)' = 0$, donc $v_j t^j \in C$, ce qui est
impossible si $v_j \ne 0$ ($t^j$ serait dans $K$). Donc $v = v_0 \in K$ et, avec
$u_i'/u_i = w_i'/w_i + e_i\eta'$,
$f = \big(v_0 + \sum_i c_i e_i \eta\big)' + \sum_i c_i w_i'/w_i$. $\square$

**Cas $t$ algébrique.** On passe à une clôture normale et on applique la trace
$\mathrm{Tr} : K(t) \to K$ : elle commute avec la dérivation et
$\mathrm{Tr}(u'/u) = N(u)'/N(u)$ ($N$ = norme). Diviser par le degré $n = [K(t) : K]$ donne
$f = (\mathrm{Tr}(v)/n)' + \sum_i (c_i/n)\,N(u_i)'/N(u_i)$, avec $\mathrm{Tr}(v)/n$,
$N(u_i) \in K$. $\square$

### 3.3 Corollaire : le critère de Liouville pour $\int f\,e^{g}$

**Corollaire.** Soient $f, g \in C(x)$, $f \ne 0$, $g$ non constant. Alors $\int f e^{g}\,dx$
est élémentaire **si et seulement si** l'équation
$$
R' + g'\,R = f
$$
a une solution $R \in C(x)$ ; dans ce cas $\int f e^g = R\,e^g$.

*Preuve.* « Si » : $(R e^g)' = (R' + g'R)e^g$. « Seulement si » : $t = e^g$ est transcendant
sur $C(x)$ quand $g$ n'est pas constant (résultat classique admis ici) et n'ajoute pas de
constante. Appliquons Liouville dans $K(t)$ avec $K = C(x)$ : $f t = v' + \sum c_i u_i'/u_i$
avec $v, u_i \in K(t)$. L'argument du cas exponentiel ci-dessus (sans la descente finale)
montre que $v \in K[t, t^{-1}]$ et $u_i \in K^* t^{e_i}$, donc que le terme logarithmique
est dans $K$. En comparant les coefficients de $t^1$ : $f = v_1' + g' v_1$. $\square$

C'est exactement l'**équation différentielle de Risch** (section 5.6) dans le cas le plus
simple.

---

## 4. Conséquences : primitives non élémentaires

Toutes les preuves suivent le même schéma : (1) placer l'intégrande dans une tour
$C(x)(t)$ ; (2) Liouville réduit l'existence d'une primitive à une équation dans $C(x)$ ; (3)
une comparaison d'**ordres de pôles** et de **degrés** montre que cette équation n'a pas de
solution rationnelle.

**Lemme d'ordre.** Pour $y \in C(x)$ et $a \in C$, si $\nu_a(y) = -m < 0$ alors
$\nu_a(y') = -m - 1$. Si $y$ est un polynôme de degré $d \ge 1$, $y'$ est de degré $d - 1$.
(Dériver $(x-a)^{-m}u$ avec $u(a) \ne 0$.)

### 4.1 $\int e^{x^2}\,dx$ n'est pas élémentaire (preuve complète)

Par le corollaire 3.3 avec $f = 1$, $g = x^2$, il suffit de montrer que
$$
R' + 2x\,R = 1 \qquad (\ast)
$$
n'a pas de solution $R \in C(x)$.

1. *Pas de pôle fini.* Si $R$ avait un pôle d'ordre $m \ge 1$ en $a \in C$, $R'$ y aurait un
   pôle d'ordre $m + 1$ et $2xR$ un pôle d'ordre au plus $m$ : le membre de gauche de $(\ast)$
   aurait un pôle d'ordre $m + 1$, alors que $1$ n'en a pas. Donc $R$ est un polynôme.
2. *Degré.* $R = 0$ ne convient pas. Si $R$ est un polynôme de degré $d \ge 0$, $2xR$ est de
   degré $d + 1 \ge 1$ et $R'$ de degré $< d$ : le membre de gauche est de degré
   $d + 1 \ge 1$, il ne peut valoir $1$.

Donc $(\ast)$ n'a pas de solution rationnelle et $\int e^{x^2}$ n'est pas élémentaire. $\square$

Le même argument vaut pour $e^{-x^2}$ ($R' - 2xR = 1$), d'où la fonction d'erreur
$\mathrm{erf}$ ; et pour $e^{P(x)}$ avec $P$ polynôme de degré $\ge 2$ (le degré de $P'R$
dépasse celui de $R'$).

**Ce qui reste élémentaire.** Pour $\int x^{k} e^{x^2}$, l'équation est $R' + 2xR = x^k$ ;
$R$ est un polynôme de degré $k - 1$.

- $k = 1$ : $R = 1/2$, $\int x e^{x^2} = \tfrac12 e^{x^2}$.
- $k = 2$ : $R = \alpha x + \beta$ donne $\alpha + 2\alpha x^2 + 2\beta x = x^2$, soit
  $\alpha = 1/2$ et $\alpha = 0$ : **pas de solution**, $\int x^2 e^{x^2}$ n'est pas
  élémentaire (par parties, $= \tfrac12 x e^{x^2} - \tfrac12\int e^{x^2}$).
- $(2x^2 + 1)e^{x^2}$ : $R = x$ convient ($1 + 2x^2$), d'où $\int (2x^2+1)e^{x^2} = x e^{x^2}$.
- $k = 3$ : $R = \tfrac12(x^2 - 1)$, $\int x^3 e^{x^2} = \tfrac12 (x^2-1) e^{x^2}$.

En général, pour $k \in \mathbb{N}$, $\int x^k e^{x^2}$ est élémentaire si et seulement si
$k$ est impair. *Preuve.* $R$ est un polynôme de degré $k - 1$ (pas de pôle : $2x$ et $x^k$
n'en ont pas). Le coefficient de $x^j$ dans $R' + 2xR$ est $(j+1)r_{j+1} + 2r_{j-1}$ ; en
partant de $r_{k-1} = 1/2$, la récurrence ne touche que les $r_j$ de la parité de $k - 1$,
tous non nuls. Pour $k$ impair, elle s'arrête sur $r_0$ et l'équation en $x^0$, $r_1 = 0$, est
compatible. Pour $k$ pair, elle impose $r_1 \ne 0$ alors que l'équation en $x^0$ exige
$r_1 = 0$ : contradiction.

### 4.2 $\int \frac{\sin x}{x}\,dx$ n'est pas élémentaire

Sur $\mathbb{C}$, avec $t = e^{ix}$ (transcendant sur $\mathbb{C}(x)$) :
$\dfrac{\sin x}{x} = \dfrac{t - t^{-1}}{2ix}$. Par Liouville et l'argument du cas exponentiel,
une primitive élémentaire imposerait $v \in K[t, t^{-1}]$, et le coefficient de $t^1$ donnerait
$$
R' + i R = \frac{1}{2ix} = -\frac{i}{2x}, \qquad R \in \mathbb{C}(x).
$$
- Un pôle de $R$ en $a \ne 0$ d'ordre $m$ donnerait un pôle d'ordre $m + 1$ à gauche, aucun à
  droite : impossible.
- Un pôle en $0$ d'ordre $m \ge 1$ donnerait un pôle d'ordre $m + 1 \ge 2$ à gauche, d'ordre 1
  à droite : impossible.
- Donc $R$ est un polynôme, $R' + iR$ aussi, et il ne peut valoir $-i/(2x)$.

Pas de primitive élémentaire (même sur $\mathbb{C}$, donc a fortiori sur $\mathbb{R}$) : c'est
le sinus intégral $\mathrm{Si}$. Même preuve pour $\cos(x)/x$ et pour $\sin(x^2)$,
$\cos(x^2)$ (équation $R' + 2ixR = \tfrac{1}{2i}$, degré : intégrales de Fresnel).

### 4.3 $\int \frac{e^x}{x}\,dx$, $\int \frac{dx}{\ln x}$ et $\int e^{e^x}dx$

*$e^x/x$.* Corollaire 3.3 avec $f = 1/x$, $g = x$ : $R' + R = 1/x$. Un pôle en $a \ne 0$ est
impossible (ordre $m + 1$ contre 0) ; un pôle en 0 d'ordre $m \ge 1$ donne $m + 1 \ge 2$ contre
1 ; un polynôme donne un polynôme. Pas de solution : $\int e^x/x$ (l'exponentielle intégrale
$\mathrm{Ei}$) n'est pas élémentaire.

*$1/\ln x$, par changement de variable.* Si $F$ élémentaire vérifiait $F'(x) = 1/\ln x$, alors
$G(u) = F(e^u)$ serait élémentaire (sur $C(u)$, en ajoutant $e^u$ à la tour) avec
$G'(u) = e^u/u$, ce qui vient d'être exclu.

*$1/\ln x$, directement par le critère des résidus.* Dans la tour $K(t)$, $K = C(x)$,
$t = \ln x$, $Dt = 1/x$, l'intégrande est $h = 1/t$ : dénominateur $d = t$ normal (lemme 1) et
sans carré, numérateur $a = 1$. Le résultant de Rothstein–Trager (section 5.4) vaut
$r(z) = \mathrm{Res}_t(t,\ 1 - z\,Dt) = 1 - z/x$, dont la racine $z = x$ **n'est pas une
constante** : $\int 1/\ln x$ (le logarithme intégral $\mathrm{li}$) n'est pas élémentaire.

*$e^{e^x}$.* Le changement $u = e^x$ donne $\int e^u/u\,du$, non élémentaire.

### 4.4 Autres classiques (énoncés seulement)

Les primitives de $\sqrt{1 - k^2\sin^2 x}$ ($0 < k < 1$, intégrales elliptiques) et de
$1/\sqrt{x^3 + 1}$ ne sont pas élémentaires ; la preuve relève du cas **algébrique** de Risch
(courbes de genre 1), hors du périmètre de ce document : énoncés admis.

---

## 5. L'algorithme de Risch, cas transcendant

### 5.1 Entrée et principe

**Entrée** : $f \in K_n = C(x)(t_1, \dots, t_n)$ où chaque $t_i$ est un monôme logarithmique,
exponentiel (ou tangent) sur $K_{i-1}$, transcendant et sans nouvelle constante (section 2.3).
**Sortie** : une primitive élémentaire de $f$, ou la **preuve** qu'il n'en existe pas.

**Principe** : on travaille dans le dernier étage $K(t)$, $K = K_{n-1}$, $t = t_n$, en
traitant $f$ comme une fraction rationnelle en $t$ à coefficients dans $K$, et on se ramène à
des problèmes dans $K$ (intégration, équation de Risch), résolus récursivement. Au bas de la
tour, $K_0 = C(x)$ : intégration des fractions rationnelles (document dédié) et équation de
Risch sur $C(x)$.

### 5.2 Décomposition de l'intégrande

Par division euclidienne en $t$ et factorisation sans carré du dénominateur, on écrit
$f = p + \dfrac{a_s}{d_s} + \dfrac{a_n}{d_n}$ où $p \in K[t]$, $d_s$ est le produit des
facteurs spéciaux du dénominateur (puissances de $t$ dans le cas exponentiel, de $1 + t^2$
dans le cas tangent, rien dans le cas logarithmique), $d_n$ le produit des facteurs normaux,
$\deg a_n < \deg d_n$. Dans le cas exponentiel, la partie spéciale $a_s/t^k$ se range avec $p$
en un **polynôme de Laurent** de $K[t, t^{-1}]$.

### 5.3 Réduction d'Hermite dans $K(t)$

La réduction d'Hermite du document sur les fractions rationnelles s'applique **telle quelle**
à la partie normale, en remplaçant $d/dx$ par $D$ : elle repose sur l'existence de Bézout
$S\,U\,DV + T\,V = 1$ pour $V$ normal sans carré, premier avec $U$, ce qui est exactement la
normalité. On obtient
$$
f = D(g) + h + r,
$$
où $g \in K(t)$, $h$ est **simple** (dénominateur normal sans carré, numérateur de degré
inférieur) et $r$ est un polynôme (de Laurent dans le cas exponentiel). Aucun calcul de
racine n'est nécessaire.

### 5.4 Partie logarithmique : critère des résidus

Pour $h = a/d$ simple, on forme
$$
r(z) = \mathrm{Res}_t\big(d,\ a - z\,D(d)\big) \in K[z].
$$

**Critère des résidus.** Si $\int h$ (plus éventuellement un polynôme de $K[t]$) est
élémentaire, alors toutes les racines de $r$ sont des **constantes** : $r = \kappa\,\tilde r$
avec $\kappa \in K$ et $\tilde r \in C[z]$. Dans ce cas
$$
h - D\Big(\sum_{\tilde r(\alpha) = 0} \alpha \log \gcd\big(d,\ a - \alpha\,D(d)\big)\Big) \in K[t]
$$
(polynôme, à cause du terme $\deg(v)\,\eta'$ que produit $D(v)/v$ dans le cas exponentiel).
Si une racine n'est pas constante, $f$ n'a **pas** de primitive élémentaire. C'est la
généralisation directe de Rothstein–Trager (le cas $K = C$, $t = x$). Le résultant se calcule
par sous-résultants comme dans Lazard–Rioboo–Trager, ce qui évite de factoriser $d$.

*Pourquoi.* Liouville impose que les coefficients des logarithmes soient des constantes ; or
le coefficient de $\log(t - \beta)$ (sur une clôture algébrique) devant le pôle simple en
$\beta$ est le résidu $a(\beta)/D(d)(\beta)$, racine de $r$. Un résidu non constant ne peut pas
être un coefficient de logarithme, et un pôle simple ne peut pas venir d'une dérivée (section
1.2).

### 5.5 Partie polynomiale, cas logarithmique

$t' = \lambda \in K$ ($t = \log u$, $\lambda = u'/u$). On cherche $q \in K[t]$ et des
logarithmes de $K$ tels que $D(q) + (\text{logs})' = p = \sum_{j=0}^{n} p_j t^j$.

Le degré de $q$ est au plus $n + 1$ et son coefficient dominant $q_{n+1}$ est une constante
(comparaison des coefficients, comme dans la preuve de Liouville). Le coefficient de $t^j$
dans $D(q)$ est $q_j' + (j+1)\,q_{j+1}\,\lambda$. On écrit chaque coefficient
$q_j = w_j + b_j$ avec $w_j \in K$ calculé à l'étage $j$ et $b_j \in C$ une constante
d'intégration qui n'est fixée qu'à l'étage suivant, et on descend :

```
IntegrerPolynomeLog(p = Σ_{j=0..n} p_j t^j, t' = λ) :
    w_{n+1} ← 0                                   // q_{n+1} = b_{n+1} est une constante
    pour j de n à 0 :
        g ← p_j − (j+1) · w_{j+1} · λ             // il reste : q_j' + (j+1) b_{j+1} λ = g
        si j > 0 :
            (w_j, c) ← IntegrationLimitee(g, λ)   // w_j ∈ K, c ∈ C, w_j' + c λ = g
        sinon :
            (w_0, c, logs) ← Integrer_dans_K(g)   // tous les logarithmes de K permis ;
                                                  // c = coefficient de log u = t obtenu
        si échec : renvoyer NON_ELEMENTAIRE
        b_{j+1} ← c / (j+1)
    renvoyer Σ_{j=0..n+1} (w_j + b_j) t^j + logs     (b_0 = 0)
```

L'**intégration limitée** « trouver $w \in K$ et $c \in C$ avec $w' + c\lambda = g$ » est un
appel récursif à l'intégrateur de $K$ où l'on autorise exactement un logarithme connu,
$t = \log u$ ; son coefficient $c$ fixe la constante $b_{j+1}$ de l'étage au-dessus. Exemples
complets en section 6.

### 5.6 Partie polynomiale, cas exponentiel : l'équation de Risch

$t' = \eta' t$ ($t = e^\eta$). Pour $p = \sum_{j=-m}^{n} p_j t^j \in K[t, t^{-1}]$, on cherche
$q = \sum q_j t^j$ avec, puisque $(q_j t^j)' = (q_j' + j\eta' q_j)\,t^j$ :

- pour $j \ne 0$ : $q_j \in K$ solution de l'**équation différentielle de Risch**
  $$y' + j\eta'\,y = p_j ;$$
- pour $j = 0$ : $q_0 = \int p_0$, appel récursif dans $K$ (logarithmes permis), corrigé des
  termes $\deg(v)\eta'$ issus de la partie logarithmique.

Si l'une des équations de Risch ($j \ne 0$, $p_j \ne 0$) n'a pas de solution dans $K$, $f$
n'a **pas** de primitive élémentaire (l'unicité de la décomposition par puissances de $t$ le
garantit, comme dans la preuve de Liouville).

### 5.7 L'équation de Risch sur $C(x)$

Problème : $f, g \in C(x)$, trouver $y \in C(x)$ avec $y' + f y = g$.

**Étape 1 : borne du dénominateur (analyse locale en chaque $a \in C$).** Soit
$\nu = \nu_a$ et supposons $\nu(y) = -m < 0$.

- $\nu(f) \ge 0$ : $\nu(y' + fy) = -m - 1$, donc $\nu(g) = -m - 1$ : $y$ n'a de pôles qu'aux
  pôles de $g$, d'ordre diminué de 1.
- $\nu(f) \le -2$ : $\nu(fy) = \nu(f) - m < -m - 1$, donc $m = \nu(f) - \nu(g)$.
- $\nu(f) = -1$, de résidu $\rho$ : le coefficient de $(x-a)^{-m-1}$ dans $y' + fy$ est
  $(\rho - m)\,y_{-m}$ ; il peut s'annuler si $\rho = m$ est un **entier positif**. Alors
  $m \le \max(\rho, -\nu(g) - 1)$ ; sinon $m = -\nu(g) - 1$.

On en déduit un polynôme $h$ multiple de tout dénominateur possible ; on pose $y = q/h$ et
l'équation devient une équation polynomiale en $q$.

**Étape 2 : borne de degré (analyse à l'infini).** Pour $f, g$ polynômes et $y$ polynôme :
- $\deg f \ge 1$ : $\deg(fy) > \deg y'$, donc $\deg y = \deg g - \deg f$ ;
- $f$ constante non nulle : $\deg y = \deg g$ ;
- $f = 0$ : intégration simple, $\deg y = \deg g + 1$.

Pour $f$ rationnel, la même comparaison à l'infini donne une borne, avec un cas de
compensation possible quand $f \sim -n/x$ à l'infini ($n$ entier positif), qui se traite comme
le cas des résidus entiers.

**Étape 3 : résolution.** Avec $y = q/h$ et $\deg q \le N$, l'équation est **linéaire** en
les coefficients de $q$ : système linéaire (ou la méthode plus économique de Rothstein,
« SPDE », qui réduit le degré pas à pas). Absence de solution = preuve de non-élémentarité
pour le coefficient concerné.

**Exemples.** $e^{x^2}$ : $f = 2x$, $g = 1$, pas de pôle, $\deg y = 0 - 1 < 0$ : impossible
(section 4.1). $x^3 e^{x^2}$ : $\deg y = 3 - 1 = 2$, solution $y = (x^2-1)/2$.

Sur une tour $K = C(x)(t_1, \dots)$, l'équation de Risch se résout de même en analysant les
ordres en chaque polynôme normal et le degré en $t$, avec des cas d'annulation plus fins
(traités par Bronstein) ; c'est la partie la plus technique de l'algorithme.

### 5.8 Pseudo-code global

```
Risch(f, tour C(x) = K_0 ⊂ ... ⊂ K_n = K(t)) :
    si n = 0 : renvoyer IntegrerFractionRationnelle(f)          // toujours élémentaire
    (g, h, r) ← HermiteReduction(f, D)        // f = D(g) + h + r, h simple, r polynôme (de Laurent)
    (L, ok)   ← ResidueReduce(h)               // résultant, critère des résidus, logs
    si non ok : renvoyer NON_ELEMENTAIRE
    p ← r + h − D(L)                            // dans K[t] (ou K[t, 1/t])
    selon le type de t :
        logarithmique : (q, ok) ← IntegrerPolynomeLog(p)    // appels récursifs dans K
        exponentiel   : (q, ok) ← pour chaque j ≠ 0 : RischDE(j η', p_j) dans K ;
                                  q_0 ← Risch(p_0, tour jusqu'à K)
        tangent       : analogue (polynôme spécial 1 + t²)
    si non ok : renvoyer NON_ELEMENTAIRE
    renvoyer g + L + q
```

### 5.9 Cas algébrique et mixte

Quand la tour contient des radicaux ou des fonctions algébriques ($\sqrt{x^3+1}$), la
réduction d'Hermite se fait dans une base intégrale du corps de fonctions algébriques et la
partie logarithmique demande de décider si un diviseur est de torsion (Trager, puis Bronstein
pour le cas mixte). Ce document **ne détaille pas** ce cas ; les implémentations complètes
sont rares.

### 5.10 Complexité, décidabilité, limites

- **Décidabilité.** Pour une tour transcendante donnée, l'algorithme décide l'existence d'une
  primitive élémentaire, **à condition de savoir tester la nullité des constantes** (les
  coefficients des résultants, les résidus). Pour des constantes construites avec $\exp$ et
  $\log$, ce test n'est connu décidable que sous des conjectures de transcendance
  (conjecture de Schanuel) ; avec $\sin$ et des valeurs absolues, la nullité des
  *fonctions* devient indécidable (théorème de Richardson, voir
  [`integration-heuristiques.md`](integration-heuristiques.md#8-vérifier-une-primitive)).
- **Coût.** Pas de borne simple : le coût est dominé par les résultants (degré en $z$ jusqu'au
  degré du dénominateur), les équations de Risch récursives et surtout la croissance des
  expressions. Chaque étage de la tour multiplie le travail.
- **Mise en œuvre.** Il faut des polynômes à coefficients dans un corps différentiel arbitraire
  (récursif), un test de zéro fiable et la construction de la tour (théorème de structure).
  C'est pourquoi les intégrateurs pratiques essaient d'abord des heuristiques (table,
  substitution, parties, Risch–Norman) et réservent Risch aux cas qu'elles ne résolvent pas,
  ou s'en servent pour **prouver** qu'une intégrale n'est pas élémentaire.

---

## 6. Exemples travaillés dans une tour

Chaque résultat a été vérifié par dérivation (voir section 7).

**(a) $\int \ln x$.** $t = \ln x$, $\lambda = 1/x$, $p = t$ ($n = 1$). Étage $j = 1$ :
$g = p_1 = 1$, intégration limitée $w_1' + c/x = 1$ : $w_1 = x$, $c = 0$, donc $b_2 = 0$.
Étage $j = 0$ : $g = p_0 - w_1\lambda = -1$, $w_0 = -x$, $c = 0$, donc $b_1 = 0$. Résultat
$\int \ln x = x\ln x - x$.

**(a') $\int \frac{\ln x}{x}$.** $p = t/x$. Étage 1 : $w_1' + c/x = 1/x$ : $w_1 = 0$, $c = 1$,
donc $b_2 = 1/2$. Étage 0 : $g = 0$, $w_0 = 0$, $b_1 = 0$. Résultat $\tfrac12 t^2 =
\tfrac12\ln^2 x$ : c'est le **logarithme lui-même** produit par l'intégration limitée qui fait
monter le degré.

**(b) $\int \ln^2 x$.** $p = t^2$. Étage 2 : $w_2' + c/x = 1$, $w_2 = x$, $b_3 = 0$. Étage 1 :
$g = 0 - 2w_2\lambda = -2$, $w_1 = -2x$, $b_2 = 0$. Étage 0 : $g = 0 - w_1\lambda = 2$,
$w_0 = 2x$, $b_1 = 0$. Résultat $x\ln^2 x - 2x\ln x + 2x$.

**(c) $\int \frac{dx}{x\ln x}$.** $h = (1/x)/t$ simple. $r(z) = \mathrm{Res}_t(t, 1/x - z/x)
= (1 - z)/x$, racine $z = 1$ constante. $\gcd(t, 1/x - 1\cdot 1/x) = \gcd(t, 0) = t$, d'où
$\int = \log t = \ln(\ln x)$ (au signe près de l'argument sur $\mathbb{R}$ :
$\ln\lvert\ln x\rvert$).

**(d) $\int \left(\frac{1}{\ln x} - \frac{1}{\ln^2 x}\right)dx$.** $f = (t - 1)/t^2$,
dénominateur $t^2$ normal. Hermite avec $V = t$, $U = 1$, $DV = 1/x$ : on cherche $B$ avec
$B\cdot DV \equiv -(t - 1) \pmod t$, soit $B/x = 1$, $B = x$. Partie rationnelle $x/t$, et
$f - D(x/t) = (t-1)/t^2 - (1/t - 1/t^2) = 0$. Résultat $x/\ln x$. **Piège** : chacun des deux
termes séparément n'est pas élémentaire (section 4.3) ; la linéarité appliquée trop tôt fait
échouer un intégrateur naïf.

**(e) $\int e^x\frac{x - 1}{x^2}dx$.** $t = e^x$, $f = \frac{x-1}{x^2}\,t$. Coefficient de
$t$ : $y' + y = (x-1)/x^2$. Borne : pôle double de $g$ en 0, $f = 1$ sans pôle, donc $y$ a au
plus un pôle simple en 0 : $y = q/x$, puis $y = 1/x$ convient
($-1/x^2 + 1/x = (x-1)/x^2$). Résultat $e^x/x$ ; même piège qu'en (d).

**(f) $\int e^x\left(\ln x + \frac1x\right)dx$** (tour $x$, $t_1 = \ln x$, $t_2 = e^x$).
Coefficient de $t_2$ : $y' + y = t_1 + 1/x$ dans $K_1 = C(x, t_1)$ ; $y = t_1$ convient.
Résultat $e^x\ln x$.

**(g) $\int x^2 e^{x^2}$** : non élémentaire (section 4.1). Un intégrateur doit soit le dire,
soit renvoyer $\tfrac12 x e^{x^2} - \tfrac12\int e^{x^2}$ (forme partiellement évaluée).

---

## 7. Pièges et cas de contrôle

**Vérification.** Chaque primitive attendue ci-dessous a été dérivée avec sympy 1.14 (script
`batterie.py` du rapport d'agent C1a) : `simplify(diff(F) - f) == 0` **et** évaluation de
$F' - f$ à 40 chiffres en 4 à 6 points rationnels (écart maximal $< 10^{-168}$). Les
non-élémentarités sont démontrées dans le texte (sections 4 et 6), pas par logiciel.

### 7.1 Primitives élémentaires (le résultat exact attendu)

| # | Intégrande $f$ | Primitive attendue $F$ | Méthode (Risch) | Remarque |
|---|---|---|---|---|
| F1 | $(2x^2+1)e^{x^2}$ | $x e^{x^2}$ | Risch DE $y'+2xy = 2x^2+1$ | $e^{x^2}$ seul n'est pas intégrable |
| F2 | $x e^{x^2}$ | $\tfrac12 e^{x^2}$ | Risch DE, $y$ constant | aussi « dérivée divise » |
| F3 | $\ln x$ | $x\ln x - x$ | polynôme log | |
| F4 | $\ln^2 x$ | $x\ln^2 x - 2x\ln x + 2x$ | polynôme log, 3 étages | |
| F5 | $\dfrac{1}{x\ln x}$ | $\ln(\ln x)$, $x > 1$ ; $\ln\lvert\ln x\rvert$ en général | résidus | |
| F6 | $\dfrac{1}{\ln x} - \dfrac{1}{\ln^2 x}$ | $\dfrac{x}{\ln x}$ | Hermite dans $K(t)$ | **ne pas** séparer les termes |
| F7 | $e^x\dfrac{x-1}{x^2}$ | $\dfrac{e^x}{x}$ | Risch DE $y'+y$ | idem |
| F8 | $\dfrac{e^x}{e^x+1}$ | $\ln(e^x + 1)$ | résidus (dénominateur normal $t+1$) | |
| F9 | $\dfrac{\ln x}{x}$ | $\tfrac12\ln^2 x$ | polynôme log | aussi substitution |
| F10 | $e^x\left(\ln x + \dfrac1x\right)$ | $e^x\ln x$ | Risch DE dans $C(x,\ln x)$ | tour à 2 étages |

### 7.2 Intégrales qui doivent rester non évaluées (jamais de résultat faux)

| Intégrande | Raison | Forme admise |
|---|---|---|
| $e^{x^2}$, $e^{-x^2}$ | §4.1 | non évaluée (ou $\tfrac{\sqrt\pi}{2}\mathrm{erfi}$, $\mathrm{erf}$ si on les ajoute) |
| $x^2 e^{x^2}$ | §4.1, $k$ pair | non évaluée, ou $\tfrac12 x e^{x^2} - \tfrac12\int e^{x^2}$ |
| $\sin(x)/x$, $\cos(x)/x$ | §4.2 | non évaluée ($\mathrm{Si}$, $\mathrm{Ci}$) |
| $\sin(x^2)$, $\cos(x^2)$ | §4.2 | non évaluée (Fresnel) |
| $e^x/x$ | §4.3 | non évaluée ($\mathrm{Ei}$) |
| $1/\ln x$ | §4.3 | non évaluée ($\mathrm{li}$) |
| $e^{e^x}$ | §4.3 | non évaluée |
| $1/\sqrt{x^3+1}$ | cas algébrique (admis) | non évaluée |

### 7.3 Pièges

1. **Linéarité prématurée** (F6, F7) : une somme de termes non intégrables peut avoir une
   primitive. Essayer la somme entière (ou réduire au même dénominateur) avant de couper.
2. **Tour illégitime** : $e^{2x}$ et $e^x$, $\ln(x^2)$ et $\ln x$, $e^{\ln x}$. Normaliser
   avant d'appliquer un critère de non-existence, sinon on « prouve » faussement qu'une
   intégrale élémentaire ne l'est pas ($\int e^{2x}/(e^x + 1)$ est élémentaire :
   $e^x - \ln(e^x + 1)$).
3. **Constantes non testables** : un résidu qui « a l'air » non constant parce qu'il n'est
   pas simplifié ($\sin^2 x + \cos^2 x$ pris pour une fonction de $x$). Le critère des résidus
   n'est une preuve que si le test de constance est fiable.
4. **Domaine réel** : les logarithmes de Liouville sont complexes ; une primitive réelle
   demande $\ln\lvert u\rvert$ et des arctangentes (voir le document sur les fractions
   rationnelles), et une primitive correcte « à une constante près » ne l'est que **par
   intervalle** de définition.
5. **Non élémentaire ≠ non intégrable** : $\int_{-\infty}^{\infty} e^{-x^2} = \sqrt\pi$ ; le
   chemin numérique de `integrerDefinie` doit rester disponible.

---

## 8. État dans le code

Référence : worktree `theorie` (branche `feature/theorie`, base `origin/main` 569e686). Les
numéros de ligne ont été vérifiés en ouvrant les fichiers.

| Élément | État | Où |
|---|---|---|
| Point d'entrée `integrer()` : constante, refus multivariable, puis `primitive()` du nœud | implémenté | `src/Noeud.cpp:265` |
| Primitive non évaluée (pas de résultat faux) ; $(\int f)' = f$ | implémenté | `src/Regles.cpp:308`–`317` ; nœud `include/ASTNode.hpp:612` |
| Table des primitives d'argument affine ($\sin$, $\cos$, $\tan$, $\exp$, $\ln$, réciproques) | implémenté | `src/Regles.cpp:226`–`256`, `:288`–`302` |
| $e^{x^2}$, $\ln(x)/x$ laissées non évaluées (testé) | implémenté | `tests/test_ast.cpp:472`–`476` |
| Corps différentiel / monômes / tour d'extensions | absent | — |
| Théorème de structure (indépendance des $\exp$/$\log$) | absent | la forme canonique fusionne seulement les produits de puissances (`src/Noeud.cpp`) |
| Réduction d'Hermite dans $K(t)$ | absent | — |
| Résultant, critère des résidus | absent (pas de résultant dans `Polynome`) | `include/Polynome.hpp:45` (API : division, PGCD, Yun, Sturm) |
| Équation différentielle de Risch | absent | — |
| Preuve de non-élémentarité | absent (l'échec rend une `IntegraleNonEvaluee`, sans preuve) | `src/Regles.cpp:198` (exemple) |
| Intégrale définie : primitive si certifiée continue, sinon Gauss–Kronrod | implémenté | `src/Integration.cpp:656` (`integraleExacte`), `:695` (`calculer`) |

Conséquence pratique : aujourd'hui, `integrer()` ne sait **rien** dire de plus que « non
évaluée » pour F1, F4–F8, F10 ; ces cas sont des objectifs, pas des régressions.

### Branche en cours (`feature/integration-symbolique`, commit 0f8b5a5)

Lu en lecture seule dans le worktree `primitives-substitution` (aucune compilation) :
`src/Primitives.cpp` organise l'intégration en une suite de méthodes
(`METHODES`, `src/Primitives.cpp:548`) essayées après la règle du nœud : table étendue
(`:335`) et changement de variable « la dérivée divise » (`:483`), avec $\ln\lvert u\rvert$
(`lnAbs`, `:141`). Le commentaire de `:546` annonce intégration par parties, fractions
rationnelles et linéarisation trigonométrique. Rien de Risch (Hermite dans $K(t)$, résidus,
équation de Risch) n'y figure à ce commit. F2 et F9 relèvent de la substitution et devraient
y passer ; F1, F3–F8, F10 demandent des méthodes à venir.

---

## 9. Références

- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, Springer — référence
  principale : corps différentiels, monômes normaux et spéciaux, théorème de Liouville (forme
  forte), réduction d'Hermite généralisée, critère des résidus, équation de Risch, théorème
  de structure.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra*, Kluwer — chapitres
  sur l'intégration (Liouville, Risch transcendant).
- J. von zur Gathen, J. Gerhard, *Modern Computer Algebra*, Cambridge University Press —
  intégration symbolique des fractions rationnelles et cas hyperexponentiel, coûts.
