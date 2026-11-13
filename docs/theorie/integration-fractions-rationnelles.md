# Intégration des fractions rationnelles

Toute fraction rationnelle à coefficients dans un corps $K$ de caractéristique 0 a une
primitive élémentaire : une fraction rationnelle plus une combinaison de logarithmes à
coefficients algébriques. Le problème n'est donc pas l'existence mais le **calcul exact sans
racines** : obtenir la partie rationnelle sans factoriser (Hermite, Horowitz–Ostrogradsky),
la partie logarithmique dans la plus petite extension algébrique possible
(Rothstein–Trager, Lazard–Rioboo–Trager), puis une forme **réelle et continue** (arctangentes
de Rioboo). C'est aussi le cas de base de l'algorithme de Risch
([`integration-fondements.md`](integration-fondements.md)).

## Table des matières

1. [Notations](#notations)
2. [Forme du résultat et méthode naïve](#1-forme-du-résultat-et-méthode-naïve)
3. [Préliminaires : division, Bézout, décomposition sans carré (Yun)](#2-préliminaires--division-bézout-décomposition-sans-carré-yun)
4. [Partie rationnelle : réduction d'Hermite](#3-partie-rationnelle--réduction-dhermite)
5. [Partie rationnelle : Horowitz–Ostrogradsky](#4-partie-rationnelle--horowitzostrogradsky)
6. [Partie logarithmique : Rothstein–Trager](#5-partie-logarithmique--rothsteintrager)
7. [Partie logarithmique : Lazard–Rioboo–Trager](#6-partie-logarithmique--lazardriobootrager)
8. [Cas réel : logarithmes et arctangentes continues](#7-cas-réel--logarithmes-et-arctangentes-continues)
9. [Exemples travaillés complets](#8-exemples-travaillés-complets)
10. [Complexité et choix d'implémentation](#9-complexité-et-choix-dimplémentation)
11. [Pièges et cas de contrôle](#10-pièges-et-cas-de-contrôle)
12. [État dans le code](#11-état-dans-le-code)
13. [Références](#12-références)

## Notations

| Symbole | Sens |
|---|---|
| $K$ | corps de caractéristique 0 (en pratique $\mathbb{Q}$, ou $\mathbb{Q}(\text{paramètres})$) ; $\overline{K}$ sa clôture algébrique |
| $A, D \in K[x]$ | numérateur et dénominateur, $D \ne 0$ ; après division, $\deg A < \deg D$ |
| $D'$ | dérivée de $D$ par rapport à $x$ |
| $\mathrm{pp}$, $\mathrm{lc}$ | partie primitive, coefficient dominant |
| $D = D_1 D_2^2 \cdots D_m^m$ | décomposition sans carré : $D_i$ sans carré, deux à deux premiers entre eux |
| $D^- = \gcd(D, D')$, $D^* = D/D^-$ | $D^* = D_1 D_2 \cdots D_m$ est la partie sans carré (« radical ») de $D$ |
| $\mathrm{Res}_x(P, Q)$ | résultant par rapport à $x$ |
| $S_j$ | sous-résultant de degré (nominal) $j$ de deux polynômes |
| $\sum_{Q(a)=0} a\,\log S(a, x)$ | somme sur les racines de $Q$ (« RootSum ») |
| $n$ | $\deg D$ |

Primitives définies à une constante près, **par intervalle** de définition sur $\mathbb{R}$.

---

## 1. Forme du résultat et méthode naïve

**Théorème.** Soient $A, D \in K[x]$, $\deg A < \deg D$. Il existe $B, C \in K[x]$, des
constantes $c_1, \dots, c_k \in \overline{K}$ et des $v_i \in \overline{K}[x]$ tels que
$$
\int \frac{A}{D} = \frac{B}{C} + \sum_{i=1}^{k} c_i \log v_i .
$$

*Preuve (par décomposition en éléments simples sur $\overline{K}$).* Avec
$D = \mathrm{lc}(D)\prod_\beta (x - \beta)^{e_\beta}$,
$$
\frac{A}{D} = \sum_\beta \sum_{j=1}^{e_\beta} \frac{a_{\beta, j}}{(x - \beta)^j},\qquad
\int \frac{a_{\beta,j}}{(x-\beta)^j} =
\begin{cases} \dfrac{-a_{\beta,j}}{(j-1)(x-\beta)^{j-1}} & j \ge 2,\\[2mm]
a_{\beta,1}\log(x - \beta) & j = 1. \end{cases}
$$
La somme des termes $j \ge 2$ est une fraction rationnelle invariante par les automorphismes
de $\overline{K}/K$, donc à coefficients dans $K$. $\square$

**Méthode naïve et pourquoi l'éviter.** Factoriser $D$ sur $\overline{K}$ exige les racines :
pour $\deg D = n$, le corps de décomposition peut être de degré $n!$ sur $K$, et la partie
rationnelle — qui est dans $K(x)$ — serait calculée en passant inutilement par ces nombres
algébriques. Les algorithmes ci-dessous ne calculent **aucune racine** pour la partie
rationnelle et **seulement les racines nécessaires** (les résidus distincts) pour la partie
logarithmique.

**Plan de calcul.**

```
IntegrerFractionRationnelle(A/D) :
    (P, A) ← division euclidienne : A/D = P + A/D, deg A < deg D     // P polynôme : ∫P immédiat
    (g, h) ← Hermite(A, D)        // A/D = g' + h, h = a/b avec b sans carré
    L      ← PartieLog(h)         // Rothstein–Trager ou Lazard–Rioboo–Trager
    si réel demandé : L ← LogVersReel(L)   // ln|.|, arctangentes de Rioboo
    renvoyer ∫P + g + L
```

---

## 2. Préliminaires : division, Bézout, décomposition sans carré (Yun)

**Bézout étendu.** Si $\gcd(a, b) = 1$, pour tout $c$ il existe un unique couple $(s, t)$ avec
$s\,a + t\,b = c$ et $\deg s < \deg b$ : on calcule $s_0, t_0$ avec $s_0 a + t_0 b = 1$ (algorithme
d'Euclide étendu), puis $s = (s_0 c) \bmod b$, $t = (c - s a)/b$ (division exacte).

**Définition.** $D$ est *sans carré* s'il n'a pas de racine multiple dans $\overline{K}$, ce qui
équivaut (caractéristique 0) à $\gcd(D, D') = 1$. La *décomposition sans carré* est
$D = \mathrm{lc}(D)\,D_1 D_2^2\cdots D_m^m$, $D_i$ unitaires sans carré deux à deux premiers
entre eux ($D_i$ = produit des $x - \beta$ de multiplicité exactement $i$).

**Algorithme de Yun.**

```
Yun(f unitaire) :
    a0 ← gcd(f, f') ;  b ← f / a0 ;  c ← f' / a0 ;  d ← c − b'
    i ← 1 ;  résultat ← []
    tant que deg b > 0 :
        a ← gcd(b, d)                   // a = f_i
        si deg a > 0 : ajouter (a, i)
        b ← b / a ;  c ← d / a ;  d ← c − b' ;  i ← i + 1
    renvoyer résultat
```

*Justification.* Si $f = \prod f_k^k$, alors $a_0 = \prod f_k^{k-1}$, $b = \prod_k f_k$ et
$c = f'/a_0 = \sum_k k f_k' \prod_{l \ne k} f_l$. Donc $d = c - b' = \sum_k (k-1) f_k'\prod_{l\ne k}
f_l$ : le terme $k = 1$ est nul. Par récurrence, à l'étape $i$, $b = \prod_{k \ge i} f_k$ et
$d = \sum_{k \ge i}(k - i) f_k' \prod_{l \ge i, l \ne k} f_l$, de sorte que
$\gcd(b, d) = f_i$ (les $f_k$, $k > i$, ne divisent pas $d$ car $\gcd(f_k, f_k') = 1$ et
$k - i \ne 0$ en caractéristique 0). $\square$

**Coût.** Des PGCD de polynômes de degrés décroissants dont la somme des degrés des $b$
successifs vaut $\sum_k k \deg f_k = n$ : $O(n^2)$ opérations dans $K$ en arithmétique
classique (estimation directe, sans compter la croissance des coefficients dans
$\mathbb{Q}$ ; von zur Gathen–Gerhard donnent les coûts avec arithmétique rapide).

**Exemple (trace).** $D = x^8 + 6x^6 + 12x^4 + 8x^2 = x^2(x^2+2)^3$.
$a_0 = x(x^2+2)^2$, $b = x(x^2+2)$, $d = 5x^2 + 2$ ;
$i = 1$ : $\gcd(b, d) = 1$ ; puis $b = x(x^2+2)$, $d = 2x^2$ ;
$i = 2$ : $\gcd = x$ ; puis $b = x^2 + 2$, $d = 0$ ; $i = 3$ : $\gcd = x^2 + 2$. Résultat
$D_1 = 1$, $D_2 = x$, $D_3 = x^2 + 2$. (Trace obtenue avec l'implémentation de référence
`algos.py`, voir section 10.)

---

## 3. Partie rationnelle : réduction d'Hermite

### 3.1 Énoncé

**Théorème (Hermite ; forme de Ostrogradsky).** Soient $\deg A < \deg D$, $D^- = \gcd(D, D')$,
$D^* = D / D^-$. Il existe un unique couple $(B, C)$ avec $\deg B < \deg D^-$,
$\deg C < \deg D^*$ et
$$
\frac{A}{D} = \left(\frac{B}{D^-}\right)' + \frac{C}{D^*}.
$$
$B/D^-$ est la **partie rationnelle** de la primitive ; $C/D^*$ a un dénominateur **sans
carré**, sa primitive est purement logarithmique.

*Unicité.* Une différence de deux décompositions donne $(B_1 - B_2)/D^-$ dont la dérivée vaut
$(C_2 - C_1)/D^*$, fraction à pôles simples. Or la dérivée d'une fraction rationnelle n'a
jamais de pôle simple (un pôle d'ordre $k \ge 1$ devient d'ordre $k + 1$) : donc
$(B_1 - B_2)/D^-$ n'a pas de pôle, c'est un polynôme ; de degré $< 0$ puisque
$\deg(B_1 - B_2) < \deg D^-$ : il est nul, puis $C_1 = C_2$. L'existence découle de
l'algorithme ci-dessous.

### 3.2 Une étape de réduction

Soit $D = U V^{j+1}$ avec $V$ sans carré, $\gcd(U, V) = 1$, $j \ge 1$. Comme $V$ est sans carré
et premier avec $U$, $\gcd(U V', V) = 1$ : par Bézout, on trouve $B, C$ avec $\deg B < \deg V$ et
$$
B\,U V' + C\,V = -\frac{A}{j}.
$$
Alors
$$
\frac{A}{U V^{j+1}} = \left(\frac{B}{V^j}\right)' + \frac{-jC - U B'}{U V^{j}} .
$$
*Vérification.* $\left(\frac{B}{V^j}\right)' = \frac{B'}{V^j} - \frac{jBV'}{V^{j+1}}$, donc
$\frac{A}{UV^{j+1}} - \left(\frac{B}{V^j}\right)' = \frac{A + jUV'B - UVB'}{UV^{j+1}}$ ; or
$A + jUV'B = A + j\left(-\frac{A}{j} - CV\right) = -jCV$, d'où le numérateur
$V(-jC - UB')$. $\square$

L'exposant de $V$ a baissé de 1 sans aucune factorisation.

### 3.3 Algorithme (version linéaire, telle que présentée par Bronstein)

```
Hermite(A, D) :                       // deg A < deg D, D unitaire
    g ← 0
    (D_1, ..., D_m) ← Yun(D)
    pour i de 2 à m, si deg D_i > 0 :
        V ← D_i ;  U ← D / V^i
        (s0, t0) ← EuclideEtendu(U V', V)          // s0 U V' + t0 V = 1, une fois par i
        pour j de i−1 à 1 :
            B ← (s0 · (−A/j)) mod V ;  C ← (−A/j − B U V') / V      // division exacte
            g ← g + B / V^j
            A ← −j C − U B'
        D ← U V
    renvoyer (g, A / D)                 // A/D initial = g' + A/D final, D final sans carré
```

**Coût.** Une paire de Bézout par $i$ (coût $O(n\deg D_i)$), puis pour chaque $j$ une
multiplication modulo $V$ et une division exacte ($O(n\deg D_i)$) : au total
$O\!\left(n\sum_i i\deg D_i\right) = O(n^2)$ opérations dans $K$ en arithmétique classique,
plus le coût de Yun (estimation directe ; la croissance des rationnels n'est pas comptée).

### 3.4 Exemples

**(a)** $\int \frac{dx}{(x^2+1)^2}$. $V = x^2 + 1$, $U = 1$, $i = 2$, $j = 1$, $A = 1$ : on
résout $B\cdot 2x + C(x^2+1) = -1$ avec $\deg B < 2$. Bézout : $(-x/2)\cdot 2x + 1\cdot(x^2+1)
= 1$, d'où $B = (x/2) \bmod V = x/2$, puis par division exacte
$C = (-1 - B\cdot 2x)/(x^2+1) = (-1 - x^2)/(x^2+1) = -1$. Nouveau numérateur
$-jC - UB' = 1 - \tfrac12 = \tfrac12$. Résultat :
$$
\frac{1}{(x^2+1)^2} = \left(\frac{x}{2(x^2+1)}\right)' + \frac{1}{2(x^2+1)},\qquad
\int \frac{dx}{(x^2+1)^2} = \frac{x}{2(x^2+1)} + \frac12\arctan x .
$$
(La paire $(B, C)$ est déterminée par la contrainte $\deg B < \deg V$ ; calculer $C$ par
division exacte plutôt qu'« à vue » évite les erreurs de signe et de facteur.)

**(b)** $\displaystyle\int \frac{x^7 - 24x^4 - 4x^2 + 8x - 8}{x^8 + 6x^6 + 12x^4 + 8x^2}\,dx$,
$D = x^2(x^2+2)^3$ (trace exacte de l'implémentation de référence) :

| étape | $V$ | $U$ | $B$ | nouveau numérateur $A$ |
|---|---|---|---|---|
| $i=2$, $j=1$ | $x$ | $(x^2+2)^3$ | $1$ | $x^6 + x^5 - 18x^3 + 8x + 8$ |
| $i=3$, $j=2$ | $x^2+2$ | $x$ | $6x$ | $x^4 + x^3 - 2x^2 - 2x + 4$ |
| $i=3$, $j=1$ | $x^2+2$ | $x$ | $3 - x$ | $x^2 + 2$ |

Partie rationnelle $g = \frac1x + \frac{6x}{(x^2+2)^2} + \frac{3-x}{x^2+2}
= \frac{3x^3 + 8x^2 + 6x + 4}{x(x^2+2)^2}$, reste $\frac{x^2+2}{x(x^2+2)} = \frac1x$. D'où
$$
\int = \frac{3x^3 + 8x^2 + 6x + 4}{x(x^2+2)^2} + \ln\lvert x\rvert .
$$

---

## 4. Partie rationnelle : Horowitz–Ostrogradsky

On écrit directement $A/D = (B/D^-)' + C/D^*$ avec $B$, $C$ à coefficients **inconnus**
($\deg B < \deg D^-$, $\deg C < \deg D^*$, soit $\deg D^- + \deg D^* = n$ inconnues). En
multipliant par $D$ :
$$
A = B'\,D^* - B\,H + C\,D^-, \qquad H = \frac{D^*\,(D^-)'}{D^-} \in K[x]
$$
($H$ est un polynôme car $D^- \mid D^* (D^-)'$ : chaque facteur $D_i^{i-1}$ de $D^-$ a sa
dérivée divisible par $D_i^{i-2}$ et $D^*$ contient $D_i$). On identifie les coefficients :
**système linéaire $n \times n$**, de solution unique (théorème 3.1).

```
Horowitz(A, D) :
    D⁻ ← gcd(D, D') ;  D* ← D / D⁻ ;  H ← D* (D⁻)' / D⁻
    B ← Σ_{k < deg D⁻} b_k x^k ;  C ← Σ_{k < deg D*} c_k x^k        // inconnues
    résoudre les n équations « coefficients de B' D* − B H + C D⁻ − A = 0 »
    renvoyer (B / D⁻, C / D*)
```

**Exemple.** Pour (b) ci-dessus : $D^- = x(x^2+2)^2$, $D^* = x(x^2+2)$, 8 inconnues ; la
solution est $B/D^- = \frac{3x^3 + 8x^2 + 6x + 4}{x(x^2+2)^2}$, $C/D^* = \frac{x^2+2}{x(x^2+2)}$ :
identique à Hermite (unicité).

**Comparaison.**

| | Hermite (linéaire) | Horowitz–Ostrogradsky |
|---|---|---|
| Préalable | Yun (décomposition sans carré complète) | un seul PGCD $\gcd(D, D')$ |
| Calcul | Bézout + divisions exactes | un système linéaire $n\times n$ |
| Coût classique | $O(n^2)$ opérations dans $K$ (+ Yun) | $O(n^3)$ par élimination de Gauss |
| Structure | itératif, résultats intermédiaires lisibles (journal d'étapes) | « tout d'un coup », simple à coder si l'on dispose d'un solveur linéaire exact |
| Paramètres symboliques | divisions dans $K(\text{paramètres})$ | idem, mais pivots à tester non nuls |

Pour SymAlgo++, Hermite est le choix naturel : `Polynome` a déjà la division et le PGCD
exacts et Yun (`sansCarre`) ; il manque seulement l'Euclide **étendu** (section 11).

---

## 5. Partie logarithmique : Rothstein–Trager

On suppose maintenant $D$ **sans carré**, $\deg A < \deg D$ et $\gcd(A, D) = 1$.

### 5.1 Théorème

**Théorème (Rothstein, Trager).** Soit
$$
R(z) = \mathrm{Res}_x\big(D,\ A - z\,D'\big) \in K[z].
$$
Alors
$$
\int \frac{A}{D} = \sum_{\substack{\alpha \in \overline{K}\\ R(\alpha) = 0}} \alpha\,\log\gcd\big(D,\ A - \alpha D'\big),
$$
et les racines de $R$ sont exactement les **résidus** de $A/D$. Le corps $K(\alpha_1, \dots,
\alpha_k)$ engendré par les racines de $R$ est la plus petite extension de $K$ sur laquelle la
primitive s'écrit $\sum c_i\log v_i$ avec des $v_i$ sans carré deux à deux premiers entre eux
et des $c_i$ distincts.

*Preuve.* Soit $D = \mathrm{lc}(D)\prod_\beta (x - \beta)$ ($n$ racines distinctes). Comme $D$
est sans carré, $D'(\beta) \ne 0$ et
$$
\frac{A}{D} = \sum_\beta \frac{\rho_\beta}{x - \beta}, \qquad \rho_\beta = \frac{A(\beta)}{D'(\beta)} .
$$
Par la formule du résultant $\mathrm{Res}_x(D, G) = \mathrm{lc}(D)^{\deg G}\prod_\beta G(\beta)$,
$$
R(z) = \kappa \prod_\beta \big(A(\beta) - z D'(\beta)\big) = \kappa' \prod_\beta (\rho_\beta - z),
\qquad \kappa' \ne 0 :
$$
les racines de $R$ sont les résidus, la multiplicité de $\alpha$ est le nombre de $\beta$ de
résidu $\alpha$. Enfin, une racine $\beta$ de $D$ annule $A - \alpha D'$ si et seulement si
$\rho_\beta = \alpha$, donc $\gcd(D, A - \alpha D') = \prod_{\rho_\beta = \alpha}(x - \beta)$, et
$\sum_\alpha \alpha\log\prod_{\rho_\beta=\alpha}(x-\beta) = \sum_\beta \rho_\beta\log(x-\beta)$.
Minimalité : dans toute écriture $\sum c_i \log v_i$ du type indiqué, le résidu en $\beta$ est
le $c_i$ du $v_i$ qui s'annule en $\beta$ ; les $c_i$ sont donc des résidus. $\square$

### 5.2 Algorithme

```
RothsteinTrager(A, D) :              // D sans carré, deg A < deg D
    R ← Res_x(D, A − z D')            // polynôme en z de degré n
    pour chaque facteur irréductible r(z) de R sur K :
        v(z, x) ← gcd(D, A − z D')  calculé dans (K[z]/(r))[x]
        ajouter  Σ_{r(α)=0} α log v(α, x)
```

**Coût et inconvénients.** Le résultant se calcule par sous-résultants (ou par
évaluation–interpolation en $z$). Il faut ensuite **factoriser $R$ sur $K$** (coûteux) et
calculer des PGCD dans des extensions algébriques $K[z]/(r)$ — c'est ce que
Lazard–Rioboo–Trager supprime.

### 5.3 Exemples

| $A/D$ | $R(z)$ (à une constante près) | résidus | primitive |
|---|---|---|---|
| $\frac{1}{x^2+1}$ | $4z^2 + 1$ | $\pm\frac{i}{2}$ | $\frac{i}{2}\log(x+i) - \frac{i}{2}\log(x-i)$ |
| $\frac{2x}{x^2-1}$ | $(z-1)^2$ | $1$ (double) | $\log(x^2 - 1)$ : **une seule** racine, pas d'extension |
| $\frac{1}{x^2-2}$ | $8z^2 - 1$ | $\pm\frac{\sqrt2}{4}$ | $\frac{\sqrt2}{4}\log\frac{x-\sqrt2}{x+\sqrt2}$ : $\sqrt2$ **nécessaire** |
| $\frac{1}{x^3+1}$ | $(3z-1)(9z^2+3z+1)$ | $\frac13$, racines de $9z^2+3z+1$ | $\frac13\log(x+1) + \sum_{9\alpha^2+3\alpha+1=0}\alpha\log(x + 3\alpha)$ |
| $\frac{1}{x^4+1}$ | $256z^4 + 1$ | 4 racines | $\sum_{256\alpha^4+1=0}\alpha\log(x - 4\alpha)$ |

Détail pour $\frac{1}{x^2+1}$ : pour $\alpha = i/2$, $\gcd(x^2+1, 1 - ix) = x + i$ (car
$1 - ix = -i(x + i)$). La forme $\frac{i}{2}\log\frac{x+i}{x-i}$ est réelle à une constante
près : $\arctan x - \frac{\pi}{2}$ (section 7).

Détail pour $\frac{1}{x^3+1}$ : $\rho_\beta = 1/(3\beta^2) = -\beta/3$ (car $\beta^3 = -1$), donc
$\beta = -3\alpha$ et $\gcd(D, A - \alpha D') = x + 3\alpha$ ; le facteur $27z^3 - 1$ du
résultant (calculé : $R = -(3z-1)(9z^2+3z+1)$) redonne bien $\alpha = -\beta/3$.

---

## 6. Partie logarithmique : Lazard–Rioboo–Trager

### 6.1 Idée

La multiplicité d'une racine $\alpha$ de $R$ est $\deg\gcd(D, A - \alpha D')$ (preuve de 5.1).
La **propriété de spécialisation des sous-résultants** dit que, si $\gcd(D, A - \alpha D')$ est
de degré $i$ (et sous une condition de non-annulation des coefficients dominants), le
sous-résultant $S_i(z, x)$ de $D$ et $A - zD'$ (calculé **une fois**, avec $z$ formel) vérifie
$S_i(\alpha, x) \propto \gcd(D, A - \alpha D')$. Il suffit donc de regrouper les racines de $R$
par multiplicité — une **décomposition sans carré** de $R$, pas une factorisation — et de
lire les PGCD dans la suite des sous-résultants : plus de calcul dans $K[z]/(r)$, plus de
factorisation obligatoire.

### 6.2 Algorithme (tel que présenté par Bronstein)

```
LazardRiobooTrager(A, D) :            // D sans carré, deg A < deg D
    (R, S_0, ..., S_{n-1}) ← suite des sous-résultants de D et A − z D' par rapport à x
                              // R = S_0 ∈ K[z] ; S_j de degré nominal j en x
    (Q_1, ..., Q_k) ← Yun(R)        // R = c · Π Q_i^i, décomposition sans carré en z
    pour i de 1 à k, si deg_z Q_i > 0 :
        si i = deg D : S ← D
        sinon :        S ← S_i
            // retirer de S le contenu en z qui s'annulerait sur les racines de Q_i :
            (A_1, ..., A_q) ← Yun(lc_x(S))
            pour j de 1 à q : S ← S / gcd(A_j, Q_i)^j
        ajouter  Σ_{Q_i(α)=0} α log S(α, x)          // S réduit modulo Q_i si l'on veut
```

On peut ensuite (facultatif) factoriser chaque $Q_i$ pour simplifier l'affichage ou pour la
conversion réelle ; si $\deg Q_i \le 2$ les racines s'écrivent par radicaux. L'énoncé exact de
la normalisation du contenu est chez Bronstein ; l'implémentation de référence du rapport
(`algos.py`, sous-résultants calculés par la définition déterminantale) donne les résultats
de la section 8 sur tous les exemples de ce document.

### 6.3 Exemple : la multiplicité 3

$\int \frac{x^4 - 3x^2 + 6}{x^6 - 5x^4 + 5x^2 + 4}$. Le dénominateur est irréductible sur
$\mathbb{Q}$, sans carré ; $\deg A < \deg D$.

- $R(z) = 45796\,(4z^2 + 1)^3$ : une seule classe, $Q_3 = 4z^2 + 1$, $i = 3$.
- Le sous-résultant de degré 3 vaut
  $800x^3z^3 - 14x^3z - 400x^2z^2 + 7x^2 - 2440xz^3 + 32xz + 792z^2 - 16$ ; réduit modulo
  $4z^2 + 1$ il vaut $107(-2z\,x^3 + x^2 + 6z\,x - 2)$, et rendu unitaire en $x$ (inverse de
  $-2z$ modulo $4z^2+1$ : $2z$) :
  $$S(z, x) = x^3 + 2z\,x^2 - 3x - 4z .$$
- Primitive (complexe) : $\displaystyle\sum_{4\alpha^2+1=0}\alpha\log\big(x^3 + 2\alpha x^2 - 3x - 4\alpha\big)$.
  Pour $\alpha = i/2$ : $x^3 + i x^2 - 3x - 2i$, qui est bien un facteur de $D$ sur
  $\mathbb{Q}(i)$ (vérifié : $D = (x^3 - ix^2 - 3x + 2i)(x^3 + ix^2 - 3x - 2i)$).

Rothstein–Trager aurait demandé le PGCD de $D$ et $A - \alpha D'$ dans $\mathbb{Q}(i)[x]$ ;
Lazard–Rioboo–Trager le lit dans un sous-résultant calculé sur $\mathbb{Q}[z]$.

---

## 7. Cas réel : logarithmes et arctangentes continues

### 7.1 Des paires de logarithmes conjugués aux arctangentes

Si $K \subset \mathbb{R}$, les racines non réelles de $Q_i$ vont par paires conjuguées
$\alpha = u \pm iv$ ($v > 0$), et $S(\alpha, x) = P + iQ$, $S(\bar\alpha, x) = P - iQ$ avec
$P, Q \in \mathbb{R}[x]$. Alors
$$
(u + iv)\log(P + iQ) + (u - iv)\log(P - iQ) = u\log(P^2 + Q^2) + v\cdot i\log\frac{P + iQ}{P - iQ}.
$$
**Lemme.** $\dfrac{d}{dx}\left[i\log\dfrac{P + iQ}{P - iQ}\right] = \dfrac{2(P'Q - PQ')}{P^2 + Q^2}
= \dfrac{d}{dx}\left[2\arctan\dfrac{P}{Q}\right]$.

*Preuve.* $\frac{d}{dx}\log\frac{P+iQ}{P-iQ} = \frac{P'+iQ'}{P+iQ} - \frac{P'-iQ'}{P-iQ}
= \frac{2i(PQ' - P'Q)}{P^2 + Q^2}$ ; multiplier par $i$. Et
$\frac{d}{dx}\arctan\frac PQ = \frac{(P'Q - PQ')/Q^2}{1 + P^2/Q^2} = \frac{P'Q - PQ'}{P^2+Q^2}$.
$\square$

D'où la contribution réelle $u\log(P^2 + Q^2) + 2v\arctan(P/Q)$.

### 7.2 Le piège : $\arctan(P/Q)$ est discontinue

$\arctan(P/Q)$ saute de $\pm\pi$ en chaque zéro réel de $Q$ où $P \ne 0$. La formule a la
bonne dérivée **sur chaque intervalle** entre deux zéros de $Q$, mais n'est pas une primitive
sur $\mathbb{R}$ : l'utiliser pour une intégrale définie qui enjambe un zéro de $Q$ donne une
valeur **fausse** (erreur d'un multiple de $2v\pi$).

Exemple : pour $\int \frac{x^4 - 3x^2 + 6}{x^6 - 5x^4 + 5x^2 + 4}$, la conversion directe
($\alpha = i/2$, $u = 0$, $v = 1/2$, $P = x^3 - 3x$, $Q = x^2 - 2$) donne
$$
F_{\text{naïf}} = \arctan\frac{x^3 - 3x}{x^2 - 2},
$$
dont la dérivée est bien l'intégrande (vérifié) mais qui saute de $\pi/2$ à $-\pi/2$ en
$x = \sqrt2$ (limites calculées : $\pi/2$ à gauche, $-\pi/2$ à droite) et de même en
$-\sqrt2$, alors que l'intégrande est continue sur $\mathbb{R}$.

### 7.3 Algorithme de Rioboo : arctangentes de polynômes

**Principe.** Remplacer $i\log\frac{A+iB}{A-iB}$ par une somme d'arctangentes **de
polynômes**, continues sur $\mathbb{R}$.

```
LogVersAtan(A, B) :      // A, B ∈ R[x], non tous deux nuls ; renvoie f avec
                          // f' = d/dx [ i log((A+iB)/(A−iB)) ], f somme d'arctan de polynômes
    si B divise A : renvoyer 2 arctan(A / B)
    si deg A < deg B : renvoyer LogVersAtan(−B, A)
    (D, C) ← Bézout : B D − A C = G = gcd(A, B)
    renvoyer 2 arctan((A D + B C) / G) + LogVersAtan(D, C)
```

*Justification.* $(A + iB)(D - iC) = (AD + BC) + i(BD - AC) = (AD + BC) + iG$, donc
$$
\frac{A + iB}{A - iB} = \frac{(AD + BC) + iG}{(AD + BC) - iG}\cdot\frac{D + iC}{D - iC},
$$
et le premier facteur donne $2\arctan\frac{AD+BC}{G}$ par le lemme 7.1, avec $G \mid AD + BC$
(donc argument polynomial). Le cas $\deg A < \deg B$ utilise $-B + iA = i(A + iB)$ et
$-B - iA = -i(A - iB)$ : $\frac{-B + iA}{-B - iA} = -\frac{A+iB}{A-iB}$, et les deux
logarithmes ne diffèrent que d'une constante. La récursion termine car les coefficients de Bézout sont de degrés strictement plus
petits. Toutes les arctangentes ont des arguments **polynomiaux** : la somme est continue sur
$\mathbb{R}$. $\square$

**Application à l'exemple.** $A = x^3 - 3x$, $B = x^2 - 2$ :

1. $B \nmid A$, $\deg A \ge \deg B$ ; Bézout $B D - A C = 1$ avec $D = \frac{x^2 - 1}{2}$,
   $C = \frac{x}{2}$ ; $AD + BC = \frac{x^5 - 3x^3 + x}{2}$ : terme
   $2\arctan\frac{x^5 - 3x^3 + x}{2}$.
2. $\mathrm{LogVersAtan}\big(\frac{x^2-1}{2}, \frac x2\big)$ : Bézout $D = 2x$, $C = 2$,
   $AD + BC = x^3$ : terme $2\arctan(x^3)$.
3. $\mathrm{LogVersAtan}(2x, 2)$ : $2 \mid 2x$ : terme $2\arctan x$.

Avec le facteur $v = 1/2$ :
$$
\boxed{\int \frac{x^4 - 3x^2 + 6}{x^6 - 5x^4 + 5x^2 + 4}\,dx = \arctan\frac{x^5 - 3x^3 + x}{2} + \arctan(x^3) + \arctan x}
$$
continue sur $\mathbb{R}$ (vérifié par dérivation ; c'est aussi la forme que renvoie sympy).

### 7.4 Conversion réelle complète (« LogVersReel »)

```
LogVersReel(Σ_{Q(α)=0} α log S(α, x)) :          // Q ∈ R[z], S ∈ R[z, x]
    écrire Q(u + iv) = P_Q(u, v) + i Q_Q(u, v) et S(u + iv, x) = A(u, v, x) + i B(u, v, x)
    pour chaque solution réelle (u, v) de P_Q = Q_Q = 0 avec v > 0 :
        ajouter  u log(A² + B²) + v · LogVersAtan(A, B)
    pour chaque racine réelle u de Q :
        ajouter  u log|S(u, x)|
```

Les solutions $(u, v)$ s'obtiennent en résolvant un système polynomial à deux inconnues
(résultant en $v$, puis racines réelles certifiées) ; pour $\deg Q \le 2$, elles s'écrivent
directement par radicaux, ce qui couvre la plupart des cas pratiques.

### 7.5 Autre exemple de discontinuité

$\int \frac{x^2+1}{x^4+x^2+1}$ : $R = (12z^2+1)^2$, $Q_2 = 12z^2 + 1$, $S = x^2 + 6z\,x - 1$
(unitaire). Pour $\alpha = \frac{i}{2\sqrt3}$ : $P = x^2 - 1$, $Q = \sqrt3\,x$, $v = \frac{1}{2\sqrt3}$.

- Forme naïve : $\frac{1}{\sqrt3}\arctan\frac{x^2-1}{\sqrt3\,x}$, de saut $-\pi/\sqrt3$ en $x = 0$
  (calculé).
- Rioboo ($\mathrm{LogVersAtan}(x^2 - 1, \sqrt3 x)$, calculé sur $\mathbb{Q}(\sqrt3)$) :
  $$\int \frac{x^2+1}{x^4+x^2+1} = \frac{1}{\sqrt3}\left(\arctan\frac{x^3+2x}{\sqrt3} + \arctan\frac{x}{\sqrt3}\right),$$
  égale (écart numérique nul) à la forme obtenue par éléments simples réels
  $\frac{1}{\sqrt3}\left(\arctan\frac{2x+1}{\sqrt3} + \arctan\frac{2x-1}{\sqrt3}\right)$.

---

## 8. Exemples travaillés complets

**(1)** $\displaystyle\int\frac{dx}{x^2+1} = \arctan x$. RT : $R = 4z^2 + 1$ ; $\alpha = i/2$,
$S = x + i$, $u = 0$, $v = 1/2$, $P = x$, $Q = 1$ : $v\cdot 2\arctan(x/1) = \arctan x$.

**(2)** $\displaystyle\int\frac{dx}{x^3+1} = \frac13\ln\lvert x+1\rvert - \frac16\ln(x^2 - x + 1)
+ \frac{1}{\sqrt3}\arctan\frac{2x-1}{\sqrt3}$.
LRT : $R \propto 27z^3 - 1$ (sans carré, $i = 1$), $S = x + 3z$. Racine réelle $u = 1/3$ :
$\frac13\ln\lvert x+1\rvert$. Racines complexes de $9z^2 + 3z + 1$ : $(u, v) = (-\frac16, \frac{\sqrt3}{6})$ ;
$S = (x - \frac12) + i\frac{\sqrt3}{2}$ : $u\log(A^2+B^2) = -\frac16\ln(x^2 - x + 1)$ et
$v\,\mathrm{LogVersAtan}(x - \frac12, \frac{\sqrt3}{2}) = \frac{\sqrt3}{6}\cdot 2\arctan\frac{2x-1}{\sqrt3}$.

**(3)** $\displaystyle\int\frac{x^4 - 3x^2 + 6}{x^6 - 5x^4 + 5x^2 + 4} = \arctan\frac{x^5-3x^3+x}{2}
+ \arctan(x^3) + \arctan x$ (sections 6.3 et 7.3).

**(4)** $\displaystyle\int\frac{dx}{(x^2+1)^2} = \frac{x}{2(x^2+1)} + \frac12\arctan x$ (Hermite 3.4 a,
puis (1)).

**(5)** $\displaystyle\int\frac{x^7 - 24x^4 - 4x^2 + 8x - 8}{x^8 + 6x^6 + 12x^4 + 8x^2} =
\frac{3x^3 + 8x^2 + 6x + 4}{x(x^2+2)^2} + \ln\lvert x\rvert$ (Hermite 3.4 b ; reste $1/x$, RT :
$R = z - 1$).

**(6)** $\displaystyle\int\frac{dx}{x^4+1} = \frac{\sqrt2}{8}\ln\frac{x^2+\sqrt2x+1}{x^2-\sqrt2x+1}
+ \frac{\sqrt2}{4}\big(\arctan(\sqrt2x+1) + \arctan(\sqrt2x-1)\big)$ ; la forme
$\frac{\sqrt2}{4}\arctan\frac{\sqrt2x}{1-x^2}$ pour la partie arctangente est discontinue en
$x = \pm1$ (saut de $\mp\frac{\sqrt2\pi}{4}$, calculé).

---

## 9. Complexité et choix d'implémentation

| Étape | Coût (opérations dans $K$, arithmétique classique) | Remarque |
|---|---|---|
| Division, partie polynomiale | $O(n^2)$ | |
| Yun | $O(n^2)$ | estimation directe (§2) |
| Hermite linéaire | $O(n^2)$ | estimation directe (§3.3) |
| Horowitz–Ostrogradsky | $O(n^3)$ | système linéaire dense |
| Résultant $R(z)$ / sous-résultants | polynomial ; $O(n^3)$ par évaluation–interpolation en $z$ ($n+1$ résultants de coût $O(n^2)$) | estimation personnelle, non tirée d'une référence |
| Rothstein–Trager | + factorisation de $R$ sur $K$ + PGCD dans $K[z]/(r)$ | le plus cher en pratique |
| Lazard–Rioboo–Trager | + Yun sur $R$ ; aucune extension algébrique | à préférer |
| Rioboo (arctangentes) | au plus $\deg A$ appels récursifs, chacun un Euclide étendu | arguments de degré $\le \deg A + \deg B$ |

Sur $\mathbb{Q}$, la **taille des coefficients** domine souvent le nombre d'opérations
(résultants : le coefficient 45796 de l'exemple 6.3 est typique) ; les implémentations
sérieuses utilisent des méthodes modulaires. Ces bornes ne prennent pas en compte cette
croissance.

**Ce qu'il faut ajouter à SymAlgo++** (dans l'ordre) : Euclide étendu sur `Polynome` ;
Hermite ; polynômes à deux variables (ou à coefficients dans $\mathbb{Q}[z]$) pour les
sous-résultants ; LRT ; un moyen de représenter $\sum_{Q(\alpha)=0}$ quand $\deg Q > 2$ (nœud
« somme sur les racines » ou refus propre) ; la conversion réelle avec Rioboo ;
$\ln\lvert u\rvert$ pour les racines réelles.

---

## 10. Pièges et cas de contrôle

**Vérification.** Chaque primitive attendue ci-dessous a été dérivée avec sympy 1.14 (script
`batterie.py`, rapport d'agent C1a) : `simplify(diff(F) - f) == 0` pour toutes, et écart
$F' - f$ évalué à 40 chiffres en 4 à 6 points ($< 10^{-168}$). Les algorithmes (Yun, Hermite
linéaire, Horowitz, LRT par sous-résultants déterminantaux, LogVersAtan) ont été réimplantés
en Python (`algos.py`) et appliqués aux exemples : leurs sorties sont celles citées.

| # | Intégrande | Primitive attendue (forme réelle) | Méthode | Piège visé |
|---|---|---|---|---|
| R1 | $\frac{1}{x^2+1}$ | $\arctan x$ | RT + conversion | forme complexe à convertir |
| R2 | $\frac{1}{x^3+1}$ | $\frac13\ln\lvert x+1\rvert - \frac16\ln(x^2-x+1) + \frac{1}{\sqrt3}\arctan\frac{2x-1}{\sqrt3}$ | LRT + LogVersReel | racine réelle et paire complexe |
| R3 | $\frac{x^4-3x^2+6}{x^6-5x^4+5x^2+4}$ | $\arctan\frac{x^5-3x^3+x}{2} + \arctan x^3 + \arctan x$ | LRT + Rioboo | **arctangente naïve discontinue** en $\pm\sqrt2$ |
| R4 | $\frac{1}{(x^2+1)^2}$ | $\frac{x}{2(x^2+1)} + \frac12\arctan x$ | Hermite | |
| R5 | $\frac{x^7-24x^4-4x^2+8x-8}{x^8+6x^6+12x^4+8x^2}$ | $\frac{3x^3+8x^2+6x+4}{x(x^2+2)^2} + \ln\lvert x\rvert$ | Hermite (3 étapes) | multiplicités 2 et 3 |
| R6 | $\frac{2x}{x^2-1}$ | $\ln\lvert x^2-1\rvert$ | LRT, $R = (z-1)^2$ | ne **pas** écrire $\ln\lvert x-1\rvert + \ln\lvert x+1\rvert$ séparément est un choix, les deux sont justes |
| R7 | $\frac{1}{x^2-2}$ | $\frac{\sqrt2}{4}\ln\left\lvert\frac{x-\sqrt2}{x+\sqrt2}\right\rvert$ | RT | extension $\mathbb{Q}(\sqrt2)$ obligatoire |
| R8 | $\frac{1}{x^4+1}$ | $\frac{\sqrt2}{8}\ln\frac{x^2+\sqrt2x+1}{x^2-\sqrt2x+1} + \frac{\sqrt2}{4}\left(\arctan(\sqrt2x+1)+\arctan(\sqrt2x-1)\right)$ | LRT + LogVersReel | degré 4, arctangente naïve discontinue en $\pm1$ |
| R9 | $\frac{x}{x^4+1}$ | $\frac12\arctan(x^2)$ | substitution $u = x^2$ (plus simple que LRT) | ordre des méthodes |
| R10 | $\frac{x^3+x}{x^2-1}$ | $\frac{x^2}{2} + \ln\lvert x^2-1\rvert$ | division puis RT | oublier la partie polynomiale |
| R11 | $\frac{1}{x(x+1)^2}$ | $\ln\lvert x\rvert - \ln\lvert x+1\rvert + \frac{1}{x+1}$ | Hermite + RT | |
| R12 | $\frac{1}{x^2+x+1}$ | $\frac{2}{\sqrt3}\arctan\frac{2x+1}{\sqrt3}$ | RT, $R = 3z^2+1$ | |
| R13 | $\frac{1}{(x-1)^3}$ | $-\frac{1}{2(x-1)^2}$ | Hermite seul, aucun logarithme | $R$ identiquement absent |
| R14 | $\frac{x^2+1}{x^4+x^2+1}$ | $\frac{1}{\sqrt3}\left(\arctan\frac{2x+1}{\sqrt3}+\arctan\frac{2x-1}{\sqrt3}\right)$ | LRT ($R=(12z^2+1)^2$) + Rioboo | naïf discontinu en 0 |

**Pièges supplémentaires.**

1. **Valeurs absolues.** Sur $\mathbb{R}$, $\int \frac{dx}{x} = \ln\lvert x\rvert$ ; l'actuelle
   règle $\ln(ax+b)/a$ de `src/Regles.cpp:186` n'est juste que là où $ax + b > 0$ (la branche
   en cours introduit `lnAbs`). Pour $u\log(A^2+B^2)$, l'argument est positif : pas de valeur
   absolue.
2. **Primitive « juste » mais fausse pour l'intégrale définie.** R3, R8, R14 sous forme
   naïve : $F(b) - F(a)$ faux dès que $[a, b]$ contient un zéro du dénominateur de
   l'arctangente. Test de contrôle : $\int_0^2$ de l'intégrande de R3 doit valoir
   $\arctan 5 + \arctan 8 + \arctan 2$ (forme de Rioboo, $F(0) = 0$), et non la valeur
   obtenue avec la forme naïve, qui diffère de $\pi$.
3. **Coefficients approchés.** Un intégrande à coefficients réels approchés (`estExact()`
   faux, `include/Polynome.hpp:65`) rend les PGCD et résultants instables : ne pas lancer
   Hermite/LRT dessus, ou convertir en rationnels exacts d'abord (choix à documenter).
4. **$\gcd(A, D) \ne 1$.** Simplifier la fraction avant RT : sinon $R(0) = 0$ et un « résidu
   nul » parasite apparaît.
5. **$\deg A \ge \deg D$.** Toujours diviser d'abord (R10).
6. **Paramètres.** Avec $K = \mathbb{Q}(a)$, les résidus dépendent de $a$ et le choix réel/complexe
   des racines dépend du signe de $a$ ($\int \frac{dx}{x^2 + a}$ : arctangente si $a > 0$,
   logarithme si $a < 0$, $-1/x$ si $a = 0$). Il faut soit discuter, soit refuser.
7. **Extension minimale.** Ne pas introduire de radicaux inutiles : R6 n'en demande aucun
   (une seule racine double de $R$), alors que la factorisation de $D$ sur $\mathbb{Q}$ donne
   deux logarithmes.

---

## 11. État dans le code

Référence : worktree `theorie` (base `origin/main` 569e686), lignes vérifiées.

| Élément | État | Où |
|---|---|---|
| Polynôme exact à une variable (`Nombre` rationnels, repli GMP) | implémenté | `include/Polynome.hpp:45` |
| Reconnaissance d'un polynôme dans une expression | implémenté | `Polynome::depuisExpression`, `include/Polynome.hpp:58` |
| Division euclidienne exacte | implémenté | `src/Polynome.cpp:252` |
| PGCD (Euclide, résultat unitaire) | implémenté | `src/Polynome.cpp:280` |
| Euclide **étendu** (Bézout) | **absent** | requis par Hermite et Rioboo |
| Décomposition sans carré (Yun) | implémenté | `src/Polynome.cpp:290` (`sansCarre`) |
| Racines réelles certifiées (Sturm) | implémenté | `src/Polynome.cpp:338` (suite), `:434` (isolement), `:520` (`racinesReelles`) |
| Factorisation | **partielle** : facteurs linéaires rationnels + reste sans carré à coefficients entiers, **non** décomposé en irréductibles | `src/Polynome.cpp:590` |
| Résultant, sous-résultants, polynômes à deux variables | absent | — |
| Primitive $(ax+b)^n$ et $1/(ax+b)$ (seuls dénominateurs gérés) | implémenté | `src/Regles.cpp:180`–`189` |
| Hermite, Horowitz | absent | — |
| Rothstein–Trager, LRT | absent | — |
| Conversion réelle, Rioboo | absent | — |
| Garde-fou de l'intégrale définie : refuse une primitive dont un dénominateur, un $\ln$, un $\tan$... s'annule ou sort du domaine sur $[a, b]$ | implémenté | `src/Integration.cpp:201` (`Analyseur`), règles `:251`–`280`, appel `:656` |

Le garde-fou de `Analyseur` protège de l'erreur de la section 7.2 **à condition** que la
primitive discontinue fasse apparaître le dénominateur $Q$ (forme $\arctan(P/Q)$ : $Q$ est un
dénominateur, donc testé). Une primitive de Rioboo, continue, passera le test et donnera la
valeur exacte. (Constat par lecture du code, non exécuté.)

Aujourd'hui, `integrer()` renvoie une intégrale non évaluée pour R1–R14 sauf R13 et les cas
$1/(ax+b)$ ($(x-1)^{-3}$ relève de la règle $(ax+b)^n$, `src/Regles.cpp:183`).

### Branche en cours (`feature/integration-symbolique`, commit 0f8b5a5)

Lu en lecture seule : la « table étendue » (`src/Primitives.cpp:305`, `trinome`) traite
$1/(ax^2+bx+c)$ à discriminant négatif ou nul ($\arctan$, $-1/(au)$) et
$1/\sqrt{ax^2+bx+c}$ ($\arcsin$, $\ln$) par mise sous forme canonique ; les trinômes à racines
réelles distinctes et tout dénominateur de degré $\ge 3$ sont explicitement laissés « aux
méthodes suivantes » (fractions rationnelles, non encore écrites à ce commit). R1, R12 et la
substitution de R9 devraient donc y passer ; R2–R8, R10, R11, R14 non.

---

## 12. Références

- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, Springer — Hermite (forme
  linéaire), Horowitz–Ostrogradsky, Rothstein–Trager, Lazard–Rioboo–Trager, conversion en
  arctangentes (algorithme de Rioboo), pseudo-code dont s'inspirent les sections 3, 6 et 7.
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra*, Kluwer —
  décomposition sans carré, Hermite, Horowitz, Rothstein–Trager.
- J. von zur Gathen, J. Gerhard, *Modern Computer Algebra*, Cambridge University Press —
  algorithme de Yun, sous-résultants, intégration des fractions rationnelles et analyses de
  coût avec arithmétique rapide.
