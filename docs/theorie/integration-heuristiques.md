# Intégration symbolique — méthodes pratiques et heuristiques

Ce qu'un intégrateur fait **avant** (ou à la place de) l'algorithme de Risch : table et
linéarité, changement de variable « la dérivée divise », intégration par parties,
intégrales trigonométriques et substitution de Weierstrass, exponentielles-polynômes,
heuristique de Risch–Norman, et vérification systématique par dérivation. Pour chaque
méthode : test de reconnaissance, algorithme, exemples vérifiés, cas où elle boucle ou se
trompe. La théorie (Liouville, Risch) est dans
[`integration-fondements.md`](integration-fondements.md), les fractions rationnelles dans
[`integration-fractions-rationnelles.md`](integration-fractions-rationnelles.md).

## Table des matières

1. [Notations](#notations)
2. [Architecture d'un intégrateur pratique](#1-architecture-dun-intégrateur-pratique)
3. [Table de primitives et linéarité](#2-table-de-primitives-et-linéarité)
4. [Changement de variable « la dérivée divise »](#3-changement-de-variable--la-dérivée-divise-)
5. [Intégration par parties](#4-intégration-par-parties)
6. [Intégrales trigonométriques et substitution de Weierstrass](#5-intégrales-trigonométriques-et-substitution-de-weierstrass)
7. [Exponentielles-polynômes](#6-exponentielles-polynômes)
8. [Heuristique de Risch–Norman (« intégrateur du pauvre »)](#7-heuristique-de-rischnorman--intégrateur-du-pauvre-)
9. [Vérifier une primitive](#8-vérifier-une-primitive)
10. [Pièges et cas de contrôle](#9-pièges-et-cas-de-contrôle)
11. [État dans le code](#10-état-dans-le-code)
12. [Références](#11-références)

## Notations

| Symbole | Sens |
|---|---|
| $f$ | intégrande, expression d'une variable $x$ |
| $F$ | primitive cherchée : $F' = f$ sur un intervalle |
| $u = g(x)$ | candidat de changement de variable ; $t$ nouvelle variable |
| $P, Q$ | polynômes en $x$ (coefficients exacts) |
| $a, b$ | constantes (rationnelles ou symboliques non nulles) |
| $\sec = 1/\cos$ | notation commode ; SymAlgo++ écrit $\cos(x)^{-1}$ |
| $I$ | l'intégrale en cours, inconnue d'une équation (parties cycliques) |

Toutes les primitives sont à une constante près, **par intervalle** de définition ; sur
$\mathbb{R}$, $\int \frac{du}{u} = \ln\lvert u\rvert$.

---

## 1. Architecture d'un intégrateur pratique

Un intégrateur pratique est une **suite de méthodes**, de la moins chère à la plus chère,
chacune devant **rejeter vite** ce qu'elle ne sait pas faire :

```
Integrer(f) :
    si f ne dépend pas de x : renvoyer f · x
    r ← RegleDuNoeud(f)               // table, argument affine, linéarité
    si r : renvoyer r
    pour M dans [TableEtendue, DeriveeDivise, ExpPolynome, Trigonometrique,
                 FractionRationnelle, Parties, RischNorman, (Risch)] :
        r ← M(f)                       // nul si M ne s'applique pas
        si r et non ContientIntegraleNonEvaluee(r) : renvoyer r
    renvoyer IntegraleNonEvaluee(f)     // jamais de résultat faux
```

Trois principes :

1. **Pas de résultat faux.** Une méthode qui ne conclut pas rend « rien » ; le résultat final
   peut être non évalué, jamais erroné. Une vérification par dérivation (section 8) est un
   filet de sécurité utile en mode test.
2. **Bornes de récursion.** Substitution et parties s'appellent récursivement : profondeur et
   nombre de candidats bornés, sinon l'intégrateur peut boucler (sections 3.5 et 4.4).
3. **Coût nul pour ce qui marche déjà.** La règle du nœud passe d'abord, seule ; les méthodes
   suivantes ne coûtent qu'aux intégrales qu'elle ne sait pas faire.

---

## 2. Table de primitives et linéarité

### 2.1 Table (argument affine $u = ax + b$, $a \ne 0$)

| $f$ | $F$ | Domaine / remarque |
|---|---|---|
| $u^n$, $n \ne -1$ | $\dfrac{u^{n+1}}{a(n+1)}$ | $n$ non entier : $u > 0$ |
| $u^{-1}$ | $\dfrac{\ln\lvert u\rvert}{a}$ | |
| $e^{u}$ | $\dfrac{e^u}{a}$ | |
| $c^{u}$, $c > 0$, $c \ne 1$ | $\dfrac{c^u}{a\ln c}$ | |
| $\ln u$ | $\dfrac{u\ln u - u}{a}$ | $u > 0$ |
| $\sin u$, $\cos u$ | $-\dfrac{\cos u}{a}$, $\dfrac{\sin u}{a}$ | |
| $\tan u$ | $-\dfrac{\ln\lvert\cos u\rvert}{a}$ | |
| $\sec^2 u$, $\csc^2 u$ | $\dfrac{\tan u}{a}$, $-\dfrac{\cot u}{a}$ | |
| $\sec u$ | $\dfrac{\ln\lvert\sec u + \tan u\rvert}{a}$ | |
| $\dfrac{1}{1+u^2}$, $\dfrac{1}{\sqrt{1-u^2}}$ | $\dfrac{\arctan u}{a}$, $\dfrac{\arcsin u}{a}$ | $\lvert u\rvert < 1$ |
| $\arctan u$ | $\dfrac{u\arctan u - \frac12\ln(1+u^2)}{a}$ | |
| $\arcsin u$ | $\dfrac{u\arcsin u + \sqrt{1-u^2}}{a}$ | |

Une table **étendue** reconnaît en plus des formes canoniques de second degré :
$\frac{1}{\alpha x^2+\beta x+\gamma}$ (discriminant $< 0$ : arctangente après mise sous forme
canonique), $\frac{1}{\sqrt{\alpha x^2+\beta x+\gamma}}$ ($\arcsin$ ou $\ln$), $\tan^2$, $\sec^2$...

### 2.2 Linéarité

$\int (\lambda f + \mu g) = \lambda\int f + \mu\int g$ pour des constantes $\lambda, \mu$, et les
facteurs constants sortent d'un produit.

**Piège : la linéarité prématurée.** Une somme peut être intégrable alors que ses termes ne
le sont pas :
$$
\int\left(\frac{1}{\ln x} - \frac{1}{\ln^2 x}\right) = \frac{x}{\ln x},\qquad
\int e^x\left(\frac1x - \frac1{x^2}\right) = \frac{e^x}{x}.
$$
Couper la somme d'abord donne deux intégrales non élémentaires
([`integration-fondements.md`](integration-fondements.md), section 4.3).
Remède : si un terme échoue, réessayer la somme entière par les méthodes globales
(réduction au même dénominateur, Risch–Norman, parties sur un des termes : voir 4.3), ou
regrouper les termes non intégrés et les traiter ensemble.

---

## 3. Changement de variable « la dérivée divise »

### 3.1 Énoncé

**Théorème (changement de variable).** Si $g$ est dérivable sur un intervalle $I$ et $G$ est
une primitive de $h$ sur $g(I)$, alors $G \circ g$ est une primitive de $(h\circ g)\,g'$ sur $I$.

Donc, si l'on sait écrire $f(x) = h(g(x))\,g'(x)$, alors $\int f = \left(\int h\right)\circ g$.

### 3.2 Test de reconnaissance

Pour un candidat $u = g(x)$ :

1. calculer $q = f / g'$ ;
2. remplacer dans $q$ toute occurrence de $g(x)$ par une nouvelle variable $t$ ;
3. si le résultat $h(t)$ **ne contient plus $x$**, alors $f = h(g)\,g'$ : intégrer $h$ en $t$,
   puis remplacer $t$ par $g(x)$.

Le test est **exact** (si $h$ ne contient plus $x$, l'identité $f = h(g)g'$ est vraie par
construction) ; ce qui est heuristique, c'est le **choix des candidats** et la capacité de
l'étape 2 à reconnaître $g$ sous toutes ses écritures.

### 3.3 Dans un AST canonique hash-consé

Avec une forme canonique automatique (produits $\prod b_i^{e_i}$ triés, sommes à
coefficients, hash-consing) :

- **Candidats** : sous-expressions de $f$ — arguments des fonctions, bases et exposants des
  puissances, fonctions elles-mêmes — obtenues par un parcours en profondeur, sans doublon
  (comparaison d'adresses), en écartant $x$ et les fonctions affines de $x$ (déjà traitées
  par la table), et en nombre borné.
- **Quotient** : $q = f \cdot (g')^{-1}$ calculé par la multiplication canonique, qui fusionne
  les puissances d'une même base : pour $f = 2x\cos(x^2)$, $g' = 2x$, la forme canonique donne
  directement $q = \cos(x^2)$.
- **Remplacement** : une substitution structurelle $g \mapsto t$ (le hash-consing rend la
  reconnaissance d'un sous-arbre identique à $g$ immédiate), **complétée** des puissances de
  $g$ écrites autrement : $x^4 = (x^2)^2$ pour $g = x^2$ ; $e^{2x} = (e^x)^2$ et plus
  généralement $e^{rv} = (e^v)^r$ ($r$ rationnel) pour $g = e^v$.
- **Si $x$ subsiste**, tenter une simplification rationnelle de $q$ (réduction par PGCD) avant
  d'abandonner le candidat : $\frac{x+1}{2x+2}$ ne devient $\frac12$ qu'après réduction.
- **Test « contient $x$ »** : avec l'infrastructure multivariable, $t$ est une variable
  distincte de $x$ ; le test est alors « l'expression ne contient aucune variable autre que
  $t$ ».

```
DeriveeDivise(f, x) :
    si profondeur ≥ PROFONDEUR_MAX : renvoyer rien
    t ← nouvelle variable
    pour u dans Candidats(f) (au plus N) :
        du ← u'
        si du = 0 : continuer                      // u constant déguisé (sin² + cos²)
        h ← Remplacer(f / du, u → t)
        si h contient x :
            h ← Remplacer(Simplifier(f / du), u → t)
            si h contient x : continuer
        H ← Integrer(h, t)                          // profondeur + 1
        si H contient une intégrale non évaluée : continuer
        renvoyer Remplacer(H, t → u)
    renvoyer rien
```

### 3.4 Exemples (vérifiés, section 9)

| $f$ | $u$ | $h(t)$ | $F$ |
|---|---|---|---|
| $2x\cos(x^2)$ | $x^2$ | $\cos t$ | $\sin(x^2)$ |
| $x e^{x^2}$ | $x^2$ | $\frac12 e^t$ | $\frac12 e^{x^2}$ |
| $\cos x\,e^{\sin x}$ | $\sin x$ | $e^t$ | $e^{\sin x}$ |
| $\frac{\ln x}{x}$ | $\ln x$ | $t$ | $\frac12\ln^2 x$ |
| $\frac{1}{x\ln x}$ | $\ln x$ | $\frac1t$ | $\ln\lvert\ln x\rvert$ |
| $\frac{x}{x^2+1}$ | $x^2+1$ | $\frac{1}{2t}$ | $\frac12\ln(x^2+1)$ |
| $x^2(x^3+1)^5$ | $x^3+1$ | $\frac13 t^5$ | $\frac{(x^3+1)^6}{18}$ |
| $\frac{e^x}{1+e^{2x}}$ | $e^x$ | $\frac{1}{1+t^2}$ (grâce à $e^{2x} = (e^x)^2$) | $\arctan(e^x)$ |
| $\frac{x}{\sqrt{1-x^4}}$ | $x^2$ | $\frac{1}{2\sqrt{1-t^2}}$ (grâce à $x^4 = (x^2)^2$) | $\frac12\arcsin(x^2)$ |
| $x^3 e^{x^2}$ | $x^2$ | $\frac12 t e^t$ (puis parties) | $\frac12(x^2-1)e^{x^2}$ |
| $\sin x\cos x$ | $\sin x$ | $t$ | $\frac12\sin^2 x$ |

### 3.5 Échecs et pièges

1. **Racines et signe.** Pour $u = x^2$, réécrire $x = \sqrt t$ est **faux** pour $x < 0$.
   Seules les réécritures sûres sont permises : $x = t^{k}$ quand $u = x^{1/k}$
   (ex. $u = \sqrt x$, $x = t^2$, valable sur le domaine $x \ge 0$ de $u$). $\int e^{x^2}$ ne
   doit **pas** « passer » avec $u = x^2$ et $q = \frac{e^t}{2\sqrt t}$.
2. **Constante multiplicative** : $\int x e^{x^2}$ demande $q = \frac12 e^t$ ; la forme
   canonique gère le $\frac12$, mais un test qui exigerait « $f$ contient exactement $g'$ »
   échouerait.
3. **Formes équivalentes non reconnues** : $\int \sin x\cos x$ a plusieurs primitives
   ($\frac12\sin^2 x$, $-\frac12\cos^2 x$, $-\frac14\cos 2x$) qui diffèrent d'une constante :
   un test de non-régression doit comparer les **dérivées**, pas les textes.
4. **Candidat dont la dérivée est nulle** sans que la forme canonique le voie
   ($\sin^2 x + \cos^2 x$) : division par zéro symbolique ; vérifier $u' \ne 0$.
5. **Boucle** : $h$ peut être plus compliquée que $f$ (substitution inverse) ; borner la
   profondeur (2 suffit en pratique) et le nombre de candidats.
6. **Succès partiel** : si $\int h$ revient non évaluée, abandonner le candidat (sinon on
   remplace une intégrale non évaluée par une autre, plus obscure).
7. **Domaine** : $\ln\lvert u\rvert$ au retour, simplifiable en $\ln u$ quand $u > 0$ est
   manifeste ($x^2 + 1$, $e^x + 1$).

---

## 4. Intégration par parties

### 4.1 Formule et choix

$\int u\,v' = u v - \int u' v$. Le choix de $u$ suit la règle **LIATE** (ou ILATE) : prendre
pour $u$ le facteur le plus haut dans l'ordre

| L | I | A | T | E |
|---|---|---|---|---|
| logarithmes | réciproques ($\arctan$, $\arcsin$) | algébriques (polynômes) | trigonométriques | exponentielles |

et pour $v'$ le reste, qu'on doit savoir intégrer. Justification : dériver $\ln$ et les
réciproques donne des fractions algébriques (simplification), dériver un polynôme baisse son
degré, alors qu'intégrer $\sin$, $\cos$, $e^{ax}$ ne les complique pas.

Cas $v' = 1$ : $\int \ln x = x\ln x - \int 1$, $\int \arctan x = x\arctan x - \int\frac{x}{1+x^2}$.

### 4.2 Terminaison

Pour $\int P(x)\,\phi(x)$, $\phi \in \{e^{ax}, \sin bx, \cos bx\}$, avec $u = P$ : chaque
étape baisse $\deg P$ de 1, terminaison en $\deg P + 1$ étapes (et formule fermée en
section 6). Pour $\int P(x)\ln^n x$ avec $u = \ln^n x$ : chaque étape baisse $n$. Pour
$\int x^m\arctan x$ : une étape, puis fraction rationnelle. Une **mesure qui décroît**
(degré, puissance du logarithme) doit être exhibée pour chaque règle ; sans elle, borner la
profondeur.

### 4.3 Cas cycliques : équation linéaire en l'intégrale

Quand deux intégrations par parties ramènent l'intégrale de départ, on résout une équation.

**$I = \int e^{ax}\sin(bx)$**, $a, b \ne 0$. Avec $u = \sin bx$ puis $u = \cos bx$
(en gardant $v' = e^{ax}$ les deux fois) :
$$
I = \frac{e^{ax}\sin bx}{a} - \frac ba\int e^{ax}\cos bx
  = \frac{e^{ax}\sin bx}{a} - \frac{b}{a}\left(\frac{e^{ax}\cos bx}{a} + \frac ba I\right),
$$
d'où $I\left(1 + \frac{b^2}{a^2}\right) = \frac{e^{ax}(a\sin bx - b\cos bx)}{a^2}$ :
$$
\int e^{ax}\sin bx = \frac{e^{ax}(a\sin bx - b\cos bx)}{a^2+b^2},\qquad
\int e^{ax}\cos bx = \frac{e^{ax}(a\cos bx + b\sin bx)}{a^2+b^2}.
$$

**$\int \sin(\ln x)$** : deux parties ($v' = 1$) donnent $I = x\sin\ln x - x\cos\ln x - I$,
d'où $I = \frac x2(\sin\ln x - \cos\ln x)$.

**$\int \sec^3 x$** : $u = \sec x$, $v' = \sec^2 x$ :
$I = \sec x\tan x - \int \sec x\tan^2 x = \sec x\tan x - \int(\sec^3 x - \sec x)$, donc
$2I = \sec x\tan x + \ln\lvert\sec x + \tan x\rvert$.

**Algorithme.** Garder la trace des intégrales rencontrées (comparaison d'adresses grâce au
hash-consing) ; si l'intégrale courante réapparaît sous la forme $\int\ldots = R + c\,I$ :

- $c \ne 1$ : $I = R/(1 - c)$ ;
- $c = 1$ : l'équation est $0 = R$, **sans information** — c'est le signe d'un mauvais choix
  (par exemple $u = e^{ax}$ à la première étape puis $u = \sin$ à la seconde : on « défait »
  la première intégration). Abandonner ce chemin.

### 4.4 Échecs et boucles

1. **Alternance de choix** ($c = 1$ ci-dessus) : boucle infinie si l'on ne détecte pas le
   retour de $I$.
2. **Parties qui compliquent** : $\int x e^{x^2}$ avec $u = x$ donne $\int e^{x^2}$, non
   élémentaire, alors que la substitution réussit : essayer la substitution **avant** les
   parties.
3. **Parties « fausses bonnes »** : $\int \frac{e^x}{x}$ avec $u = 1/x$ donne
   $\frac{e^x}{x} + \int\frac{e^x}{x^2}$, puis $\int\frac{e^x}{x^2}$ etc. : les puissances de
   $x$ croissent sans fin. Borner la profondeur.
4. **Mais** cette même étape résout $\int e^x\left(\frac1x - \frac1{x^2}\right)$ : les deux
   intégrales non élémentaires se compensent (le résultat de la méthode est
   $\frac{e^x}{x} + \int\frac{e^x}{x^2} - \int\frac{e^x}{x^2}$ ; avec le hash-consing, les deux
   nœuds non évalués sont le même nœud et la somme canonique devrait les annuler — non testé).
   Intéressant, mais fragile : Risch–Norman (section 7) le fait de façon systématique.

---

## 5. Intégrales trigonométriques et substitution de Weierstrass

### 5.1 $\int \sin^m x\cos^n x$, $m, n \in \mathbb{N}$

- **$m$ impair** : $u = \cos x$, $\sin^m x = \sin x\,(1 - \cos^2 x)^{(m-1)/2}$, d'où un polynôme
  en $u$. Ex. : $\int\sin^3 x\cos^2 x = -\frac{\cos^3 x}{3} + \frac{\cos^5 x}{5}$ ;
  $\int \sin^5 x = -\cos x + \frac23\cos^3 x - \frac15\cos^5 x$.
- **$n$ impair** : $u = \sin x$ symétriquement.
- **$m$ et $n$ pairs** : **linéariser** avec $\sin^2 = \frac{1 - \cos 2x}{2}$,
  $\cos^2 = \frac{1+\cos 2x}{2}$ (ou les formules d'Euler, $\cos^n x = 2^{-n}(e^{ix}+e^{-ix})^n$
  développé). Ex. : $\int\sin^2 x = \frac x2 - \frac{\sin 2x}{4}$ ;
  $\int\sin^2 x\cos^2 x = \frac x8 - \frac{\sin 4x}{32}$ ;
  $\int\cos^4 x = \frac{3x}{8} + \frac{\sin 2x}{4} + \frac{\sin 4x}{32}$.

Exposants négatifs : $\sec^n$, $\csc^n$ par réduction (5.3) ; mélanges quelconques : Weierstrass
(5.4).

### 5.2 Produits $\sin(ax)\cos(bx)$

$$
\sin\alpha\cos\beta = \tfrac12[\sin(\alpha+\beta) + \sin(\alpha-\beta)],\quad
\sin\alpha\sin\beta = \tfrac12[\cos(\alpha-\beta) - \cos(\alpha+\beta)],\quad
\cos\alpha\cos\beta = \tfrac12[\cos(\alpha-\beta) + \cos(\alpha+\beta)].
$$
Ex. : $\int\sin 3x\cos 5x = \frac{\cos 2x}{4} - \frac{\cos 8x}{16}$. **Piège** : $a = \pm b$
donne un terme $\sin 0 = 0$ ou $\cos 0 = 1$ (donc un terme en $x$) : ne pas diviser par $a - b$.

### 5.3 Puissances de $\tan$ et de $\sec$ (réductions)

$$
\int\tan^n x = \frac{\tan^{n-1}x}{n-1} - \int\tan^{n-2}x,\qquad
\int\sec^n x = \frac{\sec^{n-2}x\tan x}{n-1} + \frac{n-2}{n-1}\int\sec^{n-2}x\qquad(n \ge 2).
$$
(Vérifiées par dérivation pour $n = 2, \dots, 8$.) La première vient de
$\tan^n = \tan^{n-2}(\sec^2 - 1)$, la seconde d'une intégration par parties. Ex. :
$\int\tan^3 x = \frac{\tan^2 x}{2} + \ln\lvert\cos x\rvert$ ;
$\int\tan^4 x = \frac{\tan^3 x}{3} - \tan x + x$. Mesure décroissante : $n$.

### 5.4 Substitution de Weierstrass $t = \tan(x/2)$

Pour une fraction rationnelle $R(\sin x, \cos x)$ :
$$
\sin x = \frac{2t}{1+t^2},\quad \cos x = \frac{1-t^2}{1+t^2},\quad \tan x = \frac{2t}{1-t^2},\quad
dx = \frac{2\,dt}{1+t^2},
$$
et $\int R(\sin x, \cos x)\,dx = \int R\!\left(\frac{2t}{1+t^2}, \frac{1-t^2}{1+t^2}\right)\frac{2\,dt}{1+t^2}$
est une **fraction rationnelle en $t$** : toujours élémentaire, traitée par le document dédié.
C'est la méthode **complète** pour cette classe, mais à garder en dernier recours : elle
produit des expressions lourdes.

**Piège : discontinuité en $x = \pi + 2k\pi$.** $t = \tan(x/2)$ n'est défini que sur
$]-\pi + 2k\pi, \pi + 2k\pi[$ ; si $G(t)$ est la primitive en $t$, $G(\tan(x/2))$ est une
primitive **sur chaque intervalle**, mais elle saute en $x = \pi + 2k\pi$ de
$G(-\infty) - G(+\infty)$ (limites de $G$). Correction (quand ces limites sont finies) :
$$
F(x) = G\!\left(\tan\frac x2\right) + \big(G(+\infty) - G(-\infty)\big)\left\lfloor\frac{x+\pi}{2\pi}\right\rfloor,
$$
prolongée par continuité aux points $\pi + 2k\pi$.

**Exemple.** $\int\frac{dx}{2 + \cos x}$ : $G(t) = \int\frac{2\,dt}{3 + t^2} =
\frac{2}{\sqrt3}\arctan\frac{t}{\sqrt3}$, d'où
$G(\tan\frac x2) = \frac{2}{\sqrt3}\arctan\frac{\tan(x/2)}{\sqrt3}$, qui saute de
$-\frac{2\pi}{\sqrt3}$ en $x = \pi$ (limite calculée). Primitive continue :
$$
F(x) = \frac{2}{\sqrt3}\arctan\frac{\tan(x/2)}{\sqrt3} + \frac{2\pi}{\sqrt3}\left\lfloor\frac{x+\pi}{2\pi}\right\rfloor ,
$$
ou, **sans partie entière**, la forme continue
$$
\int\frac{dx}{A + B\cos x} = \frac{1}{\sqrt{A^2-B^2}}\left(x - 2\arctan\frac{B\sin x}{A + \sqrt{A^2-B^2} + B\cos x}\right)\qquad (A > \lvert B\rvert),
$$
vérifiée symboliquement (dérivée égale à l'intégrande pour $A, B$ formels) ; son dénominateur
$A + \sqrt{A^2-B^2} + B\cos x \ge A - \lvert B\rvert + \sqrt{A^2-B^2} > 0$, donc pas de saut.
Pour $A = 2$, $B = 1$ : $F(x) = \frac{1}{\sqrt3}\left(x - 2\arctan\frac{\sin x}{2+\sqrt3+\cos x}\right)$,
qui coïncide avec la forme corrigée (écart nul calculé de part et d'autre de $\pi$ et $3\pi$).

**Autres formes utiles** (plus simples que Weierstrass) : $\int\frac{dx}{\sin x} =
\ln\lvert\tan\frac x2\rvert$ ; $\int\frac{dx}{1+\sin x} = \tan x - \sec x$ (aussi
$-\frac{2}{1+\tan(x/2)}$, à une constante près **par intervalle**).

---

## 6. Exponentielles-polynômes

### 6.1 $\int P(x)\,e^{ax}$, $a \ne 0$

**Ansatz.** La primitive est de la forme $Q(x)e^{ax}$ avec $\deg Q = \deg P$, et
$(Qe^{ax})' = (Q' + aQ)e^{ax}$ : il faut résoudre $Q' + aQ = P$. C'est l'**équation de Risch**
dans son cas le plus simple ([`integration-fondements.md`](integration-fondements.md), section 5.7),
toujours résoluble ici. Formule fermée :
$$
\int P(x)e^{ax}\,dx = e^{ax}\sum_{k=0}^{\deg P}\frac{(-1)^k P^{(k)}(x)}{a^{k+1}}
$$
(vérifiée symboliquement pour $P$ générique de degré 5). Ex. :
$\int (x^2+1)e^{2x} = e^{2x}\left(\frac{x^2}{2} - \frac x2 + \frac34\right)$ ;
$\int x e^x = (x - 1)e^x$.

### 6.2 $\int P(x)\cos(bx)$ et $\int P(x)\sin(bx)$, $b \ne 0$

$$
\int P\cos bx = \sin bx\sum_{k\ge0}\frac{(-1)^k P^{(2k)}}{b^{2k+1}} + \cos bx\sum_{k\ge0}\frac{(-1)^k P^{(2k+1)}}{b^{2k+2}},
$$
$$
\int P\sin bx = -\cos bx\sum_{k\ge0}\frac{(-1)^k P^{(2k)}}{b^{2k+1}} + \sin bx\sum_{k\ge0}\frac{(-1)^k P^{(2k+1)}}{b^{2k+2}}
$$
(vérifiées symboliquement, $P$ générique de degré 5). Ex. :
$\int x^3\cos x = (x^3 - 6x)\sin x + (3x^2 - 6)\cos x$ ;
$\int x^2\sin x = -x^2\cos x + 2x\sin x + 2\cos x$.

### 6.3 $\int P(x)\,e^{ax}\cos(bx)$, $\int P(x)\,e^{ax}\sin(bx)$

Passer par $e^{(a+ib)x}$ : $\int P e^{(a+ib)x} = e^{(a+ib)x}\sum_k\frac{(-1)^kP^{(k)}}{(a+ib)^{k+1}}$,
puis prendre parties réelle et imaginaire ; ou, sans complexes, l'ansatz
$e^{ax}(Q_1\cos bx + Q_2\sin bx)$ avec $\deg Q_i \le \deg P$ et un système linéaire. Ex. :
$\int x e^x\sin x = \frac{e^x}{2}(x\sin x - x\cos x + \cos x)$.

**Reconnaissance.** $f$ est un produit dont un facteur est polynomial en $x$ et les autres
de la forme $e^{ax+c}$, $\sin(bx+c)$, $\cos(bx+c)$ à arguments affines ; tout autre facteur
fait rejeter. Coût : $O(\deg P)$ dérivées.

---

## 7. Heuristique de Risch–Norman (« intégrateur du pauvre »)

### 7.1 Principe

L'algorithme de Risch traite la tour **étage par étage** ; l'approche « parallèle »
(Risch–Norman) traite **tous les générateurs à la fois** :

1. Écrire $f \in C(x, \theta_1, \dots, \theta_k)$ où les $\theta_i$ sont les fonctions présentes
   ($e^{u}$, $\ln u$, $\tan u$...) et forment un corps stable par dérivation :
   $\theta_i' \in C(x, \theta)$. La dérivation est
   $D = \frac{\partial}{\partial x} + \sum_i \theta_i'\,\frac{\partial}{\partial \theta_i}$.
2. Écrire $f = p/q$ avec $p, q \in C[x, \theta]$.
3. **Ansatz** : $F = \dfrac{P(x, \theta)}{d} + \sum_j c_j\log v_j$, avec
   - un dénominateur candidat $d$, typiquement $\gcd(q, Dq)$ (les facteurs normaux de $q$ avec
     leur multiplicité diminuée de 1, plus les facteurs spéciaux) ;
   - des logarithmes candidats $v_j$ : les facteurs irréductibles de $q$ et les polynômes
     spéciaux ($1 + \theta^2$ pour $\theta = \tan u$) ;
   - $P$ à coefficients **inconnus**, de degré total borné (borne heuristique, par exemple
     $1 + \max(\deg p, \deg q)$).
4. Identifier les coefficients de $D(F) = f$ (après multiplication par le dénominateur
   commun) : **système linéaire** en les coefficients de $P$ et les $c_j$.
5. Solution trouvée : $F$ (à vérifier par dérivation). Pas de solution : **échec**, pas une
   preuve.

Bronstein a diffusé une implémentation compacte de cette idée sous le nom de
« poor man's integrator » (`pmint`) ; le principe ci-dessus est le sien, les détails de ses
bornes ne sont pas reproduits ici.

### 7.2 Exemples

**$\int\tan x$.** $\theta = \tan x$, $\theta' = 1 + \theta^2$ (spécial). $f = \theta$, $q = 1$.
Ansatz $F = P(x, \theta) + c\log(1 + \theta^2)$ :
$D(F) = P_x + (1+\theta^2)P_\theta + c\,\frac{2\theta(1+\theta^2)}{1+\theta^2}
= P_x + (1+\theta^2)P_\theta + 2c\,\theta$. Avec $P = 0$, $c = \frac12$ :
$F = \frac12\ln(1 + \tan^2 x) = -\ln\lvert\cos x\rvert$.

**$\int e^x\frac{x-1}{x^2}$.** $\theta = e^x$, $D = \partial_x + \theta\,\partial_\theta$,
$q = x^2$, $Dq = 2x$, $d = \gcd(x^2, 2x) = x$. Ansatz $F = P(x,\theta)/x$ :
$D(F) = \frac{x(P_x + \theta P_\theta) - P}{x^2}$, à égaler à $\frac{\theta(x-1)}{x^2}$, soit
$x(P_x + \theta P_\theta) - P = \theta x - \theta$. $P = \theta$ convient (et $P = \lambda x$
donne la constante d'intégration). $F = e^x/x$ — là où linéarité et substitution échouent.

### 7.3 Forces et échecs

- **Forces** : un seul système linéaire, pas de récursion sur la tour, traite d'un coup des
  intégrandes mêlant $\exp$, $\ln$, $\tan$ ; résout F6, F7 de
  [`integration-fondements.md`](integration-fondements.md#7-pièges-et-cas-de-contrôle).
- **Borne de degré trop petite** : solution manquée ; trop grande : système énorme (le nombre
  d'inconnues croît comme $\binom{N + k + 1}{k + 1}$ pour un degré total $N$ et $k$
  générateurs).
- **Logarithmes manquants** : $\int\frac{dx}{x^2-2} = \frac{\sqrt2}{4}\ln\left\lvert\frac{x-\sqrt2}{x+\sqrt2}\right\rvert$
  demande des logarithmes de facteurs **sur $\mathbb{Q}(\sqrt2)$** ; avec les seuls facteurs
  sur $\mathbb{Q}$ ($v = x^2 - 2$), l'ansatz $c\ln(x^2-2)$ donne $\frac{2cx}{x^2-2} \ne
  \frac{1}{x^2-2}$ : échec. Remède : traiter la partie logarithmique par Rothstein–Trager /
  LRT.
- **Nouveaux générateurs** : si la primitive demande une fonction absente de $f$ (un
  $\arctan$ pour une fraction réelle, un $\ln$ d'une expression qui ne divise pas $q$),
  l'ansatz ne peut pas la trouver.
- **Non-élémentarité** : l'échec ne prouve rien ($\int e^{x^2}$ échoue, mais pour une bonne
  raison que seul Risch démontre).

---

## 8. Vérifier une primitive

### 8.1 Principe

$F$ est une primitive de $f$ si et seulement si $F' - f = 0$. La dérivation est
**algorithmique et sûre** ; tout le problème est le **test de zéro** de $F' - f$.

```
Verifier(F, f) :
    e ← Simplifier(F' − f)             // forme canonique + réduction rationnelle
    si e = 0 : renvoyer VRAI           // preuve
    // sinon : peut-être nul mais non reconnu (sin² + cos² − 1, identités de log...)
    pour k points aléatoires x_k du domaine :
        si |e(x_k)| > tolérance : renvoyer FAUX       // preuve de fausseté (aux arrondis près)
    renvoyer PROBABLEMENT_VRAI
```

- Pour une **fraction rationnelle** $e$, la forme normale (numérateur/dénominateur premiers
  entre eux) décide : $e = 0$ si et seulement si son numérateur est le polynôme nul.
- Pour un **polynôme**, l'évaluation en un point aléatoire d'un grand ensemble fini est un test
  probabiliste rigoureux (lemme de Schwartz–Zippel : un polynôme non nul de degré $d$ s'annule
  en un point tiré uniformément dans $S$ avec probabilité $\le d/\lvert S\rvert$).
- Pour des expressions **transcendantes**, l'évaluation numérique est seulement une
  heuristique (arrondis, points mal choisis, branches de $\ln$ et $\sqrt{\ }$).

### 8.2 Indécidabilité (énoncé seulement)

**Théorème de Richardson (1968, forme usuelle).** Pour la classe des expressions en $x$
construites à partir des rationnels, de $\pi$, de $\ln 2$, de $x$, par addition,
multiplication et composition avec $\sin$, $\exp$ et la valeur absolue, il n'existe **pas**
d'algorithme décidant si une expression est identiquement nulle sur $\mathbb{R}$.

Conséquence : aucun intégrateur ne peut garantir la vérification de toutes ses primitives ;
il faut se restreindre à des classes où le test de zéro est décidable (fractions
rationnelles, polynômes en exponentielles et logarithmes **algébriquement indépendants** —
d'où l'importance du théorème de structure de Risch — avec un test de zéro des constantes
qui, pour les constantes exp-log, n'est connu décidable que sous la conjecture de Schanuel).

### 8.3 Pièges de la vérification

1. **Branches** : $\frac{d}{dx}\ln\lvert u\rvert = \frac{u'}{u}$ partout où $u \ne 0$, mais
   $\ln u$ n'est défini que pour $u > 0$ : une vérification numérique en un point où $u < 0$
   doit utiliser $\ln\lvert u\rvert$.
2. **$\sqrt{x^2} = \lvert x\rvert$**, pas $x$ : une simplification qui l'oublie « vérifie » des
   primitives fausses pour $x < 0$.
3. **Primitives discontinues** : $F' = f$ en tout point n'empêche pas un saut de $F$
   (arctangente naïve, Weierstrass) ; la vérification par dérivation ne le voit **pas**. Pour
   l'intégrale définie, il faut en plus la continuité de $F$ sur $[a, b]$ (ce que fait
   l'`Analyseur` de `src/Integration.cpp`).
4. **Formes différentes** : deux primitives justes diffèrent d'une constante ; comparer des
   textes de sortie dans les tests est fragile, comparer les dérivées (ou $F_1 - F_2$
   constante) est robuste.

---

## 9. Pièges et cas de contrôle

**Vérification.** Chaque primitive ci-dessous a été dérivée avec sympy 1.14 (script
`batterie.py` du rapport d'agent C1a) : `simplify(diff(F) - f) == 0` pour toutes sauf H24
(où `simplify` ne conclut pas seul ; `expand_trig` et la réécriture en exponentielles donnent
0), **et** $F' - f$ évalué à 40 chiffres en 3 à 6 points du domaine ($< 10^{-168}$). Les
formules générales (6.1, 6.2, 4.3, réductions 5.3, forme continue 5.4) ont été vérifiées
symboliquement (`verif_formules.py`).

### 9.1 Table et linéarité

| # | $f$ | $F$ attendue | Remarque |
|---|---|---|---|
| H1 | $3x^2 - 2x + 5$ | $x^3 - x^2 + 5x$ | |
| H2 | $\frac{1}{2x+3}$ | $\frac12\ln\lvert 2x+3\rvert$ | valeur absolue sur $\mathbb{R}$ |
| H3 | $e^{3x-1}$ | $\frac13 e^{3x-1}$ | |
| H4 | $2^x$ | $\frac{2^x}{\ln 2}$ | |
| H33 | $\sec x$ | $\ln\lvert\sec x + \tan x\rvert$ | |

### 9.2 La dérivée divise

| # | $f$ | $F$ attendue | Candidat |
|---|---|---|---|
| H5 | $2x\cos(x^2)$ | $\sin(x^2)$ | $x^2$ |
| H6 | $\cos x\,e^{\sin x}$ | $e^{\sin x}$ | $\sin x$ |
| H7 | $\frac{x}{x^2+1}$ | $\frac12\ln(x^2+1)$ | $x^2+1$ |
| H8 | $x^2(x^3+1)^5$ | $\frac{(x^3+1)^6}{18}$ | $x^3+1$ (ne pas développer !) |
| H9 | $\sin x\cos x$ | $\frac12\sin^2 x$ | $\sin x$ ; autres formes à une constante près |
| H10 | $\tan x$ | $-\ln\lvert\cos x\rvert$ | $\cos x$ |
| H11 | $\frac{e^x}{1+e^{2x}}$ | $\arctan(e^x)$ | $e^x$, avec $e^{2x} = (e^x)^2$ |
| H12 | $\frac{x}{\sqrt{1-x^4}}$ | $\frac12\arcsin(x^2)$ | $x^2$, avec $x^4 = (x^2)^2$ |
| H13 | $x^3e^{x^2}$ | $\frac12(x^2-1)e^{x^2}$ | $x^2$, puis parties en $t$ |
| H35 | $x\sqrt{x+1}$ | $\frac25(x+1)^{5/2} - \frac23(x+1)^{3/2}$ | $\sqrt{x+1}$ ou $x+1$ |
| H37 | $\frac{1}{x\sqrt{x^2-1}}$ | $\arctan\sqrt{x^2-1}$ | $\sqrt{x^2-1}$ ; valable pour $\lvert x\rvert > 1$ (vérifié en $x = -2$) |

### 9.3 Parties (y compris cycliques)

| # | $f$ | $F$ attendue | Remarque |
|---|---|---|---|
| H14 | $x e^x$ | $(x-1)e^x$ | |
| H15 | $x^2\sin x$ | $-x^2\cos x + 2x\sin x + 2\cos x$ | 3 étapes |
| H16 | $x\ln x$ | $\frac{x^2}{2}\ln x - \frac{x^2}{4}$ | $u = \ln x$ (L avant A) |
| H17 | $\arctan x$ | $x\arctan x - \frac12\ln(1+x^2)$ | $v' = 1$ |
| H18 | $e^{2x}\sin 3x$ | $\frac{e^{2x}}{13}(2\sin 3x - 3\cos 3x)$ | cyclique, $c = -9/4 \ne 1$ |
| H36 | $e^{ax}\cos bx$ | $\frac{e^{ax}(a\cos bx + b\sin bx)}{a^2+b^2}$ | paramètres $a, b > 0$ |
| H19 | $\sin(\ln x)$ | $\frac x2(\sin\ln x - \cos\ln x)$ | cyclique, $x > 0$ |
| H20 | $\sec^3 x$ | $\frac12(\sec x\tan x + \ln\lvert\sec x+\tan x\rvert)$ | cyclique |

### 9.4 Trigonométrie

| # | $f$ | $F$ attendue | Méthode |
|---|---|---|---|
| H21 | $\sin^3 x\cos^2 x$ | $-\frac{\cos^3 x}{3} + \frac{\cos^5 x}{5}$ | $m$ impair |
| H38 | $\sin^5 x$ | $-\cos x + \frac23\cos^3 x - \frac15\cos^5 x$ | $m$ impair |
| H22 | $\sin^2 x$ | $\frac x2 - \frac{\sin 2x}{4}$ | linéarisation |
| H23 | $\sin^2 x\cos^2 x$ | $\frac x8 - \frac{\sin 4x}{32}$ | linéarisation |
| H34 | $\cos^4 x$ | $\frac{3x}{8} + \frac{\sin 2x}{4} + \frac{\sin 4x}{32}$ | linéarisation |
| H24 | $\sin 3x\cos 5x$ | $\frac{\cos 2x}{4} - \frac{\cos 8x}{16}$ | produit → somme |
| H25 | $\tan^3 x$ | $\frac{\tan^2 x}{2} + \ln\lvert\cos x\rvert$ | réduction |
| H26 | $\tan^4 x$ | $\frac{\tan^3 x}{3} - \tan x + x$ | réduction |
| H27 | $\frac{1}{2+\cos x}$ | $\frac{1}{\sqrt3}\left(x - 2\arctan\frac{\sin x}{2+\sqrt3+\cos x}\right)$ | **continue** ; la forme Weierstrass brute saute de $-\frac{2\pi}{\sqrt3}$ en $\pi$ |
| H28 | $\frac{1}{\sin x}$ | $\ln\lvert\tan\frac x2\rvert$ | |
| H32 | $\frac{1}{1+\sin x}$ | $\tan x - \sec x$ | |

### 9.5 Exponentielles-polynômes

| # | $f$ | $F$ attendue |
|---|---|---|
| H29 | $(x^2+1)e^{2x}$ | $e^{2x}\left(\frac{x^2}{2} - \frac x2 + \frac34\right)$ |
| H30 | $x^3\cos x$ | $(x^3-6x)\sin x + (3x^2-6)\cos x$ |
| H31 | $x e^x\sin x$ | $\frac{e^x}{2}(x\sin x - x\cos x + \cos x)$ |

### 9.6 Risch–Norman et sommes à ne pas couper

| $f$ | $F$ attendue | Remarque |
|---|---|---|
| $e^x\frac{x-1}{x^2}$ | $\frac{e^x}{x}$ | F7 de [`integration-fondements.md`](integration-fondements.md#71-primitives-élémentaires-le-résultat-exact-attendu) |
| $\frac{1}{\ln x} - \frac{1}{\ln^2 x}$ | $\frac{x}{\ln x}$ | F6 ; la linéarité seule échoue |
| $\tan x$ | $-\ln\lvert\cos x\rvert$ | ansatz avec $\log(1+\tan^2 x)$ |

### 9.7 Doivent rester non évaluées (et ne pas boucler)

$e^{x^2}$ (piège : $u = x^2$ et $x = \sqrt t$ interdit), $\frac{\sin x}{x}$, $\frac{e^x}{x}$
(piège : parties sans fin), $\frac{1}{\ln x}$, $\sin(x^2)$, $x^2 e^{x^2}$ (ou forme partielle
$\frac12 xe^{x^2} - \frac12\int e^{x^2}$). Contrôle : la réponse est une `IntegraleNonEvaluee`
(ou contient une), en un temps borné.

---

## 10. État dans le code

Référence : worktree `theorie` (base `origin/main` 569e686), lignes vérifiées.

| Méthode | État | Où |
|---|---|---|
| Constante, $x$, paramètre, $\pi$ | implémenté | `src/Regles.cpp:58`, `:65`, `:73`, `:78` ; point d'entrée `src/Noeud.cpp:265` |
| Linéarité (somme ; un terme non intégré reste non évalué, le reste est calculé) | implémenté | `src/Regles.cpp:103`–`109` |
| Facteurs constants d'un produit | implémenté | `src/Regles.cpp:144`–`159` |
| Produit de deux facteurs dépendant de $x$ | **non évalué** (ni parties ni substitution) | `src/Regles.cpp:158` |
| Table d'argument affine : $(ax+b)^n$, $\ln$, $c^{ax+b}$, $\sin$, $\cos$, $\tan$, $\exp$, $\ln$, $\arcsin$, $\arccos$, $\arctan$ | implémenté (détection : `coefficientLineaire`, `src/Regles.cpp:23`) | `src/Regles.cpp:180`–`199`, `:226`–`256`, `:288`–`302` |
| $\ln\lvert u\rvert$ | **absent** : $\int\frac{1}{ax+b}$ et $\int\tan$ rendent $\ln(\ldots)$ sans valeur absolue | `src/Regles.cpp:186`, `:241` |
| Dérivée divise, parties, trigonométrie, Weierstrass, exp-polynômes, Risch–Norman | absent | — |
| Vérification par dérivation | dans les tests seulement (dérivée numérique centrée) | `tests/test_ast.cpp:31`–`41` (`verifierPrimitive`) |
| Garde-fou de continuité pour l'intégrale définie | implémenté | `src/Integration.cpp:201` (`Analyseur`), `:251`–`280` |
| Test existant : $x\sin x$ non évalué, $e^{x^2}$ et $\ln(x)/x$ non évalués | implémenté | `tests/test_ast.cpp:459`–`477` |

Constat : sur `main`, $\ln(x)/x$ (F9 du document fondements) reste non évaluée alors que
la substitution la résout ; c'est attendu à ce stade (pas de substitution sur `main`).

### Branche en cours (`feature/integration-symbolique`, commit 0f8b5a5)

Lu en lecture seule dans le worktree `primitives-substitution`, sans compilation :

| Élément | État à ce commit | Où (`src/Primitives.cpp`) |
|---|---|---|
| Suite de méthodes après la règle du nœud | implémenté | `METHODES` `:548`, point d'entrée `ASTNode::integrer` `:559` |
| Table étendue : trinômes ($\arctan$, $\arcsin$, $\ln$), $\sec^2$, $\csc^2$, $\tan^2$, hyperboliques, $c^u$ symbolique | implémenté | `trinome` `:305`, `tableEtendue` `:335` |
| $\ln\lvert u\rvert$ simplifié en $\ln u$ si $u$ manifestement positif | implémenté | `estPositifOuNul` `:108`, `lnAbs` `:141` |
| Dérivée divise : candidats bornés (8 candidats, 64 nœuds), profondeur 2, remplacement des puissances $x^{kn} = (x^n)^k$ et $e^{rv} = (e^v)^r$, réduction rationnelle si $x$ subsiste, abandon si résultat partiel | implémenté | constantes `:380`–`382`, `Remplacement` `:190`, `Candidats` `:405`, `parSubstitution` `:483` |
| Réécriture $x = t^{k}$ seulement quand la cible est $x^{1/k}$ (sûre, cf. 3.5.1) | implémenté | `Remplacement::appliquerA`, cas `Variable` (`:230`–`235`) |
| Parties, fractions rationnelles, linéarisation trigonométrique | annoncé, non écrit | commentaire `:546` |

Par lecture du code (non exécuté), H5–H12 et H35 relèvent des méthodes déjà présentes ; H13
demande aussi les parties. Les autres cas de la section 9 sont des objectifs.

---

## 11. Références

- M. Bronstein, *Symbolic Integration I: Transcendental Functions*, Springer — cadre général,
  fractions rationnelles, Risch (la référence de l'heuristique de Risch–Norman n'est pas
  donnée ici faute de certitude sur l'ouvrage qui la traite).
- K. O. Geddes, S. R. Czapor, G. Labahn, *Algorithms for Computer Algebra*, Kluwer —
  intégration (méthodes heuristiques et Risch).
- J. von zur Gathen, J. Gerhard, *Modern Computer Algebra*, Cambridge University Press —
  tests de zéro probabilistes (Schwartz–Zippel), intégration symbolique.
- Le théorème de Richardson est cité d'après son énoncé usuel ; la publication originale n'est
  pas référencée ici faute de pouvoir en vérifier les détails.
