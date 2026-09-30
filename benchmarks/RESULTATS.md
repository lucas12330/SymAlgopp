# Résultats des benchmarks SymAlgo++ / GiNaC

Mesures réalisées avec `make bench` (`benchmarks/bench_suite.cpp`, Google Benchmark,
compilation `-O3`), sur un Intel Celeron N4120 @ 1.10 GHz (Arch Linux, GCC 16).
Les compteurs mémoire proviennent d'un remplacement global de `operator new/delete`,
appliqué de la même façon aux deux bibliothèques :

* `octets_alloues` : octets alloués par itération (pression sur l'allocateur) ;
* `octets_resultat` : mémoire encore occupée par le résultat (empreinte).

Lancer un scénario précis :

```bash
./bin/bench_suite --benchmark_filter='Derivee10'
```

## Référence : avant les optimisations de fond

Branche `feature/mesure-performances` (commit de ce fichier).

| Scénario | SymAlgo++ | GiNaC | Rapport |
| :--- | ---: | ---: | :--- |
| Évaluation en un point | 137 ns | 30 966 ns | SymAlgo++ ×226 |
| Évaluation sur 10 000 points | 1,31 ms | 309 ms | SymAlgo++ ×236 |
| Dérivée première (simplifiée) | 7,1 µs | 28,2 µs | SymAlgo++ ×4 |
| **Dérivée 10e de exp(sin x)·x²** | **4 406 ms** (901 Mo alloués) | **7,4 ms** (0,8 Mo) | **GiNaC ×600** |
| Collecte de 100 termes semblables | 145 µs | 250 µs | voir note |
| Série de exp(sin x), ordre 10 | 14 µs | 4 758 µs | SymAlgo++ ×340 |
| **Mémoire du résultat (dérivée 6e)** | **375 Ko** | **4,9 Ko** | **GiNaC ×77** |

**Note sur la collecte** : le résultat de SymAlgo++ n'est pas réellement collecté
(la somme simplifiée garde 100 termes au lieu de 8 : le simplificateur ne regroupe
que des termes adjacents d'un arbre binaire). GiNaC produit la forme réduite.

**Diagnostic** : l'évaluation numérique et les séries sont déjà nettement en tête.
Les points faibles sont structurels : sans forme canonique ni partage des
sous-expressions identiques, les dérivées successives explosent en taille
(chaque dérivation double les termes, rien ne se regroupe). C'est la cible des
chantiers « hash-consing » et « forme canonique n-aire ».

## Après le chantier « hash-consing »

Branche `feature/hash-consing` : compteur de références intrusif (`Ref`), type de
noeud par `enum` (fin des `dynamic_cast`), hash-consing (chaque expression est unique
en mémoire), mémorisation de la simplification, dérivation sur graphe partagé.
Médiane de 3 répétitions.

| Scénario | SymAlgo++ (référence) | SymAlgo++ (maintenant) | GiNaC | Rapport actuel |
| :--- | ---: | ---: | ---: | :--- |
| Évaluation en un point | 137 ns | 137 ns | 31 244 ns | SymAlgo++ ×228 |
| Évaluation sur 10 000 points | 1,31 ms | 1,32 ms | 303 ms | SymAlgo++ ×229 |
| Dérivée première | 7,1 µs · 3,0 Ko | **3,2 µs · 1,6 Ko** | 28,5 µs · 4,0 Ko | SymAlgo++ ×9 |
| Dérivée 10e | 4 406 ms · 901 Mo | **0,40 ms · 107 Ko** | 7,39 ms · 833 Ko | **SymAlgo++ ×18** |
| Collecte (100 termes) | 145 µs | 65 µs | 251 µs | voir note |
| Série d'ordre 10 | 14 µs | 6,7 µs | 4 794 µs | SymAlgo++ ×716 |
| Mémoire du résultat (dérivée 6e) | 375 Ko | **4,6 Ko** | 4,9 Ko | **SymAlgo++ −6 %** |

Toutes les dérivées successives sont validées contre les coefficients de Taylor
calculés indépendamment (test `derivees_successives_contre_series`).

**Collecte** : toujours pas de regroupement réel des termes semblables (la somme
garde 100 termes) ; c'est l'objet du chantier « forme canonique n-aire ».

## Après le chantier « forme canonique n-aire »

Branche `feature/forme-canonique` : sommes et produits n-aires triés, coefficients
rationnels exacts (`Nombre`, repli GMP), construction canonique automatique.
Médiane de 3 répétitions.

| Scénario | Après hash-consing | Forme canonique | GiNaC | Rapport actuel |
| :--- | ---: | ---: | ---: | :--- |
| Évaluation en un point | 137 ns | 112 ns | 31 409 ns | SymAlgo++ ×280 |
| Évaluation sur 10 000 points | 1,32 ms | 1,07 ms | 307 ms | SymAlgo++ ×287 |
| Dérivée première | 3,2 µs | 8,8 µs | 28,3 µs | SymAlgo++ ×3,2 |
| Dérivée 10e | 0,40 ms | 1,36 ms | 7,23 ms | SymAlgo++ ×5,3 |
| **Collecte (100 termes)** | 65 µs, *non collectés* | **143 µs, 8 termes** | 238 µs | **SymAlgo++ ×1,7** |
| Série d'ordre 10 | 6,7 µs | 4,0 µs | 4 757 µs | SymAlgo++ ×1 190 |
| **Mémoire du résultat (dérivée 6e)** | 4,6 Ko | **4,2 Ko** | 4,9 Ko | **SymAlgo++ −14 %** |

**Lecture** : les dérivées sont plus lentes qu'après le hash-consing seul, car chaque
construction fait désormais un vrai travail (regroupement des termes, fusion des
exposants, arithmétique exacte), mais le résultat est une forme réduite : la collecte
est enfin effective et plus rapide que GiNaC, et deux limites jusque-là non déterminées
(`x·ln x` et `x^x` en 0) sont résolues grâce à la simplification des quotients.
Tous les scénarios sont désormais à l'avantage de SymAlgo++.

## Après le chantier « évaluation compilée »

Branche `feature/evaluation-compilee` : compilation en programme linéaire (sous-expressions
partagées calculées une fois, registres réutilisés), évaluation par blocs vectorisés,
compilation automatique quand elle est rentable. Médiane de 3 répétitions.

| Scénario | Avant | Maintenant | GiNaC | Rapport actuel |
| :--- | ---: | ---: | ---: | :--- |
| Évaluation en un point | 112 ns | 123 ns (arbre, pas de partage) | 30 820 ns | SymAlgo++ ×250 |
| **Évaluation d'une grande expression** (dérivée 8e) | 3 640 ns | **1 109 ns** | 1 607 µs | **SymAlgo++ ×1 450** |
| Évaluation sur 10 000 points (point par point) | 1,07 ms | 1,07 ms | 302 ms | SymAlgo++ ×283 |
| **Évaluation sur 10 000 points (par blocs)** | — | **0,61 ms** | 302 ms | **SymAlgo++ ×494** |
| Dérivée première | 8,8 µs | 8,7 µs | 28,1 µs | SymAlgo++ ×3,2 |
| Dérivée 10e | 1,36 ms | 1,38 ms | 7,45 ms | SymAlgo++ ×5,4 |
| Collecte (100 termes) | 143 µs | 146 µs | 255 µs | SymAlgo++ ×1,7 |
| Série d'ordre 10 | 4,0 µs | 4,1 µs | 4 784 µs | SymAlgo++ ×1 170 |
| Mémoire du résultat (dérivée 6e) | 4,2 Ko | 4,2 Ko | 4,9 Ko | SymAlgo++ −14 % |

**Choix mesurés** : le programme compilé n'est utilisé point par point que si le partage
des sous-expressions le rend rentable (sans partage, son interprétation coûte 5 à 15 %
de plus que l'arbre). Les fonctions sin/cos/exp ne sont pas vectorisées (libmvec) car
leur erreur peut atteindre 4 ULP, contre moins d'1 ULP pour la libm scalaire.

## Bilan des quatre chantiers

| Scénario | Référence initiale | Final | GiNaC |
| :--- | ---: | ---: | ---: |
| Dérivée 10e | 4 406 ms · 901 Mo | 1,38 ms · 0,56 Mo | 7,45 ms · 0,83 Mo |
| Mémoire du résultat (dérivée 6e) | 375 Ko | 4,2 Ko | 4,9 Ko |
| Collecte de 100 termes | non collectés | 8 termes, 146 µs | 255 µs |
| Grande expression, un point | 3 640 ns (après forme canonique) | 1 109 ns | 1 607 µs |
| 10 000 points | 1,31 ms | 0,61 ms | 302 ms |

SymAlgo++ est désormais devant GiNaC sur tous les scénarios mesurés.
