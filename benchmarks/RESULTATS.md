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
