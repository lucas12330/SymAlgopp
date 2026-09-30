/**
 * @file bench_suite.cpp
 * @brief Benchmarks comparatifs SymAlgo++ / GiNaC (Google Benchmark).
 *
 * Chaque scénario est mesuré sur les deux bibliothèques avec la même expression.
 * La mémoire est mesurée en remplaçant operator new/delete : les compteurs
 * « octets_alloues » (allocations par itération) et « octets_resultat » (mémoire
 * encore occupée par le résultat) s'appliquent de la même façon aux deux.
 */

#include <benchmark/benchmark.h>
#include <ginac/ginac.h>

#include <atomic>
#include <cstdlib>
#include <new>
#include <vector>

#include "ASTNode.hpp"
#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"

using namespace symalgo;

// ============================================================================
// Comptage des allocations (remplacement global de operator new/delete)
// ============================================================================

namespace {

std::atomic<long long> g_octetsAlloues{0}; // cumul des allocations
std::atomic<long long> g_octetsVivants{0}; // mémoire actuellement occupée

// En-tête placé avant chaque bloc pour retrouver sa taille à la libération
constexpr std::size_t EN_TETE = alignof(std::max_align_t);

void* allouer(std::size_t taille) {
    void* brut = std::malloc(taille + EN_TETE);
    if (!brut) throw std::bad_alloc();
    *static_cast<std::size_t*>(brut) = taille;
    g_octetsAlloues += static_cast<long long>(taille);
    g_octetsVivants += static_cast<long long>(taille);
    return static_cast<char*>(brut) + EN_TETE;
}

void liberer(void* p) noexcept {
    if (!p) return;
    void* brut = static_cast<char*>(p) - EN_TETE;
    g_octetsVivants -= static_cast<long long>(*static_cast<std::size_t*>(brut));
    std::free(brut);
}

} // namespace

void* operator new(std::size_t n) { return allouer(n); }
void* operator new[](std::size_t n) { return allouer(n); }
void operator delete(void* p) noexcept { liberer(p); }
void operator delete[](void* p) noexcept { liberer(p); }
void operator delete(void* p, std::size_t) noexcept { liberer(p); }
void operator delete[](void* p, std::size_t) noexcept { liberer(p); }

namespace {

// Octets alloués par itération, à appeler après la boucle de mesure
void rapporterAllocations(benchmark::State& state, long long debut) {
    state.counters["octets_alloues"] = benchmark::Counter(
        static_cast<double>(g_octetsAlloues - debut) / static_cast<double>(state.iterations()));
}

// ============================================================================
// Expressions communes
// ============================================================================

// f(x) = sin(x)*cos(x) + x^3 - 2x
ExprPtr fSymAlgo() {
    const ExprPtr x = var("x");
    return (ast_sin(x) * ast_cos(x)) + ast_pow(x, 3.0) - (cst(2.0) * x);
}

GiNaC::ex fGinac(const GiNaC::symbol& x) {
    return GiNaC::sin(x) * GiNaC::cos(x) + GiNaC::pow(x, 3) - 2 * x;
}

constexpr int TERMES_COLLECTE = 50;
constexpr int POINTS_TRACE = 10000;

} // namespace

// ============================================================================
// 1. Évaluation en un point
// ============================================================================

static void BM_SymAlgo_Eval(benchmark::State& state) {
    const EquationClassique eq(fSymAlgo());
    for (auto _ : state) {
        double res = eq.eval(5.0);
        benchmark::DoNotOptimize(res);
    }
}
BENCHMARK(BM_SymAlgo_Eval);

static void BM_GiNaC_Eval(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex eq = fGinac(x);
    for (auto _ : state) {
        double res = GiNaC::ex_to<GiNaC::numeric>(eq.subs(x == 5.0).evalf()).to_double();
        benchmark::DoNotOptimize(res);
    }
}
BENCHMARK(BM_GiNaC_Eval);

// ============================================================================
// 2. Évaluation sur 10 000 points (tableau de valeurs, tracé)
// ============================================================================

static void BM_SymAlgo_EvalTableau(benchmark::State& state) {
    const EquationClassique eq(fSymAlgo());
    std::vector<double> y(POINTS_TRACE);
    for (auto _ : state) {
        for (int i = 0; i < POINTS_TRACE; ++i) y[i] = eq.eval(i * 1e-3);
        benchmark::DoNotOptimize(y.data());
    }
    state.SetItemsProcessed(state.iterations() * POINTS_TRACE);
}
BENCHMARK(BM_SymAlgo_EvalTableau);

// Même tableau, par l'évaluation compilée et vectorisée par blocs
static void BM_SymAlgo_EvalTableauVectorise(benchmark::State& state) {
    const EquationClassique eq(fSymAlgo());
    std::vector<double> xs(POINTS_TRACE);
    for (int i = 0; i < POINTS_TRACE; ++i) xs[i] = i * 1e-3;
    for (auto _ : state) {
        std::vector<double> y = eq.eval(xs);
        benchmark::DoNotOptimize(y.data());
    }
    state.SetItemsProcessed(state.iterations() * POINTS_TRACE);
}
BENCHMARK(BM_SymAlgo_EvalTableauVectorise);

static void BM_GiNaC_EvalTableau(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex eq = fGinac(x);
    std::vector<double> y(POINTS_TRACE);
    for (auto _ : state) {
        for (int i = 0; i < POINTS_TRACE; ++i) {
            y[i] = GiNaC::ex_to<GiNaC::numeric>(eq.subs(x == i * 1e-3).evalf()).to_double();
        }
        benchmark::DoNotOptimize(y.data());
    }
    state.SetItemsProcessed(state.iterations() * POINTS_TRACE);
}
BENCHMARK(BM_GiNaC_EvalTableau)->Unit(benchmark::kMillisecond);

// Évaluation d'une grande expression à sous-expressions partagées :
// dérivée 8e de exp(sin x) * x^2 (compilée automatiquement par EquationClassique)
static void BM_SymAlgo_EvalGrandeExpression(benchmark::State& state) {
    const ExprPtr x = var("x");
    EquationClassique d(ast_exp(ast_sin(x)) * ast_pow(x, 2.0));
    for (int k = 0; k < 8; ++k) d = d.derivee();
    double v = 0.3;
    for (auto _ : state) {
        double res = d.eval(v);
        benchmark::DoNotOptimize(res);
        v += 1e-9;
    }
}
BENCHMARK(BM_SymAlgo_EvalGrandeExpression);

static void BM_GiNaC_EvalGrandeExpression(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex d = (GiNaC::exp(GiNaC::sin(x)) * GiNaC::pow(x, 2)).diff(x, 8);
    double v = 0.3;
    for (auto _ : state) {
        double res = GiNaC::ex_to<GiNaC::numeric>(d.subs(x == v).evalf()).to_double();
        benchmark::DoNotOptimize(res);
        v += 1e-9;
    }
}
BENCHMARK(BM_GiNaC_EvalGrandeExpression)->Unit(benchmark::kMicrosecond);

// ============================================================================
// 3. Dérivation (première dérivée, puis dérivée 10e)
// ============================================================================

static void BM_SymAlgo_Derivee(benchmark::State& state) {
    const EquationClassique eq(fSymAlgo());
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        EquationClassique derivee = eq.derivee();
        benchmark::DoNotOptimize(derivee);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_SymAlgo_Derivee);

static void BM_GiNaC_Derivee(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex eq = fGinac(x);
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        GiNaC::ex d = eq.diff(x);
        benchmark::DoNotOptimize(d);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_GiNaC_Derivee);

// d^10/dx^10 [exp(sin x) * x^2] : mesure la croissance des expressions
static void BM_SymAlgo_Derivee10(benchmark::State& state) {
    const ExprPtr x = var("x");
    const EquationClassique eq(ast_exp(ast_sin(x)) * ast_pow(x, 2.0));
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        EquationClassique d = eq;
        for (int k = 0; k < 10; ++k) d = d.derivee();
        benchmark::DoNotOptimize(d);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_SymAlgo_Derivee10)->Unit(benchmark::kMillisecond);

static void BM_GiNaC_Derivee10(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex eq = GiNaC::exp(GiNaC::sin(x)) * GiNaC::pow(x, 2);
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        GiNaC::ex d = eq.diff(x, 10);
        benchmark::DoNotOptimize(d);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_GiNaC_Derivee10)->Unit(benchmark::kMillisecond);

// ============================================================================
// 4. Collecte des termes semblables : sum_{i=1..50} (i x^(i mod 7) + i sin x)
// ============================================================================

static void BM_SymAlgo_Collecte(benchmark::State& state) {
    const ExprPtr x = var("x");
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        ExprPtr somme = cst(0.0);
        for (int i = 1; i <= TERMES_COLLECTE; ++i) {
            somme = somme + cst(i) * ast_pow(x, i % 7) + cst(i) * ast_sin(x);
        }
        ExprPtr resultat = somme->simplifier();
        benchmark::DoNotOptimize(resultat);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_SymAlgo_Collecte);

static void BM_GiNaC_Collecte(benchmark::State& state) {
    GiNaC::symbol x("x");
    const long long debut = g_octetsAlloues;
    for (auto _ : state) {
        GiNaC::ex somme = 0;
        for (int i = 1; i <= TERMES_COLLECTE; ++i) {
            somme = somme + i * GiNaC::pow(x, i % 7) + i * GiNaC::sin(x);
        }
        benchmark::DoNotOptimize(somme);
    }
    rapporterAllocations(state, debut);
}
BENCHMARK(BM_GiNaC_Collecte);

// ============================================================================
// 5. Développement en série de exp(sin x) à l'ordre 10
// ============================================================================

static void BM_SymAlgo_Serie(benchmark::State& state) {
    const ExprPtr x = var("x");
    const ExprPtr f = ast_exp(ast_sin(x));
    for (auto _ : state) {
        ExprPtr dl = f->DL(0.0, 10);
        benchmark::DoNotOptimize(dl);
    }
}
BENCHMARK(BM_SymAlgo_Serie);

static void BM_GiNaC_Serie(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex f = GiNaC::exp(GiNaC::sin(x));
    for (auto _ : state) {
        GiNaC::ex s = f.series(x == 0, 11);
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(BM_GiNaC_Serie);

// ============================================================================
// 6. Empreinte mémoire du résultat : dérivée 6e de exp(sin x) * x^2
// ============================================================================

static void BM_SymAlgo_Memoire(benchmark::State& state) {
    const ExprPtr x = var("x");
    const EquationClassique eq(ast_exp(ast_sin(x)) * ast_pow(x, 2.0));
    double octets = 0.0;
    for (auto _ : state) {
        const long long avant = g_octetsVivants;
        EquationClassique d = eq;
        for (int k = 0; k < 6; ++k) d = d.derivee();
        octets = static_cast<double>(g_octetsVivants - avant);
        benchmark::DoNotOptimize(d);
    }
    state.counters["octets_resultat"] = octets;
}
BENCHMARK(BM_SymAlgo_Memoire)->Unit(benchmark::kMicrosecond);

static void BM_GiNaC_Memoire(benchmark::State& state) {
    GiNaC::symbol x("x");
    const GiNaC::ex eq = GiNaC::exp(GiNaC::sin(x)) * GiNaC::pow(x, 2);
    double octets = 0.0;
    for (auto _ : state) {
        const long long avant = g_octetsVivants;
        GiNaC::ex d = eq.diff(x, 6);
        octets = static_cast<double>(g_octetsVivants - avant);
        benchmark::DoNotOptimize(d);
    }
    state.counters["octets_resultat"] = octets;
}
BENCHMARK(BM_GiNaC_Memoire)->Unit(benchmark::kMicrosecond);

// ============================================================================
// 7. Équations différentielles : matrice compagnon (SymAlgo++ seul)
// ============================================================================

static void BM_SymAlgo_EqDiff_MatriceCompagnon(benchmark::State& state) {
    const int ordre = static_cast<int>(state.range(0));
    // y^(n) + 2*y^(n-1) + ... + (n+1)*y = 0
    EquationDifferentielle eq;
    for (int i = 0; i <= ordre; ++i) eq.ajouterTerme(i, static_cast<double>(ordre - i + 1));
    for (auto _ : state) {
        Eigen::MatrixXd A = eq.getMatriceCompagnon();
        benchmark::DoNotOptimize(A);
    }
}
BENCHMARK(BM_SymAlgo_EqDiff_MatriceCompagnon)->DenseRange(1, 30, 2);

BENCHMARK_MAIN();
