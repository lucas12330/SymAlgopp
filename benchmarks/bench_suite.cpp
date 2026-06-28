#include <benchmark/benchmark.h>
#include <ginac/ginac.h>
#include "EquationClassique.hpp"
#include "EquationDifferentielle.hpp"
#include "ASTNode.hpp"

// ---------------------------------------------------------
// 1. Benchmark: Evaluation d'Equations Classiques
// ---------------------------------------------------------

// SymAlgo++ : f(x) = sin(x)*cos(x) + x^3 - 2*x
static void BM_SymAlgo_Eval(benchmark::State& state) {
    ExprPtr expr = (ast_sin(var("x")) * ast_cos(var("x"))) + ast_pow(var("x"), 3.0) - (cst(2.0) * var("x"));
    EquationClassique eq(expr);

    for (auto _ : state) {
        double res = eq.eval(5.0);
        benchmark::DoNotOptimize(res);
    }
}
BENCHMARK(BM_SymAlgo_Eval);

// GiNaC : f(x) = sin(x)*cos(x) + x^3 - 2*x
static void BM_GiNaC_Eval(benchmark::State& state) {
    GiNaC::symbol x("x");
    GiNaC::ex eq = GiNaC::sin(x) * GiNaC::cos(x) + GiNaC::pow(x, 3) - 2 * x;

    for (auto _ : state) {
        GiNaC::ex res = eq.subs(x == 5.0).evalf();
        double d_res = GiNaC::ex_to<GiNaC::numeric>(res).to_double();
        benchmark::DoNotOptimize(d_res);
    }
}
BENCHMARK(BM_GiNaC_Eval);

// ---------------------------------------------------------
// 2. Benchmark: Dérivation d'Equations Classiques
// ---------------------------------------------------------

static void BM_SymAlgo_Derivee(benchmark::State& state) {
    ExprPtr expr = (ast_sin(var("x")) * ast_cos(var("x"))) + ast_pow(var("x"), 3.0) - (cst(2.0) * var("x"));
    EquationClassique eq(expr);

    for (auto _ : state) {
        EquationClassique* derivee = eq.derivee();
        benchmark::DoNotOptimize(derivee);
        delete derivee; // On libère la mémoire pour éviter les fuites pendant le bench
    }
}
BENCHMARK(BM_SymAlgo_Derivee);

static void BM_GiNaC_Derivee(benchmark::State& state) {
    GiNaC::symbol x("x");
    GiNaC::ex eq = GiNaC::sin(x) * GiNaC::cos(x) + GiNaC::pow(x, 3) - 2 * x;

    for (auto _ : state) {
        GiNaC::ex d = eq.diff(x);
        benchmark::DoNotOptimize(d);
    }
}
BENCHMARK(BM_GiNaC_Derivee);

// ---------------------------------------------------------
// 3. Benchmark: Equations Différentielles (Matrice Compagnon)
// ---------------------------------------------------------
// Evalue le temps de création et de structuration du système
// matriciel pour des équations différentielles d'ordre N
static void BM_SymAlgo_EqDiff_MatriceCompagnon(benchmark::State& state) {
    int ordre = state.range(0);
    
    // Création de l'équation : y^(n) + 2*y^(n-1) + ... + n*y = 0
    EquationDifferentielle eq;
    for(int i = 0; i <= ordre; ++i) {
        eq.ajouterTerme(i, (double)(ordre - i + 1));
    }
    
    for (auto _ : state) {
        Eigen::MatrixXd A = eq.getMatriceCompagnon();
        benchmark::DoNotOptimize(A);
    }
}
// Tester pour les ordres 1, 3, 5, 7, 9, 11, 13, 15
BENCHMARK(BM_SymAlgo_EqDiff_MatriceCompagnon)->DenseRange(1, 15, 2);

// ---------------------------------------------------------
// 4. Benchmark: Limite en x = 0 (sin(x)/x)
// ---------------------------------------------------------
static void BM_SymAlgo_Limite(benchmark::State& state) {
    ExprPtr expr = ast_sin(var("x")) / var("x");
    EquationClassique eq(expr);

    for (auto _ : state) {
        EquationClassique* lim = eq.limite(0.0);
        benchmark::DoNotOptimize(lim);
        delete lim;
    }
}
BENCHMARK(BM_SymAlgo_Limite);

static void BM_GiNaC_Limite(benchmark::State& state) {
    GiNaC::symbol x("x");
    GiNaC::ex eq = GiNaC::sin(x) / x;

    for (auto _ : state) {
        GiNaC::ex s = eq.series(x==0, 1);
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(BM_GiNaC_Limite);

// ---------------------------------------------------------
// 5. Benchmark: Développement Limité (Ordre 5)
// ---------------------------------------------------------
static void BM_SymAlgo_DL(benchmark::State& state) {
    ExprPtr expr = ast_exp(var("x")) * ast_sin(var("x"));
    EquationClassique eq(expr);

    for (auto _ : state) {
        EquationClassique* dl = eq.DL(0.0, 5);
        benchmark::DoNotOptimize(dl);
        delete dl;
    }
}
BENCHMARK(BM_SymAlgo_DL);

static void BM_GiNaC_DL(benchmark::State& state) {
    GiNaC::symbol x("x");
    GiNaC::ex eq = GiNaC::exp(x) * GiNaC::sin(x);

    for (auto _ : state) {
        GiNaC::ex s = eq.series(x==0, 6); // order 5 implies terms up to O(x^6)
        benchmark::DoNotOptimize(s);
    }
}
BENCHMARK(BM_GiNaC_DL);

// ---------------------------------------------------------
// 6. Benchmark: Intégrale Définie sur [0, 5]
// ---------------------------------------------------------
static void BM_SymAlgo_IntegraleDefinie(benchmark::State& state) {
    ExprPtr expr = ast_pow(var("x"), 3.0) + (cst(2.0) * var("x"));
    EquationClassique eq(expr);

    for (auto _ : state) {
        double val = eq.integraleDefinie(0.0, 5.0);
        benchmark::DoNotOptimize(val);
    }
}
BENCHMARK(BM_SymAlgo_IntegraleDefinie);

static void BM_GiNaC_IntegraleDefinie(benchmark::State& state) {
    GiNaC::symbol x("x");
    GiNaC::ex eq = GiNaC::pow(x, 3) + 2 * x;

    for (auto _ : state) {
        GiNaC::ex int_eq = GiNaC::integral(x, 0, 5, eq);
        GiNaC::ex val = int_eq.eval_integ();
        benchmark::DoNotOptimize(val);
    }
}
BENCHMARK(BM_GiNaC_IntegraleDefinie);

BENCHMARK_MAIN();
