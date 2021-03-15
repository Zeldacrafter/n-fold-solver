#include <benchmark/benchmark.h>
#include <iostream>
#include <fstream>

#include "nfold_class.cpp"
#include "solver_class.cpp"

static void BM_Solution_Small(benchmark::State& state) {
    static bool setupDone = false;
    static NFold<int> nfold;
    if(!setupDone) {
        std::string path = "input/max-small";

        std::ifstream inputFile(path + ".in");
        int n, r, s, t;
        inputFile >> n >> r >> s >> t;
        nfold = NFold<int>(n, r, s, t);
        inputFile >> nfold;

        setupDone = true;
    }

    for (auto _ : state)
        solver::Solver(nfold).solve();
}
BENCHMARK(BM_Solution_Small)->Unit(benchmark::kMillisecond);

static void BM_Solution_Large(benchmark::State& state) {
    static bool setupDone = false;
    static NFold<int> nfold;
    if(!setupDone) {
        std::string path = "input/max-large";

        std::ifstream inputFile(path + ".in");
        int n, r, s, t;
        inputFile >> n >> r >> s >> t;
        nfold = NFold<int>(n, r, s, t);
        inputFile >> nfold;

        setupDone = true;
    }

    for (auto _ : state)
        solver::Solver(nfold).solve();
}
BENCHMARK(BM_Solution_Large)->Unit(benchmark::kMillisecond);

static void BM_Solution_Largest(benchmark::State& state) {
    static bool setupDone = false;
    static NFold<int> nfold;
    if(!setupDone) {
        std::string path = "input/input-huge";

        std::ifstream inputFile(path + ".in");
        int n, r, s, t;
        inputFile >> n >> r >> s >> t;
        nfold = NFold<int>(n, r, s, t);
        inputFile >> nfold;

        setupDone = true;
    }

    for (auto _ : state)
        solver::Solver(nfold).solve();
}
BENCHMARK(BM_Solution_Largest)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();