#ifdef USING_BOOST
#include <signal.h>
#include <boost/stacktrace.hpp>
#endif

#include "src/template.h"
#include "src/bruteforce.cpp"
#include "src/static_solver_class.cpp"

#ifdef USING_BOOST
void handler(int) {
    std::cerr << boost::stacktrace::stacktrace() << std::endl;
    _exit(1);
}
#endif

int main(int argc, char* argv[]) {
#if not defined(N_NFOLD) || not defined(R_NFOLD) || not defined(S_NFOLD) || not defined(T_NFOLD)
    std::cout << "Variables N_NFOLD, R_NFOLD, S_NFOLD and T_NFOLD must be set at compile-time." << std::endl;
    exit(1);
#else
#ifdef USING_BOOST
    ::signal(SIGABRT, handler);
    ::signal(SIGSEGV, handler);
#endif
    using namespace std;

    // No args or first arg starts with 'n':
    //   Normal mode: Input parsing as specified in the task
    // Otherwise:
    //   Input parsing for input in https://github.com/benjamincarman/cpp-nfold-ILP
    bool normal = argc == 1 || argv[1][0] == 'n';

    int n, r, s, t;
    cin >> n >> r >> s >> t;
    assert(N_NFOLD == n);
    assert(R_NFOLD == r);
    assert(S_NFOLD == s);
    assert(T_NFOLD == t);

    StaticNFold<int, N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD> nfold;
    Vec<int> initSol(n*t);
    if(normal) {
        cin >> nfold;
    } else {
        cin >> nfold.c >> nfold.l >> nfold.u >> nfold.as[0] >> nfold.bs[0] >> nfold.b >> initSol;
        FOR(i, 1, n) {
            nfold.as[i] = nfold.as[0];
            nfold.bs[i] = nfold.bs[0];
        }
    }

    cout << "Input:\n" << nfold << endl;

    auto [res, resWgt] =
    normal ? *static_solver::StaticSolver<int, N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD>(nfold).solve()
           :  static_solver::StaticSolver<int, N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD>(nfold).solve(initSol);

    cout << "Solution found:\n" << dvar(pp(res), resWgt) << endl;
    return 0;
#endif
}