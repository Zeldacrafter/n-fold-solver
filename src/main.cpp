#ifdef USING_BOOST
#include <signal.h>
#include <boost/stacktrace.hpp>
#endif

#include <Eigen/Dense>

#include "utils.hh"
#include "solver_class.hh"

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

    StaticNFold<int> nfold(n, r, s, t);
    Vec<int> initSol(n*t);
    if(normal) {
        cin >> nfold;
    } else {
        cin >> nfold.c >> nfold.l >> nfold.u >> nfold.as[0] >> nfold.bs[0] >> nfold.b >> initSol;
        for(int i = 1; i < n; ++i) {
            nfold.as[i] = nfold.as[0];
            nfold.bs[i] = nfold.bs[0];
        }
    }

    //cout << "Input:\n" << nfold << endl;

    if(normal) {
        auto maybeRes = StaticSolver<int>(nfold).solve();
        if(maybeRes) {
            cout << (*maybeRes).second << std::endl << (*maybeRes).first;
        } else {
            cout << "No solution exists";
            return 1;
        }
    } else {
        auto res = StaticSolver<int>(nfold).solve(initSol);
        cout << res.second << std::endl << res.first;
    }
    return 0;
#endif
}
