#include <signal.h>
#include <boost/stacktrace.hpp>

#include "src/template.h"
#include "src/bruteforce.cpp"
#include "src/solver_class.cpp"

void handler(int) {
    std::cerr << boost::stacktrace::stacktrace() << std::endl;
    _exit(1);
}

int main(int argc, char* argv[]) {
    ::signal(SIGABRT, handler);
    ::signal(SIGSEGV, handler);
    using namespace std;

    // No args or first arg starts with 'n':
    //   Normal mode: Input parsing as specified in the task
    // Otherwise:
    //   Input parsing for input in https://github.com/benjamincarman/cpp-nfold-ILP
    bool normal = argc == 1 || argv[1][0] == 'n';

    int n, r, s, t;
    cin >> n >> r >> s >> t;

    NFold<int> nfold(n, r, s, t);
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
        normal ? *solver::Solver(nfold).solve()
               :  solver::Solver(nfold).solve(initSol);

    cout << "Solution found:\n" << dvar(pp(res), resWgt) << endl;
    return 0;
}
