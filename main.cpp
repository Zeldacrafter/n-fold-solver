#include "template.h"
#include "bruteforce.cpp"
#include "Solver.cpp"

#include <signal.h>
#include <boost/stacktrace.hpp>

void handler(int) {
    std::cerr << boost::stacktrace::stacktrace() << std::endl;
    _exit(1);
}

int main() {
    ::signal(SIGABRT, handler);
    ::signal(SIGSEGV, handler);
    using namespace std;

    int n, r, s, t;
    cin >> n >> r >> s >> t;

    NFold<int> nfold(n, r, s, t);
    cin >> nfold;

    //cout << constructAInit(nfold) << endl;
    //cout << nfold << endl;
    //cout << dvar(pp(nfold.l), pp(nfold.u)) << endl;
    //cout << dvar(pp(nfold.b), pp(nfold.c)) << endl;

    //cout << "Best: " << pp(bruteForceBest(nfold)) << endl;
    //cout << "Worst: " << pp(bruteForceWorst(nfold)) << endl;

    //cout << Solver(nfold).solve() << endl;
    dout << nfold << endl
         << dvar(pp(nfold.l), pp(nfold.u)) << endl
         << dvar(pp(nfold.b), pp(nfold.c)) << endl;

    dout << pp(Solver(nfold).solve()) << endl;
    /*
    auto [aInit, sol] = constructAInit(nfold);

    dout << aInit << endl
         << dvar(pp(aInit.l), pp(aInit.u)) << endl
         << dvar(pp(aInit.b), pp(aInit.c)) << endl
         << dvar(pp(sol)) << endl;

    assert(aInit*sol == aInit.b);

    Solver solver(aInit);
    auto [initSol, wgt] = solver.solve(sol);
    cout << dvar(pp(initSol), wgt) << endl;
    assert(aInit*initSol == aInit.b);
    cout << Solver(nfold).solve(initSol) << endl;
     */
    return 0;
}