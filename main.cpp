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

    cout << nfold << endl;
    cout << dvar(pp(nfold.l), pp(nfold.u)) << endl;
    cout << dvar(pp(nfold.b), pp(nfold.c)) << endl;

    cout << "Best: " << pp(bruteForceBest(nfold)) << endl;
    cout << "Worst: " << pp(bruteForceWorst(nfold)) << endl;

    cout << Solver(nfold).solve() << endl;
    return 0;
}