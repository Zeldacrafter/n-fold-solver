#include "src/template.h"
#include "src/bruteforce.cpp"
#include "src/Solver.cpp"

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

    dout << nfold << endl;

    dout << pp(*solver::Solver(nfold).solve()) << endl;
    return 0;
}
