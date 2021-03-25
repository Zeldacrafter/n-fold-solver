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
#ifdef USING_BOOST
    ::signal(SIGABRT, handler);
    ::signal(SIGSEGV, handler);
#endif
    using namespace std;

    int n, r, s, t;
    cin >> n >> r >> s >> t;

    n_fold<int> nfold(n, r, s, t);
    cin >> nfold;

    auto maybeRes = n_fold_solver<int>(nfold).solve();
    if(maybeRes) {
        cout << (*maybeRes).second << std::endl << (*maybeRes).first;
    } else {
        cout << "No solution exists";
        return 1;
    }
    return 0;
}
