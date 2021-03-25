#include <Eigen/Dense>

#include "utils.hh"
#include "solver_class.hh"

int main() {
#if not defined(N_NFOLD) || not defined(R_NFOLD) || not defined(S_NFOLD) || not defined(T_NFOLD)
    std::cout << "Variables N_NFOLD, R_NFOLD, S_NFOLD and T_NFOLD must be set at compile-time." << std::endl;
    exit(1);
#else
    using namespace std;

    int n, r, s, t;
    cin >> n >> r >> s >> t;
    assert(N_NFOLD == n);
    assert(R_NFOLD == r);
    assert(S_NFOLD == s);
    assert(T_NFOLD == t);

    n_fold<int, N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD> nfold;
    cin >> nfold;

    auto maybeRes = n_fold_solver<int, N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD>(nfold).solve();
    if(maybeRes) {
        cout << (*maybeRes).second << std::endl << (*maybeRes).first;
    } else {
        cout << "No solution exists";
        return 1;
    }
    return 0;
#endif
}
