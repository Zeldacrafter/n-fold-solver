#include <Eigen/Dense>

#include "utils.hh"
#include "solver_class.hh"

int main() {
#if not defined(N_NFOLD) || not defined(R_NFOLD) || not defined(S_NFOLD) || not defined(T_NFOLD)
    std::cout << "Variables N_NFOLD, R_NFOLD, S_NFOLD and T_NFOLD must be set at compile-time." << std::endl;
    return 1;
#else
    using namespace std;

    int n, r, s, t;
    cin >> n >> r >> s >> t;
    assertm(n == N_NFOLD, "The size n of the input must match N_NFOLD, the size set at compile-time");
    assertm(r == R_NFOLD, "The size r of the input must match R_NFOLD, the size set at compile-time");
    assertm(s == S_NFOLD, "The size s of the input must match S_NFOLD, the size set at compile-time");
    assertm(t == T_NFOLD, "The size t of the input must match T_NFOLD, the size set at compile-time");

    n_fold<N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD> nfold;
    cin >> nfold;

    auto maybeRes = n_fold_solver<N_NFOLD, R_NFOLD, S_NFOLD, T_NFOLD>(nfold).solve();
    if(maybeRes) {
        cout << "Maximum cost: " << (*maybeRes).second << endl
             << "Result vector: " << (*maybeRes).first << endl;
    } else {
        cout << "No solution exists";
        return 1;
    }
    return 0;
#endif
}
