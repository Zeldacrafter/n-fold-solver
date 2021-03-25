#include <Eigen/Dense>

#include "utils.hh"
#include "solver_class.hh"

int main(int argc, char* argv[]) {
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
