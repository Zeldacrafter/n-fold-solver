#include <Eigen/Dense>

#include "template.hh"
#include "bruteforce.cc"

int main() {
    int n, r, s, t;
    std::cin >> n >> r >> s >> t;

    NFold<int> nfold(n, r, s, t);
    std::cin >> nfold;

    std::cout << nfold;

    std::cout << bruteForce(nfold) << std::endl;
    return 0;
}
