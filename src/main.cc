#include <bits/stdc++.h>
#include <Eigen/Dense>

#include "bruteforce.cc"

template<typename T>
using mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>

int main() {
    std::ifstream inp("input/max-small.in");
    int content;
    cout << "Solution: " << bruteForce(inp) << endl;
    return 0;
}
