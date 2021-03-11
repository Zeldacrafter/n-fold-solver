#include <bits/stdc++.h>
#include <Eigen/Dense>

#include "template.hh"
#include "bruteforce.cc"


int main() {
    NFold<int> nfold(std::cin);
    std::cout << bruteForce(nfold) << std::endl;
    return 0;
}
