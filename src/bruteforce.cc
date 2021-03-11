#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <iostream>

#include "template.hh"
#include "nFold.cc"

int bruteForce(const NFold<int>& nfold) {
    int best = std::numeric_limits<int>::min();

    Vec<int> x = nfold.l;
    bool ok = true;
    while (ok) {
        if (nfold * x == nfold.b)
            ckmax(best, x.dot(nfold.c));
        F0R (i, nfold.n * nfold.t + 1) {
            if (i == nfold.n * nfold.t) {
                ok = false;
                break;
            }
            if (++x(i) > nfold.u(i)) {
                x(i) = nfold.l(i);
            } else {
                break;
            }
        }
    }

    return best;
}

