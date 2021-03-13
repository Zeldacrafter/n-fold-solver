#pragma once

#include <Eigen/Dense>
#include <iostream>

#include "template.h"
#include "NFold.cpp"

std::pair<Vec<int>, int> bruteForceBest(const NFold<int>& nfold) {
    int best = std::numeric_limits<int>::min();
    Vec<int> bestVec(1);

    Vec<int> x = nfold.l;
    bool ok = true;
    while (ok) {
        if (nfold * x == nfold.b && ckmax(best, x.dot(nfold.c))) {
            // new best found
            bestVec = x;
        }
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

    return std::make_pair(bestVec, best);
}

std::pair<Vec<int>, int> bruteForceWorst(const NFold<int>& nfold) {
    int worst = std::numeric_limits<int>::max();
    Vec<int> worstVec(1);

    Vec<int> x = nfold.l;
    bool ok = true;
    while (ok) {
        if (nfold * x == nfold.b && ckmin(worst, x.dot(nfold.c))) {
            // new worst found
            worstVec = x;
        }
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

    return std::make_pair(worstVec, worst);
}