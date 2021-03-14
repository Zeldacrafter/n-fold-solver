#pragma once

#include <Eigen/Dense>
#include <iostream>

#include "template.h"
#include "NFold.cpp"

template <typename T>
std::pair<Vec<T>, T> bruteForceBest(const NFold<T>& nfold) {
    int best = std::numeric_limits<T>::min();
    Vec<T> bestVec(1);

    Vec<T> x = nfold.l;
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

template <typename T>
std::pair<Vec<T>, T> bruteForceWorst(const NFold<T>& nfold) {
    T worst = std::numeric_limits<T>::max();
    Vec<T> worstVec(1);

    Vec<T> x = nfold.l;
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