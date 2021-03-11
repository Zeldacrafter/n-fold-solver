#include <bits/stdc++.h>
#include <Eigen/Dense>

#include "template.cc"

using namespace Eigen;

using vec = VectorXi;
using mat = MatrixXi;

void readVec(vec& v, istream& inp) {
    F0R (i, SZ(v)) inp >> v(i);
}

void readMatr(mat& m, istream& inp) {
    F0R (r, m.rows()) {
        F0R (c, m.cols()) {
            inp >> m(r, c);
        }
    }
}

int bruteForce(istream& inp) {

    int n, r, s, t; 
    inp >> n >> r >> s >> t;

    vec l(n * t), u(n * t), b(r + n * s), c(n * t);
    readVec(l, inp); readVec(u, inp); readVec(b, inp); readVec(c, inp);
    vector<mat> as(n, mat(r, t)), bs(n, mat(s, t));
    for (auto& m : as) readMatr(m, inp);
    for (auto& m : bs) readMatr(m, inp);

    mat A(r + n * s, t * n);
    F0R (i, n) {
        F0R (rr, r)
            F0R (cc, t)
                A(rr, i * t + cc) = as[i](rr, cc);
        F0R (rr, s)
            F0R (cc, t)
                A(r + i * s + rr, i * t + cc) = bs[i](rr, cc);
    }

    cout << A << endl;
    int best = numeric_limits<int>::min();

    vec x = l;
    bool ok = true;
    while (ok) {
        if (A * x == b)
            ckmax(best, c.dot(x));
        F0R (i, n * t + 1) {
            if (i == n * t) {
                ok = false;
                break;
            }
            if (++x(i) > u(i)) {
                x(i) = l(i);
            } else {
                break;
            }
        }
    }

    return best;
}

