#pragma once

#include <Eigen/StdVector>
#include <tsl/hopscotch_hash.h>

#include "template.h"


template<typename U, int S>
std::istream& operator>>(std::istream& inp, sVec<U, S>& x) {
    F0R(i, SZ(x)) inp >> x(i);
    return inp;
}

template<typename T, int S1, int S2>
std::istream& operator>>(std::istream& inp, sMat<T, S1, S2>& x) {
    F0R(r, x.rows())
        F0R(c, x.cols())
            inp >> x(r, c);
    return inp;
}

template<typename U, int N, int R, int S, int T>
class StaticNFold {
public:
    sVec<U, N*T> l, u, c;
    sVec<U, R + N*S> b;
    std::array<sMat<U, R, T>, N> as;
    std::array<sMat<U, S, T>, N> bs;

    StaticNFold() = default;

    // Function call operator accesses the matrix element in that position.
    U operator()(size_t row, size_t col) {
        assert(row < R + N*S);
        assert(col < N*T);

        if(row < R) {
            return as[col/T](row, col%T); // A block
        } else {
            size_t blockRow = (row - R)/S;
            if(col >= blockRow*T && col < (blockRow + 1)*T) {
                return bs[blockRow]( (row-R)%S, col%T ); // B block
            } else {
                return 0; // 0 block
            }
        }
    }

    // multiplication with a column vector
    sVec<U, R + N*S> operator*(const sVec<U, N*T>& x) const {
        sVec<U, R + N*S> res = sVec<U, R + N*S>::Zero();
        F0R(i, N) res.head(R) += as[i] * x.segment(i * T, T);
        F0R(i, N) res.segment(R + i*S, S) = bs[i] * x.segment(i * T, T);
        return res;
    }

    friend std::ostream& operator<<(std::ostream& outp, StaticNFold<U, N, R, S, T> x) {
        outp << dvar(N, R, S, T) << std::endl
             << dvar(pp(x.l), pp(x.u)) << std::endl
             << dvar(pp(x.c), pp(x.b)) << std::endl
             << "A:" << std::endl;
        F0R(rr, R + N*S) {
            F0R(cc, N*T)
                outp << x(rr, cc) << ' ';
            outp << std::endl;
        }
        return outp;
    }

    friend std::istream& operator>>(std::istream& inp, StaticNFold<U, N, R, S, T>& x) {
        return inp >> x.l >> x.u >> x.b >> x.c >> x.as >> x.bs;
    }
};

