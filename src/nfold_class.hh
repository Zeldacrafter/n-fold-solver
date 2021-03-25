#ifndef N_FOLD_NFOLD_CLASS_HH
#define N_FOLD_NFOLD_CLASS_HH

#include <Eigen/StdVector>
#include <tsl/hopscotch_hash.h>

#include "utils.hh"

template<typename U, int N, int R, int S, int T>
class n_fold {
public:
    sVec<U, N*T> l, u, c;
    sVec<U, R + N*S> b;
    std::array<sMat<U, R, T>, N> as;
    std::array<sMat<U, S, T>, N> bs;

    n_fold() = default;

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
        for(int i = 0; i < N; ++i) {
            res.head(R) += as[i] * x.segment(i * T, T);
            res.segment(R + i*S, S) = bs[i] * x.segment(i * T, T);
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& outp, n_fold<U, N, R, S, T> x) {
        outp << "n: " << N << ", r: " << R << ", s: " << S << ", t: " << T << std::endl
             << "l: " << x.l << std::endl
             << "u: " << x.u << std::endl
             << "c: " << x.c << std::endl
             << "b: " << x.b << std::endl
             << "A:" << std::endl;
        for(int rr = 0; rr < R + N*S; ++rr) {
            for(int cc = 0; cc < N*T; ++cc) {
                outp << x(rr, cc) << ' ';
            }
            outp << std::endl;
        }
        return outp;
    }

    friend std::istream& operator>>(std::istream& inp, n_fold<U, N, R, S, T>& x) {
        inp >> x.l >> x.u >> x.b >> x.c;
        for(int i = 0; i < x.as.size(); ++i) {
            inp >> x.as[i];
        }
        for(int i = 0; i < x.bs.size(); ++i) {
            inp >> x.bs[i];
        }
        return inp;
    }
};

#endif //N_FOLD_NFOLD_CLASS_HH