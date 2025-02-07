#ifndef N_FOLD_NFOLD_CLASS_HH
#define N_FOLD_NFOLD_CLASS_HH

#include <Eigen/StdVector>
#include <tsl/hopscotch_hash.h>

#include "utils.hh"

/**
 * Class that represents an n-fold ILP instance.
 * @tparam N Size parameter for the n-fold.
 * @tparam R Size parameter for the n-fold.
 * @tparam S Size parameter for the n-fold.
 * @tparam T Size parameter for the n-fold.
 * @tparam U Value type in the n-fold.
 */
template<int N, int R, int S, int T, typename U = int>
class n_fold {
public:
    sVec<U, N*T> l, u, c;
    sVec<U, R + N*S> b;
    std::array<sMat<U, R, T>, N> as;
    std::array<sMat<U, S, T>, N> bs;

    n_fold() = default;

    /**
     * Return the largest entry in the matrix.
     * @return The largest entry in the matrix.
     */
    U getDelta() {
        U delta = std::numeric_limits<U>::min();
        for(auto& x : as) {
            delta = std::max(delta, x.template lpNorm<Eigen::Infinity>());
        }
        for(auto& x : bs) {
            delta = std::max(delta, x.template lpNorm<Eigen::Infinity>());
        }
        return delta;
    }

    /**
     * Access the element at a position in that matrix.
     * @param row The row of the wanted entry.
     * @param col The column of the wanted entry.
     * @return The value of the wanted entry.
     */
    U operator()(size_t row, size_t col) {
        assertm(row >= 0 && row < R + N*S, "The specified row of the n-fold is invalid.");
        assertm(col >= 0 && col < N*T, "The specified column of the n-fold is invalid.");

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

    /**
     * Multiplication of a vector with a column vector.
     * @param x The column vector.
     * @return The result of the multiplication.
     */
    sVec<U, R + N*S> operator*(const sVec<U, N*T>& x) const {
        sVec<U, R + N*S> res = sVec<U, R + N*S>::Zero();
        for(int i = 0; i < N; ++i) {
            res.head(R) += as[i] * x.segment(i * T, T);
            res.segment(R + i*S, S) = bs[i] * x.segment(i * T, T);
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& outp, n_fold<N, R, S, T, U> x) {
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

    friend std::istream& operator>>(std::istream& inp, n_fold<N, R, S, T, U>& x) {
        inp >> x.l >> x.u >> x.b >> x.c;
        for(auto& mat : x.as) {
            inp >> mat;
        }
        for(auto& mat : x.bs) {
            inp >> mat;
        }
        return inp;
    }
};

#endif //N_FOLD_NFOLD_CLASS_HH