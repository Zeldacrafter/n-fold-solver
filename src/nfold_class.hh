#ifndef N_FOLD_NFOLD_CLASS_HH
#define N_FOLD_NFOLD_CLASS_HH

#include <Eigen/StdVector>
#include <tsl/hopscotch_hash.h>

#include "utils.hh"

template<typename U>
class n_fold {
public:
    size_t N, R, S, T;
    Vec<U> l, u, c;
    Vec<U> b;
    std::vector<Mat<U>> as;
    std::vector<Mat<U>> bs;

    n_fold(size_t n, size_t r, size_t s, size_t t) :
        N{n}, R{r}, S{s}, T{t},
        l(N*T), u(N*T), c(N*T), b(R + N*S),
        as(N, Mat<U>(R, T)), bs(N, Mat<U>(S, T)) {};

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
    Vec<U> operator*(const Vec<U>& x) const {
        Vec<U> res = Vec<U>::Zero(R + N*S);
        for(int i = 0; i < N; ++i) {
            res.head(R) += as[i] * x.segment(i * T, T);
            res.segment(R + i*S, S) = bs[i] * x.segment(i * T, T);
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& outp, n_fold<U> x) {
        outp << "n: " << x.N << ", r: " << x.R << ", s: " << x.S << ", t: " << x.T << std::endl
             << "l: " << x.l << std::endl
             << "u: " << x.u << std::endl
             << "c: " << x.c << std::endl
             << "b: " << x.b << std::endl
             << "A:" << std::endl;
        for(int rr = 0; rr < x.R + x.N*x.S; ++rr) {
            for(int cc = 0; cc < x.N*x.T; ++cc) {
                outp << x(rr, cc) << ' ';
            }
            outp << std::endl;
        }
        return outp;
    }

    friend std::istream& operator>>(std::istream& inp, n_fold<U>& x) {
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