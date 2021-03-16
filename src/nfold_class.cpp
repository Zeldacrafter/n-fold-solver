#pragma once

#include <Eigen/Dense>
#include <tsl/hopscotch_hash.h>

#include "template.h"


template<typename T>
std::istream& operator>>(std::istream& inp, Vec<T>& x) {
    F0R(i, SZ(x)) inp >> x(i);
    return inp;
}

template<typename T>
std::istream& operator>>(std::istream& inp, Mat<T>& x) {
    F0R(r, x.rows())
    F0R(c, x.cols())
    inp >> x(r, c);
    return inp;
}

template<typename T>
class NFold {
public:
    size_t n, r, s, t;
    Vec<T> l, u, b, c;
    std::vector<Mat<T>> as, bs;

    NFold() = default;

    NFold(size_t _n, size_t _r, size_t _s, size_t _t)
            : n{_n}, r{_r}, s{_s}, t{_t},
              l(n*t), u(n*t), b(r + n*s), c(n*t),
              as(n, Mat<T>(r, t)), bs(n, Mat<T>(s, t)) { }

    // Function call operator accesses the matrix element in that position.
    T operator()(size_t row, size_t col) {
        assert(static_cast<size_t>(row) < r + n*s);
        assert(static_cast<size_t>(col) < n*t);

        if(row < r) {
            return as[col/t](row, col%t); // A block
        } else {
            size_t blockRow = (row - r)/s;
            if(col >= blockRow*t && col < (blockRow + 1)*t) {
                return bs[blockRow]( (row-r)%s, col%t ); // B block
            } else {
                return 0; // 0 block
            }
        }
    }

    // multiplication with a column vector
    Vec<T> operator*(const Vec<T>& x) const {
        assert(static_cast<size_t>(SZ(x)) == n*t);

        Vec<T> res(r + n*s);
        res.fill(0);
        F0R(i, n) res.head(r) += as[i] * x.segment(i * t, t);
        F0R(i, n) res.segment(r + i*s, s) = bs[i] * x.segment(i * t, t);
        return res;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& outp, NFold<U> x) {
        outp << dvar(x.n, x.r, x.s, x.t) << std::endl
             << dvar(pp(x.l), pp(x.u)) << std::endl
             << dvar(pp(x.c), pp(x.b)) << std::endl
             << "A:" << std::endl;
        F0R(rr, x.r + x.n*x.s) {
            F0R(cc, x.n*x.t)
                outp << x(rr, cc) << ' ';
            outp << std::endl;
        }
        return outp;
    }

    template <typename U>
    friend std::istream& operator>>(std::istream& inp, NFold<U>& x) {
        return inp >> x.l >> x.u >> x.b >> x.c >> x.as >> x.bs;
    }
};


