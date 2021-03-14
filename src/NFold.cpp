#pragma once

#include <Eigen/Dense>
#include "template.h"
#include "../third-party/hopscotch-map/include/tsl/hopscotch_map.h"


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

    NFold() {};

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
             << "A: " << std::endl;
        F0R(rr, x.r + x.n*x.s) {
            F0R(cc, x.n*x.t)
                outp << x(rr, cc) << ' ';
            outp << std::endl;
        }
        return outp;
    }

    T getDelta() {
        T delta = std::numeric_limits<T>::min();
        for(Mat<T>& x : as)
            ckmax(delta, x.maxCoeff());
        for(Mat<T>& x : bs)
            ckmax(delta, x.maxCoeff());
        return delta;
    }

    template <typename U>
    friend std::istream& operator>>(std::istream& inp, NFold<U>& x) {
        return inp >> x.l >> x.u >> x.b >> x.c >> x.as >> x.bs;
    }
};


// Constructs an nFold as described by Jansens paper in chapter 4.
// This is used to find an initial solution for the original input nfold.
template <typename T>
std::pair<NFold<T>, Vec<T>> constructAInit(const NFold<T>& x) {
    const int n = x.n, r = x.r, s = x.s, t = x.t;

    NFold<T> res(n, r, s, t + r + s);

    //Construct new matrix
    F0R(i, n) {
        res.as[i].block(0, 0,     r, t) = x.as[i];
        if(!i) res.as[i].block(0, t, r, r).setIdentity();
        else   res.as[i].block(0, t, r, r).setZero();
        res.as[i].block(0, t + r, r, s).setZero();

        res.bs[i].block(0, 0,     s, t) = x.bs[i];
        res.bs[i].block(0, t,     s, r).setZero();
        res.bs[i].block(0, t + r, s, s).setIdentity();
    }

    // Construct new righthand side
    res.b = x.b - x*x.l;

    // Construct upper and lower bound
    F0R(i, n) {
        res.l.segment(i*(t + r + s), t) = (x.u - x.l).segment(i*t, t);
        res.u.segment(i*(t + r + s), t) = (x.u - x.l).segment(i*t, t);
        if(!i) {
            res.l.segment(i*(t + r + s) + t, r) = res.b.segment(0, r);
            res.u.segment(i*(t + r + s) + t, r) = res.b.segment(0, r);
        } else {
            res.l.segment(i*(t + r + s) + t, r).setZero();
            res.u.segment(i*(t + r + s) + t, r).setZero();
        }
        res.l.segment(i*(t + r + s) + t + r, s) = res.b.segment(r + i*s, s);
        res.u.segment(i*(t + r + s) + t + r, s) = res.b.segment(r + i*s, s);
    }
    res.l = res.l.array().min(0);
    res.u = res.u.array().max(0);

    // Construct cost vector
    res.c.setZero();
    res.c.segment(t, r).setOnes();
    F0R(j, r) // Identity matrix for A_1
       if(res.b(j) >= 0)
           res.c(t + j) = -1;
    F0R(i, n) { // Identity matrix for B_i
        res.c.segment(t + r + i*(t + r + s), s).setOnes();
        F0R(j, s)
            if(res.b(r + i*s + j) >= 0)
                res.c(t + r + i*(t + r + s) + j) = -1;
    }

    // Construct init sol
    Vec<T> initSol = Vec<T>::Zero((t + r + s)*n);
    initSol.segment(t, r) = res.b.segment(0, r);
    F0R(i, n) {
        initSol.segment((t + r + s)*i + t + r, s) = res.b.segment(r + i*s, s);
    }

    return std::make_pair(res, initSol);
}
