#include "template.hh"

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

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
    int n, r, s, t;
    Vec<T> l, u, b, c;
    std::vector<Mat<T>> as, bs;

    NFold(int _n, int _r, int _s, int _t) : n{_n}, r{_r}, s{_s}, t{_t} {}

    NFold(int _n, int _r, int _s, int _t, 
            Vec<T> _l, Vec<T> _u,
            std::vector<Mat<T>> _as, std::vector<Mat<T>> _bs)
            : n{_n}, r{_r}, s{_s}, t{_t}, l{_l}, u{_u}, as{_as}, bs{_bs} { }

    // Function call operator accesses the matrix element in that position.
    T operator()(int row, int col) {
        assert(row >= 0 && row < r + n*s);
        assert(col >= 0 && col < n*t);

        if(row < r) {
            return as[col/t](row, col%t); // A block
        } else {
            int blockRow = (row - r)/s;
            if( col >= blockRow*t && col < (blockRow + 1)*t ) {
                return bs[blockRow]( (row-r)%s, col%t ); // B block
            } else {
                return 0; // 0 block
            }
        }
    }

    // multiplication with a column vector
    Vec<T> operator*(const Vec<T>& x) const {
        assert(SZ(x) == n*t);

        Vec<T> res(r + n*s);
        F0R(i, r + n*s) {
            res[i] = 0;
            if(i < r) {
                // In A block
                int col = 0;
                for(const Mat<T>& m : as) {
                    F0R(_, m.cols()) {
                        res[i] += m(i, col % t) * x(col);
                        ++col;
                    }
                }
            } else {
                // In B block
                int blockRow = (i - r)/s;
                const Mat<T>& m = bs[blockRow];
                F0R(col, t) {
                    res[i] += m((i - r) % s, col) * x(blockRow*t + col);
                }
            }
        }
        return res;
    }

    template <typename U>
    friend std::ostream& operator<<(std::ostream& outp, NFold<U>& x) {
        outp << dvar(x.n, x.r, x.s, x.t) << std::endl
             << "A: " << std::endl;
        F0R(rr, x.r + x.n*x.s) {
            F0R(cc, x.n*x.t)
                outp << x(rr, cc) << ' ';
            outp << std::endl;
        }
        return outp;
    }

    template <typename U>
    friend std::istream& operator>>(std::istream& inp, NFold<U>& x) {
        x.l = Vec<T>(x.n*x.t);
        x.u = Vec<T>(x.n*x.t);
        x.b = Vec<T>(x.r + x.n*x.s);
        x.c = Vec<T>(x.n*x.t);
        x.as = std::vector<Mat<T>>(x.n, Mat<T>(x.r, x.t)); 
        x.bs = std::vector<Mat<T>>(x.n, Mat<T>(x.s, x.t));
        return inp >> x.l >> x.u >> x.b >> x.c >> x.as >> x.bs;
    }
};
