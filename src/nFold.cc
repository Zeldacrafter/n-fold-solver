template<typename T>
class NFold {
  public:
    using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  public:
    int n, r, s, t;
    Vec l, u, b, c;
    Mat A;

    NFold(int _n, int _r, int _s, int _t) : n{_n}, r{_r}, s{_s}, t{_t} {}

    NFold(int _n, int _r, int _s, int _t, 
            Vec _l, Vec _u,
            std::vector<Mat> as, std::vector<Mat> bs)
            : n{_n}, r{_r}, s{_s}, t{_t}, l{_l}, u{_u}, A(r + n*s, t*n) {

        F0R (i, n) {
            F0R (rr, r)
                F0R (cc, t)
                    A(rr, i*t + cc) = as[i](rr, cc);
            F0R (rr, s)
                F0R (cc, t)
                    A(r + i*s + rr, i * t + cc) = bs[i](rr, cc);
        }
    }

    template <typename U>
    friend std::istream& operator>>(std::istream& inp, NFold<U>& x) {
        assert(x.n && x.r && x.s && x.t);

        x.l = Vec(x.n*x.t);
        F0R (i, SZ(x.l)) inp >> x.l(i);

        x.u = Vec(x.n*x.t);
        F0R (i, SZ(x.u)) inp >> x.u(i);

        x.b = Vec(x.r + x.n*x.s);
        F0R (i, SZ(x.b)) inp >> x.b(i);

        x.c = Vec(x.n*x.t);
        F0R (i, SZ(x.c)) inp >> x.c(i);
        
        std::vector<Mat> as(x.n, Mat(x.r, x.t)); 
        for (auto& m : as) 
            F0R (row, m.rows())
                F0R (col, m.cols())
                    inp >> m(row, col);

        std::vector<Mat> bs(x.n, Mat(x.s, x.t));
        for (auto& m : bs)
            F0R (row, m.rows())
                F0R (col, m.cols())
                    inp >> m(row, col);

        x.A = Mat(x.r + x.n*x.s, x.t*x.n);
        F0R (i, x.n) {
            F0R (rr, x.r)
                F0R (cc, x.t)
                    x.A(rr, i * x.t + cc) = as[i](rr, cc);
            F0R (rr, x.s)
                F0R (cc, x.t)
                    x.A(x.r + i * x.s + rr, i * x.t + cc) = bs[i](rr, cc);
        }

        return inp;
    }
};
