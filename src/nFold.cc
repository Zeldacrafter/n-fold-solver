template<typename T>
class NFold {
  public:
    using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

  private:
    void readMatr(Mat& m, std::istream& inp) {
        F0R (r, m.rows()) {
            F0R (c, m.cols()) {
                inp >> m(r, c);
            }
        }
    }

  public:
    int n, r, s, t;
    Vec l, u, b, c;
    Mat A;

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

    NFold(std::istream& inp) {
        inp >> n >> r >> s >> t;

        l = Vec(n*t);
        F0R (i, SZ(l)) inp >> l(i);

        u = Vec(n*t);
        F0R (i, SZ(u)) inp >> u(i);

        b = Vec(r + n*s);
        F0R (i, SZ(b)) inp >> b(i);

        c = Vec(n*t);
        F0R (i, SZ(c)) inp >> c(i);
        
        std::vector<Mat> as(n, Mat(r, t)); 
        for (auto& m : as) readMatr(m, inp);

        std::vector<Mat> bs(n, Mat(s, t));
        for (auto& m : bs) readMatr(m, inp);

        A = Mat(r + n*s, t*n);
        F0R (i, n) {
            F0R (rr, r)
                F0R (cc, t)
                    A(rr, i * t + cc) = as[i](rr, cc);
            F0R (rr, s)
                F0R (cc, t)
                    A(r + i * s + rr, i * t + cc) = bs[i](rr, cc);
        }
    }
};
