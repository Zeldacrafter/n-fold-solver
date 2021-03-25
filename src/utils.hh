#ifndef N_FOLD_SOLVER_UTILS_HH
#define N_FOLD_SOLVER_UTILS_HH

#include <iostream>
#include <boost/container_hash/hash.hpp>
#include <Eigen/Dense>

///////////////////////////////////////////////
/// Typedefs for Eigen matrices and vectors ///
///////////////////////////////////////////////

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T, int S1, int S2>
using sMat = Eigen::Matrix<T, S1, S2>;
template<typename T, int S1>
using sVec = Eigen::Matrix<T, S1, 1>;

////////////////////////////////////////////////
/// Input operators for vectors and matrices ///
////////////////////////////////////////////////

template<typename U>
std::istream& operator>>(std::istream& inp, Vec<U>& x) {
    for(int i = 0; i < x.size(); ++i) {
        inp >> x(i);
    }
    return inp;
}

template<typename U, int S>
std::istream& operator>>(std::istream& inp, sVec<U, S>& x) {
    for(int i = 0; i < x.size(); ++i) {
        inp >> x(i);
    }
    return inp;
}

template<typename T>
std::istream& operator>>(std::istream& inp, Mat<T>& x) {
    for(int r = 0; r < x.rows(); ++r) {
        for(int c = 0; c < x.cols(); ++c) {
            inp >> x(r, c);
        }
    }
    return inp;
}

template<typename T, int S1, int S2>
std::istream& operator>>(std::istream& inp, sMat<T, S1, S2>& x) {
    for(int r = 0; r < x.rows(); ++r) {
        for(int c = 0; c < x.cols(); ++c) {
            inp >> x(r, c);
        }
    }
    return inp;
}

////////////////////////////////////
/// Output operators for vectors ///
////////////////////////////////////

template <typename K>
std::ostream& operator<<(std::ostream& o, Vec<K>& v) {
    o << '<';
    for(int i = 0; i < v.size(); ++i) {
        if (i) o << ", ";
        o << v(i);
    }
    return o << '>';
}

template <typename K, int S>
std::ostream& operator<<(std::ostream& o, sVec<K, S> v) {
    o << '<';
    for(int i = 0; i < v.size(); ++i) {
        if (i) o << ", ";
        o << v(i);
    }
    return o << '>';
}

namespace utils {

    /**
     * Hash function for static size Eigen vectors.
     * @tparam K The type in the vector.
     * @tparam S The size of the vector.
     */
    template <typename K>
    struct staticVectorHash {
        size_t operator()(const Vec<K>& v) const {
            return boost::hash_range(v.data(), v.data() + v.size());
        }
    };

    /**
     * Implementation of the mathematical signum function.
     * @tparam T The type of the value to get the sign off.
     * @param val The value to get the sign of.
     * @return {@code 1} if {@code val > 0},
     *         {@code -1} if {@code val < 0} and
     *         {@code 0} if {@code val == 0}.
     */
    template <typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }
}

#endif //N_FOLD_SOLVER_UTILS_HH
