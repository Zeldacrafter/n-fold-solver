#ifndef NFOLD_TEMPLATE_H
#define NFOLD_TEMPLATE_H

#ifdef USING_BOOST
#include <boost/stacktrace.hpp>
#endif
#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <utility>
#include <memory>

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T, int S1, int S2>
using sMat = Eigen::Matrix<T, S1, S2>;
template<typename T, int S1>
using sVec = Eigen::Matrix<T, S1, 1>;

#define DEBUG 1

// https://google.github.io/styleguide/cppguide.html#Template_metaprogramming
// everything from the template metaprogramming section was
// successfully ignored
#define ALL(x) (x).begin(), (x).end()
#define RALL(x) (x).rbegin(), (x).rend()
#define SZ(x) (int)(x).size()
#define FOR(a, b, c) for (auto a = (b); (a) < (c); ++(a))
#define F0R(a, b) FOR (a, static_cast<decltype(b)>(0), (b))
template <typename T>
bool ckmin(T& a, const T& b) { return a > b && (a = b, true); }
template <typename T>
bool ckmax(T& a, const T& b) { return a < b && (a = b, true); }
// Output to 'cerr' if 'DEBUG' flag is set. Do nothing otherwise.
#ifndef DEBUG
#define DEBUG 0
#endif
#define dout if (DEBUG) std::cerr
// Output all passed variables with their corresponding name and value.
#define dvarimpl(...) mkDB(#__VA_ARGS__, __VA_ARGS__) << " "
#define dvar(...) " " << dvarimpl(__VA_ARGS__)
#define dvarr(...) " \033[31m" << dvarimpl(__VA_ARGS__)
#define dvarb(...) " \033[34m" << dvarimpl(__VA_ARGS__)
#define dvarg(...) " \033[32m" << dvarimpl(__VA_ARGS__)
#define dvary(...) " \033[33m" << dvarimpl(__VA_ARGS__)
#define dvarc(...) " \033[36m" << dvarimpl(__VA_ARGS__)
#define dvari(...) " \033[7m" << dvarimpl(__VA_ARGS__)

///////////////////////////////////////////////////////////////
// Utility functions.
///////////////////////////////////////////////////////////////

namespace impl {
    template <typename T, typename F, size_t... Is>
    F for_each(T& t, F f, std::index_sequence<Is...>) {
        auto l = { (f(std::get<Is>(t), Is), 0)... };
        (void) l;
        return f;
    }
}

template <typename... Ts, typename F>
F for_each(std::tuple<Ts...>& t, F f) {
    return impl::for_each(t, f, std::make_index_sequence<sizeof...(Ts)>{});
}

template <typename... Ts, typename F>
F for_each(const std::tuple<Ts...>& t, F f) {
    return impl::for_each(t, f, std::make_index_sequence<sizeof...(Ts)>{});
}


// IsC indicates whether a type defines a 'const_iterator'.
// IsC::value is true if 'const_iterator' exists and false otherwise.
template <typename T> std::true_type const_iterator_check(typename T::const_iterator*);
template <typename T> std::false_type const_iterator_check(...);
template <typename T> struct IsC : decltype(const_iterator_check<T>(nullptr)) {};
// No new input/output for string as those already exist.
template <> struct IsC<std::string> : std::false_type {};
template <typename U, int S> struct IsC<sVec<U, S>> : std::false_type {};
#ifdef USING_BOOST
template <> struct IsC<boost::stacktrace::stacktrace> : std::false_type {};
#endif

///////////////////////////////////////////////////////////////
// Begin Output
///////////////////////////////////////////////////////////////

// Forward declarations.
template <typename T>
std::enable_if_t<IsC<T>::value, std::ostream&> operator<<(std::ostream&, const T&);
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream&, const std::pair<T1, T2>&);

// Print each tuple element.
template <typename... Ts>
std::ostream& operator<<(std::ostream& o, const std::tuple<Ts...>& t) {
    o << '(';
    for_each(t, [&](auto& x, size_t i) { if(i) o << ", "; o << x; });
    return o << ')';
}

// Special case for 1-tuple to avoid confusing parentheses
template <typename T>
std::ostream& operator<<(std::ostream& o, const std::tuple<T>& t) {
    return o << std::get<0>(t);
}

// Output for pairs via above defined tuple output routine.
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& o, const std::pair<T1, T2>& p) {
    return o << '(' << p.first << ", " << p.second << ')';
}

// Output every element in a container with 'begin' and 'end' iterators.
template <typename T>
std::enable_if_t<IsC<T>::value, std::ostream&> operator<<(std::ostream& o, const T& c) {
    o << '[';
    for (auto it = c.cbegin(); it != c.cend(); ++it)
        o << *it << (next(it) != c.cend() ? ", " : "");
    return o << ']';
}

///////////////////////////////////////////////////////////////
// Debug output
///////////////////////////////////////////////////////////////

template <typename... Ts>
struct DB {
    std::string n;
    std::tuple<Ts...> d;
    DB(const std::string& ns, Ts... ds) : n{ns}, d{ds...} {}
    friend std::ostream& operator<<(std::ostream& o, const DB& db) {
        int i = 0;
        for_each(db.d, [&](const auto& e, int idx) {
            (idx ? o << " " : o) << "[";
            while (i < SZ(db.n) and std::isspace(db.n[i])) ++i;
            int br = 0, str = 0, chr = 0, esc = 0;
            while (i < SZ(db.n)) {
                if (db.n[i] == '\\') esc = not esc;
                if (not chr and not esc and db.n[i] == '\"') str = not str;
                if (not str and not esc and db.n[i] == '\'') chr = not chr;
                if (not str and not chr) {
                    br += brt(db.n[i]);
                    if (db.n[i] == ',' and br == 0) {
                        ++i;
                        break;
                    }
                }
                if (db.n[i] != '\\') esc = false;
                o << db.n[i++];
            }
            o << ": " << e << "]";
        });
        return o;
    }
    static inline int brt(char c) {
        switch (c) {
            case '(': case '[': case '{': return 1;
            case ')': case ']': case '}': return -1;
            default: return 0;
        }
    }
};
template <typename... Ts>
DB<Ts...> mkDB(const std::string& n, Ts... d) { return DB<Ts...>(n, d...); }

///////////////////////////////////////////////////////////////
// Pretty output
///////////////////////////////////////////////////////////////

// PrettyPrint struct that contains a value to be printed and
// a list of seperators which indicate how different dimensions
// of multidimensional values should be seperated.
template <typename T, size_t N>
struct PP {
    // Value to print.
    const T& v;
    // Pointer to seperator list.
    std::shared_ptr<std::array<std::string, N>> se;
    // Index of next seperator.
    size_t idx;
    PP(const T& value, std::shared_ptr<std::array<std::string, N>> p, size_t i = 0)
            : v{value}, se{p}, idx{i} { }
};

// If a value is not a pair, tuple or std-library-continer just print it.
// Pairs and tuples are implemented via template specialization further down.
template <typename T, size_t M>
std::enable_if_t<not IsC<T>::value, std::ostream&>
operator<<(std::ostream& o, const PP<T, M>& p) {
    return o << p.v;
}

// Prints every tuple element.
template <size_t M, typename... Ts>
std::ostream& operator<<(std::ostream& o, const PP<std::tuple<Ts...>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '(';
    for_each(p.v, [&](auto& x, size_t i) {
        if(i) o << sep;
        o << PP<std::decay_t<decltype(x)>, M>(x, p.second, p.idx + 1);
    });
    return o << ')';
}

// Print pairs with the specified seperator for that level.
template <typename T1, typename T2, size_t M>
std::ostream& operator<<(std::ostream& o, const PP<std::pair<T1, T2>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '(';
    return o << PP<T1, M>(p.v.first, p.se, p.idx + 1) << sep
             << PP<T2, M>(p.v.second, p.se, p.idx + 1) << ')';
}

template <typename K, size_t M>
std::ostream& operator<<(std::ostream& o, const PP<const Vec<K>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '<';
    F0R (i, static_cast<size_t>(p.v.size())) {
        if (i) o << sep;
        o << p.v(i);
    }
    return o << '>';
}

template <typename K, size_t M>
std::ostream& operator<<(std::ostream& o, const PP<Vec<K>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '<';
    F0R (i, static_cast<size_t>(p.v.size())) {
        if (i) o << sep;
        o << p.v(i);
    }
    return o << '>';
}

template <typename K, size_t M, int S>
std::ostream& operator<<(std::ostream& o, const PP<const sVec<K, S>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '<';
    F0R (i, static_cast<size_t>(p.v.size())) {
        if (i) o << sep;
        o << p.v(i);
    }
    return o << '>';
}

template <typename K, size_t M, int S>
std::ostream& operator<<(std::ostream& o, const PP<sVec<K, S>, M>& p) {
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    o << '<';
    F0R (i, S) {
        if (i) o << sep;
        o << p.v(i);
    }
    return o << '>';
}

// Print std-library-container with the specified seperator.
template <typename T, size_t M>
std::enable_if_t<IsC<T>::value, std::ostream&>
operator<<(std::ostream& o, const PP<T, M>& p) {
    // Seperator for the current layer (or default)
    const std::string& sep = p.idx < M ? (*p.se)[p.idx] : " ";
    // Print every container element
    o << '[';
    for (auto it = p.v.cbegin(); it != p.v.cend(); ++it)
        o << PP<typename T::value_type, M>(*it, p.se, p.idx + 1)
          << (next(it) != p.v.cend() ? sep : "");
    return o << ']';
}

// This function is the main way for a user to interface with the PrettyPrinter.
template <typename T, typename... Ts, size_t N = sizeof...(Ts)>
PP<T, N> pp(const T& value, Ts... seps) {
    return PP<T, N>(value, std::make_shared<std::array<std::string, N>>(std::array<std::string, N>{seps...}));
}

///////////////////////////////////////////////////////////////
// Begin Input
///////////////////////////////////////////////////////////////

// Forward declarations.
template <typename T>
std::enable_if_t<IsC<T>::value, std::istream&> operator>>(std::istream&, T&);
template <typename T1, typename T2>
std::istream& operator>>(std::istream&, std::pair<T1, T2>&);

// Read a tuple.
template <typename... Ts>
std::istream& operator>>(std::istream& i, std::tuple<Ts...>& t) {
    for_each(t, [&](auto& x, int) { i >> x; });
    return i;
}

// Read the contents of a 'pair' object.
template <typename T1, typename T2>
std::istream& operator>>(std::istream& i, std::pair<T1, T2>& p) {
    return i >> p.fi >> p.se;
}

// Read containers with 'begin' and 'end' iterators.
template <typename T>
std::enable_if_t<IsC<T>::value, std::istream&> operator>>(std::istream& i, T& v) {
    for (auto& x : v) i >> x;
    return i;
}

#endif //NFOLD_TEMPLATE_H
