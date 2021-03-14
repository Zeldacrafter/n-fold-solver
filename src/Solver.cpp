#pragma once

#include <map>
#include <boost/container_hash/hash.hpp>

#include "NFold.cpp"

typedef int T;

template <typename K>
struct std::hash<Vec<K>> {
    size_t operator()(const Vec<K>& v) const {
        return boost::hash_range(v.data(), v.data() + v.size());
    }
};

namespace solver {
    using std::pair, std::make_pair;
    using std::vector;
    using std::optional, std::nullopt;

    typedef vector<std::pair<long long, Vec<T>>> stateVec;
    typedef tsl::hopscotch_map<Vec<T>, std::pair<long long, Vec<T>>> graphMap;

    class Solver {

    public:
        NFold<T> x;
        T l_a; // TODO: Do we want to handle this?

        explicit Solver(NFold<T>& _x) : x{_x} { }

        optional<pair<Vec<T>, T>> solve() {
            optional<Vec<T>> initSolution = findInitSol();
            if (initSolution) {
                assert(x * *initSolution == x.b);
                return optional(solve(*initSolution));
            } else {
                return nullopt;
            }
        }

        pair<Vec<T>, T> solve(const Vec<T> &initSolution, bool findZero = false) {
            assert(x * initSolution == x.b);

            Vec<T> z0 = initSolution;
            // improve on z0 via fixpoint algorithm.
            for (bool changed = true; changed;) {
                changed = false;

                Vec<T> l = x.l - z0;
                Vec<T> u = x.u - z0;

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step;
                Vec<T> augRes = solveAugIp(l, u);
                assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                // The resulting vector has to be an integer solution to the nfold.
                Vec<T> nextCandidate = z0 + augRes;
                assert(x * nextCandidate == x.b);

                T currWeight = nextCandidate.dot(x.c);
                if (currWeight > z0.dot(x.c)) {
                    z0 = nextCandidate;
                    changed = true;

                    // If we are only looking for a zero weight solution we can just return here.
                    // This is the case for finding an initial solution.
                    if (findZero && !currWeight) {
                        return make_pair(nextCandidate, 0);
                    }
                }
/*
                // Number of bits of lambda we have to guess is bounded
                // by ceil(log(gamma)) + 1
                T gamma = (x.u - x.l).maxCoeff();
                T maxIters = ceil(log2(gamma)) + 1; // TODO: Cleaner expression please

                F0R(i, maxIters + 1) {
                    T lambda = 1 << i;

                    Vec<T> l = (((x.l - z0).array() + lambda - 1) / lambda);
                    Vec<T> u = (((x.u - z0).array()) / lambda);

                    // After ever augmentation step we should have Ay = 0 for the result of the augmentation step;
                    Vec<T> augRes = solveAugIp(l, u);
                    assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                    // The resulting vector has to be an integer solution to the nfold.
                    Vec<T> nextCandidate = z0 + lambda * augRes;
                    assert(x * nextCandidate == x.b);

                    T currWeight = nextCandidate.dot(x.c);
                    if (currWeight > z0.dot(x.c)) {
                        z0 = nextCandidate;
                        changed = true;

                        // If we are only looking for a zero weight solution we can just return here.
                        // This is the case for finding an initial solution.
                        if (findZero && !currWeight) {
                            return make_pair(nextCandidate, 0);
                        }
                    }
                }
*/
            }

            return make_pair(z0, z0.dot(x.c));
        }

    private:
        Vec<T> solveAugIp(const Vec<T> &l, const Vec<T> &u) {
            Vec<T> zero = Vec<T>::Zero(x.r + x.s);

            // Track the weight and path to the current position.
            graphMap curr;
            curr[zero] = make_pair(0, Vec<T>(0));

            F0R(block, x.n) {
                // Move to the next block and save the result in curr.
                processBlock(block, curr, l, u);
            }

            return curr[zero].second;
        }

        void processBlock(size_t block, graphMap &curr, const Vec<T> &l, const Vec<T> &u) {
            Vec<T> zero = Vec<T>::Zero(x.r + x.s);

            Mat<T> M(x.r + x.s, x.t);
            M << x.as[block], x.bs[block];

            F0R(col, x.t) {
                size_t yPos = block * x.t + col;

                graphMap next;
                for (const auto& [oldPos, wgtAndVec] : curr) {
                    auto [wgt, oldVec] = wgtAndVec;
                    FOR(y, l(yPos), u(yPos) + 1) {

                        Vec<T> candidate = y * M.col(col) + oldPos;
                        T candidateWeight = wgt + x.c(yPos) * y;
                        // In the last step we only want to add elements if they are an possible solution
                        // We also only add an edge if it is the first or the one with the highest weight on its
                        // path to that note.
                        if((col != x.t - 1 || candidate.tail(x.s).isZero()) &&
                               (!next.count(candidate) || next[candidate].first < candidateWeight)) {
                            Vec<T> newVec(SZ(oldVec) + 1);
                            newVec << oldVec, y;
                            next[candidate] = make_pair(wgt + x.c(yPos) * y, newVec);
                        }
                    }
                }

                curr = move(next);
            }

            assert(SZ(curr));
        }

        optional<Vec<T>> findInitSol() {
            auto[aInit, initSol] = constructAInit(x);
            auto[sol, weight] = Solver(aInit).solve(initSol, true);

            if(!weight) {
                // Solution found.
                Vec<T> res(x.n * x.t);
                F0R(i, x.n) {
                    res.segment(i * x.t, x.t) = sol.segment(i * (x.t + x.s + x.r), x.t);
                }
                res = res + x.l; // The constructed solution is offset by l. We need to adjust for that here.

                assert(x * res == x.b);
                return optional{res};
            } else {
                // No solution exists.
                return nullopt;
            }
        }
    };
}
