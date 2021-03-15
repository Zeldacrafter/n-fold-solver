#pragma once

#include <map>
#include <boost/container_hash/hash.hpp>

#include "NFold.cpp"

template <typename K>
struct std::hash<Vec<K>> {
    size_t operator()(const Vec<K>& v) const {
        return boost::hash_range(v.data(), v.data() + v.size());
    }
};

namespace solver {
    /* TODO: What do we do about L_A? It is supposed to be a upper bound on our solution.
     *       But since L_A is very large any input that would benefit from it would not
     *       terminate before the heat death of the universe.
     */
    using std::pair, std::make_pair;
    using std::vector;
    using std::optional, std::nullopt;


    template <typename T>
    class Solver {
        typedef tsl::hopscotch_map<Vec<T>, std::pair<long long, Vec<T>>> graphLayer;

    public:
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

            Vec<T>& z0 = initSolution;
            // As long as a better solution exists it is guaranteed that we find one.
            // This means we can improve our initial solution with a fix-point algorithm.
            for (bool changed = true; changed;) {
                changed = false;

                Vec<T> l = x.l - z0;
                Vec<T> u = x.u - z0;

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step.
                // This means that we still have A(z0 + y) = b.
                Vec<T> augRes = solveAugIp(l, u);
                assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                // The resulting vector has to be an integer solution to the nfold.
                Vec<T> nextCandidate = z0 + augRes;
                assert(x * nextCandidate == x.b);

                T currWeight = nextCandidate.dot(x.c);
                if (currWeight > z0.dot(x.c)) {
                    z0 = nextCandidate;
                    changed = true;

                    // Heuristic:
                    // If we know that our optimal solution is at best 0 we can just return here.
                    // This is the case for finding an initial solution
                    // and cuts down on execution time by quite a bit.
                    if (findZero && !currWeight) {
                        return make_pair(nextCandidate, 0);
                    }
                }
                 /*
                // TODO: Why would this thing with lambda be an improvement?
                //       A lambda > 1 only makes the search space a strict subspace
                //       compared to the one with lambda = 1. Why would we find more solutions like this?
                //       Maybe read up on that in the paper again?
                // Number of bits of lambda we have to guess is bounded
                // by ceil(log(gamma)) + 1
                T gamma = (x.u - x.l).maxCoeff();
                T maxIters = ceil(log2(gamma)) + 1; // TODO: Cleaner expression please

                for(int i = maxIters; i >= 0; --i) {
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
        NFold<T> x;

        Vec<T> solveAugIp(const Vec<T> &l, const Vec<T> &u) {
            Vec<T> zero = Vec<T>::Zero(x.r + x.s);

            // Track the weight and path to all nodes in the current layer of the graph.
            graphLayer curr;
            curr[zero] = make_pair(0, Vec<T>(0));

            F0R(block, x.n) {
                // Move to the next block and save the result in curr.
                processBlock(block, curr, l, u);
                assert(SZ(curr));
            }

            return curr[zero].second;
        }

        void processBlock(size_t block, graphLayer &curr, const Vec<T> &l, const Vec<T> &u) {
            Vec<T> zero = Vec<T>::Zero(x.r + x.s);

            Mat<T> M(x.r + x.s, x.t);
            M << x.as[block], x.bs[block];

            F0R(col, x.t) {
                size_t yPos = block * x.t + col;

                graphLayer next;
                for (const auto& [oldPos, wgtAndVec] : curr) {
                    auto [wgt, oldVec] = wgtAndVec;
                    FOR(y, l(yPos), u(yPos) + 1) {

                        Vec<T> candidate = y * M.col(col) + oldPos;
                        T candidateWeight = wgt + x.c(yPos) * y;
                        // We only add an edge if there is no other edge to that node (yet)
                        // or the new path to that node has a higher total weight. (longest path)

                        // In the last step we only want to add elements if they are an possible solution.
                        // This means that they need to solve the corresponding B block fully (last s elements are 0)
                        // since the B block is fixed after this.
                        if((!next.count(candidate) || next[candidate].first < candidateWeight)
                                && (col != x.t - 1 || candidate.tail(x.s).isZero())) {
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
