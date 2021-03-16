#pragma once

#include <map>
#include <boost/container_hash/hash.hpp>
#include <tsl/hopscotch_map.h>

#include "nfold_class.cpp"

// TODO: https://google.github.io/styleguide/cppguide.html#std_hash
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
    using ::std::pair, ::std::make_pair;
    using ::std::optional, ::std::nullopt;

    template <typename T>
    class Solver {
      public:
        explicit Solver(NFold<T>& _x) : x{_x} { }

        optional<pair<Vec<T>, T>> solve() {
            optional<Vec<T>> initSolution = findInitSol(x);
            if (initSolution) {
                assert(x * *initSolution == x.b);
                return optional(solve(*initSolution));
            } else {
                return nullopt;
            }
        }

        pair<Vec<T>, T> solve(const Vec<T> &initSolution, const optional<T> knownBest = nullopt) {
            assert(x * initSolution == x.b);

            Vec<T> z0 = initSolution;
            // As long as a better solution exists it is guaranteed that we find one.
            // This means we can improve our initial solution with a fix-point algorithm.
            for (bool changed = true; changed;) {
                changed = false;
                dout << "Starting augmentation step" << std::endl;

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step.
                // This means that we still have A(z0 + y) = b.
                Vec<T> augRes = solveAugIp(x.l - z0, x.u - z0);
                assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                // The resulting vector has to be an integer solution to the nfold.
                Vec<T> nextCandidate = z0 + augRes;
                assert(x * nextCandidate == x.b);

                if (T currWeight = nextCandidate.dot(x.c); currWeight > z0.dot(x.c)) {
                    dout << dvar(currWeight) << std::endl;
                    z0 = nextCandidate;
                    changed = true;

                    // Heuristic:
                    // If we know that our optimal solution is if one exists and are simply
                    // concerned with finding one/determining if one exists we can return preemptively
                    // This is the case for finding an initial solution
                    // and cuts down on execution time by quite a bit.
                    if (knownBest && *knownBest == currWeight) {
                        return make_pair(nextCandidate, currWeight);
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
        typedef tsl::hopscotch_map<Vec<T>, pair<T, Vec<T>>> graphLayer;

        NFold<T> x;

        Vec<T> solveAugIp(const Vec<T> &l, const Vec<T> &u) {
            Vec<T> zero = Vec<T>::Zero(x.r + x.s);

            // Track the weight and path to all nodes in the current layer of the graph.
            graphLayer curr;
            curr[zero] = make_pair(0, Vec<T>(0));

            F0R(block, x.n) {
                // Move to the next block and save the result in curr.
                curr = processBlock(block, std::move(curr), l, u);
                assert(SZ(curr));
            }

            return curr[zero].second;
        }

        graphLayer processBlock(size_t block, graphLayer curr, const Vec<T> &l, const Vec<T> &u) {
            dout << dvar(block, SZ(curr)) << std::endl;
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
            return curr;
        }

        // Finds an initial solution to an NFold instance or reports
        // that none exists.
        static optional<Vec<T>> findInitSol(NFold<T>& x) {
            auto[aInit, initSol] = constructAInit(x);
            auto[sol, weight] = Solver<T>(aInit).solve(initSol, 0);

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

        // Constructs an nFold as described by Jansens paper in chapter 4.
        // This is used to find an initial solution for the original input nfold.
        static std::pair<NFold<T>, Vec<T>> constructAInit(const NFold<T>& x) {
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
    };
} // namespace solver
