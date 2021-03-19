#pragma once

#include <map>
#include <boost/container_hash/hash.hpp>
#include <tsl/hopscotch_map.h>

#include "static_nfold_class.cpp"

// TODO: https://google.github.io/styleguide/cppguide.html#std_hash
template <typename K, int S>
struct std::hash<sVec<K, S>> {
    size_t operator()(const sVec<K, S>& v) const {
        return boost::hash_range(v.data(), v.data() + v.size());
    }
};

namespace static_solver {
    /* TODO: What do we do about L_A? It is supposed to be a upper bound on our solution.
     *       But since L_A is very large any input that would benefit from it would Not
     *       terminate before the heat death of the universe.
     */
    using ::std::pair, ::std::make_pair;
    using ::std::optional, ::std::nullopt;

    template <typename U, int N, int R, int S, int T>
    class StaticSolver {
    public:
        explicit StaticSolver(StaticNFold<U, N, R, S, T>& _x) : x{_x} { }

        optional<pair<sVec<U, N*T>, U>> solve() {
            optional<sVec<U, N*T>> initSolution = findInitSol(x);
            if (initSolution) {
                assert(x * *initSolution == x.b);
                return optional(solve(*initSolution));
            } else {
                return nullopt;
            }
        }

        pair<sVec<U, N*T>, U> solve(const sVec<U, N*T> &initSolution, const optional<U> knownBest = nullopt) {
            assert(x * initSolution == x.b);

            sVec<U, N*T> z0 = initSolution;
            // As long as a better solution eSts it is guaranteed that we find one.
            // This means we can improve our initial solution with a fix-point algorithm.
            for (bool changed = true; changed;) {
                changed = false;
                dout << "Starting augmentation step" << std::endl;

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step.
                // This means that we still have A(z0 + y) = b.
                sVec<U, N*T> augRes = solveAugIp(x.l - z0, x.u - z0, startLayer);
                assert(x * augRes == (sVec<U, R + N*S>::Zero()));

                // The resulting vector has to be an integer solution to the Nfold.
                sVec<U, N*T> nextCandidate = z0 + augRes;
                assert(x * nextCandidate == x.b);

                if (U currWeight = nextCandidate.dot(x.c); currWeight > z0.dot(x.c)) {
                    dout << dvar(currWeight) << std::endl;
                    z0 = nextCandidate;
                    changed = true;

                    // Heuristic:
                    // If we know that our optimal solution is if one eSts and are simply
                    // concerned with finding one/determining if one eSts we can return preemptively
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
               T maTers = ceil(log2(gamma)) + 1; // TODO: Cleaner eRession please

               for(int i = maTers; i >= 0; --i) {
                   T lambda = 1 << i;

                   Vec<U> l = (((x.l - z0).array() + lambda - 1) / lambda);
                   Vec<U> u = (((x.u - z0).array()) / lambda);

                   // After ever augmentation step we should have Ay = 0 for the result of the augmentation step;
                   Vec<U> augRes = solveAugIp(l, u);
                   assert(x * augRes == Vec<U>::Zero(SZ(x.b)));

                   // The resulting vector has to be an integer solution to the Nfold.
                   Vec<U> NextCandidate = z0 + lambda * augRes;
                   assert(x * NextCandidate == x.b);

                   T currWeight = NextCandidate.dot(x.c);
                   if (currWeight > z0.dot(x.c)) {
                       z0 = NextCandidate;
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
        template<int X>
        using graphLayer = tsl::hopscotch_map<
                                 sVec<U, R + S>,
                                 pair<U, sVec<U, X>>,
                                 std::hash<sVec<U, R + S>>,
                                 std::equal_to<sVec<U, R + S>>,
                                 Eigen::aligned_allocator<pair<sVec<U, R + S>, pair<U, sVec<U, X>>>>>;

        StaticNFold<U, N, R, S, T> x;
        const graphLayer<0> startLayer = { {sVec<U, R + S>::Zero(), make_pair(0, sVec<U, 0>())} };

        ///////////////////////////////////////////////////////////////
        /// SolveAugIp                                              ///
        /// This routine had to be implemented recursively          ///
        /// to make it possible to check all types at compile time. ///
        /// This makes it possible to use the statically sized      ///
        /// vectors and matrices of Eigen which is a significant    ///
        /// performance boost overall and well worth                ///
        /// the added complexity in the implementation.             ///
        ///////////////////////////////////////////////////////////////

        template<int BLOCK = 0>
        std::enable_if_t<BLOCK < N, sVec<U, N*T>>
        solveAugIp(const sVec<U, N*T> &l, const sVec<U, N*T> &u, graphLayer<BLOCK*T> curr) {
            // We are not yet at the last block.
            // This means that there is at least one more block to process now:

            sMat<U, R + S, T> M;
            M << x.as[BLOCK], x.bs[BLOCK];
            // Process the next block.
            graphLayer<BLOCK*T + T> next = processSubBlocks<BLOCK>(std::move(curr), l, u, M);
            // Move to processing the following block (if one exists).
            return solveAugIp<BLOCK + 1>(l, u, std::move(next));
        }

        template<int BLOCK>
        std::enable_if_t<BLOCK == N, sVec<U, N*T>>
        solveAugIp(const sVec<U, N*T>&, const sVec<U, N*T>&, graphLayer<BLOCK*T> curr) {
            // If we are at block N we are done.
            // There are no more blocks to process and we can just return.
            assert(curr.count(sVec<U, R + S>::Zero()));
            return curr[sVec<U, R + S>::Zero()].second;
        }

        ////////////////////////
        /// ProcessSubBlocks ///
        ////////////////////////

        template<int BLOCK, int COL = 0>
        std::enable_if_t<COL < T, graphLayer<BLOCK*T + T>>
        processSubBlocks(graphLayer<BLOCK*T + COL> curr,
                         const sVec<U, N*T> &l, const sVec<U, N*T> &u,
                         const sMat<U, R + S, T> M) {

            size_t yPos = BLOCK * T + COL;
            graphLayer<BLOCK*T + COL + 1> next;

            for (const auto& [oldPos, wgtAndVec] : curr) {
                auto [wgt, oldVec] = wgtAndVec;
                FOR(y, l(yPos), u(yPos) + 1) {

                    sVec<U, R + S> candidate = y * M.col(COL) + oldPos;
                    U candidateWeight = wgt + x.c(yPos) * y;
                    // We only add an edge if there is No other edge to that Node (yet)
                    // or the New path to that Node has a higher total weight. (longest path)

                    // In the last step we only want to add elements if they are an possible solution.
                    // This means that they Need to solve the corresponding B block fully (last s elements are 0)
                    // since the B block is fixed after this.
                    if((!next.count(candidate) || next[candidate].first < candidateWeight)
                       && (COL != T - 1 || candidate.tail(S).isZero())) {
                        sVec<U, BLOCK*T + COL + 1> NewVec;
                        NewVec << oldVec, y;
                        next[candidate] = make_pair(wgt + x.c(yPos) * y, NewVec);
                    }
                }
            }

            assert(SZ(next));
            return processSubBlocks<BLOCK, COL + 1>(std::move(next), l, u, M);
        }

        template<int BLOCK, int COL>
        std::enable_if_t<COL == T, graphLayer<BLOCK*T + T>>
        processSubBlocks(graphLayer<BLOCK*T + COL> curr,
                         const sVec<U, N*T>&, const sVec<U, N*T>&,
                         const sMat<U, R + S, T>&) {
            return curr;
        }

        /////////////////////////////
        /// Find Initial Solution ///
        ////////////////////////////

        // Finds an initial solution to an NFold instance or reports
        // that none exists.
        static optional<sVec<U, N*T>> findInitSol(StaticNFold<U, N, R, S, T>& x) {
            auto [aInit, initSol] = constructAInit(x);
            auto [sol, weight] = StaticSolver<U, N, R, S, T + R + S>(aInit).solve(initSol, 0);

            if(!weight) {
                // Solution found.
                sVec<U, N*T> res(N * T);
                F0R(i, N) {
                    res.segment(i * T, T) = sol.segment(i * (T + S + R), T);
                }
                res = res + x.l; // The constructed solution is offset by l. We Need to adjust for that here.

                assert(x * res == x.b);
                return optional{res};
            } else {
                // No solution eSts.
                return nullopt;
            }
        }

        // Constructs an NFold as described by Jansens paper in chapter 4.
        // This is used to find an initial solution for the original input Nfold.
        static std::pair<StaticNFold<U, N, R, S, T + R + S>, sVec<U, N * (T + R + S)>> constructAInit(const StaticNFold<U, N, R, S, T>& x) {
            StaticNFold<U, N, R, S, T + R + S> res;

            //Construct New matrix
            F0R(i, N) {
                res.as[i].block(0, 0,     R, T) = x.as[i];
                if(!i) res.as[i].block(0, T, R, R).setIdentity();
                else   res.as[i].block(0, T, R, R).setZero();
                res.as[i].block(0, T + R, R, S).setZero();

                res.bs[i].block(0, 0,     S, T) = x.bs[i];
                res.bs[i].block(0, T,     S, R).setZero();
                res.bs[i].block(0, T + R, S, S).setIdentity();
            }

            // Construct New righthand side
            res.b = x.b - x*x.l;

            // Construct upper and lower bound
            F0R(i, N) {
                res.l.segment(i*(T + R + S), T) = (x.u - x.l).segment(i*T, T);
                res.u.segment(i*(T + R + S), T) = (x.u - x.l).segment(i*T, T);
                if(!i) {
                    res.l.segment(i*(T + R + S) + T, R) = res.b.segment(0, R);
                    res.u.segment(i*(T + R + S) + T, R) = res.b.segment(0, R);
                } else {
                    res.l.segment(i*(T + R + S) + T, R).setZero();
                    res.u.segment(i*(T + R + S) + T, R).setZero();
                }
                res.l.segment(i*(T + R + S) + T + R, S) = res.b.segment(R + i*S, S);
                res.u.segment(i*(T + R + S) + T + R, S) = res.b.segment(R + i*S, S);
            }
            res.l = res.l.array().min(0);
            res.u = res.u.array().max(0);

            // Construct cost vector
            res.c.setZero();
            res.c.segment(T, R).setOnes();
            F0R(j, R) // Identity matrix for A_1
                if(res.b(j) >= 0)
                    res.c(T + j) = -1;
            F0R(i, N) { // Identity matrix for B_i
                res.c.segment(T + R + i*(T + R + S), S).setOnes();
                F0R(j, S)
                    if(res.b(R + i*S + j) >= 0)
                        res.c(T + R + i*(T + R + S) + j) = -1;
            }

            // Construct init sol
            sVec<U, N*(T + R + S)> initSol = sVec<U, N*(T + R + S)>::Zero();
            initSol.segment(T, R) = res.b.segment(0, R);
            F0R(i, N) {
                initSol.segment((T + R + S)*i + T + R, S) = res.b.segment(R + i*S, S);
            }

            return std::make_pair(res, initSol);
        }
    };
} // Namespace solver
