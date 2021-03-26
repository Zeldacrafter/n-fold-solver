#ifndef N_FOLD_SOLVER_CLASS_HH
#define N_FOLD_SOLVER_CLASS_HH

#include <map>
#include <stack>
#pragma GCC diagnostic push // Ignore shadow warnings in external library.
#pragma GCC diagnostic ignored "-Wshadow"
#include "../third-party/hopscotch-map/include/tsl/hopscotch_map.h"
#pragma GCC diagnostic pop

#include "utils.hh"
#include "nfold_class.hh"
#include "prefix_tree_class.hh"

template <typename U, int N, int R, int S, int T>
class n_fold_solver {
public:
    explicit n_fold_solver(n_fold<U, N, R, S, T>& _x) : x{_x} { }

    std::optional<std::pair<sVec<U, N*T>, U>> solve() {
        std::optional<sVec<U, N*T>> initSolution = findInitSol(x);
        if (initSolution) {
            assert(x * *initSolution == x.b);
            return std::optional(solve(*initSolution));
        } else {
            return std::nullopt;
        }
    }

    std::pair<sVec<U, N*T>, U> solve(const sVec<U, N*T> &initSolution,
                                     const std::optional<U> knownBest = std::nullopt) {
        assert(x * initSolution == x.b);

        sVec<U, N*T> z0 = initSolution;
        // As long as a better solution eSts it is guaranteed that we find one.
        // This means we can improve our initial solution with a fix-point algorithm.
        for (bool changed = true; changed;) {
            changed = false;

            // After ever augmentation step we should have Ay = 0 for the result of the augmentation step.
            // This means that we still have A(z0 + y) = b.
            //sVec<U, N*T> augRes = solveAugIp(x.l - z0, x.u - z0, startLayer);
            std::optional<sVec<U, N*T>> augRes = solveAugIp(x.l - z0, x.u - z0);
            if(!augRes) {
                // No solution was found. We cannot improve our result further;
                break;
            }
            assert(x * *augRes == (sVec<U, R + N*S>::Zero()));

            // The resulting vector has to be an integer solution to the Nfold.
            sVec<U, N*T> nextCandidate = z0 + *augRes;
            assert(x * nextCandidate == x.b);

            if (U currWeight = nextCandidate.dot(x.c); currWeight > z0.dot(x.c)) {
                z0 = nextCandidate;
                changed = true;

                // Heuristic:
                // If we know what our optimal solution is and we are simply
                // concerned with finding one/determining if one exists we can return preemptively.
                // This is the case for finding an initial solution
                // and cuts down on execution time by quite a bit.
                if (knownBest && *knownBest == currWeight) {
                    break;
                }
            }
        }

        return std::make_pair(z0, z0.dot(x.c));
    }

private:

    n_fold<U, N, R, S, T> x;
    prefix_tree<U> nodes;

    using graphLayer = tsl::hopscotch_map<
            sVec<U, R + S>,
            std::pair<U, size_t>,
            utils::staticVectorHash<U, R + S>,
            std::equal_to<sVec<U, R + S>>,
            Eigen::aligned_allocator<std::pair<sVec<U, R + S>, std::pair<U, size_t>>>>;

    std::optional<sVec<U, N*T>> solveAugIp(const sVec<U, N*T> &l, const sVec<U, N*T> &u) {
        nodes.clear();
        sVec<U, R + S> zero = sVec<U, R + S>::Zero();

        graphLayer curr;
        int startIndex = nodes.add(U(0), prefix_tree<U>::NO_PARENT);
        assert(startIndex == 0);
        curr[zero] = std::make_pair(U(0), startIndex);

        for(int block = 0; block < N; ++block) {

            sMat<U, R + S, T> M;
            M << x.as[block], x.bs[block];

            for(int col = 0; col < T; ++col) {
                size_t yPos = block*T + col;

                graphLayer next;
                for (const auto& [oldPos, wgtAndIdx] : curr) {
                    auto& [wgt, idx] = wgtAndIdx;
                    int amtFound = 0;
                    for(int y = l(yPos); y <= u(yPos); ++y) {
                        sVec<U, R + S> candidate = y * M.col(col) + oldPos;
                        U candidateWeight = wgt + x.c(yPos) * y;

                        // In the last step we only want to add elements if they are an possible solution.
                        // This means that they need to solve the corresponding B block fully (last s elements are 0)
                        // since the B block is fixed after this.
                        if(col != T - 1 || candidate.tail(S).isZero()) {
                            // We only add an edge if there is no other edge to that node (yet)
                            // or the new path to that node has a higher total weight. (longest path)

                            if(auto ptr = next.find(candidate);
                                    ptr == next.end() || ptr->second.first < candidateWeight) {
                                if (ptr != next.end()) {
                                    // We already have an element at that position.
                                    nodes.remove(ptr->second.second, false);
                                }
                                int insertionIndex = nodes.add(y, idx);
                                next[candidate] = std::make_pair(candidateWeight, insertionIndex);
                                amtFound++;
                            }
                        }
                    }

                    // If this is a dead end we can just remove the node and free up space.
                    if(!amtFound) {
                        nodes.remove(idx);
                    }
                }

                curr = std::move(next);
            }

        }

        if(curr.size()) {
            std::vector<U> path = nodes.constructPath(curr[zero].second);
            assert(path.size() == N*T);
            return sVec<U, N*T>(path.data());
        } else {
            return std::nullopt;
        }
    }

    /////////////////////////////
    /// Find Initial Solution ///
    /////////////////////////////

    // Finds an initial solution to an NFold instance or reports
    // that none exists.
    static std::optional<sVec<U, N*T>> findInitSol(n_fold<U, N, R, S, T>& x) {
        auto [aInit, initSol] = constructAInit(x);
        auto [sol, weight] = n_fold_solver<U, N, R, S, T + R + S>(aInit).solve(initSol, 0);

        if(!weight) {
            // Solution found.
            sVec<U, N*T> res(N * T);
            for(int i = 0; i < N; ++i) {
                res.segment(i * T, T) = sol.segment(i * (T + S + R), T);
            }
            res = res + x.l; // The constructed solution is offset by l. We Need to adjust for that here.

            assert(x * res == x.b);
            return res;
        } else {
            // No solution exists.
            return std::nullopt;
        }
    }

    // Constructs an NFold as described by Jansens paper in chapter 4.
    // This is used to find an initial solution for the original input Nfold.
    static std::pair<n_fold<U, N, R, S, T + R + S>, sVec<U, N * (T + R + S)>>
    constructAInit(const n_fold<U, N, R, S, T>& x) {
        n_fold<U, N, R, S, T + R + S> res;

        //Construct New matrix
        for(int i = 0; i < N; ++i) {
            res.as[i].block(0, 0, R, T) = x.as[i];
            if(!i) {
                res.as[i].block(0, T, R, R).setIdentity();
            } else {
                res.as[i].block(0, T, R, R).setZero();
            }
            res.as[i].block(0, T + R, R, S).setZero();

            res.bs[i].block(0, 0,     S, T) = x.bs[i];
            res.bs[i].block(0, T,     S, R).setZero();
            res.bs[i].block(0, T + R, S, S).setIdentity();
        }

        // Construct New righthand side
        res.b = x.b - x*x.l;

        // Construct upper and lower bound
        for(int i = 0; i < N; ++i) {
            res.l.segment(i*(T + R + S), T) = (x.u - x.l).segment(i*T, T);
            if(!i) res.l.segment(i*(T + R + S) + T, R) = res.b.segment(0, R);
            else   res.l.segment(i*(T + R + S) + T, R).setZero();
            res.l.segment(i*(T + R + S) + T + R, S) = res.b.segment(R + i*S, S);
        }
        res.u = res.l;
        res.l = res.l.array().min(0);
        res.u = res.u.array().max(0);

        // Construct cost vector
        res.c.setZero();
        res.c.segment(T, R).setOnes();
        for(int j = 0; j < R; ++j) {
            res.c(T + j) = -utils::sgn(res.b(j));
        }
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < S; ++j) {
                res.c(T + R + i*(T + R + S) + j) = -utils::sgn(res.b(R + i*S + j));
            }
        }

        // Construct init sol
        sVec<U, N*(T + R + S)> initSol = sVec<U, N*(T + R + S)>::Zero();
        initSol.segment(T, R) = res.b.segment(0, R);
        for(int i = 0; i < N; ++i) {
            initSol.segment((T + R + S)*i + T + R, S) = res.b.segment(R + i*S, S);
        }

        return std::make_pair(res, initSol);
    }
};

#endif //N_FOLD_SOLVER_CLASS_HH