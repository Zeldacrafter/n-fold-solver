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

template <typename U>
class n_fold_solver {
public:
    size_t N, R, S, T;
    explicit n_fold_solver(n_fold<U>& _x)
        : x{_x}, N(_x.N), R(_x.R), S(_x.S), T(_x.T)  {
    }

    std::optional<std::pair<Vec<U>, U>> solve() {
        std::optional<Vec<U>> initSolution = findInitSol(x);
        if (initSolution) {
            assert(x * *initSolution == x.b);
            return std::optional(solve(*initSolution));
        } else {
            return std::nullopt;
        }
    }

    std::pair<Vec<U>, U> solve(const Vec<U> &initSolution,
                                     const std::optional<U> knownBest = std::nullopt) {
        assert(x * initSolution == x.b);

        Vec<U> z0 = initSolution;
        // As long as a better solution eSts it is guaranteed that we find one.
        // This means we can improve our initial solution with a fix-point algorithm.
        for (bool changed = true; changed;) {
            changed = false;

            // After ever augmentation step we should have Ay = 0 for the result of the augmentation step.
            // This means that we still have A(z0 + y) = b.
            //sVec<U, N*T> augRes = solveAugIp(x.l - z0, x.u - z0, startLayer);
            std::optional<Vec<U>> augRes = solveAugIp(x.l - z0, x.u - z0);
            if(!augRes) {
                // No solution was found. We cannot improve our result further;
                break;
            }
            assert(x * *augRes == (Vec<U>::Zero(R + N*S)));

            // The resulting vector has to be an integer solution to the Nfold.
            Vec<U> nextCandidate = z0 + *augRes;
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

    n_fold<U> x;
    prefix_tree<U> nodes;

    using graphLayer = tsl::hopscotch_map<
            Vec<U>,
            std::pair<U, size_t>,
            utils::staticVectorHash<U>,
            std::equal_to<Vec<U>>,
            Eigen::aligned_allocator<std::pair<Vec<U>, std::pair<U, size_t>>>>;

    std::optional<Vec<U>> solveAugIp(const Vec<U> &l, const Vec<U> &u) {
        nodes.clear();
        Vec<U> zero = Vec<U>::Zero(R + S);

        graphLayer curr;
        int startIndex = nodes.add(U(0), prefix_tree<U>::NO_PARENT);
        assert(startIndex == 0);
        curr[zero] = std::make_pair(U(0), startIndex);

        for(int block = 0; block < N; ++block) {

            Mat<U> M(R + S, T);
            M << x.as[block], x.bs[block];

            for(int col = 0; col < T; ++col) {
                size_t yPos = block*T + col;

                graphLayer next;
                for (const auto& [oldPos, wgtAndIdx] : curr) {
                    auto& [wgt, idx] = wgtAndIdx;
                    int amtFound = 0;
                    for(int y = l(yPos); y <= u(yPos); ++y) {
                        Vec<U> candidate = y * M.col(col) + oldPos;
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
                                    nodes.remove(ptr->second.second);
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
            Vec<U> res(path.size());
            for(int i = 0; i < path.size(); ++i)
                res(i) = path[i];
            return res;
        } else {
            return std::nullopt;
        }
    }

    /////////////////////////////
    /// Find Initial Solution ///
    /////////////////////////////

    // Finds an initial solution to an NFold instance or reports
    // that none exists.
    static std::optional<Vec<U>> findInitSol(n_fold<U>& x) {
        size_t N = x.N, R = x.R, S = x.S, T = x.T;

        auto [aInit, initSol] = constructAInit(x);
        auto [sol, weight] = n_fold_solver<U>(aInit).solve(initSol, 0);

        if(!weight) {
            // Solution found.
            Vec<U> res(N * T);
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
    static std::pair<n_fold<U>, Vec<U>>
    constructAInit(const n_fold<U>& x) {
        size_t N = x.N, R = x.R, S = x.S, T = x.T;
        n_fold<U> res(N, R, S, T + R + S);

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
        Vec<U> initSol = Vec<U>::Zero(N*(T + R + S));
        initSol.segment(T, R) = res.b.segment(0, R);
        for(int i = 0; i < N; ++i) {
            initSol.segment((T + R + S)*i + T + R, S) = res.b.segment(R + i*S, S);
        }

        return std::make_pair(res, initSol);
    }
};

#endif //N_FOLD_SOLVER_CLASS_HH