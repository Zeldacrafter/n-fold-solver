#pragma once

#include <map>
#include "NFold.cpp"
#include "bruteforce.cpp"


typedef int T;

struct cmpVecs {
    bool operator()(const Vec<T>& a, const Vec<T>& b) const {
        assert(SZ(a) == SZ(b));
        F0R(i, SZ(a))
            if(a(i) != b(i))
                return a(i) < b(i);
        return false;
    }
};

class Solver {
public:
    NFold<T> x;
    T l_a;

    Solver(NFold<T>& _x) : x{_x} {
        T delta = x.getDelta();

        T part1 = 1;
        F0R(i, x.s) part1 *= 2*x.r*delta + 1;

        T part2 = 1;
        F0R(i, x.r) part2 *= 2*x.r*x.getDelta()*part1 + 1;

        l_a = part1 * part2;
    }

    T solve() {
        using std::endl;
        dout << endl << dvar(l_a) << endl;

        Vec<T> initSolution = findInitSol();
        assert(x * initSolution == x.b);
        dout << "Found init solution " << pp(initSolution)
             << " with weight " << initSolution.dot(x.c) << endl;

        Vec<T> z0 = initSolution;
        bool changed = true; // TODO: Actually calc iterations on N.
        while(changed) {
            changed = false;

            // Number of bits of lambda we have to guess is bounded
            // by ceil(log(gamma)) + 1
            T gamma = (x.u - x.l).maxCoeff();
            T maxIters = ceil(log2(gamma)) + 1; // TODO: Cleaner expression please

            Vec<T> currBest = z0;
            assert(x*currBest == x.b);

            F0R(i, maxIters + 1) {
                T lambda = 1 << i;

                Vec<T> l = (((x.l - z0).array() + lambda - 1)/lambda)
                        .max(-l_a)
                        .matrix();
                Vec<T> u = (((x.u - z0).array())/lambda)
                        .min(l_a)
                        .matrix();
                dout << dvar(i, lambda, pp(l), pp(u), pp(z0)) << endl;

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step;
                Vec<T> augRes = solveAugIp(l, u);
                dout << "Augmentation result: " << dvar(pp(augRes), pp(x * augRes)) << endl;
                assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                // The resulting vector has to be an integer solution to the nfold.
                Vec<T> nextCandidate = z0 + lambda*augRes;
                assert(x * nextCandidate == x.b);

                dout << "Found valid candidate with cost " << dvar(nextCandidate.dot(x.c)) << endl;
                if(nextCandidate.dot(x.c) > currBest.dot(x.c)) {
                    dout << "New best found with weight " << nextCandidate.dot(x.c) << endl;
                    currBest = nextCandidate;
                    changed = true;
                }
            }
            z0 = currBest; // TODO: Why not replace it earlier?
        }

        return z0.dot(x.c);
    }

    Vec<T> solveAugIp(const Vec<T>& l, const Vec<T>& u) {
        // Track the weight and path to the current position.
        std::pair<long long, Vec<T>> state = std::make_pair(0, Vec<T>(0));

        F0R(block, x.n) {
            // Move to the next block
            state = processBlock(block, state, l, u);
        }

        return state.second;
    }

    std::pair<long long, Vec<T>> processBlock(size_t block, const std::pair<long long, Vec<T>>& state,
                                              const Vec<T>& l, const Vec<T>& u) {
        T delta = x.getDelta();
        std::map<Vec<T>, std::pair<long long, Vec<T>>, cmpVecs> curr;
        curr[Vec<T>::Zero(x.r + x.s)] = state;

        Mat<T> M(x.r + x.s, x.t);
        M << x.as[block], x.bs[block];

        F0R(col, x.t) {
            // We are moving to U^(block)_col now.

            std::map<Vec<T>, std::pair<long long, Vec<T>>, cmpVecs> next;

            for(auto& [oldPos, wgtAndVec] : curr) {
                auto& [wgt, oldVec] = wgtAndVec;

                size_t yPos = block*x.t + col;
                FOR(y, l(yPos), u(yPos) + 1) {
                    Vec<T> candidate = y*M.col(col) + oldPos;

                    if(candidate.lpNorm<Eigen::Infinity>() <= delta*l_a) {
                        long long newWgt = wgt + x.c(yPos)*y;
                        // cannot just call ckmax because we need to replace 0 entries.
                        if(!next.count(candidate) || newWgt > next[candidate].first) {

                            Vec<T> newBest(SZ(oldVec) + 1);
                            newBest << oldVec, y;
                            next[candidate] = std::make_pair(wgt + x.c(yPos)*y, newBest);
                        }
                    }
                }
            }
            curr = std::move(next);
        }

        assert(curr.count(Vec<T>::Zero(x.r + x.s)));
        return curr[Vec<T>::Zero(x.r + x.s)];
    }

    Vec<T> findInitSol() {
        // TODO: Actual implementation here
        return bruteForceWorst(x).first;
    }

};
