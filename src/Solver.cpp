#pragma once

#include <map>
#include <boost/container_hash/hash.hpp>

#include "NFold.cpp"
#include "bruteforce.cpp"


typedef int T;

struct cmpVecs {
    bool operator()(const Vec<T>& a, const Vec<T>& b) const {
        assert(SZ(a) == SZ(b));
        return std::lexicographical_compare(a.data(), a.data() + SZ(a), b.data(), b.data() + SZ(b));
    }
};

template <typename K>
struct std::hash<Vec<K>> {
    size_t operator()(const Vec<K>& v) const {
        return boost::hash_range(v.data(), v.data() + v.size());
    }
};

class Solver {
    typedef std::vector<std::pair<long long, Vec<T>>> stateVec;
    typedef tsl::hopscotch_map<Vec<T>, stateVec> graphMap;
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

    std::pair<Vec<T>, T> solve() {
        Vec<T> initSolution = findInitSol();
        assert(x * initSolution == x.b);
        return solve(initSolution);
    }

    std::pair<Vec<T>, T> solve(Vec<T>& initSolution, bool findZero = false) {
        assert(x*initSolution == x.b);

        Vec<T> z0 = initSolution;
        for(bool changed = true; changed;) {
            changed = false;

            // Number of bits of lambda we have to guess is bounded
            // by ceil(log(gamma)) + 1
            T gamma = (x.u - x.l).maxCoeff();
            T maxIters = ceil(log2(gamma)) + 1; // TODO: Cleaner expression please

            Vec<T> currBest = z0;
            assert(x*currBest == x.b);

            F0R(i, maxIters + 1) {
                T lambda = 1 << i;

                Vec<T> l = (((x.l - z0).array() + lambda - 1)/lambda);   //.max(-l_a).matrix();
                Vec<T> u = (((x.u - z0).array())/lambda);                //.min(l_a).matrix();

                // After ever augmentation step we should have Ay = 0 for the result of the augmentation step;
                Vec<T> augRes = solveAugIp(l, u);
                assert(x * augRes == Vec<T>::Zero(SZ(x.b)));

                // The resulting vector has to be an integer solution to the nfold.
                Vec<T> nextCandidate = z0 + lambda*augRes;
                assert(x * nextCandidate == x.b);

                if(nextCandidate.dot(x.c) > currBest.dot(x.c)) {
                    currBest = nextCandidate;
                    changed = true;

                    // If we are only looking for a zero weight solution we can just return here.
                    if(findZero && nextCandidate.dot(x.c) == 0) {
                        return std::make_pair(nextCandidate, 0);
                    }
                }
            }
            z0 = currBest; // TODO: Why not replace it earlier?
        }

        return std::make_pair(z0, z0.dot(x.c));
    }

    Vec<T> solveAugIp(const Vec<T>& l, const Vec<T>& u) {
        Vec<T> zero = Vec<T>::Zero(x.r + x.s);

        // Track the weight and path to the current position.
        graphMap curr;
        curr[zero].emplace_back(0, Vec<T>(0));

        F0R(block, x.n) {
            // Move to the next block
            processBlock(block, curr, l, u);
        }

        // TODO: with stl
        Vec<T> best;
        long long bestWgt = std::numeric_limits<T>::min();

        for(auto [wgt, v] : curr[zero])
            if(ckmax(bestWgt, wgt))
                best = v;
        return best;
    }

    void processBlock(size_t block, graphMap& curr, const Vec<T>& l, const Vec<T>& u) {
        Vec<T> zero = Vec<T>::Zero(x.r + x.s);

        Mat<T> M(x.r + x.s, x.t);
        M << x.as[block], x.bs[block];

        F0R(col, x.t) {
            size_t yPos = block * x.t + col;

            graphMap next;
            for(auto& [oldPos, wgtAndVecs] : curr) {
                for(auto& [wgt, oldVec] : wgtAndVecs) {

                    FOR(y, l(yPos), u(yPos) + 1) {

                        Vec<T> candidate = y * M.col(col) + oldPos;
                        Vec<T> newVec(SZ(oldVec) + 1);
                        newVec << oldVec, y;
                        next[candidate].emplace_back(wgt + x.c(yPos) * y, newVec);

                        // TODO: if (candidate.lpNorm<Eigen::Infinity>() <= delta * l_a) {
                    }
                }
            }

            curr = std::move(next);
        }

        for (auto it = curr.begin(); it != curr.end();) {
            if (!it->first.tail(x.s).isZero()) it = curr.erase(it);
            else ++it;
        }

        assert(SZ(curr));
    }

    Vec<T> findInitSol() {
        auto [aInit, initSol] = constructAInit(x);
        Solver solver(aInit);
        auto [sol, wgt] = solver.solve(initSol, true);

        // TODO: Handle non zero wgt -> No solution

        Vec<T> res(x.n*x.t);
        F0R(i, x.n)
            res.segment(i*x.t, x.t) = sol.segment(i*(x.t + x.s + x.r), x.t);
        res = res + x.l; // The original solution is offset by l

        assert(x*res == x.b);
        return res;
    }

};
