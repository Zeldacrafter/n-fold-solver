#include <gtest/gtest.h>
#include <fstream>

#include "bruteforce.cpp"
#include "NFold.cpp"
#include "Solver.cpp"

class SolverFixture : public ::testing::TestWithParam<std::string> {
protected:
    NFold<int> nfold;
    long long wantedOutput;

    SolverFixture() : nfold(0, 0, 0, 0) {
        std::string path = GetParam();

        std::ifstream inputFile(path + ".in");
        int n, r, s, t;
        inputFile >> n >> r >> s >> t;
        nfold = NFold<int>(n, r, s, t);
        inputFile >> nfold;

        std::ifstream outputFile(path + ".out");
        outputFile >> wantedOutput;
    }
};

TEST_P(SolverFixture, BruteForceTest) {
    ASSERT_EQ(wantedOutput, bruteForceBest(nfold).second);
}

TEST_P(SolverFixture, AlgorithmTest) {
    auto solution = solver::Solver(nfold).solve();
    ASSERT_TRUE(solution);
    ASSERT_EQ(wantedOutput, solution->second);
}

INSTANTIATE_TEST_CASE_P(
        BruteForceTests,
        SolverFixture,
        ::testing::Values(
                "input/max-small",
                "input/max-large"
        )
);