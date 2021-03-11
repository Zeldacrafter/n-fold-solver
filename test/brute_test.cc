#include <gtest/gtest.h>
#include <fstream>

#include "../src/bruteforce.cc"

class BruteForceFixture : public ::testing::TestWithParam<std::string> {
protected:
    int bruteForceOutput;
    int wantedOutput;

    BruteForceFixture() {
        std::string path = GetParam();

        std::ifstream inputFile(path + ".in");
        int n, r, s, t;
        inputFile >> n >> r >> s >> t;
        NFold<int> nFold(n, r, s, t);
        inputFile >> nFold;
        bruteForceOutput = bruteForce(nFold);

        std::ifstream outputFile(path + ".out");
        outputFile >> wantedOutput;
    }
};

TEST_P(BruteForceFixture, BruteForceTest) {
    ASSERT_EQ(wantedOutput, bruteForceOutput);
}

INSTANTIATE_TEST_CASE_P(
    BruteForceTests,
    BruteForceFixture,
    ::testing::Values(
        "input/max-small",
        "input/max-large"
    )
);
