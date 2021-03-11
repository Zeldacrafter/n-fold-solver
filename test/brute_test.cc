#include <gtest/gtest.h>
#include <fstream>

#include "../src/bruteforce.cc"

class BruteForceFixture : public ::testing::TestWithParam<string> {
protected:
    string path;

    BruteForceFixture() {
        path = GetParam();
    }
};

TEST_P(BruteForceFixture, BruteForceTest) {
    string path = GetParam();

    ifstream outputFile(path + ".out");
    int solution;
    outputFile >> solution;
    cout << "For path " << path << ".out: " << solution << endl;

    ifstream inputFile(path + ".in");
    ASSERT_EQ(solution, bruteForce(inputFile));
}

INSTANTIATE_TEST_CASE_P(
    BruteForceTests,
    BruteForceFixture,
    ::testing::Values(
        "input/max-small",
        "input/max-large"
    )
);
