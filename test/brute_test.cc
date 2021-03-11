#include <gtest/gtest.h>
#include <fstream>

#include "../src/bruteforce.cc"

class BruteForceFixture : public ::testing::TestWithParam<string> {
protected:
    int bruteForceOutput;
    int wantedOutput;

    BruteForceFixture() {
        string path = GetParam();

        ifstream inputFile(path + ".in");
        bruteForceOutput = bruteForce(inputFile);

        ifstream outputFile(path + ".out");
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
