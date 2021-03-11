#include <gtest/gtest.h>
#include <fstream>

#include "../src/bruteforce.cc"

void testFile(string path) {
    ifstream outputFile(path + ".out");
    int solution;
    outputFile >> solution;

    ifstream inputFile(path + ".in");

    ASSERT_EQ(solution, bruteForce(inputFile));
}

TEST(GoogleTestCi, BruteForceTest) {
    testFile("input/max-small");
    testFile("input/max-large");
}
