/// \file  func.h
/// \brief Declaration of random number generators (LCG, Xorshift32, MiddleSquare) 
///        and statistical test functions for randomness analysis.

#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdint>


/**
 * @brief Linear Congruential Generator (LCG) random number generator
 */
class LCG {
private:
    uint32_t state; ///< Current state of the generator
public:
    /**
     * @brief Generates the next random number in sequence
     * @return Next 32-bit unsigned random number
     */
    LCG(uint32_t seed);

    /**
     * @brief Generates the next random number in sequence
     * @return Next 32-bit unsigned random number
     */
    uint32_t next();
};

/**
 * @brief Xorshift32 random number generator
 */
class Xorshift32 {
private:
    uint32_t state; ///< Current state of the generator

public:
    /**
     * @brief Constructor that initializes the generator with a seed
     * @param seed Initial seed value
     */
    Xorshift32(uint32_t seed);

    /**
     * @brief Generates the next random number using xorshift algorithm
     * @return Next 32-bit unsigned random number
     */
    uint32_t next();
};

/**
 * @brief Middle Square random number generator
 */
class MiddleSquare {
private:
    uint32_t state; ///< Current state of the generator
public:

    /**
     * @brief Constructor that initializes the generator with a seed
     * @param seed Initial seed value
     */
    MiddleSquare(uint32_t seed);

    /**
     * @brief Generates the next random number using middle-square method
     * @return Next 32-bit unsigned random number
     */
    uint32_t next();
};

/**
 * @brief Calculates incomplete gamma function using series expansion
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Value of the incomplete gamma function
 */
double incompleteGammaSeries(double a, double x);

/**
 * @brief Calculates incomplete gamma function using continued fraction
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Value of the incomplete gamma function
 */
double incompleteGammaContinuedFrac(double a, double x);

/**
 * @brief Calculates the regularized upper incomplete gamma function Q(a,x)
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Value of Q(a,x)
 */
double gammaQ(double a, double x);

/**
 * @brief Calculates the standard normal cumulative distribution function
 * @param x Input value
 * @return Probability that a random variable is less than or equal to x
 */
double Phi(double x);

/**
 * @brief Performs the monobit test on a binary sequence
 * @param bits Vector of bits (0s and 1s)
 * @return P-value from the test
 */
double monobitTestP(const std::vector<int>& bits);

/**
 * @brief Performs the block frequency test on a binary sequence
 * @param bits Vector of bits (0s and 1s)
 * @param M Block size
 * @return P-value from the test
 */
double blockFrequencyTestP(const std::vector<int>& bits, int M);

/**
 * @brief Performs the runs test on a binary sequence
 * @param bits Vector of bits (0s and 1s)
 * @return P-value from the test
 */
double runsTestP(const std::vector<int>& bits);

/**
 * @brief Performs the longest run of ones test on a binary sequence
 * @param bits Vector of bits (0s and 1s)
 * @return P-value from the test
 */
double longestRunOnesTestP(const std::vector<int>& bits);

/**
 * @brief Performs the cumulative sums test on a binary sequence
 * @param bits Vector of bits (0s and 1s)
 * @return P-value from the test
 */
double cumulativeSumsTestP(const std::vector<int>& bits);

/**
 * @brief Calculates basic statistics for a set of values
 * @param values Input values
 * @param[out] mean Calculated mean
 * @param[out] stddev Calculated standard deviation
 * @param[out] cov Calculated coefficient of variation
 */
void calculateStats(const std::vector<uint32_t>& values, double& mean, double& stddev, double& cov);

/**
 * @brief Performs chi-square goodness-of-fit test
 * @param values Input values to test
 * @return P-value from the test
 */
double chiSquareTestP(const std::vector<uint32_t>& values);

/**
 * @brief Converts 32-bit values to a sequence of bits
 * @param values Input values
 * @return Vector of bits (0s and 1s)
 */
std::vector<int> generateBits(const std::vector<uint32_t>& values);

/**
 * @brief Performs a battery of NIST statistical tests on a bit sequence
 * @param bits Input bit sequence to test
 * @param[out] p1 P-value from monobit test
 * @param[out] p2 P-value from block frequency test
 * @param[out] p3 P-value from runs test
 * @param[out] p4 P-value from longest run of ones test
 * @param[out] p5 P-value from cumulative sums test
 */
void performNISTTests(const std::vector<int>& bits, double& p1, double& p2, double& p3, double& p4, double& p5);

#endif // FUNC_H