/// \file func.cpp
/// \brief Implementation of random number generators and statistical tests for randomness analysis

#include "func.h"

/**
 * @brief Constructs a Linear Congruential Generator with given seed
 * @param seed Initial seed value
 */
LCG::LCG(uint32_t seed) {
    state = seed;
}

/**
 * @brief Generates next random number using LCG algorithm
 * @return 32-bit unsigned random number
 */
uint32_t LCG::next() {
    uint64_t res = uint64_t(1664525) * state + 1013904223;
    state = static_cast<uint32_t>(res);
    return state;
}

/**
 * @brief Constructs a Xorshift32 generator with given seed
 * @param seed Initial seed value
 */
Xorshift32::Xorshift32(uint32_t seed) {
    state = seed;
}

/**
 * @brief Generates next random number using Xorshift algorithm
 * @return 32-bit unsigned random number
 */
uint32_t Xorshift32::next() {
    uint32_t x = state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    state = x;
    if (state == 0) state = 0xBAADF00D;
    return state;
}

/**
 * @brief Constructs a Middle Square generator with given seed
 * @param seed Initial seed value
 */
MiddleSquare::MiddleSquare(uint32_t seed) {
    state = seed;
}

/**
 * @brief Generates next random number using Middle Square method
 * @return 32-bit unsigned random number
 */
uint32_t MiddleSquare::next() {
    uint64_t squared = static_cast<uint64_t>(state) * state;
    state = static_cast<uint32_t>((squared >> 16) & 0xFFFFFFFF);
    if (state == 0) state = 0xBAADF00D;
    return state;
}

/**
 * @brief Calculates incomplete gamma function using series expansion
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Computed value of incomplete gamma function
 */
double incompleteGammaSeries(double a, double x) {
    const int ITMAX = 1000;
    const double EPS = 1e-9;
    double sum = 1.0 / a;
    double term = sum;
    for (int n = 1; n <= ITMAX; ++n) {
        term *= x / (a + n);
        sum += term;
        if (fabs(term) < fabs(sum) * EPS) break;
    }
    return sum * exp(-x + a * log(x));
}

/**
 * @brief Calculates incomplete gamma function using continued fraction
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Computed value of incomplete gamma function
 */
double incompleteGammaContinuedFrac(double a, double x) {
    const int ITMAX = 1000;
    const double EPS = 1e-9;
    const double FPMIN = 1e-30;
    double b = x + 1.0 - a;
    double c = 1.0 / FPMIN;
    double d = 1.0 / b;
    double f = d;
    for (int i = 1; i <= ITMAX; ++i) {
        double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        double delta = d * c;
        f *= delta;
        if (fabs(delta - 1.0) < EPS) break;
    }
    return exp(-x + a * log(x)) * f;
}

/**
 * @brief Computes regularized upper incomplete gamma function Q(a,x)
 * @param a Shape parameter
 * @param x Upper limit of integration
 * @return Value of Q(a,x)
 */
double gammaQ(double a, double x) {
    if (x < 0.0 || a <= 0.0) {
        return 0.0;
    }
    if (x == 0.0) {
        return 1.0;
    }
    if (x < a + 1.0) {
        double gam = incompleteGammaSeries(a, x);
        double log_gamma_a = lgamma(a);
        double P = gam / exp(log_gamma_a);
        return 1.0 - P;
    } else {
        double Q = incompleteGammaContinuedFrac(a, x);
        double log_gamma_a = lgamma(a);
        Q /= exp(log_gamma_a);
        return Q;
    }
}

/**
 * @brief Computes standard normal cumulative distribution function
 * @param x Input value
 * @return Probability that a random variable is less than or equal to x
 */
double Phi(double x) {
    return 0.5 * (1.0 + std::erf(x / M_SQRT2));
}

/**
 * @brief Performs the monobit test on a binary sequence (NIST SP 800-22)
 * @param bits Vector of bits (0s and 1s) to test
 * @return P-value indicating the probability of randomness (values >= 0.01 suggest randomness)
 * @details Tests the proportion of ones and zeros in the sequence to verify they are approximately equal
 */
double monobitTestP(const std::vector<int>& bits) {
    int n = bits.size();
    long sum = 0;
    for (int b : bits) sum += b;
    long S = 2 * sum - n;
    double s_obs = fabs(S) / sqrt(n);
    double p = std::erfc(s_obs / M_SQRT2);
    return p;
}

/**
 * @brief Performs the block frequency test on a binary sequence (NIST SP 800-22)
 * @param bits Vector of bits (0s and 1s) to test
 * @param M Block size (default = 128)
 * @return P-value indicating the probability of randomness
 * @details Tests the proportion of ones within M-bit blocks to verify they are approximately 50%
 */
double blockFrequencyTestP(const std::vector<int>& bits, int M = 128) {
    int n = bits.size();
    if (n < M) {
        return 0.0;
    }
    int N = n / M;
    if (N <= 0) {
        return 0.0;
    }
    double chi_sq = 0.0;
    for (int i = 0; i < N; ++i) {
        int count1 = 0;
        for (int j = 0; j < M; ++j) {
            if (bits[i * M + j] == 1) count1++;
        }
        double pi = (double)count1 / M;
        chi_sq += (pi - 0.5) * (pi - 0.5);
    }
    chi_sq *= 4.0 * M;
    double p = gammaQ(N / 2.0, chi_sq / 2.0);
    return p;
}

/**
 * @brief Performs the runs test on a binary sequence (NIST SP 800-22)
 * @param bits Vector of bits (0s and 1s) to test
 * @return P-value indicating the probability of randomness
 * @details Tests the total number of runs (both 0-runs and 1-runs) in the sequence
 */
double runsTestP(const std::vector<int>& bits) {
    int n = bits.size();
    int count1 = 0;
    for (int b : bits) {
        if (b == 1) count1++;
    }
    double pi = (double)count1 / n;
    if (fabs(pi - 0.5) >= (2.0 / sqrt(n))) {
        return 0.0;
    }
    int runs = 1;
    for (int i = 1; i < n; ++i) {
        if (bits[i] != bits[i-1]) {
            runs++;
        }
    }
    double runs_expected = 2.0 * n * pi * (1 - pi);
    double sigma = 2 * sqrt(n) * pi * (1 - pi);
    if (sigma <= 0) {
        return 0.0;
    }
    double Z = fabs(runs - runs_expected) / sigma;
    double p = std::erfc(Z / M_SQRT2);
    return p;
}

/**
 * @brief Performs the longest run of ones test on a binary sequence (NIST SP 800-22)
 * @param bits Vector of bits (0s and 1s) to test
 * @return P-value indicating the probability of randomness
 * @details Tests the longest run of ones within M-bit blocks of the sequence
 */
double longestRunOnesTestP(const std::vector<int>& bits) {
    int n = bits.size();
    int M;
    if (n < 128) {
        return 0.0;
    } else if (n < 6272) {
        M = 8;
    } else if (n < 750000) {
        M = 128;
    } else {
        M = 10000;
    }
    int N = n / M;
    if (N <= 0) {
        return 0.0;
    }
    int K;
    std::vector<int> v;
    std::vector<double> p;
    if (M == 8) {
        K = 3;
        v.assign(4, 0);
        p = {0.2148, 0.3672, 0.2305, 0.1875};
    } else if (M == 128) {
        K = 5;
        v.assign(6, 0);
        p = {0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124};
    } else if (M == 10000) {
        K = 6;
        v.assign(7, 0);
        p = {0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727};
    }
    for (int i = 0; i < N; ++i) {
        int maxRun = 0;
        int currentRun = 0;
        for (int j = 0; j < M; ++j) {
            if (bits[i * M + j] == 1) {
                currentRun++;
                if (currentRun > maxRun) {
                    maxRun = currentRun;
                }
            } else {
                currentRun = 0;
            }
        }
        if (M == 8) {
            if (maxRun <= 1) v[0]++;
            else if (maxRun == 2) v[1]++;
            else if (maxRun == 3) v[2]++;
            else if (maxRun >= 4) v[3]++;
        } else if (M == 128) {
            if (maxRun <= 4) v[0]++;
            else if (maxRun == 5) v[1]++;
            else if (maxRun == 6) v[2]++;
            else if (maxRun == 7) v[3]++;
            else if (maxRun == 8) v[4]++;
            else if (maxRun >= 9) v[5]++;
        } else if (M == 10000) {
            if (maxRun <= 10) v[0]++;
            else if (maxRun == 11) v[1]++;
            else if (maxRun == 12) v[2]++;
            else if (maxRun == 13) v[3]++;
            else if (maxRun == 14) v[4]++;
            else if (maxRun == 15) v[5]++;
            else if (maxRun >= 16) v[6]++;
        }
    }
    double chi_sq = 0.0;
    for (size_t i = 0; i < v.size(); ++i) {
        double expected = N * p[i];
        double diff = v[i] - expected;
        chi_sq += (diff * diff) / expected;
    }
    double p_value = gammaQ(K / 2.0, chi_sq / 2.0);
    return p_value;
}

/**
 * @brief Performs the cumulative sums test on a binary sequence (NIST SP 800-22)
 * @param bits Vector of bits (0s and 1s) to test
 * @return P-value indicating the probability of randomness
 * @details Tests the maximum deviation of the cumulative sum of the sequence from zero
 */
double cumulativeSumsTestP(const std::vector<int>& bits) {
    int n = bits.size();
    long cumSum = 0;
    long maxDep = 0;
    long minDep = 0;
    for (int b : bits) {
        cumSum += (b == 1 ? 1 : -1);
        if (cumSum > maxDep) maxDep = cumSum;
        if (cumSum < minDep) minDep = cumSum;
    }
    long maxAbs = std::max(std::labs(maxDep), std::labs(minDep));
    double zs = maxAbs / sqrt(n);
    double z = maxAbs;
    double p = 1.0;
    int k_start = (int)ceil(( - (double)n / z + 1.0) / 4.0);
    int k_end   = (int)floor((   (double)n / z - 1.0) / 4.0);
    for (int k = k_start; k <= k_end; ++k) {
        double term = Phi((4*k + 1) * zs) - Phi((4*k - 1) * zs);
        p -= term;
    }
    k_start = (int)ceil(( - (double)n / z - 3.0) / 4.0);
    k_end   = (int)floor((   (double)n / z - 1.0) / 4.0);
    for (int k = k_start; k <= k_end; ++k) {
        double term = Phi((4*k + 3) * zs) - Phi((4*k + 1) * zs);
        p += term;
    }
    return p;
}

/**
 * @brief Calculates basic statistics for a set of values
 * @param values Input vector of 32-bit unsigned integers
 * @param[out] mean Calculated mean value
 * @param[out] stddev Calculated standard deviation
 * @param[out] cov Calculated coefficient of variation (stddev/mean)
 */
void calculateStats(const std::vector<uint32_t>& values, double& mean, double& stddev, double& cov) {
    long double sum = 0.0;
    for (uint32_t v : values) {
        sum += v;
    }
    mean = sum / values.size();

    long double variance = 0.0;
    for (uint32_t v : values) {
        variance += (v - mean) * (v - mean);
    }
    variance /= values.size();

    stddev = std::sqrt(variance);
    cov = (mean != 0.0) ? stddev / mean : 0;
}

/**
 * @brief Performs chi-square goodness-of-fit test on random numbers
 * @param values Vector of 32-bit unsigned integers to test
 * @return P-value indicating the probability of uniform distribution
 * @details Tests if the numbers are uniformly distributed across the range
 */
double chiSquareTestP(const std::vector<uint32_t>& values) {
    int n = values.size();
    int k = 1 + floor(log2(n));

    std::vector<int> counts(k, 0); 

    const double divisor = (double(UINT32_MAX) + 1.0);

    for (long i = 0; i < n; i++) {
        double uptd_elem = double(values[i]) / divisor;
        int interval_num = int(floor(uptd_elem * k));
        counts[interval_num]++;
    }

    double chi2 = 0.0;
    double pN = double(n) / double(k);

    for (int j = 0; j < k; ++j) {
        double diff = double(counts[j]) - pN;
        chi2 += (diff * diff) / pN;
    }

    return gammaQ((k - 1) / 2.0, chi2 / 2.0);
}

/**
 * @brief Converts 32-bit values to a sequence of bits
 * @param values Vector of 32-bit unsigned integers
 * @return Vector of bits (0s and 1s) extracted from the input values
 * @details Extracts all 32 bits from each input value (MSB first)
 */
std::vector<int> generateBits(const std::vector<uint32_t>& values) {
    std::vector<int> bits;
    bits.reserve(values.size() * 32);
    for (uint32_t v : values) {
        for (int bit = 31; bit >= 0; --bit) {
            bits.push_back((v >> bit) & 1);
        }
    }
    return bits;
}

/**
 * @brief Performs battery of NIST statistical tests on bit sequence
 * @param bits Input bit sequence to test
 * @param[out] p1 P-value from monobit test
 * @param[out] p2 P-value from block frequency test
 * @param[out] p3 P-value from runs test
 * @param[out] p4 P-value from longest run of ones test
 * @param[out] p5 P-value from cumulative sums test
 */
void performNISTTests(const std::vector<int>& bits, double& p1, double& p2, double& p3, double& p4, double& p5) {
    p1 = monobitTestP(bits);
    p2 = blockFrequencyTestP(bits, 128);
    p3 = runsTestP(bits);
    p4 = longestRunOnesTestP(bits);
    p5 = cumulativeSumsTestP(bits);
}