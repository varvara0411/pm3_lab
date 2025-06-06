/*#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <array>


class LCG {
private:
    uint32_t state; 

public:
    
    LCG(uint32_t seed) : state(seed) {}

    
    uint32_t next() {
        uint64_t res = uint64_t(1664525) * state + 1013904223u;
        state = static_cast<uint32_t>(res);
        return state;
    }
};


class Xorshift32 {
private:
    uint32_t state; 

public:
    
    Xorshift32(uint32_t seed) : state(seed) {
        if (state == 0) {
            state = 0xBAADF00D;  
        }
    }

    
    uint32_t next() {
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return state;
    }
};


class Xorshift128 {
private:
    uint64_t s[2]; 

public:
    
    Xorshift128(uint64_t seed) {
        s[0] = seed;
        s[1] = seed ^ 0x123456789ABCDEF0ULL;
        if ((s[0] | s[1]) == 0) {
            s[0] = 0xBAADF00DCAFEBABEULL;
            s[1] = 0xDEADBEEF01234567ULL;
        }
        // Прогрев генератора
        for (int i = 0; i < 20; ++i) {
            next();
        }
    }

    
    uint64_t next32() {
        uint64_t s1 = s[0];
        uint64_t s0 = s[1];
        s[0] = s0;
        s1 ^= s1 << 23;
        s[1] = s1 ^ s0 ^ (s0 >> 18) ^ (s1 >> 5);
        return s[1];
    }

    
    uint32_t next() {
        return static_cast<uint32_t>(next32());
    }
};


static double gammaQ(double a, double x);


static double incompleteGammaSeries(double a, double x) {
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


static double incompleteGammaContinuedFrac(double a, double x) {
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

static double gammaQ(double a, double x) {
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


static double Phi(double x) {
    return 0.5 * (1.0 + std::erf(x / M_SQRT2));
}


double monobitTestP(const std::vector<int>& bits) {
    int n = bits.size();
    long sum = 0;
    for (int b : bits) sum += b;
    long S = 2 * sum - n;
    double s_obs = fabs(S) / sqrt(n);
    double p = std::erfc(s_obs / M_SQRT2);
    return p;
}


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


int main() {
    const int SAMPLE_COUNT = 20;      
    const int SAMPLE_SIZE = 1000;     
    const uint32_t SEED_BASE1 = 1000; 
    const uint32_t SEED_BASE2 = 2000; 
    const uint32_t SEED_BASE3 = 3000; 
    const uint32_t SEED_BASE4 = 5000; 

    std::cout << std::fixed << std::setprecision(4);
    for (int method = 1; method <= 4; ++method) {
        std::string methodName;
        if (method == 1) methodName = "Custom PRNG 1 (LCG)";
        if (method == 2) methodName = "Custom PRNG 2 (Xorshift)";
        if (method == 3) methodName = "Custom PRNG 3 (Xorshift128)";
        if (method == 4) methodName = "std::mt19937";
        std::cout << "Results for " << methodName << ":\n";
        for (int s = 1; s <= SAMPLE_COUNT; ++s) {
            uint32_t seed = (method == 1 ? SEED_BASE1 : method == 2 ? SEED_BASE2 : method == 3 ? SEED_BASE3 : SEED_BASE4) + s;
            LCG lcg(0);
            Xorshift32 xorGen(0);
            Xorshift128 xor128Gen(0);
            std::mt19937 mt;
            if (method == 1) {
                lcg = LCG(seed);
            } else if (method == 2) {
                xorGen = Xorshift32(seed);
            } else if (method == 3) {
                xor128Gen = Xorshift128(seed);
            } else {
                mt.seed(seed);
            }
            std::vector<uint32_t> values;
            values.reserve(SAMPLE_SIZE);
            for (int i = 0; i < SAMPLE_SIZE; ++i) {
                uint32_t rv;
                if (method == 1) {
                    rv = lcg.next();
                } else if (method == 2) {
                    rv = xorGen.next();
                } else if (method == 3) {
                    rv = xor128Gen.next();
                } else {
                    rv = mt();
                }
                values.push_back(rv);
            }
            long double sum = 0.0;
            long double sumSq = 0.0;
            for (uint32_t v : values) {
                sum += v;
                sumSq += (long double)v * v;
            }
            long double mean = sum / SAMPLE_SIZE;
            long double var = (sumSq / SAMPLE_SIZE) - mean * mean;
            if (var < 0) var = 0;
            long double stddev = sqrt(var);
            long double cov = (mean != 0.0) ? stddev / mean : 0.0;
            std::vector<int> bits;
            bits.reserve(SAMPLE_SIZE * 32);
            for (uint32_t v : values) {
                for (int bit = 31; bit >= 0; --bit) {
                    bits.push_back((v >> bit) & 1U);
                }
            }
            double p_chi = 0.0;
            std::array<int, 16> freq{};
            for (uint32_t v : values) {
                unsigned index = v >> 28;
                freq[index]++;
            }
            double chi_stat = 0.0;
            double expected = (double)SAMPLE_SIZE / 16.0;
            for (int j = 0; j < 16; ++j) {
                double diff = freq[j] - expected;
                chi_stat += diff * diff / expected;
            }
            p_chi = gammaQ(15.0 / 2.0, chi_stat / 2.0);
            bool chi_pass = (p_chi > 0.05);

            //double p_chi = 0.0;
            //std::vector<int> freq;
            //int num_bins = 1 + static_cast<int>(std::log2(SAMPLE_SIZE)); // Правило Стерджеса
            //freq.resize(num_bins, 0);
            //double bin_size = (double)UINT32_MAX / num_bins;
            //for (uint32_t v : values) {
            //    unsigned index = static_cast<unsigned>(v / bin_size);
            //    if (index >= num_bins) index = num_bins - 1; // На случай округления
            //    freq[index]++;
            //}
            //double chi_stat = 0.0;
            //double expected = (double)SAMPLE_SIZE / num_bins;
            //for (int j = 0; j < num_bins; ++j) {
            //    double diff = freq[j] - expected;
            //    chi_stat += diff * diff / expected;
            //}
            //p_chi = gammaQ((num_bins - 1.0) / 2.0, chi_stat / 2.0);
            //bool chi_pass = (p_chi > 0.05);

            double p1 = monobitTestP(bits);
            double p2 = blockFrequencyTestP(bits, 128);
            double p3 = runsTestP(bits);
            double p4 = longestRunOnesTestP(bits);
            double p5 = cumulativeSumsTestP(bits);
            bool pass1 = (p1 > 0.01);
            bool pass2 = (p2 > 0.01);
            bool pass3 = (p3 > 0.01);
            bool pass4 = (p4 > 0.01);
            bool pass5 = (p5 > 0.01);
            std::cout << "Sample " << s << " - Mean: " << (double)mean
                      << ", StdDev: " << (double)stddev
                      << ", CoV: " << (double)cov << "\n";
            std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
                      << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
            std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
                      << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
                      << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
                      << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
                      << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
        }
        std::cout << std::string(80, '-') << "\n";
    }

    std::ofstream csv("prng_timing.csv");
    if (!csv) {
        std::cerr << "Error opening CSV file for writing.\n";
        return 1;
    }
    csv << "N,LCG,Xorshift,Xorshift128,MT19937\n";
    static volatile unsigned long long dummy_sum = 0ULL;
    for (int N = 1000; N <= 1000000; N += 100000) {
        LCG lcg(12345);
        Xorshift32 xorGen(12345);
        Xorshift128 acgGen(12335);
        std::mt19937 mt(12345);
        auto t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += lcg.next();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto lcg_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += xorGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto xor_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += acgGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto acg_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += mt();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto mt_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        double lcg_us = lcg_ns / 1000.0;
        double xor_us = xor_ns / 1000.0;
        double acg_us = acg_ns / 1000.0;
        double mt_us = mt_ns / 1000.0;
        csv << N << "," << std::fixed << std::setprecision(1)
            << lcg_us << "," << xor_us << "," << acg_us << "," << mt_us << "\n";
    }
    csv.close();
    return 0;
}*/

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <array>

/**
 * @brief Линейный конгруэнтный генератор псевдослучайных чисел (LCG)
 */
class LCG {
private:
    uint32_t state; /**< Внутреннее состояние генератора */

public:
    /**
     * @brief Конструктор LCG
     * @param seed Начальное значение (семя) генератора
     */
    LCG(uint32_t seed) : state(seed) {}

    /**
     * @brief Сгенерировать следующее псевдослучайное число
     * @return Следующее состояние (32-битное число)
     */
    uint32_t next() {
        uint64_t res = uint64_t(1664525) * state + 1013904223u;
        state = static_cast<uint32_t>(res);
        return state;
    }
};

/**
 * @brief Генератор псевдослучайных чисел Xorshift32
 */
class Xorshift32 {
private:
    uint32_t state; /**< Внутреннее состояние генератора */

public:
    /**
     * @brief Конструктор Xorshift32
     * @param seed Начальное значение (семя) генератора
     */
    Xorshift32(uint32_t seed) : state(seed) {
        if (state == 0) {
            state = 0xBAADF00D;  ///< Предотвращение состояния 0
        }
    }

    /**
     * @brief Сгенерировать следующее псевдослучайное число
     * @return Следующее состояние (32-битное число)
     */
    uint32_t next() {
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return state;
    }
};

/**
 * @brief Генератор псевдослучайных чисел методом серединных квадратов 
 */
class MiddleSquare {
private:
    uint32_t state;  // Текущее состояние (32 бита)

public:
    /**
     * @brief Конструктор MiddleSquare
     * @param seed Начальное значение (семя) генератора
     */
    MiddleSquare(uint32_t seed) {
        //if (state < 10000) {  // Защита от слишком маленького seed
        //    state = 0xBAADF00D;  // Используем 32-битное "магическое" число
        //}
        state = seed;
    }

    /**
     * @brief Генерация следующего псевдослучайного числа
     * @return 32-битное число
     */
    uint32_t next() {
        uint64_t squared = static_cast<uint64_t>(state) * state;
        
        // Берём середину 64-битного квадрата (32 бита из центра)
        state = static_cast<uint32_t>((squared >> 16) & 0xFFFFFFFF);
        
        return state;
    }
};

/**
 * @brief Вычисление верхней регуляризованной гамма-функции Q(a, x)
 * @param a Параметр a > 0
 * @param x Переменная x >= 0
 * @return Значение Q(a, x)
 */
static double gammaQ(double a, double x);

/**
 * @brief Вычисление несобственного интеграла гамма-функции через ряд
 * @param a Параметр a
 * @param x Переменная x
 * @return Нерегуляризованная функция нижней неполной гаммы
 */
static double incompleteGammaSeries(double a, double x) {
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
 * @brief Вычисление несобственного интеграла гамма-функции через непрерывную дробь
 * @param a Параметр a
 * @param x Переменная x
 * @return Нерегуляризованная функция верхней неполной гаммы
 */
static double incompleteGammaContinuedFrac(double a, double x) {
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

static double gammaQ(double a, double x) {
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
 * @brief Стандартная нормальная функция распределения Φ(x)
 * @param x Аргумент
 * @return Значение CDF нормального распределения
 */
static double Phi(double x) {
    return 0.5 * (1.0 + std::erf(x / M_SQRT2));
}

/**
 * @brief Тест однобитовой (monobit) равномерности
 * @param bits Вектор бит (0 или 1)
 * @return p-значение теста
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
 * @brief Тест блочной частотности (block frequency)
 * @param bits Вектор бит
 * @param M Размер блока (по умолчанию 128)
 * @return p-значение теста
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
 * @brief Тест числа серий (runs test)
 * @param bits Вектор бит
 * @return p-значение теста
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
 * @brief Тест самой длинной серии из единиц (longest run of ones)
 * @param bits Вектор бит
 * @return p-значение теста
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
 * @brief Тест кумулятивных сумм (cumulative sums test)
 * @param bits Вектор бит
 * @return p-значение теста
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
 * @brief Вычисление статистик (среднее, стандартное отклонение, коэффициент вариации)
 * @param values Вектор значений
 * @param mean Возвращаемое среднее значение
 * @param stddev Возвращаемое стандартное отклонение
 * @param cov Возвращаемый коэффициент вариации
 */
void calculateStats(const std::vector<uint32_t>& values, double& mean, double& stddev, double& cov) {
    long double sum = 0.0;
    long double sumSq = 0.0;
    for (uint32_t v : values) {
        sum += v;
        sumSq += (long double)v * v;
    }
    mean = sum / values.size();
    double var = (sumSq / values.size()) - mean * mean;
    if (var < 0) var = 0;
    stddev = sqrt(var);
    cov = (mean != 0.0) ? stddev / mean : 0.0;
}

/**
 * @brief Тест хи-квадрат на равномерность распределения
 * @param values Вектор значений
 * @return p-значение теста
 */
double chiSquareTestP(const std::vector<uint32_t>& values) {
    //std::array<int, 16> freq{};
    //for (uint32_t v : values) {
    //    unsigned index = v >> 28;
    //    freq[index]++;
    //}
    //double chi_stat = 0.0;
    //double expected = (double)values.size() / 16.0;
    //for (int j = 0; j < 16; ++j) {
    //    double diff = freq[j] - expected;
    //    chi_stat += diff * diff / expected;
    //}
    //return gammaQ(15.0 / 2.0, chi_stat / 2.0);
    
    
    if (values.empty()) return 0.0; // обработка пустого ввода
    
    // Вычисление числа интервалов по правилу Старджесса: k = 1 + 3.322 * log10(n)
    size_t n = values.size();
    size_t k = 1 + static_cast<size_t>(std::round(3.322 * std::log10(n)));
    k = std::max<size_t>(1, std::min<size_t>(k, 256)); // ограничиваем разумными пределами
    
    // Создание массива частот
    std::vector<int> freq(k, 0);
    
    // Распределение значений по корзинам
    for (uint32_t v : values) {
        // Масштабируем значение до диапазона [0, k-1]
        size_t index = static_cast<size_t>((static_cast<double>(v) / UINT32_MAX) * k);
        index = std::min(index, k - 1); // гарантируем, что не выйдем за границы
        freq[index]++;
    }
    
    // Вычисление статистики хи-квадрат
    double chi_stat = 0.0;
    double expected = static_cast<double>(n) / k;
    for (int count : freq) {
        double diff = count - expected;
        chi_stat += diff * diff / expected;
    }
    
    // Вычисление p-value (степени свободы = k - 1)
    return gammaQ((k - 1) / 2.0, chi_stat / 2.0);
}

/**
 * @brief Генерация битовой последовательности из значений
 * @param values Вектор значений
 * @return Вектор битов (0 или 1)
 */
std::vector<int> generateBits(const std::vector<uint32_t>& values) {
    std::vector<int> bits;
    bits.reserve(values.size() * 32);
    for (uint32_t v : values) {
        for (int bit = 31; bit >= 0; --bit) {
            bits.push_back((v >> bit) & 1U);
        }
    }
    return bits;
}

/**
 * @brief Выполнение статистических тестов NIST
 * @param bits Вектор битов
 * @param p1 Возвращаемое p-значение monobit теста
 * @param p2 Возвращаемое p-значение block frequency теста
 * @param p3 Возвращаемое p-значение runs теста
 * @param p4 Возвращаемое p-значение longest run теста
 * @param p5 Возвращаемое p-значение cumulative sums теста
 */
void performNISTTests(const std::vector<int>& bits, double& p1, double& p2, double& p3, double& p4, double& p5) {
    p1 = monobitTestP(bits);
    p2 = blockFrequencyTestP(bits, 128);
    p3 = runsTestP(bits);
    p4 = longestRunOnesTestP(bits);
    p5 = cumulativeSumsTestP(bits);
}

/**
 * @brief Точка входа: измерение характеристик генераторов и статистические тесты
 * @return Код возврата (0 при успешном завершении)
 */
int main() {
    //const int SAMPLE_COUNT = 20;      /**< Количество выборок на метод */
    //const int SAMPLE_SIZE = 1000;     /**< Размер одной выборки */
    //const uint32_t SEED_BASE1 = 1000; /**< Базовое семя для LCG */
    //const uint32_t SEED_BASE2 = 2000; /**< Базовое семя для Xorshift32 */
    //const uint32_t SEED_BASE3 = 3000; /**< Базовое семя для Xorshift128 */
    //const uint32_t SEED_BASE4 = 5000; /**< Базовое семя для std::mt19937 */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(10000, UINT32_MAX);

    const int SAMPLE_COUNT = 20;
    const int SAMPLE_SIZE = 1000;

    std::cout << std::fixed << std::setprecision(4);

    // Тестирование LCG
    std::cout << "Results for Custom PRNG 1 (LCG):\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        //uint32_t seed = SEED_BASE1 + s;
        //LCG lcg(seed);

        uint32_t seed = dist(gen);  // Случайный seed
        LCG lcg(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);
        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            values.push_back(lcg.next());
        }

        double mean, stddev, cov;
        calculateStats(values, mean, stddev, cov);
        double p_chi = chiSquareTestP(values);
        bool chi_pass = (p_chi > 0.01);

        std::vector<int> bits = generateBits(values);
        double p1, p2, p3, p4, p5;
        performNISTTests(bits, p1, p2, p3, p4, p5);
        
        bool pass1 = (p1 > 0.01);
        bool pass2 = (p2 > 0.01);
        bool pass3 = (p3 > 0.01);
        bool pass4 = (p4 > 0.01);
        bool pass5 = (p5 > 0.01);

        std::cout << "Sample " << s << " - Mean: " << mean
                  << ", StdDev: " << stddev
                  << ", CoV: " << cov << "\n";
        std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
                  << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
                  << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
                  << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
                  << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

    // Тестирование Xorshift32
    std::cout << "Results for Custom PRNG 2 (Xorshift):\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        //uint32_t seed = SEED_BASE2 + s;
        //Xorshift32 xorGen(seed);

        uint32_t seed = dist(gen);  // Случайный seed
        Xorshift32 xorGen(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);
        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            values.push_back(xorGen.next());
        }

        double mean, stddev, cov;
        calculateStats(values, mean, stddev, cov);
        double p_chi = chiSquareTestP(values);
        bool chi_pass = (p_chi > 0.01);

        std::vector<int> bits = generateBits(values);
        double p1, p2, p3, p4, p5;
        performNISTTests(bits, p1, p2, p3, p4, p5);
        
        bool pass1 = (p1 > 0.01);
        bool pass2 = (p2 > 0.01);
        bool pass3 = (p3 > 0.01);
        bool pass4 = (p4 > 0.01);
        bool pass5 = (p5 > 0.01);

        std::cout << "Sample " << s << " - Mean: " << mean
                  << ", StdDev: " << stddev
                  << ", CoV: " << cov << "\n";
        std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
                  << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
                  << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
                  << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
                  << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

    // Тестирование MiddleSquare (вместо Xorshift128)
    std::cout << "Results for MiddleSquare PRNG:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        //uint64_t seed = SEED_BASE3 + s;
        //MiddleSquare msGen(seed); // <- Заменяем Xorshift128 на MiddleSquare

        uint32_t seed = dist(gen);  // 32-битный seed
        MiddleSquare msGen(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);
        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            values.push_back(msGen.next()); // Генерация значений
        }

        double mean, stddev, cov;
        calculateStats(values, mean, stddev, cov);
        double p_chi = chiSquareTestP(values);
        bool chi_pass = (p_chi > 0.01);

        std::vector<int> bits = generateBits(values);
        double p1, p2, p3, p4, p5;
        performNISTTests(bits, p1, p2, p3, p4, p5);
    
        bool pass1 = (p1 > 0.01);
        bool pass2 = (p2 > 0.01);
        bool pass3 = (p3 > 0.01);
        bool pass4 = (p4 > 0.01);
        bool pass5 = (p5 > 0.01);

        std::cout << "Sample " << s << " - Mean: " << mean
                << ", StdDev: " << stddev
                << ", CoV: " << cov << "\n";
        std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
                << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
                << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
                << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
                << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
                << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

    // Тестирование Xorshift128
    //std::cout << "Results for Custom PRNG 3 (Xorshift128):\n";
    //for (int s = 1; s <= SAMPLE_COUNT; ++s) {
    //    uint32_t seed = SEED_BASE3 + s;
    //    Xorshift128 xor128Gen(seed);
    //    std::vector<uint32_t> values;
    //    values.reserve(SAMPLE_SIZE);
    //    for (int i = 0; i < SAMPLE_SIZE; ++i) {
    //        values.push_back(xor128Gen.next());
    //    }

    //    double mean, stddev, cov;
    //    calculateStats(values, mean, stddev, cov);
    //    double p_chi = chiSquareTestP(values);
     //   bool chi_pass = (p_chi > 0.05);

     //   std::vector<int> bits = generateBits(values);
    //    double p1, p2, p3, p4, p5;
     //   performNISTTests(bits, p1, p2, p3, p4, p5);
        
     //   bool pass1 = (p1 > 0.01);
     //   bool pass2 = (p2 > 0.01);
     //   bool pass3 = (p3 > 0.01);
     //   bool pass4 = (p4 > 0.01);
      //  bool pass5 = (p5 > 0.01);

        //std::cout << "Sample " << s << " - Mean: " << mean
                //  << ", StdDev: " << stddev
                //  << ", CoV: " << cov << "\n";
        //std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
         //         << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
        //std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
         //         << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
        //          << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
         //         << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
         //         << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
    //}
    //std::cout << std::string(80, '-') << "\n";

    // Тестирование std::mt19937
    std::cout << "Results for std::mt19937:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        //uint32_t seed = SEED_BASE4 + s;
        //std::mt19937 mt(seed);

        uint32_t seed = dist(gen);  // Случайный seed
        std::mt19937 mt(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);
        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            values.push_back(mt());
        }

        double mean, stddev, cov;
        calculateStats(values, mean, stddev, cov);
        double p_chi = chiSquareTestP(values);
        bool chi_pass = (p_chi > 0.01);

        std::vector<int> bits = generateBits(values);
        double p1, p2, p3, p4, p5;
        performNISTTests(bits, p1, p2, p3, p4, p5);
        
        bool pass1 = (p1 > 0.01);
        bool pass2 = (p2 > 0.01);
        bool pass3 = (p3 > 0.01);
        bool pass4 = (p4 > 0.01);
        bool pass5 = (p5 > 0.01);

        std::cout << "Sample " << s << " - Mean: " << mean
                  << ", StdDev: " << stddev
                  << ", CoV: " << cov << "\n";
        std::cout << "  Chi-square Uniform: " << (chi_pass ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p_chi << std::setprecision(4) << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED") << ", "
                  << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED") << ", "
                  << "Runs: " << (pass3 ? "PASSED" : "FAILED") << ", "
                  << "LongestRun: " << (pass4 ? "PASSED" : "FAILED") << "\n"
                  << "CumSum: " << (pass5 ? "PASSED" : "FAILED") << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

    // Тестирование производительности
    std::ofstream csv("times.csv");
    if (!csv) {
        std::cerr << "Error opening CSV file for writing.\n";
        return 1;
    }
    csv << "N,LCG,Xorshift,MiddleSquare,MT19937\n";
    static volatile unsigned long long dummy_sum = 0ULL;
    for (int N = 1000; N <= 1000000; N += 100000) {
        LCG lcg(12345);
        Xorshift32 xorGen(12345);
        //Xorshift128 acgGen(12335);
        MiddleSquare acgGen(12335);
        std::mt19937 mt(12345);
        
        auto t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += lcg.next();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto lcg_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();
        
        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += xorGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto xor_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += acgGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto acg_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        dummy_sum = 0ULL;
        for (int i = 0; i < N; ++i) {
            dummy_sum += mt();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto mt_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();
        
        csv << N << "," << std::fixed << std::setprecision(1)
            << lcg_ns << "," << xor_ns << "," << acg_ns << "," << mt_ns << "\n";
    }
    csv.close();
    return 0;
}