/// \file main.cpp
/// \brief Random Number Generator Testing Application
/// 
/// This program evaluates and compares four pseudo-random number generators:
/// 1. Linear Congruential Generator (LCG)
/// 2. Xorshift32
/// 3. Middle Square Method
/// 4. Mersenne Twister (std::mt19937)
///
/// The evaluation includes statistical tests and performance benchmarks.

#include "func.h"

/**
 * @brief Main function for random number generator testing
 * @return Exit code (0 for success)
 * 
 * @details The program performs the following operations:
 * 1. Generates test sequences for each PRNG algorithm
 * 2. Saves sequences to files for further analysis
 * 3. Conducts statistical analysis:
 *    - Basic statistics (mean, standard deviation, coefficient of variation)
 *    - Uniformity test (chi-square)
 *    - Randomness tests (NIST test suite)
 * 4. Measures and compares generation performance
 * 5. Outputs results in both human-readable and CSV formats
 */
int main() {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(10000, UINT32_MAX);

    const int SAMPLE_COUNT = 20; 
    const int SAMPLE_SIZE = 1000; 

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Results for LCG:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        
        uint32_t seed = dist(gen);  
        LCG lcg(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);

        std::string filename = "/home/install/pm3_lab/pm3_lab/samples/lcg_sample_" + std::to_string(s) + ".txt";
        std::ofstream outfile(filename);

        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            uint32_t num = lcg.next();
            values.push_back(num);
            outfile << num << "\n";  
        }

        outfile.close();  

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
                  << " (p=" << std::setprecision(3) << p_chi << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p1 << "), ";
        std::cout << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p2 << "), ";
        std::cout << "Runs: " << (pass3 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p3 << "), \n";
        std::cout << "LongestRun: " << (pass4 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p4 << "), ";
        std::cout << "CumSum: " << (pass5 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p5 << std::setprecision(4) << ")\n";
    }
    std::cout << std::string(80, '-') << "\n";

    
    std::cout << "Results for Xorshift32:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        

        uint32_t seed = dist(gen);  
        Xorshift32 xorGen(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);

        std::string filename = "/home/install/pm3_lab/pm3_lab/samples/xorshift32_sample_" + std::to_string(s) + ".txt";
        std::ofstream outfile(filename);

        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            uint32_t num = xorGen.next();
            values.push_back(num);
            outfile << num << "\n";  
        }

        outfile.close();  

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
                  << " (p=" << std::setprecision(3) << p_chi << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p1 << "), ";
        std::cout << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p2 << "), ";
        std::cout << "Runs: " << (pass3 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p3 << "), \n";
        std::cout << "LongestRun: " << (pass4 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p4 << "), ";
        std::cout << "CumSum: " << (pass5 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p5 << std::setprecision(4) << ")\n";
    }
    std::cout << std::string(80, '-') << "\n";

    
    std::cout << "Results for MiddleSquare:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        

        uint32_t seed = dist(gen);  
        MiddleSquare msGen(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);

        std::string filename = "/home/install/pm3_lab/pm3_lab/samples/middlesquare_sample_" + std::to_string(s) + ".txt";
        std::ofstream outfile(filename);

        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            uint32_t num = msGen.next();
            values.push_back(num);
            outfile << num << "\n";  
        }

        outfile.close();  

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
                  << " (p=" << std::setprecision(3) << p_chi << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p1 << "), ";
        std::cout << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p2 << "), ";
        std::cout << "Runs: " << (pass3 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p3 << "), \n";
        std::cout << "LongestRun: " << (pass4 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p4 << "), ";
        std::cout << "CumSum: " << (pass5 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p5 << std::setprecision(4) << ")\n";
    }
    std::cout << std::string(80, '-') << "\n";


    std::cout << "Results for std::mt19937:\n";
    for (int s = 1; s <= SAMPLE_COUNT; ++s) {
        

        uint32_t seed = dist(gen);  
        std::mt19937 mt(seed);

        std::vector<uint32_t> values;
        values.reserve(SAMPLE_SIZE);

        std::string filename = "/home/install/pm3_lab/pm3_lab/samples/mt19937_sample_" + std::to_string(s) + ".txt";
        std::ofstream outfile(filename);

        for (int i = 0; i < SAMPLE_SIZE; ++i) {
            uint32_t num = mt();
            values.push_back(num);
            outfile << num << "\n"; 
        }

        outfile.close(); 

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
                  << " (p=" << std::setprecision(3) << p_chi << "), ";
        std::cout << "Monobit: " << (pass1 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p1 << "), ";
        std::cout << "BlockFreq: " << (pass2 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p2 << "), ";
        std::cout << "Runs: " << (pass3 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p3 << "), \n";
        std::cout << "LongestRun: " << (pass4 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p4 << "), ";
        std::cout << "CumSum: " << (pass5 ? "PASSED" : "FAILED")
                  << " (p=" << std::setprecision(3) << p5 << std::setprecision(4) << ")\n";
    }
    std::cout << std::string(80, '-') << "\n";

    std::ofstream csv("times.csv");
    csv << "Size,LCG,Xorshift32,MiddleSquare,MT19937\n";

    const std::vector<int> time_test_sizes = {1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    for (int size : time_test_sizes) {
        LCG lcg(12345);
        Xorshift32 xorGen(12345);
        
        MiddleSquare acgGen(12335);
        std::mt19937 mt(12345);
        
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < size; ++i) {
            lcg.next();
        }
        auto t2 = std::chrono::high_resolution_clock::now();
        auto lcg_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();
        
        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < size; ++i) {
            xorGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto xor_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < size; ++i) {
            acgGen.next();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto acg_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();

        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < size; ++i) {
            mt();
        }
        t2 = std::chrono::high_resolution_clock::now();
        auto mt_ns = std::chrono::duration<double, std::micro>(t2 - t1).count();
        
        csv << size << "," << std::fixed << std::setprecision(4)
            << lcg_ns << "," << xor_ns << "," << acg_ns << "," << mt_ns << "\n";
    }
    csv.close();
    return 0;
}

