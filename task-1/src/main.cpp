#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <cstdlib>
#include <algorithm>

#include "grid.h"
#include "funcs.h"
#include "algo.h"

struct Function {
    std::string name;
    double (*func_f)(double, double);
    double (*func_g)(double, double);
};

class Benchmark {
private:
    std::vector<Function> functions_;
    std::vector<int> nthreads_;
    std::vector<int> grid_sizes_;
    std::vector<int> block_sizes_;
    double eps_;

    static int check_time(Grid& grid, int num_threads = 1) {
        omp_set_num_threads(num_threads);
        auto start = std::chrono::high_resolution_clock::now();
        int iters = process(grid);
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    static int mean_time(std::vector<int> times) {
        int sum = 0;
        for (const auto time: times) {
            sum += time;
        }
        return sum / times.size();
    }

    static int max_time(std::vector<int> times) {
        return *std::max_element(times.begin(), times.end());
    }

    static int min_time(std::vector<int> times) {
        return *std::min_element(times.begin(), times.end());
    }

public:

    Benchmark(
        std::vector<Function>& functions,
        std::vector<int> nthreads,
        std::vector<int> grid_sizes,
        std::vector<int> block_sizes,
        double eps
    ) : functions_(functions), nthreads_(nthreads), grid_sizes_(grid_sizes), block_sizes_(block_sizes), eps_(eps) {}

    void run(const char* fname, int nruns = 5) {
        std::ofstream fout(fname);
        if (!fout.is_open()) {
            throw std::runtime_error("Failed to open output file");
        }
        for (const auto& [name, func_f, func_g] : functions_) {
            std::cout << "Running bench for " << name << std::endl;
            fout << "Function: " << name << std::endl;
            fout 
            << "  " 
            << std::setw(10) << "Nthreads" 
            << std::setw(10) << "N" 
            << std::setw(18) << "Block_size" 
            << std::setw(20) << "Time (mean), ms" 
            << std::setw(20) << "Time (min), ms" 
            << std::setw(20) << "Time (max), ms" 
            << std::endl;

            for (const auto nthread: nthreads_) {
                for (const auto N: grid_sizes_) {
                    for (const auto block_size: block_sizes_) {
                        std::vector<int> times;
                        for (int i = 0; i < nruns; ++i) {
                            auto grid = Grid(func_f, func_g, N, block_size, eps_);
                            times.push_back(check_time(grid));
                        }
                        fout 
            << "  " 
            << std::setw(10) << nthread
            << std::setw(10) << N
            << std::setw(18) << block_size
            << std::setw(20) << mean_time(times) 
            << std::setw(20) << min_time(times) 
            << std::setw(20) << max_time(times) 
            << std::endl;
                    }
                }
            }
            fout << std::endl;
        }
    }
};

int main() {
    std::srand(42);

    std::vector<Function> functions = {{"Book function", book::f, book::g}, {"Model function 1", model1::f, model1::g}};
    std::vector<int> nthreads = {1, 4, 8};
    std::vector<int> grid_sizes = {500};
    std::vector<int> block_sizes = {64};
    double eps = 0.1;

    auto bench = Benchmark(
        functions,
        nthreads,
        grid_sizes,
        block_sizes,
        eps
    );

    bench.run("../bench.txt", 5);
}
