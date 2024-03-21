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
    std::vector<double> epss_;

    struct BenchResult {
        int32_t iters;
        int64_t time;
    };

    static BenchResult check_time(Grid& grid, int num_threads) {
        omp_set_num_threads(num_threads);
        auto start = std::chrono::high_resolution_clock::now();
        int iters = process(grid);
        auto end = std::chrono::high_resolution_clock::now();
        return {iters, std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()};
    }

public:

    Benchmark(
        std::vector<Function>& functions,
        std::vector<int>& nthreads,
        std::vector<int>& grid_sizes,
        std::vector<int>& block_sizes,
        std::vector<double>& epss
    ) : 
    functions_(functions), 
    nthreads_(nthreads), 
    grid_sizes_(grid_sizes), 
    block_sizes_(block_sizes), 
    epss_(epss) {}

    void run(int nruns = 5) {
        std::ofstream fout;
        for (const auto& [name, func_f, func_g] : functions_) {
            std::string fname = name + ".csv";
            fout.open(fname);
            std::cout << "Running bench for " << name << std::endl;
            fout 
            << "iterations,"
            << "nthreads," 
            << "n,"
            << "block_size,"
            << "eps,";
            for (int i = 0; i < nruns; ++i) {
                fout << "time_" << i + 1;
                if (i < nruns - 1) {
                    fout << ",";
                }
            }
            fout << std::endl;

            for (const auto eps: epss_) {
                for (const auto nthread: nthreads_) {
                    for (const auto N: grid_sizes_) {
                        for (const auto block_size: block_sizes_) {
                            std::vector<BenchResult> results;
                            for (int i = 0; i < nruns; ++i) {
                                auto grid = Grid(func_f, func_g, N, block_size, eps);
                                results.push_back(check_time(grid, nthread));
                            }
                            fout 
                            << results[0].iters << ","
                            << nthread << "," 
                            << N << ","
                            << block_size << ","
                            << eps << ",";
                            for (int i = 0; i < nruns; ++i) {
                                fout << results[i].time;
                                if (i < nruns - 1) {
                                    fout << ",";
                                }
                            }
                            fout << std::endl;
                        }
                    }
                }
            }
            fout.close();
        }
    }
};

int main() {
    std::srand(42);

    std::vector<Function> functions = {
        {"Book", book::f, book::g}, 
        {"Model1", model1::f, model1::g}
    };
    std::vector<int> nthreads = {1, 2, 4, 8};
    std::vector<int> grid_sizes = {1000};
    std::vector<int> block_sizes = {2, 4, 8, 16, 32, 64};
    std::vector<double> eps = {0.1};

    auto bench = Benchmark(
        functions,
        nthreads,
        grid_sizes,
        block_sizes,
        eps
    );

    bench.run(5);
}
