#include <iostream>
#include <chrono>
#include <omp.h>
#include "grid.h"
#include "funcs.h"
#include "algo.h"

namespace tests {
    void check_time(Grid& grid, int num_threads = 1) {
        omp_set_num_threads(num_threads);
        auto start = std::chrono::high_resolution_clock::now();
        int iters = process(grid);
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "- Num threads: " << num_threads 
                  << ", iters: " << iters 
                  << ", execution time: " 
                  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms.\n";
    }

    void run_book_tests(int N, int block_size, double eps) {
        std::cout << "Executing book tests...\n";

        {
            Grid grid(book::f, book::g, N, block_size, eps);
            check_time(grid, 1);
        }

        {
            Grid grid(book::f, book::g, N, block_size, eps);
            check_time(grid, 8);
        }
    }

    void run_all_tests(int N = 1024, int block_size = 64, double eps = 0.1) {
        run_book_tests(N, block_size, eps);
    }
}

int main() {
    int N = 1024; // Grid size
    int block_size = 64; // Number of blocks
    double eps = 0.1; // Epsilon for error calculation

    tests::run_all_tests(N, block_size, eps);
}
