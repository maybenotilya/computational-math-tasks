#include "grid.h"
#include <fstream>
#include <cstdlib>

using std::vector;

vector<vector<double>> Grid::init_u(int N, double (*g)(double, double)) {
    vector<vector<double>> u(N + 2, vector<double>(N + 2, 0.0));
    double h = 1.0 / (N + 1);
    for (int i = 0; i <= N + 1; ++i) {
        u[i][0] = g(h * i, 0);
        u[i][N + 1] = g(h * i, 1);
    }

    for (int j = 0; j <= N + 1; ++j) {
        u[0][j] = g(0, h * j);
        u[N + 1][j] = g(1, h * j);
    }

    // for (int i = 1; i < N + 1; ++i) {
    //     for (int j = 1; j < N + 1; ++j) {
    //         u[i][j] = -100 + std::rand() % 201;
    //     }
    // }

    return u;
}

vector<vector<double>> Grid::init_f(int N, double (*f)(double, double)) {
    vector<vector<double>> fmatrix(N + 2, vector<double>(N + 2, 0.0));
    double h = 1.0 / (N + 1);
    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            fmatrix[i][j] = f(h * i, h * j);
        }
    }
    return fmatrix;
}

Grid::Grid(
    double (*func_f)(double, double), 
    double (*func_g)(double, double),
    int N,
    int block_size, 
    double eps
) : N(N), block_size(block_size), eps(eps), u(init_u(N, func_g)), f(init_f(N, func_f)), h(1.0 / (N + 1)) {}

void Grid::save(const char* file_name) {
    std::ofstream fout(file_name);
    if (!fout.is_open()) {
        throw std::runtime_error("Failed to save grid.");
    }
    for (const auto& row : u) {
        for (const auto& el: row) {
            fout << el << ' ';
        }
        fout << std::endl;
    }
    fout.close();
}