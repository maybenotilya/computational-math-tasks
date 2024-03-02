#include "grid.h"

double** Grid::create_matrix(int N) {
    double** matrix = new double*[N];
    for (int i = 0; i < N; ++i) {
        matrix[i] = new double[N];
    }
    return matrix;
}

double** Grid::init_u(int N, double (*g)(double, double)) {
    double** u = create_matrix(N + 2);
    double h = 1.0 / (N + 1);
    for (int i = 0; i <= N + 1; ++i) {
        u[i][0] = g(h * i, 0);
        u[i][N + 1] = g(h * i, 1);
    }

    for (int j = 0; j <= N + 1; ++j) {
        u[0][j] = g(0, h * j);
        u[N + 1][j] = g(1, h * j);
    }

    return u;
}

double** Grid::init_f(int N, double (*f)(double, double)) {
    double** fmatrix = create_matrix(N + 2);
    double h = 1.0 / (N + 1);
    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            fmatrix[i][j] = f(h * i, h * j);
        }
    }
    return fmatrix;
}

void Grid::free_matrix(double** matrix, int N) {
    for (int i = 0; i < N; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

Grid::Grid(
    double (*func_f)(double, double), 
    double (*func_g)(double, double),
    int N,
    int block_size, 
    double eps
) : N(N), block_size(block_size), eps(eps), u(init_u(N, func_g)), f(init_f(N, func_f)), h(1.0 / (N + 1)) {}

Grid::~Grid() {
    free_matrix(u, N);
    free_matrix(f, N);
}