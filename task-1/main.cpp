#include <iostream>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <cstring>

// Algorithm configuration
struct Config {
    int N;
    int NB;
    int BLOCK_SIZE;
    double eps;
    double h;
};

// To compare floats
const double EPSILON = 1e-19;

double** create_matrix(int N) {
    double** matrix = new double*[N];
    for (int i = 0; i < N; ++i) {
        matrix[i] = new double[N];
    }
    return matrix;
}

void free_matrix(double** matrix, int N) {
    for (int i = 0; i < N; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

double** init_u(int N, double (*g)(double, double)) {
    double** u = create_matrix(N + 2);
    double h = 1 / (N + 1);
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

double** init_f(int N, double (*f)(double, double)) {
    double** fmatrix = create_matrix(N + 2);
    double h = 1.0 / (N + 1);
    for (int i = 1; i < N + 1; ++i) {
        for (int j = 1; j < N + 1; ++j) {
            fmatrix[i][j] = f(h * i, h * j);
        }
    }
    return fmatrix;
}

double processBlock(double** u, double** f, Config& cfg, int block_i, int block_j) {
    double dmax = 0;
    for (int i = cfg.BLOCK_SIZE * block_i; i < cfg.BLOCK_SIZE * (block_i + 1) && i < cfg.N + 2; ++i) {
        // Avoid calculations on border
        if (i == 0 || i == cfg.N + 1) {
            continue;
        }
        for (int j = cfg.BLOCK_SIZE * block_j; j < cfg.BLOCK_SIZE * (block_j + 1) && j < cfg.N + 2; ++j) {
            // Avoid calculations on border
            if (j == 0 || j == cfg.N + 1) {
                continue;
            }
            double temp = u[i][j];
            u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - cfg.h * cfg.h * f[i][j]);
            double dm = fabs(temp - u[i][j]);
            if (dmax < dm) {dmax = dm;}
        }
    }
    return dmax;
}

int process(double** u, double** f, Config& cfg) {
    int iters = 0;
    double dmax = 0;
    double* dm = new double[cfg.NB];
    memset(dm, 0, cfg.NB * sizeof(cfg.NB));

    do {
        ++iters;
        dmax = 0;
        for (int nx = 0; nx < cfg.NB; ++nx) {
            dm[nx] = 0;
            int i;
            int j;
            double d;
#pragma omp parallel for shared(cfg, u, f, nx, dm) private(i, j, d) 
            for (i = 0; i < nx + 1; ++i) {
                j = nx - i;
                d = processBlock(u, f, cfg, i, j);
                if (dm[i] < d) {dm[i] = d;}
            }
        }

        for (int nx = cfg.NB - 2; nx > -1; --nx) {
            dm[nx] = 0;
            int i;
            int j;
            double d;
#pragma omp parallel for shared(cfg, u, f, nx, dm) private(i, j, d) 
            for (i = 0; i < nx + 1; ++i) {
                j = 2 * (cfg.NB - 1) - nx - i;
                d = processBlock(u, f, cfg, i, j);
                if (dm[i] < d) {dm[i] = d;}
            }
        }
        for (int i = 0; i < cfg.NB; ++i) {
            if (dmax < dm[i]) {dmax = dm[i];}
        }
    } while (dmax > cfg.eps); 
    delete[] dm;
    return iters;
}

void run_time_test(Config& cfg, double (*func_g)(double, double), double (*func_f)(double, double), int num_threads = 1) {
    double** u = init_u(cfg.N, func_g);
    double** f = init_f(cfg.N, func_f);

    omp_set_num_threads(num_threads);
    auto start = std::chrono::high_resolution_clock::now();
    int iters = process(u, f, cfg);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Num threads: " << num_threads 
              << ", iters: " << iters 
              << ", execution time: " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms.\n";

    free_matrix(u, cfg.N);
    free_matrix(f, cfg.N);
}

int main() {
    int N = 2000; // Grid size
    int NB = 40; // Number of blocks
    double eps = 0.001; // Epsilon for error calculation
    double h = 1.0 / (N + 1); // Grid step

    auto func_f = [](double x, double y){ return 0.0; }; 
    auto func_g = [](double x, double y){
        if (fabs(y) < EPSILON) {
            return 100 - 200 * x;
        }
        if (fabs(x) < EPSILON) {
            return 100 - 200 * y;
        }
        if (fabs(y - 1) < EPSILON) {
            return -100 + 200 * x;
        }
        if (fabs(x - 1) < EPSILON) {
            return -100 + 200 * y;
        }
        return 0.0;
    }; // Boundary values

    Config cfg = Config{N, NB, N / NB, eps, h};
    run_time_test(cfg, func_g, func_f, 1);
    run_time_test(cfg, func_g, func_f, 8);
}
