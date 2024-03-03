#include "grid.h"
#include <cmath>
#include <omp.h>

double process_block(Grid& grid, int block_i, int block_j) {
    double dmax = 0;
    auto& u = grid.u;
    auto& f = grid.f;
    for (int i = grid.block_size * block_i; i < grid.block_size * (block_i + 1) && i < grid.N + 2; ++i) {
        // Avoid calculations on border
        if (i == 0 || i == grid.N + 1) {
            continue;
        }
        for (int j = grid.block_size * block_j; j < grid.block_size * (block_j + 1) && j < grid.N + 2; ++j) {
            // Avoid calculations on border
            if (j == 0 || j == grid.N + 1) {
                continue;
            }
            double temp = u[i][j];
            u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - grid.h * grid.h * f[i][j]);
            double dm = fabs(temp - u[i][j]);
            if (dmax < dm) {dmax = dm;}
        }
    }
    return dmax;
}

int process(Grid& grid) {
    int iters = 0;
    double dmax = 0;
    int NB = grid.N / grid.block_size + (grid.N % grid.block_size == 0 ? 0 : 1);
    double* dm = new double[NB]{};

    do {
        ++iters;
        dmax = 0;
        for (int nx = 0; nx < NB; ++nx) {
            dm[nx] = 0;
            int i;
            int j;
            double d;
#pragma omp parallel for shared(grid, nx, dm) private(i, j, d) 
            for (i = 0; i < nx + 1; ++i) {
                j = nx - i;
                d = process_block(grid, i, j);
                if (dm[i] < d) {dm[i] = d;}
            }
        }

        for (int nx = NB - 2; nx > -1; --nx) {
            dm[nx] = 0;
            int i;
            int j;
            double d;
#pragma omp parallel for shared(grid, nx, dm) private(i, j, d) 
            for (i = 0; i < nx + 1; ++i) {
                j = 2 * (NB - 1) - nx - i;
                d = process_block(grid, i, j);
                if (dm[i] < d) {dm[i] = d;}
            }
        }
        for (int i = 0; i < NB; ++i) {
            if (dmax < dm[i]) {dmax = dm[i];}
        }
    } while (dmax > grid.eps); 
    delete[] dm;
    return iters;
}
