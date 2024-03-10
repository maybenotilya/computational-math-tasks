#include "grid.h"
#include <cmath>
#include <iostream>
#include <omp.h>

double process_block(Grid& grid, int block_i, int block_j) {
    double dmax = 0;
    int i0 = 1 + block_i * grid.block_size;
    int j0 = 1 + block_j * grid.block_size;
    int i1 = std::min(i0 + grid.block_size, grid.N + 1);
    int j1 = std::min(j0 + grid.block_size, grid.N + 1);

    for (int i = i0; i < i1; ++i) {
        for (int j = j0; j < j1; ++j) {
            double temp = grid.u[i][j];
            grid.u[i][j] = 0.25 * (grid.u[i - 1][j] + grid.u[i + 1][j] + grid.u[i][j - 1] + grid.u[i][j + 1] - grid.h * grid.h * grid.f[i][j]);
            double dm = fabs(temp - grid.u[i][j]);
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

        for (int nx = NB - 2; nx >= 0; --nx) {
            int i;
            int j;
            double d;
#pragma omp parallel for shared(grid, nx, dm) private(i, j, d) 
            for (i = NB - nx - 1; i < NB; ++i) {
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
