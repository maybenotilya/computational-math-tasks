#ifndef GRID_HEADER
#define GRID_HEADER


class Grid {
private:
    static double** create_matrix(int N);

    static double** init_u(int N, double (*g)(double, double));

    static double** init_f(int N, double (*f)(double, double));

    static void free_matrix(double** matrix, int N);

public:
    int N;
    int block_size;
    double eps;
    double** u;
    double** f;
    double h;

    Grid(
        double (*func_f)(double, double), 
        double (*func_g)(double, double),
        int N,
        int block_size, 
        double eps
    );
    
    ~Grid();
};


#endif