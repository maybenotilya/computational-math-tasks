#ifndef GRID_HEADER
#define GRID_HEADER

#include <vector>

using std::vector;

class Grid {
private:
    static vector<vector<double>> init_u(int N, double (*g)(double, double));

    static vector<vector<double>> init_f(int N, double (*f)(double, double));

public:
    int N;
    int block_size;
    double eps;
    vector<vector<double>> u;
    vector<vector<double>> f;
    double h;

    Grid(
        double (*func_f)(double, double), 
        double (*func_g)(double, double),
        int N,
        int block_size, 
        double eps
    );
    
    ~Grid() {};
};


#endif