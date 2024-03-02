#include <cmath>

const double EPSILON = 1e-19;

namespace book {
    double f(double x, double y) {
        return 0.0;
    }

    double g(double x, double y) {
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
    }
}