#include <cmath>

const double EPSILON = 1e-19;

namespace book {
    double f(double x, double y) {
        return 0.0;
    }

    double g(double x, double y) {
        if (fabs(y) < EPSILON) {
            return 100.0 - 200.0 * x;
        }
        if (fabs(x) < EPSILON) {
            return 100.0 - 200.0 * y;
        }
        if (fabs(y - 1) < EPSILON) {
            return -100.0 + 200.0 * x;
        }
        if (fabs(x - 1) < EPSILON) {
            return -100.0 + 200.0 * y;
        }
        return 0.0;
    }
}

namespace model1 {
    double f(double x, double y) {
        return 5000.0 / (x * y);
    }

    double g(double x, double y) {
        return 100.0*x*x + 200.0*y*y*y;
    }
}