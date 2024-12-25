#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0]*6;
    double x2 = x[1]*6;

    // Bound constrained chained Mifflin 2 function computation
    return -x1 + 2 * (x1 * x1 + x2 * x2 - 1) + 1.75 * fabs(x1 * x1 + x2 * x2 - 1) + 1;
}