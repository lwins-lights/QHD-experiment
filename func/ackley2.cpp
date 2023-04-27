#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    double x1 = x[0] * 64;
    double x2 = x[1] * 64;
    return -200 / 64 * exp(-0.2*sqrt(x1*x1+x2*x2));
}