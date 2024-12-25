#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    return 4 * (x[0] * x[0] + x[1] * x[1]);
}

void get_subgradient(const double *x, double *ret) {
    ret[0] = 8 * x[0];
    ret[1] = 8 * x[1];
}