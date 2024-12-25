#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 1;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    double x1 = x[0] * 7.6;
    return floor(x1) + 4;
}

void get_subgradient(const double *x, double *ret) {
    ret[0] = 0;
}