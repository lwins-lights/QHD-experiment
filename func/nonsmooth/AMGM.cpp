#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

// B = 50
double get_potential(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0] * 20;
    double x2 = x[1] * 20;
    
    return pow((0.5*(fabs(x1) + fabs(x2)) - sqrt(2.0)*(fabs(x1)*fabs(x2))), 2);
}

void get_subgradient(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0] * 20;
    double x2 = x[1] * 20;

    // Calculate subgradient based on the cases
    if (x1 >= 0 && x2 >= 0) {
        ret[0] = x2*(0.5 - sqrt(2.0) * x2) + x1*(4*x2*x2 - 2*sqrt(2.0)*x2 + 0.5);
        ret[1] = x1*(0.5 - sqrt(2.0) * x1) + x2*(4*x1*x1 - 2*sqrt(2.0)*x1 + 0.5);
    } else if (x1 >= 0 && x2 < 0) {
        ret[0] = x2*(-0.5 - sqrt(2.0) * x2) + x1*(4*x2*x2 + 2*sqrt(2.0)*x2 + 0.5);
        ret[1] = x1*(-0.5 + sqrt(2.0) * x1) + x2*(4*x1*x1 - 2*sqrt(2.0)*x1 + 0.5);
    } else if (x1 < 0 && x2 >= 0) {
        ret[0] = x2*(-0.5 + sqrt(2.0) * x2) + x1*(4*x2*x2 - 2*sqrt(2.0)*x2 + 0.5);
        ret[1] = x1*(-0.5 - sqrt(2.0) * x1) + x2*(4*x1*x1 + 2*sqrt(2.0)*x1 + 0.5);
    } else if (x1 < 0 && x2 < 0) {
        ret[0] = x2*(0.5 + sqrt(2.0)*x2) + x1*(4*x2*x2 + 2*sqrt(2.0)*x2 + 0.5);
        ret[1] = x1*(0.5 + sqrt(2.0)*x1) + x2*(4*x1*x1 + 2*sqrt(2.0)*x1 + 0.5);
    }
}