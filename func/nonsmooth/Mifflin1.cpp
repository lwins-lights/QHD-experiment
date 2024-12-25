#include <cmath>
#include <algorithm>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

// G = 200
double get_potential(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0] * 10;
    double x2 = x[1] * 10;
    
    // Compute the two components of the function
    double term1 = x1*x1 + x2*x2 - 1;
    double term2 = 0;

    // Return the maximum of the two components
    return 20*max({term1, term2}) - x1;
}

void get_subgradient(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0] * 10;
    double x2 = x[1] * 10;
    
    // Compute the two components of the function
    double term1 = x1*x1 + x2*x2 - 1;
    double term2 = 0;

    if (term1 == term2) {
        // gradient of term1
        ret[0] = 40*x1 - 1;
        ret[1] = 40*x2;
    } else if (term1 > term2) {
        // gradient of term1
        ret[0] = 40*x1 - 1;
        ret[1] = 40*x2;
    } else {
        ret[0] = -1;
        ret[1] = 0;
    }
}