#include <cmath>
#include <algorithm>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0] + 1.5;
    double x2 = x[1] - 2;
    
    // Compute the two components of the function
    double term1 = x1 * x1 + (x2 - 1) * (x2 - 1) + x2 - 1;
    double term2 = -x1 * x1 - (x2 - 1) * (x2 - 1) + x2 + 1;

    // Return the maximum of the two components
    return max(term1, term2);
}