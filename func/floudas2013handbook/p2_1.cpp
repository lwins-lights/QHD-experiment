/*
 *  Problem 1 in Section 2 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5)^5 while compressing the Lipschitz constant by a factor
 *  global minimum adjusted to 0 by displacement
 * 
 *  the 1st, 2nd and 4th component of vector x is reversed such that the QHD discretization 
 *  will not lose the global minimum
 */

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 5;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    const double compress = 10;
    const double inf = 250;
    const double x1 = -x[0] + 0.5; 
    const double x2 = -x[1] + 0.5;
    const double x3 = x[2] + 0.5;
    const double x4 = -x[3] + 0.5;
    const double x5 = x[4] + 0.5;
    double f = 42 * x1
             + 44 * x2
             + 45 * x3
             + 47 * x4
             + 47.5 * x5
             - 50 * (x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5);
    double cons = 20 * x1 + 12 * x2 + 11 * x3 + 7 * x4 + 4 * x5;
    if (cons > 40) {
        f = inf;
    }
    return (f + 17) / compress;
}