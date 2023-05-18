/*
 *  Problem 4 in Section 4 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5) while compressing the Lipschitz constant by a factor
 *  global minimum adjusted to 0 by displacement
*/

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 1;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    const double zooming = 10;
    const double compress = 10;
    double z = x[0] * zooming;
    double f = pow(z, 2) * 4
             - pow(z, 3) * 4
             + pow(z, 4);
    return f / zooming / compress;
}