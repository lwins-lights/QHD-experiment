/*
 *  Problem 6 in Section 4 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5) while keeping the Lipschitz constant
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
    double z = x[0] * zooming;
    double f = pow(z, 6)
             - 15 * pow(z, 4)
             + 27 * pow(z, 2)
             + 250;
    return (f - 7) / zooming;
}