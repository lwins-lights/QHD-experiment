/*
 *  Problem 5 in Section 4 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5) while keeping the Lipschitz constant
 *  global minimum adjusted to 0 by displacement
*/

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    const double zooming = 10;
    double x1 = x[0] * zooming;
    double x2 = x[1] * zooming;
    double f = 2 * pow(x1, 2)
             - 1.05 * pow(x1, 4)
             + pow(x1, 6) / 6
             - x1 * x2
             + pow(x2, 2);
    return f / zooming;
}