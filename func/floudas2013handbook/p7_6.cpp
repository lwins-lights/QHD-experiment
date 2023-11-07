/*
 *  Problem 6 in Section 7 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5)^3 while compressing the Lipschitz constant by a factor
 *  global minimum adjusted to 0 by displacement
 * 
 */

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 3;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double l2_penalty(const double x, const double c) {
    if (x < 0) {
        return 0;
    } else {
        return c * x * x;
    }
}

double get_potential(const double *x) {
    const double zooming = 99;
    const double compress = 10;
    const double t1 = x[0] * zooming + 50.5;
    const double t2 = x[1] * zooming + 50.5;
    const double t3 = x[2] * zooming + 50.5;
    double cons = 0.01 * t2 / t3 + 0.01 * t1 + 0.0005 * t1 * t3;
    double f = 0.5 * t1 / t2 - t1 - 5.0 / t2 / t2
             + l2_penalty(cons - 1, 10);
    return (f + 83.254) / zooming / compress;
}