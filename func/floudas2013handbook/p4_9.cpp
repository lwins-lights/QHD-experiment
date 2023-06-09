/*
 *  Problem 9 in Section 4 of Handbook of Test Problems in Local and Global Optimization
 *  
 *  domain zoomed to [-0.5, 0.5)^2 while compressing the Lipschitz constant by a factor
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

double l2_penalty(const double x, const double c) {
    if (x < 0) {
        return 0;
    } else {
        return c * x * x;
    }
}

double get_potential(const double *x) {
    const double compress = 10;
    double x1 = (x[0] + 0.5) * 3;
    double x2 = (x[1] + 0.5) * 4;
    double e1 = 2 + 2 * pow(x1, 4) - 8 * pow(x1, 3) + 8 * pow(x1, 2);
    double e2 = 4 * pow(x1, 4) - 32 * pow(x1, 3) + 88 * pow(x1, 2) - 96 * x1 + 36;
    double f = - x1 - x2 
               + l2_penalty(x2 - e1, 1) 
               + l2_penalty(x2 - e2, 1);
    return (f + 5.50796) / compress;
}