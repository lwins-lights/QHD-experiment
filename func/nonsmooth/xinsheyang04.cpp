#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-10, -10};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0, 0};     // *** ERRORNEOUS: THIS NEEDS TO BE FIXED LATER ***

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];
    
    return (pow(sin(x1), 2) + pow(sin(x2), 2) - exp(-x1*x1-x2*x2)) * exp(-pow(sin(sqrt(abs(x1))), 2) - pow(sin(sqrt(abs(x2))), 2)) + 1;
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];

    // Calculate subgradient based on the cases
    if (x1 >= 0 && x2 >= 0) {
        ret[0] = exp(-pow(sin(sqrt(x1)), 2) - pow(sin(sqrt(x2)), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) - 
                 sin(sqrt(x1))*cos(sqrt(x1))*exp(-pow(sin(sqrt(x1)), 2) - pow(sin(sqrt(x2)), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt(x1);
        ret[1] = exp(-pow(sin(sqrt(x2)), 2) - pow(sin(sqrt(x1)), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) - 
                 sin(sqrt(x2))*cos(sqrt(x2))*exp(-pow(sin(sqrt(x2)), 2) - pow(sin(sqrt(x1)), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt(x2);
    } else if (x1 >= 0 && x2 < 0) {
        ret[0] = exp(-pow(sin(sqrt(x1)), 2) - pow(sin(sqrt(-x2)), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) - 
                 sin(sqrt(x1))*cos(sqrt(x1))*exp(-pow(sin(sqrt(x1)), 2) - pow(sin(sqrt(-x2)), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt(x1);
        ret[1] = exp(-pow(sin(sqrt(-x2)), 2) - pow(sin(sqrt(x1)), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) + 
                 sin(sqrt(-x2))*cos(sqrt(-x2))*exp(-pow(sin(sqrt(-x2)), 2) - pow(sin(sqrt(x1)), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt(-x2);
    } else if (x1 < 0 && x2 >= 0) {
        ret[0] = exp(-pow(sin(sqrt(-x1)), 2) - pow(sin(sqrt(x2)), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) + 
                 sin(sqrt(-x1))*cos(sqrt(-x1))*exp(-pow(sin(sqrt(-x1)), 2) - pow(sin(sqrt(x2)), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt(-x1);
        ret[1] = exp(-pow(sin(sqrt(x2)), 2) - pow(sin(sqrt(-x1)), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) - 
                 sin(sqrt(x2))*cos(sqrt(x2))*exp(-pow(sin(sqrt(x2)), 2) - pow(sin(sqrt(-x1)), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt(x2);
    } else {
        ret[0] = exp(-pow(sin(sqrt(-x1)), 2) - pow(sin(sqrt(-x2)), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) + 
                 sin(sqrt(-x1))*cos(sqrt(-x1))*exp(-pow(sin(sqrt(-x1)), 2) - pow(sin(sqrt(-x2)), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt(-x1);
        ret[1] = exp(-pow(sin(sqrt(-x2)), 2) - pow(sin(sqrt(-x1)), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) + 
                 sin(sqrt(-x2))*cos(sqrt(-x2))*exp(-pow(sin(sqrt(-x2)), 2) - pow(sin(sqrt(-x1)), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt(-x2);
    }
}