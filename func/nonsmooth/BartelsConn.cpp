#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-500, -500};   // lower bound of the actual function range for each dimension
const double ub[] = {500, 500};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"

double get_obj(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    
    return abs(x1*x1 + x2*x2 + x1*x2) + abs(sin(x1)) + abs(cos(x2)) - 1;
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];

    double abs1 = x1*x1 + x2*x2 + x1*x2;
    double abs2 = sin(x1);
    double abs3 = cos(x2);

    // Calculate subgradient based on the cases
    if (abs1 >= 0 && abs2 >= 0 && abs3 >= 0) {
        ret[0] = 2*x1 + x2 + cos(x1);
        ret[1] = 2*x2 + x1 - sin(x2);
    } else if (abs1 >= 0 && abs2 >= 0 && abs3 < 0) {
        ret[0] = 2*x1 + x2 + cos(x1);
        ret[1] = 2*x2 + x1 + sin(x2);
    } else if (abs1 >= 0 && abs2 < 0 && abs3 >= 0) {
        ret[0] = 2*x1 + x2 - cos(x1);
        ret[1] = 2*x2 + x1 - sin(x2);
    } else if (abs1 >= 0 && abs2 < 0 && abs3 < 0) {
        ret[0] = 2*x1 + x2 - cos(x1);
        ret[1] = 2*x2 + x1 + sin(x2);
    } else if (abs1 < 0 && abs2 >= 0 && abs3 >= 0) {
        ret[0] = -2*x1 - x2 + cos(x1);
        ret[1] = -2*x2 - x1 - sin(x2);
    } else if (abs1 < 0 && abs2 >= 0 && abs3 < 0) {
        ret[0] = -2*x1 - x2 + cos(x1);
        ret[1] = -2*x2 - x1 + sin(x2);
    } else if (abs1 < 0 && abs2 < 0 && abs3 >= 0) {
        ret[0] = -2*x1 - x2 - cos(x1);
        ret[1] = -2*x2 - x1 - sin(x2);
    } else if (abs1 < 0 && abs2 < 0 && abs3 < 0) {
        ret[0] = -2*x1 - x2 - cos(x1);
        ret[1] = -2*x2 - x1 + sin(x2);
    }
}

