#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 3;

/* parameters for encapsulation */
const double lb[] = {-10, -10, -10};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10, 10};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0.0, 2.0*acos(0.0), 0.0}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];

    return log(fabs(x1*x2) + 0.001) + cos(x1+x2) + log(fabs(x2*x3) + 0.001) + cos(x2+x3) - (log(0.001) - 1.0)*(dim - 1.0);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];

    double abs1 = x1*x2;
    double abs2 = x2*x3;

    bool zero1 = x1*x2 + 0.001 == 0;
    bool zero2 = x2*x3 + 0.001 == 0;

    if (abs1 >= 0 && abs2 >= 0) {
        ret[0] = zero1 ? 0 : x2 / (x1*x2 + 0.001) - sin(x1 + x2);
        ret[1] = (zero1 || zero2) ? 0 : x1 / (x1*x2 + 0.001) - sin(x1 + x2) + x3 / (x2*x3 + 0.001) - sin(x2 + x3);
        ret[2] = zero2 ? 0 : x2 / (x2*x3 + 0.001) - sin(x2 + x3);
    } else if (abs1 >= 0 && abs2 < 0) {
        ret[0] = zero1 ? 0 : x2 / (x1*x2 + 0.001) - sin(x1 + x2);
        ret[1] = (zero1 || zero2) ? 0 : x1 / (x1*x2 + 0.001) - sin(x1 + x2) - x3 / (-x2*x3 + 0.001) - sin(x2 + x3);
        ret[2] = zero2 ? 0 : -x2 / (-x2*x3 + 0.001) - sin(x2 + x3);
    } else if (abs1 < 0 && abs2 >= 0) {
        ret[0] = zero1 ? 0 : -x2 / (-x1*x2 + 0.001) - sin(x1 + x2);
        ret[1] = (zero1 || zero2) ? 0 : -x1 / (-x1*x2 + 0.001) - sin(x1 + x2) + x3 / (x2*x3 + 0.001) - sin(x2 + x3);
        ret[2] = zero2 ? 0 : x2 / (x2*x3 + 0.001) - sin(x2 + x3);
    } else {
        ret[0] = zero1 ? 0 : -x2 / (-x1*x2 + 0.001) - sin(x1 + x2);
        ret[1] = (zero1 || zero2) ? 0 : -x1 / (-x1*x2 + 0.001) - sin(x1 + x2) - x3 / (-x2*x3 + 0.001) - sin(x2 + x3);
        ret[2] = zero2 ? 0 : -x2 / (-x2*x3 + 0.001) - sin(x2 + x3);
    }
}