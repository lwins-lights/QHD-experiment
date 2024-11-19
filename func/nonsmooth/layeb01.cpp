#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {-5};   // lower bound of the actual function range for each dimension
const double ub[] = {5};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10000000;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {1.0}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    double x1 = x[0];
    return 100.0*sqrt(fabs(exp(pow(x1-1.0, 2)) - 1.0));
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];

    double abs1 = exp(pow(x1-1.0, 2)) - 1.0;

    if (abs1 == 0) {
        ret[0] = 0;
    } else if (abs1 > 0) {
        ret[0] = 100*exp(pow(x1-1, 2))*(x1-1) / sqrt(exp(pow(x1-1, 2)) - 1);
    } else {
        ret[0] = -100*exp(pow(x1-1, 2))*(x1-1) / sqrt(-1*exp(pow(x1-1, 2)) + 1);
    }
}