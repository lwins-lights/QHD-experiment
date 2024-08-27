#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-1, -1};   // lower bound of the actual function range for each dimension
const double ub[] = {1, 1};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0, 0}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    return pow(abs(x1), 2) + pow(abs(x2), 3);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];


    if (x1 >= 0 && x2 >= 0) {
        ret[0] = 2*x1;
        ret[1] = 3*x2*x2;
    } else if (x1 >= 0 && x2 < 0) {
        ret[0] = 2*x1;
        ret[1] = -3*x2*x2;
    } else if (x1 < 0 && x2 >= 0) {
        ret[0] = 2*x1;
        ret[1] = 3*x2*x2;
    } else {
        ret[0] = 2*x1;
        ret[1] = -3*x2*x2;
    }
}

