#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-15, -3};   // lower bound of the actual function range for each dimension
const double ub[] = {-5, 3};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {-10.0, 0.0}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    return 100*x2*x2 + 0.01*abs(x1+10);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    double abs1 = x1+10;

    if (abs1 >= 0) {
        ret[0] = 0.01;
        ret[1] = 200*x2;
    } else {
        ret[0] = -0.01;
        ret[1] = 200*x2;
    }
}