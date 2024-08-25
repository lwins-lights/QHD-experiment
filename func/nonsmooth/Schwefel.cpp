#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {-500, -500};   // lower bound of the actual function range for each dimension
const double ub[] = {500, 500};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"


double get_obj(const double *x) {
    double x1 = x[0];
    return 418.9828872724336 - x1*sin(sqrt(fabs(x1)));
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];

    if (x1 >= 0) {
        ret[0] = -1*sin(sqrt(x1)) - 0.5*sqrt(x1)*cos(sqrt(x1));
    } else {
        ret[0] = -1*sin(sqrt(-x1)) - 0.5*sqrt(-x1)*cos(sqrt(-x1));
    }
}