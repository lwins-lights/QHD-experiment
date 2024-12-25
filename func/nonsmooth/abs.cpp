/*
    f(x) = |x|
*/

#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {-1};        // lower bound of the actual function range for each dimension
const double ub[] = {1};         // upper bound
const double compress_coef = 1;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1000;       // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0.0};   // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    double x1 = x[0];
    return fabs(x1);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];  
    if (x1 >= 0) {
        ret[0] = 1;
    } else {
        ret[0] = -1;
    }
}