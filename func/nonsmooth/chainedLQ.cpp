#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-5, -5};   // lower bound of the actual function range for each dimension
const double ub[] = {5, 5};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {sqrt(2.0)/2.0, sqrt(2.0)/2.0}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = -x1 - x2;
    double term2 = -x1 - x2 + x1*x1 + x2*x2 - 1;
    return max({term1, term2}) + sqrt(2.0);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    
    double term1 = -x1 - x2;
    double term2 = -x1 - x2 + x1*x1 + x2*x2 - 1;

    if (term1 >= term2) {
        ret[0] = -1;
        ret[1] = -1;
    } else {
        ret[0] = -1 + 2*x1;
        ret[1] = -1 + 2*x2;
    }
}