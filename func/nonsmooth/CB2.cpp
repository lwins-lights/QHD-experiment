#include <cmath>
#include <algorithm>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-5, -5};   // lower bound of the actual function range for each dimension
const double ub[] = {5, 5};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {1.139037652, 0.8995599}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    
    double x1 = x[0];
    double x2 = x[1];
    
    // Compute the two components of the function
    double term1 = x1 * x1 + pow(x2, 4);
    double term2 = pow(2 - x1, 2) + pow(2 - x2, 2);
    double term3 = 2 * exp(x2 - x1);

    // Return the maximum of the two components
    return max({term1, term2, term3}) - 1.9522245;;
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    // Compute the three components of the function
    double term1 = x1 * x1 + pow(x2, 4);
    double term2 = pow(2 - x1, 2) + pow(2 - x2, 2);
    double term3 = 2 * exp(x2 - x1);

    if (term1 >= term2 && term1 >= term3) {
        // gradient for term1
        ret[0] = 2*x1;
        ret[1] = 4*x2*x2*x2;
    } else if (term2 >= term1 && term2 >= term3) {
        // gradient for term2
        ret[0] = -2*(2-x1);
        ret[1] = -2*(2-x2);
    } else if (term3 >= term2 && term3 >= term1) {
        // gradient for term3
        ret[0] = -2 * exp(x2 - x1);
        ret[1] = 2 * exp(x2 - x1);
    }
}