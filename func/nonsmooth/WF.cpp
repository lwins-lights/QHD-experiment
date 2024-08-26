#include <cmath>
#include <algorithm>
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
    // Set bound 
    double x1 = x[0];
    double x2 = x[1]; 
    
    double term1 = 0.5 * (x1 + 10*x1/(x1+0.1) + 2*x2*x2);
    double term2 = 0.5 * (-1*x1 + 10*x1/(x1+0.1) + 2*x2*x2);
    double term3 = 0.5 * (x1 - 10*x1/(x1+0.1) - 2*x2*x2);


    // Return the maximum of the components
    return max({term1, term2, term3});
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = 0.5 * (x1 + 10*x1/(x1+0.1) + 2*x2*x2);
    double term2 = 0.5 * (-1*x1 + 10*x1/(x1+0.1) + 2*x2*x2);
    double term3 = 0.5 * (x1 - 10*x1/(x1+0.1) + 2*x2*x2);

    if (term1 >= term2 && term1 >= term3) {
        // gradient for term1
        ret[0] = (0.5*x1*x1 + 0.1*x1 + 0.505) / pow(x1 + 0.1, 2);
        ret[1] = 2*x2;
    } else if (term2 >= term1 && term2 >= term3) {
        // gradient for term2
        ret[0] = (-0.5*x1*x1 - 0.1*x1 + 0.495) / pow(x1 + 0.1, 2);
        ret[1] = 2*x2;
    } else if (term3 >= term2 && term3 >= term1) {
        // gradient for term3
        ret[0] = (0.5*x1*x1 + 0.1*x1 - 0.495) / pow(x1 + 0.1, 2);
        ret[1] = 2*x2;
    }
}