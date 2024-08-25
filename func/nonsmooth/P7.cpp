#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-5, -5};   // lower bound of the actual function range for each dimension
const double ub[] = {5, 5};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = (x1+4)*(x1+4)+ x2*x2;
    double term2 = (x1-4)*(x1-4)+ x2*x2;
    double term3 = x1*x1 + (x2+2)*(x2+2);

    return max({term1, term2, term3}) - 16;
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = (x1+4)*(x1+4)+ x2*x2;
    double term2 = (x1-4)*(x1-4)+ x2*x2;
    double term3 = x1*x1 + (x2+2)*(x2+2);

    if (term1 >= term2 && term1 >= term3) {
        ret[0] = 2*(x1 + 4);
        ret[1] = 2*x2;
    } else if (term2 >= term1 && term2 >= term3) {
        ret[0] = 2*(x1 - 4);
        ret[1] = 2*x2;
    } else {
        ret[0] = 2*x1;
        ret[1] = 2*(x2+2);
    }
}