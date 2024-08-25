#include <cmath>
#include "potential.hpp"
#include <algorithm>

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-100, -100};   // lower bound of the actual function range for each dimension
const double ub[] = {100, 100};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = abs(x1);
    double term2 = abs(x2);

    return max({term1, term2});
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    double term1 = abs(x1);
    double term2 = abs(x2);

    if (term1 >= term2) {
        if (x1 >= 0) {
            ret[0] = 1;
        } else {
            ret[0] = -1;
        }
        ret[1] = 0;
    } else {
        ret[0] = 0;
        if (x2 >= 0) {
            ret[1] = 1;
        } else {
            ret[1] = -1;
        }
    }
}