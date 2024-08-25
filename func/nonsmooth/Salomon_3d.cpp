#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 3;

/* parameters for encapsulation */
const double lb[] = {-100, -100, -100};   // lower bound of the actual function range for each dimension
const double ub[] = {100, 100, 100};       // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0, 0, 0}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double pi = 2*acos(0.0);

    return 1 - cos(2*pi*sqrt(x1*x1 + x2*x2 + x3*x3)) + 0.1*sqrt(x1*x1 + x2*x2 + x3*x3);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double pi = 2*acos(0.0);

    ret[0] = x1*(2*pi*sin(2*pi*sqrt(x1*x1 + x2*x2 + x3*x3)) + 0.1) / sqrt(x1*x1 + x2*x2 + x3*x3);
    ret[1] = x2*(2*pi*sin(2*pi*sqrt(x1*x1 + x2*x2 + x3*x3)) + 0.1) / sqrt(x1*x1 + x2*x2 + x3*x3);
    ret[2] = x3*(2*pi*sin(2*pi*sqrt(x1*x1 + x2*x2 + x3*x3)) + 0.1) / sqrt(x1*x1 + x2*x2 + x3*x3);
}