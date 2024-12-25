#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-15, -15};   // lower bound of the actual function range for each dimension
const double ub[] = {30, 30};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0, 0};     

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);
    
    return -20*exp(-0.2*(sqrt(0.5*x1*x1 + 0.5*x2*x2))) - exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2))) + 20 + exp(1.0);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);

    if (x1 == 0 && x2 == 0) {
        ret[0] = 0;
        ret[1] = 0;
    } else {
        ret[0] = (2*sqrt(2)*x1*exp(-0.1*sqrt(2)*sqrt(x1*x1+x2*x2))) / (sqrt(x1*x1+x2*x2)) + pi*sin(2*pi*x1)*exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2)));
        ret[1] = (2*sqrt(2)*x2*exp(-0.1*sqrt(2)*sqrt(x1*x1+x2*x2))) / (sqrt(x1*x1+x2*x2)) + pi*sin(2*pi*x2)*exp(0.5*(cos(2*pi*x1)+cos(2*pi*x2)));
    }
}