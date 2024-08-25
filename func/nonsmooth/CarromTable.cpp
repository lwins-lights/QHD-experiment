#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-10, -10};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10};     // upper bound
const double compress_coef = 1.00;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {9.646157266349, 9.646157266349}; // the pinned point will be guaranteed to be picked by the QHD discretization


double get_obj(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);
    
    return -exp(2*abs(1 - sqrt(x1*x1 + x2*x2)/pi))*pow(cos(x1), 2)*pow(cos(x2), 2)/30.0 + 24.1568155165;
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);

    double abs1 = 1 - sqrt(x1*x1 + x2*x2)/pi;

    // Calculate subgradient based on the cases
    if (abs1 >= 0) {
        ret[0] = x1*exp(2*(1-sqrt(x1*x1 + x2*x2)/pi))*cos(x1)*cos(x1)*cos(x2)*cos(x2) / (15*pi*sqrt(x1*x1 + x2*x2)) + 
                    (1/15)*exp(2*(1-sqrt(x1*x1 + x2*x2)/pi))*sin(x1)*cos(x1)*cos(x2)*cos(x2);
        ret[1] = x2*exp(2*(1-sqrt(x2*x2 + x1*x1)/pi))*cos(x2)*cos(x2)*cos(x1)*cos(x1) / (15*pi*sqrt(x2*x2 + x1*x1)) + 
                    (1/15)*exp(2*(1-sqrt(x2*x2 + x1*x1)/pi))*sin(x2)*cos(x2)*cos(x1)*cos(x1);
    } else {
        ret[0] = -x1*exp(2*(-1+sqrt(x1*x1 + x2*x2)/pi))*cos(x1)*cos(x1)*cos(x2)*cos(x2) / (15*pi*sqrt(x1*x1 + x2*x2)) + 
                    (1/15)*exp(2*(-1+sqrt(x1*x1 + x2*x2)/pi))*sin(x1)*cos(x1)*cos(x2)*cos(x2);
        ret[1] = -x2*exp(2*(-1+sqrt(x2*x2 + x1*x1)/pi))*cos(x2)*cos(x2)*cos(x1)*cos(x1) / (15*pi*sqrt(x2*x2 + x1*x1)) + 
                    (1/15)*exp(2*(-1+sqrt(x2*x2 + x1*x1)/pi))*sin(x2)*cos(x2)*cos(x1)*cos(x1);
    } 
}