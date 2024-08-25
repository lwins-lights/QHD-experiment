#include <cmath>
#include <algorithm>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-15, -15};   // lower bound of the actual function range for each dimension
const double ub[] = {15, 15};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 5;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {-10, 15}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    
    return 0.0001*pow((fabs(sin(x1)*sin(x2)*exp(100 - sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI)) + 1), (0.1)) - 0.0001;
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];

    double abs1 = sin(x1)*sin(x2);

    // Calculate subgradient based on the cases
    if (abs1 >= 0) {
        ret[0] = 0.00001*(exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*cos(x1)*sin(x2) - x1*exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*sin(x1)*sin(x2)/(M_PI*sqrt(x1*x1 + x2*x2))) / 
                pow((sin(x1)*sin(x2)*exp(100 - sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI) + 1), 0.9);
        ret[1] = 0.00001*(exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*cos(x2)*sin(x1) - x2*exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*sin(x1)*sin(x2)/(M_PI*sqrt(x1*x1 + x2*x2))) / 
                pow((sin(x1)*sin(x2)*exp(100 - sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI) + 1), 0.9);
    } else {
        ret[0] = 0.00001*(-exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*cos(x1)*sin(x2) + x1*exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*sin(x1)*sin(x2)/(M_PI*sqrt(x1*x1 + x2*x2))) / 
                pow((-sin(x1)*sin(x2)*exp(100 - sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI) + 1), 0.9);;
        ret[1] = 0.00001*(-exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*cos(x2)*sin(x1) + x2*exp(100 - sqrt(x1*x1 + x2*x2)/M_PI)*sin(x1)*sin(x2)/(M_PI*sqrt(x1*x1 + x2*x2))) / 
                pow((-sin(x1)*sin(x2)*exp(100 - sqrt(pow(x1, 2) + pow(x2, 2)) / M_PI) + 1), 0.9);;
    } 
}
