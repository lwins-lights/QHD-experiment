#include <cmath>
#include "potential.hpp"


using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {0, 0};   // lower bound of the actual function range for each dimension
const double ub[] = {14, 14};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {2.00, 2.00}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);
    double t1 = sin(pi * (x1 - 2)) / pi / (x1 - 2);
    double t2 = sin(pi * (x2 - 2)) / pi / (x2 - 2);
    
    if (x1 == 2) {
        t1 = 1;
    }
    if (x2 == 2) {
        t2 = 1;
    }
    
    return (1 - 2*pow(abs(t1 * t2), 5))*(2 + (x1 - 7)*(x1 - 7) + 2*(x2 - 7)*(x2 - 7)) + 77;
}

void get_obj_subg(const double *x, double *ret) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    double pi = 2*acos(0.0);

    double abs1 = (sin(pi*(x1 - 2)))*(sin(pi*(x2 - 2)))/(pi*pi*(x1 -2)*(x2 - 2));

    // Calculate subgradient based on the cases
    if (abs1 >= 0) {
        ret[0] = 2*(x1-7)*(1-2*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 5))) + 
                 (pow(x1-7, 2) + 2*pow(x2 - 7, 2) + 2) * 
                 (10*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 6)*pow(x2 - 2, 5)) - 
                  10*pow(sin(pi*(x1-2)), 4)*cos(pi*(x1 - 2))*pow(sin(pi*(x2-2)), 5) / (pow(pi, 9)*pow(x1 - 2, 5)*pow(x2 - 2, 5)));
        ret[1] = 4*(x2-7)*(1-2*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 5))) + 
                 (pow(x1-7, 2) + 2*pow(x2 - 7, 2) + 2) * 
                 (10*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 6)) - 
                  10*pow(sin(pi*(x1-2)), 5)*cos(pi*(x2 - 2))*pow(sin(pi*(x2-2)), 4) / (pow(pi, 9)*pow(x1 - 2, 5)*pow(x2 - 2, 5)));
    } else {
        ret[0] = 2*(x1-7)*(1+2*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 5))) + 
                 (pow(x1-7, 2) + 2*pow(x2 - 7, 2) + 2) * 
                 (-10*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 6)*pow(x2 - 2, 5)) + 
                  10*pow(sin(pi*(x1-2)), 4)*cos(pi*(x1 - 2))*pow(sin(pi*(x2-2)), 5) / (pow(pi, 9)*pow(x1 - 2, 5)*pow(x2 - 2, 5)));
        ret[1] = 4*(x2-7)*(1+2*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 5))) + 
                 (pow(x1-7, 2) + 2*pow(x2 - 7, 2) + 2) * 
                 (-10*pow(sin(pi*(x1-2)), 5)*pow(sin(pi*(x2-2)), 5) / (pow(pi, 10)*pow(x1 - 2, 5)*pow(x2 - 2, 6)) + 
                  10*pow(sin(pi*(x1-2)), 5)*cos(pi*(x2 - 2))*pow(sin(pi*(x2-2)), 4) / (pow(pi, 9)*pow(x1 - 2, 5)*pow(x2 - 2, 5)));
    } 
}