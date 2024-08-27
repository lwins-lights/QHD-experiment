#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-500, -500};   // lower bound of the actual function range for each dimension
const double ub[] = {500, 500};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {-300.3376328023, 500}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double t1 = sqrt(fabs(x2 + x1 + 1));
    double t2 = sqrt(fabs(x2 - x1 + 1));
    return x1*sin(t2)*cos(t1) + (x2 + 1)*sin(t1)*cos(t2) + 500.802160296664;
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    
    double abs1 = x2 + x1 + 1;
    double abs2 = x2 - x1 + 1;

    if (abs1 >= 0 && abs2 >= 0) {
        ret[0] = (x2+1)*sin(sqrt(abs2))*sin(sqrt(abs1)) / (2*sqrt(abs2)) - 
                (x1*sin(sqrt(abs2))*sin(sqrt(abs1))) / (2*sqrt(abs1)) - 
                (x1*cos(sqrt(abs2))*cos(sqrt(abs1))) / (2*sqrt(abs2)) + 
                (x2+1)*cos(sqrt(abs2))*cos(sqrt(abs1)) / (2*sqrt(abs1)) + 
                sin(sqrt(abs2))*cos(sqrt(abs1));
        ret[1] = -(x2+1)*sin(sqrt(abs2))*sin(sqrt(abs1)) / (2*sqrt(abs2)) - 
                (x1*sin(sqrt(abs2))*sin(sqrt(abs1))) / (2*sqrt(abs1)) +
                (x1*cos(sqrt(abs2))*cos(sqrt(abs1))) / (2*sqrt(abs2)) + 
                (x2+1)*cos(sqrt(abs2))*cos(sqrt(abs1)) / (2*sqrt(abs1)) + 
                sin(sqrt(abs1))*cos(sqrt(abs2));
    } else if (abs1 >= 0 && abs2 < 0) {
        ret[0] = -(x2+1)*sin(sqrt(-abs2))*sin(sqrt(abs1)) / (2*sqrt(-abs2)) - 
                (x1*sin(sqrt(-abs2))*sin(sqrt(abs1))) / (2*sqrt(abs1)) +
                (x1*cos(sqrt(-abs2))*cos(sqrt(abs1))) / (2*sqrt(-abs2)) + 
                (x2+1)*cos(sqrt(-abs2))*cos(sqrt(abs1)) / (2*sqrt(abs1)) + 
                sin(sqrt(-abs2))*cos(sqrt(abs1));
        ret[1] = (x2+1)*sin(sqrt(-abs2))*sin(sqrt(abs1)) / (2*sqrt(-abs2)) - 
                (x1*sin(sqrt(-abs2))*sin(sqrt(abs1))) / (2*sqrt(abs1)) -
                (x1*cos(sqrt(-abs2))*cos(sqrt(abs1))) / (2*sqrt(-abs2)) + 
                (x2+1)*cos(sqrt(-abs2))*cos(sqrt(abs1)) / (2*sqrt(abs1)) + 
                sin(sqrt(abs1))*cos(sqrt(-abs2));
    } else if (abs1 < 0 && abs2 >= 0) {
        ret[0] = (x2+1)*sin(sqrt(abs2))*sin(sqrt(-abs1)) / (2*sqrt(abs2)) + 
                (x1*sin(sqrt(abs2))*sin(sqrt(-abs1))) / (2*sqrt(-abs1)) - 
                (x1*cos(sqrt(abs2))*cos(sqrt(-abs1))) / (2*sqrt(abs2)) - 
                (x2+1)*cos(sqrt(abs2))*cos(sqrt(-abs1)) / (2*sqrt(-abs1)) + 
                sin(sqrt(abs2))*cos(sqrt(-abs1));
        ret[1] = -(x2+1)*sin(sqrt(abs2))*sin(sqrt(-abs1)) / (2*sqrt(abs2)) + 
                (x1*sin(sqrt(abs2))*sin(sqrt(-abs1))) / (2*sqrt(-abs1)) +
                (x1*cos(sqrt(abs2))*cos(sqrt(-abs1))) / (2*sqrt(abs2)) - 
                (x2+1)*cos(sqrt(abs2))*cos(sqrt(-abs1)) / (2*sqrt(-abs1)) + 
                sin(sqrt(-abs1))*cos(sqrt(abs2));
    } else {
        ret[0] = -(x2+1)*sin(sqrt(-abs2))*sin(sqrt(-abs1)) / (2*sqrt(-abs2)) + 
                (x1*sin(sqrt(-abs2))*sin(sqrt(-abs1))) / (2*sqrt(-abs1)) + 
                (x1*cos(sqrt(-abs2))*cos(sqrt(-abs1))) / (2*sqrt(-abs2)) - 
                (x2+1)*cos(sqrt(-abs2))*cos(sqrt(-abs1)) / (2*sqrt(-abs1)) + 
                sin(sqrt(-abs2))*cos(sqrt(-abs1));
        ret[1] = (x2+1)*sin(sqrt(-abs2))*sin(sqrt(-abs1)) / (2*sqrt(-abs2)) + 
                (x1*sin(sqrt(-abs2))*sin(sqrt(-abs1))) / (2*sqrt(-abs1)) -
                (x1*cos(sqrt(-abs2))*cos(sqrt(-abs1))) / (2*sqrt(-abs2)) - 
                (x2+1)*cos(sqrt(-abs2))*cos(sqrt(-abs1)) / (2*sqrt(-abs1)) + 
                sin(sqrt(-abs1))*cos(sqrt(-abs2));
    }
}