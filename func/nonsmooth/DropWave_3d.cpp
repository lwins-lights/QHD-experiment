#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 3;

/* parameters for encapsulation */
const double lb[] = {-5.12, -5.12, -5.12};   // lower bound of the actual function range for each dimension
const double ub[] = {5.12, 5.12, 5.12};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0.00, 0.00, 0.00}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    
    double SumOfPower = pow(x1, 2) + pow(x2, 2) + pow(x3, 2);

    return -(1.00 + cos(12*sqrt(SumOfPower))) / (2 + 0.5*SumOfPower) + 1.00;
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];

    double SumOfPower = pow(x1, 2) + pow(x2, 2) + pow(x3, 2);

    if ((0.5*SumOfPower + 2) == 0 || sqrt(SumOfPower) * (0.5*SumOfPower + 2) == 0) {
        ret[0] = 0.00;
        ret[1] = 0.00;
        ret[2] = 0.00;
    } else {
        ret[0] = 12*x1*sin(12*sqrt(SumOfPower)) / (sqrt(SumOfPower) * (0.5*SumOfPower + 2)) 
            + x1*(cos(12*sqrt(SumOfPower)) + 1.00) / pow((0.5*SumOfPower + 2), 2);
        ret[1] = 12*x2*sin(12*sqrt(SumOfPower)) / (sqrt(SumOfPower) * (0.5*SumOfPower + 2)) 
                + x2*(cos(12*sqrt(SumOfPower)) + 1.00) / pow((0.5*SumOfPower + 2), 2);
        ret[2] = 12*x3*sin(12*sqrt(SumOfPower)) / (sqrt(SumOfPower) * (0.5*SumOfPower + 2)) 
                + x3*(cos(12*sqrt(SumOfPower)) + 1.00) / pow((0.5*SumOfPower + 2), 2);
    }
}