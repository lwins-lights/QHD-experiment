#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {1e-8, 1e-8};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10};       // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {1.60086, 0.468498}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double part0 = pow(cos(x1), 2)*pow(cos(x2), 2);
    double part1 = fabs(pow(cos(x1), 4)+pow(cos(x2), 4) - 2.0*part0);
    double part2 = sqrt(pow(x1, 2) + 2*pow(x2, 2));
    return -part1/part2 + 0.672438;
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];

    double part0 = pow(cos(x1), 2)*pow(cos(x2), 2);
    double abs1 = pow(cos(x1), 4)+pow(cos(x2), 4) - 2.0*part0;

    if (abs1 >= 0) {
        ret[0] = ((pow(cos(x1), 2) - pow(cos(x2), 2)) * (4*(x1*x1+2*x2*x2)*sin(x1)*cos(x1) - x1*pow(cos(x2), 2) + x1*pow(cos(x1), 2))) / pow(x1*x1+2*x2*x2, 1.5);
        ret[1] = -(2*(pow(cos(x1), 2) - pow(cos(x2), 2)) * (cos(x2)*(2*(x1*x1+2*x2*x2)*sin(x2)+x2*cos(x2)) - x2*pow(cos(x1), 2))) / pow(x1*x1+2*x2*x2, 1.5);
    } else {
        ret[0] = -((pow(cos(x1), 2) - pow(cos(x2), 2)) * (4*(x1*x1+2*x2*x2)*sin(x1)*cos(x1) - x1*pow(cos(x2), 2) + x1*pow(cos(x1), 2))) / pow(x1*x1+2*x2*x2, 1.5);
        ret[1] = (2*(pow(cos(x1), 2) - pow(cos(x2), 2)) * (cos(x2)*(2*(x1*x1+2*x2*x2)*sin(x2)+x2*cos(x2)) - x2*pow(cos(x1), 2))) / pow(x1*x1+2*x2*x2, 1.5);
    }
}