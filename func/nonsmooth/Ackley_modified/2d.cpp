#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {0, 0};   // Lower bound of the actual function range for each dimension
const double ub[] = {14, 14};     // Upper bound
const double compress_coef = 0.95;  // The actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // Specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {7, 7};     // Dynamically adjustable pinned values

// Objective function for any dimension
double get_obj(const double *x) {
    double pi = 2*acos(0.0);
    double result = 0.0;
    for (int i = 0; i < dim - 1; ++i) {
        double x1 = x[i];
        double x2 = x[i + 1];
        result += (-20)*exp(-0.1*abs(x[i]))*exp(-0.1*abs(x[i+1])) - exp(cos(2*pi*x[i]) / 2) * exp(cos(2*pi*x[i+1]) / 2);
    }
    return result;
}

// Subgradient function for any dimension
void get_obj_subg(const double *x, double *ret) {
    double pi = 2*acos(0.0);
    for (int i = 0; i < dim; ++i) {
        ret[i] = 0.0; // Initialize gradient
    }

    for (int i = 0; i < dim - 1; ++i) {
        double x1 = x[i];
        double x2 = x[i + 1];

        if (x1 > 0) {
            ret[i] += 2*exp(-0.1*abs(x2)-0.1*x1) + pi*sin(2*pi*x1)*exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2)));
        } else if (x1 < 0) {
            ret[i] += -2*exp(0.1*x1 - 0.1*abs(x2)) + pi*sin(2*pi*x1)*exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2)));
        }

        if (x2 > 0) {
            ret[i + 1] += 2*exp(-0.1*abs(x1)-0.1*x2) + pi*sin(2*pi*x1)*exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2)));
        } else if (x2 < 0) {
            ret[i + 1] += -2*exp(0.1*x2 - 0.1*abs(x1)) + pi*sin(2*pi*x1)*exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2)));
        }
    }
}