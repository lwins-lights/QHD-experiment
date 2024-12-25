#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 3;

/* parameters for encapsulation */
const double lb[] = {-10, -10};   // Lower bound of the actual function range for each dimension
const double ub[] = {10, 10};     // Upper bound
const double compress_coef = 0.95;  // The actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // Specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {0, 0};     // Dynamically adjustable pinned values

// Objective function for any dimension
double get_obj(const double *x) {
    double result = 0.0;
    for (int i = 0; i < dim - 1; ++i) {
        double x1 = x[i];
        double x2 = x[i + 1];
        result += (pow(sin(x1), 2) + pow(sin(x2), 2) - exp(-x1*x1-x2*x2)) * exp(-pow(sin(sqrt(abs(x1))), 2) - pow(sin(sqrt(abs(x2))), 2));
    }
    return result;
}

// Subgradient function for any dimension
void get_obj_subg(const double *x, double *ret) {
    for (int i = 0; i < dim; ++i) {
        ret[i] = 0.0; // Initialize gradient
    }

    for (int i = 0; i < dim - 1; ++i) {
        double x1 = x[i];
        double x2 = x[i + 1];

        double sqrt_x1 = (x1 >= 0) ? sqrt(x1) : sqrt(-x1);
        double sqrt_x2 = (x2 >= 0) ? sqrt(x2) : sqrt(-x2);

        if (x1 > 0) {
            ret[i] += exp(-pow(sin(sqrt_x1), 2) - pow(sin(sqrt_x2), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) - 
                 sin(sqrt_x1)*cos(sqrt_x1)*exp(-pow(sin(sqrt_x1), 2) - pow(sin(sqrt_x2), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt_x1;
        } else if (x1 < 0) {
            ret[i] += exp(-pow(sin(sqrt_x1), 2) - pow(sin(sqrt_x2), 2))*(2*x1*exp(-x1*x1-x2*x2) + 2*sin(x1)*cos(x1)) + 
                 sin(sqrt_x1)*cos(sqrt_x1)*exp(-pow(sin(sqrt_x1), 2) - pow(sin(sqrt_x2), 2))*(-exp(-x1*x1-x2*x2) + pow(sin(x1), 2) + pow(sin(x2), 2)) / sqrt_x1;
        }

        if (x2 > 0) {
            ret[i + 1] += exp(-pow(sin(sqrt_x2), 2) - pow(sin(sqrt_x1), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) - 
                 sin(sqrt_x2)*cos(sqrt_x2)*exp(-pow(sin(sqrt_x2), 2) - pow(sin(sqrt_x1), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt_x2;
        } else if (x2 < 0) {
            ret[i + 1] += exp(-pow(sin(sqrt_x2), 2) - pow(sin(sqrt_x1), 2))*(2*x2*exp(-x2*x2-x1*x1) + 2*sin(x2)*cos(x2)) + 
                 sin(sqrt_x2)*cos(sqrt_x2)*exp(-pow(sin(sqrt_x2), 2) - pow(sin(sqrt_x1), 2))*(-exp(-x2*x2-x1*x1) + pow(sin(x2), 2) + pow(sin(x1), 2)) / sqrt_x2;
        }
    }
}