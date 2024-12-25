#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 2;

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
        result += abs(x[i] * sin(x[i]) + 0.1 * x[i]) * (x[i+1] * sin(x[i+1]) + 0.1 * x[i+1]);
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

        double abs1 = x1 * sin(x1) + 0.1 * x1;
        double sign = (abs1 >= 0) ? 1.0 : -1.0;

        double term1 = x2 * (sin(x2) + 0.1) * (sin(x1) + x1 * cos(x1) + 0.1);
        double term2 = x1 * (sin(x1) + 0.1) * (sin(x2) + x2 * cos(x2) + 0.1);

        ret[i] += sign * term1;
        ret[i + 1] += sign * term2;
    }
}