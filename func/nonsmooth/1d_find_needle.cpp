#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {-0.5, -0.5};        // lower bound of the actual function range for each dimension
const double ub[] = {0.5, 0.5};        // upper bound
const double compress_coef = 0.95;   // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;              // specifies how fast the encapsulated function will grow out of the "bound"


double get_obj(const double *x) {
    const int N = 2;
    /* map [-L, L) to [-\pi, \pi) */
    double z = x[0] / 0.5 * M_PI;

    /* make N valleys */
    double f = -cos(N * z);

    /* make the central one deeper */
    if (z < M_PI / N && z >= -M_PI / N) {
        f *= 2;
    } else {
        f += 1;
    }

    /* make minimum 0 and scale */
    return (f + 2) / N;
}

void get_obj_subg(const double *x, double *ret) {
    const int N = 2;
    double z = x[0] / 0.5 * M_PI;

    double f = sin(N * z) * N / 0.5 * M_PI;

    if (z < M_PI / N && z >= -M_PI / N) {
        f *= 2;
    }

    ret[0] = f / N;
}