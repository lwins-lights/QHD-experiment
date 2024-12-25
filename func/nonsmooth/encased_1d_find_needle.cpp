#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {-L, -L};        // lower bound of the actual function range for each dimension
const double ub[] = {L, L};        // upper bound
const double compress_coef = 0.95;   // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_obj(const double *x) {
    const int N = 2;
    /* map [-L, L) to [-\pi, \pi) */
    double z = x[0] / L * M_PI;

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
    double z = x[0] / L * M_PI;

    double f = sin(N * z) * N / L * M_PI;

    if (z < M_PI / N && z >= -M_PI / N) {
        f *= 2;
    }

    ret[0] = f / N;
}

double get_potential(const double *x) {
    double z[dim];                  // true coordinates to be mapped
    double e[dim];                  // how far it is from outside the bound
    double scale[dim];              // the scaling factor from x to z
    double offset[dim];             // the offset factor from x to z
    double func_val;

    for (int i = 0; i < dim; i++) {
        scale[i] = (ub[i] - lb[i]) / (2 * L * compress_coef);
        offset[i] = (ub[i] + lb[i]) / 2;
        z[i] = x[i] * scale[i] + offset[i];
        if (z[i] > ub[i]) {
            e[i] = z[i] - ub[i];
            z[i] = ub[i];
        } else if (z[i] < lb[i]) {
            e[i] = z[i] - lb[i];
            z[i] = lb[i];
        } else {
            e[i] = 0;
        }
    }

    func_val = get_obj(z);
    for (int i = 0; i < dim; i++) {
        func_val += abs(e[i]) * slope;
    }

    return func_val;
}

int sgn(const double x) {
    return (x > 0) - (x < 0);
}

void get_subgradient(const double *x, double *ret) {
    double z[dim];                  // true coordinates to be mapped
    double e[dim];                  // how far it is from outside the bound
    double scale[dim];              // the scaling factor from x to z
    double offset[dim];             // the offset factor from x to z
    double func_subg[dim];

    for (int i = 0; i < dim; i++) {
        scale[i] = (ub[i] - lb[i]) / (2 * L * compress_coef);
        offset[i] = (ub[i] + lb[i]) / 2;
        z[i] = x[i] * scale[i] + offset[i];
        if (z[i] > ub[i]) {
            e[i] = z[i] - ub[i];
            z[i] = ub[i];
        } else if (z[i] < lb[i]) {
            e[i] = z[i] - lb[i];
            z[i] = lb[i];
        } else {
            e[i] = 0;
        }
    }

    get_obj_subg(z, func_subg);
    for (int i = 0; i < dim; i++) {
        if (abs(e[i]) != 0) {
            func_subg[i] = sgn(e[i]) * slope;
        }
        ret[i] = func_subg[i] * scale[i];
    }
}