#include <cmath>
#include <cstdio>
#include "potential.hpp"

using namespace std;

void get_potential_params(int &var_dim) {
    var_dim = dim;
}

void get_pinned_point(double *ret, double L) {
    double scale[dim];              // the scaling factor from x to z
    double offset[dim];             // the offset factor from x to z

    for (int i = 0; i < dim; i++) {
        scale[i] = (ub[i] - lb[i]) / (2 * L * compress_coef);
        offset[i] = (ub[i] + lb[i]) / 2;
        ret[i] = (pinned[i] - offset[i]) / scale[i];
    }
}

double get_potential(const double *x, double L) {
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
    if (func_val != func_val) {
        printf("FATAL ERROR: Objective = NaN at (%f", z[0]);
        for (int i = 1; i < dim; i++) {
            printf(" %f", z[i]);
        }
        printf(")\n");
    }

    for (int i = 0; i < dim; i++) {
        func_val += abs(e[i]) * slope;
    }

    return func_val;
}

int sgn(const double x) {
    return (x > 0) - (x < 0);
}

void get_subgradient(const double *x, double *ret, double L) {
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
        if (func_subg[i] != func_subg[i]) {
            printf("FATAL ERROR: Subgradient = (%f", func_subg[0]);
            for (int j = 1; j < dim; j++) {
                printf(" %f", func_subg[j]);
            }
            printf(") at (%f", z[0]);
            for (int j = 1; j < dim; j++) {
                printf(" %f", z[j]);
            }
            printf(")\n");
            break;
        }
    }

    for (int i = 0; i < dim; i++) {
        if (abs(e[i]) != 0) {
            func_subg[i] = sgn(e[i]) * slope;
        }
        ret[i] = func_subg[i] * scale[i];
    }
}