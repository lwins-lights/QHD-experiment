#include <cmath>

using namespace std;

const double L = 10;
const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {1e-6, 1e-6};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 1;             // specifies how fast the encapsulated function will grow out of the "bound"

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

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