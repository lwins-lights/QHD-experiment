#include <cmath>

using namespace std;

const double L = 3;
const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {-500, -500};   // lower bound of the actual function range for each dimension
const double ub[] = {500, 500};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;             // specifies how fast the encapsulated function will grow out of the "bound"

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_obj(const double *x) {
    double x1 = x[0];
    double x2 = x[1];

    double t1 = sqrt(fabs(x2 + x1 + 1));
    double t2 = sqrt(fabs(x2 - x1 + 1));
    return x1*sin(t2)*cos(t1) + (x2 + 1)*sin(t1)*cos(t2) + 500.736340;
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
    } else if (abs1 < 0 && abs2 < 0) {
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