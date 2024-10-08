#include <cmath>
#include <algorithm>

using namespace std;

const double L = 100;
const int dim = 2;

/* parameters for encapsulation */
const double lb[] = {0, 0};   // lower bound of the actual function range for each dimension
const double ub[] = {100, 100};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

int sgn(const double x) {
    return (x > 0) - (x < 0);
}

double get_obj(const double *x) {
    
    double x1 = x[0];
    double x2 = x[1];
    
    double zh1 = 9 - x1 - x2;
    double zh2 = (x1 - 3)*(x1 - 3) + (x2 - 2)*(x2 - 2) - 16;
    double zh3 = x1*x2 - 14;
    
    double term1 = zh1;
    double term2 = 100*(1 + zh2)*sgn(zh2);
    double term3 = 100*(1 + zh3)*sgn(zh3);
    double term4 = 100*(1 - x1)*sgn(x1);
    double term5 = 100*(1 - x2)*sgn(x2);
    
    
    return max({term1, term2, term3, term4, term5});
}

void get_obj_subg(const double *x, double *ret) {
    
    double x1 = x[0];
    double x2 = x[1];

    double zh1 = 9 - x1 - x2;
    double zh2 = (x1 - 3)*(x1 - 3) + (x2 - 2)*(x2 - 2) - 16;    // (x1 - 3)^2 + (x2 - 2)^2 - 16
    double zh3 = x1*x2 - 14;
    
    double term1 = zh1;
    double term2 = 100*(1 + zh2)*sgn(zh2);
    double term3 = 100*(1 + zh3)*sgn(zh3);
    double term4 = 100*(1 - x1)*sgn(x1);
    double term5 = 100*(1 - x2)*sgn(x2);

    double maxVal = max({term1, term2, term3, term4, term5});

    // Calculate subgradient based on the cases
    if (fabs(maxVal - term1) < 1e-8) {
        ret[0] = -1;
        ret[1] = -1;
    } else if (fabs(maxVal - term2) < 1e-8) {
        ret[0] = 2*(x1 - 3)*sgn(zh2);
        ret[1] = 2*(x2 - 2)*sgn(zh2);
    } else if (fabs(maxVal - term3) < 1e-8) {
        ret[0] = 100*x2*sgn(zh3);
        ret[1] = 100*x1*sgn(zh3);
    } else if (fabs(maxVal - term4) < 1e-8) {
        ret[0] = -100*sgn(x1);
        ret[1] = 0;
    } else if (fabs(maxVal - term5) < 1e-8) {
        ret[0] = 0;
        ret[1] = 100*sgn(x2);
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