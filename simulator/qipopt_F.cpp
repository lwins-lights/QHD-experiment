#include <cmath>

using namespace std;

const double L = 1;
const int dim = 5;

void get_F_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

void get_F(const double *x, double *F) {
    double z[dim], s[dim];
    
    for (int i = 0; i < dim; i++) {
        z[i] = x[i] + 1;
    }
    s[0] = z[2] - z[3] + z[4];
    s[1] = z[3];
    s[2] = -z[0] + 2 * z[3];
    s[3] = z[0] - z[1] - 2 * z[2] + 3 * z[4];
    s[4] = -z[0] - 3 * z[3] + 5;
    for (int i = 0; i < dim; i++) {
        F[i] = z[i] * s[i];
    }
}