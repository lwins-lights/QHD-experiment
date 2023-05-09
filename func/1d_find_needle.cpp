/*
 *  one deepest valley among N consecutive sine valleys
 */

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 1;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    const int N = 8;
    /* map [-L, L) to [-\pi, \pi) */
    double z = x[0] / L * M_PI;

    /* make N valleys */
    double f = -cos(N * z);

    /* make the central one deeper */
    if (z < M_PI / N && z >= -M_PI / N) {
        f *= 2;
    }

    /* make minimum 0 and scale */
    return (f + 2) / N;
}