/*
 *  one deepest valley among N consecutive sine valleys
 */

#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 4;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    const int N = 2;
    double ret = 0;

    for (int i = 0; i < dim; i++) {
        /* map [-L, L) to [-\pi, \pi) */
        double z = x[i] / L * M_PI;

        /* make N valleys */
        double f = -cos(N * z);

        /* make the central one deeper */
        if (z < M_PI / N && z >= -M_PI / N) {
            f *= 2;
        } else {
            f += 1;
        }

        /* make minimum 0 and scale */
        ret += (f + 2) / N;
    }

    return ret;
}