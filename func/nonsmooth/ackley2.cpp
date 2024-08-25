#include <cmath>

using namespace std;

const double L = 0.5;
const int dim = 1;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    double x1 = x[0] * 1000;
    return 418.9829*x1*sin(sqrt(fabs(x1))) - 2.545567497236334e-05;
}