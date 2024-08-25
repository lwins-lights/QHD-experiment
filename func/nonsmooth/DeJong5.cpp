#include <cmath>
#include <array>

using namespace std;

const double L = 0.5;
const int dim = 2;

void get_potential_params(double &var_L, int &var_dim) {
    var_L = L;
    var_dim = dim;
}

double get_potential(const double *x) {
    double x1 = x[0] * 131.072;
    double x2 = x[1] * 131.072;

    std::array<double, 25> a1 = {-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32};
    std::array<double, 25> a2 = {-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32};

    double sum = 0.0;

    for (int ii = 0; ii < 25; ++ii) {
        double a1i = a1[ii];
        double a2i = a2[ii];
        double term1 = static_cast<double>(ii + 1);
        double term2 = std::pow(x1 - a1i, 6);
        double term3 = std::pow(x2 - a2i, 6);
        double new_term = 1.0 / (term1 + term2 + term3);
        sum += new_term;
    }

    return 1.0 / (0.002 + sum) - 0.99800383779444934440;
}

void get_subgradient(const double *x, double *ret) {
    double x1 = x[0] * 131.072;
    double x2 = x[1] * 131.072;

    std::array<double, 25> a1 = {-32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32, -32, -16, 0, 16, 32};
    std::array<double, 25> a2 = {-32, -32, -32, -32, -32, -16, -16, -16, -16, -16, 0, 0, 0, 0, 0, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32};

    double sum = 0.0;

    for (int ii = 0; ii < 25; ++ii) {
        double a1i = a1[ii];
        double a2i = a2[ii];
        double term1 = static_cast<double>(ii + 1);
        double term2 = std::pow(x1 - a1i, 6);
        double term3 = std::pow(x2 - a2i, 6);
        double new_term = 1.0 / (term1 + term2 + term3);
        sum += new_term;
    }
    double denomenator = -1 / pow(0.02 + sum, 2);
    double xnum = 0;
    double ynum = 0;
    for (int ii = 0; ii < 25; ++ii) {
        double a1i = a1[ii];
        double a2i = a2[ii];
        double term1 = static_cast<double>(ii + 1);
        double term2 = std::pow(x1 - a1i, 6);
        double term3 = std::pow(x2 - a2i, 6);
        xnum += (-1 / pow(term1 + term2 + term3, 2)) * 6*pow(x1 - a1i, 5);
        ynum += (-1 / pow(term1 + term2 + term3, 2)) * 6*pow(x2 - a2i, 5);
    }

    ret[0] = denomenator * xnum;
    ret[1] = denomenator * ynum;
}

