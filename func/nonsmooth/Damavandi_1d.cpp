#include <cmath>
#include <cstdio>
#include "potential.hpp"

using namespace std;

const int dim = 1;

/* parameters for encapsulation */
const double lb[] = {0};                // lower bound of the actual function range for each dimension
const double ub[] = {14};               // upper bound
const double compress_coef = 1.00;      // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 10;                // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {2.00};         // the pinned point will be guaranteed to be picked by the QHD discretization

const double width = 2;
const int sinc_taylor_cutoff = 10;
const double sinc_eps = 1e-6;

double sinc(double x) {
    if (abs(x) > sinc_eps) {
        return sin(x) / x;
    }
    double ret = 1;
    double fact = 1;    // fact = (2k+1)!
    for (int k = 1; k < sinc_taylor_cutoff; k++) {
        fact *= (2 * k) * (2 * k + 1);
        ret += pow(-1, k) * pow(x, 2 * k) / fact; 
    }
    return ret;
}

double Dsinc(double x) {
    if (abs(x) > sinc_eps) {
        return (cos(x) - sinc(x)) / x;
    }
    double ret = 0;
    double fact = 1;    // fact = (2k+1)!
    for (int k = 1; k < sinc_taylor_cutoff; k++) {
        fact *= (2 * k) * (2 * k + 1);
        ret += pow(-1, k) * (2 * k) * pow(x, 2 * k - 1) / fact; 
    }
    return ret;
}

double get_obj(const double *x) {
    double z = x[0];
    double pi = 2 * acos(0.0);
    double t = sinc(pi * (z - 2) / width);
    return (1 - pow(abs(t), 5)) * (2 + (z - 7) * (z - 7));
}

void get_obj_subg(const double *x, double *ret) {
    double z = x[0];
    double pi = 2 * acos(0.0);
    double t = sinc(pi * (z - 2) / width);
    double a = 1 - pow(abs(t), 5);
    double b = 2 + (z - 7) * (z - 7);
    double Db = 2 * (z - 7);
    double Da;
    if (t > 0) {
        Da = pi / width * Dsinc(pi * (z - 2) / width);
    } else {
        Da = - pi / width * Dsinc(pi * (z - 2) / width);
    }
    Da *= - 5 * pow(abs(t), 4);

    ret[0] = Da * b + a * Db;
}