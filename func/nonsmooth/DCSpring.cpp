#include <cmath>
#include "potential.hpp"

using namespace std;

const int dim = 3;

/* parameters for encapsulation */
const double lb[] = {0, 0, 0};   // lower bound of the actual function range for each dimension
const double ub[] = {10, 10, 10};     // upper bound
const double compress_coef = 0.95;  // the actual function will be compressed into a hypercube "bound": (compress_coef * [-L,L])^dim
const double slope = 100;             // specifies how fast the encapsulated function will grow out of the "bound"
const double pinned[] = {5.00, 5.00, 5.00}; // the pinned point will be guaranteed to be picked by the QHD discretization

double get_obj(const double *x) {
    // Shift the bounds to be around the origin
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double K = 5;
    double alpha = 5;
    
    double sqrtSum = sqrt(pow(x1 - alpha, 2) + pow(x2 - alpha, 2) + pow(x3 - alpha, 2));
    double term1 = pow(x1 - alpha, 2) - cos(K * sqrtSum);
    double term2 = pow(x2 - alpha, 2) - cos(K * sqrtSum);
    double term3 = pow(x3 - alpha, 2) - cos(K * sqrtSum);

    return 0.1*(term1 + term2 + term3);
}

void get_obj_subg(const double *x, double *ret) {
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double K = 5;
    double alpha = 5;
    
    double sqrtSum = sqrt(pow(x1 - alpha, 2) + pow(x2 - alpha, 2) + pow(x3 - alpha, 2));

    ret[0] = 0.1*(15*(x1-alpha)*sin(K*sqrtSum)/sqrtSum + 2*(x1 - alpha));
    ret[1] = 0.1*(15*(x2-alpha)*sin(K*sqrtSum)/sqrtSum + 2*(x2 - alpha));
    ret[2] = 0.1*(15*(x3-alpha)*sin(K*sqrtSum)/sqrtSum + 2*(x3 - alpha));


}