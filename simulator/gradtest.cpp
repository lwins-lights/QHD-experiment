#include <iostream>
#include <execution>
#include <algorithm>
#include <complex>
#include <math.h>
#include <omp.h>
#include "potential.hpp"
#include "config.hpp"
#include <cnpy.h>
#include <random>
#include <string.h>

using namespace std;

void add_delta(const double *x, const int i, const double d,
               double *x_new, const double L, const int size) {
    /*
        add d*\delta_i to vector x
    */
    
    for (int j = 0; j < size; j++) {
        x_new[j] = x[j];
    }
    x_new[i] += d;
}

void grad(const int dim, const double stepsize, const double L, 
          const double *x, double *ret) {
    /*
        save \nabla V to ret[]
        V is given by get_potential()
    */

    const double V = get_potential(x, L);
    double x_new[dim];
    double V_new;
    
    for (int i = 0; i < dim; i++) {
        add_delta(x, i, stepsize, x_new, L, dim);
        V_new = get_potential(x_new, L);
        ret[i] = (V_new - V) / stepsize;
    }
}

void gradtest(const double L, const double eps, const int num) {

    int dim;
    get_potential_params(dim);

    double x[dim], xd_disc[dim], xd_subg[dim];
    int err_cnt = 0;
    double temp, maxtemp, tot_lip, max_lip, lip;

    default_random_engine gen;
    uniform_real_distribution<double> uniform(0, 1);

    maxtemp = tot_lip = max_lip = 0;

    for (int k = 0; k < num; k++) {

        /* randomize x */
        for (int j = 0; j < dim; j++) {
            x[j] = (uniform(gen) - 0.5) * 2 * L;
        }

        /* get gradient by the oracle and the finite difference */
        get_subgradient(x, xd_subg, L);
        grad(dim, eps, L, x, xd_disc);

        /* lip = the norm of the gradient at x */
        lip = 0;

        /* comapre */
        for (int j = 0; j < dim; j++) {
            lip += xd_subg[j] * xd_subg[j];
            temp = abs(xd_subg[j] - xd_disc[j]);
            maxtemp = max(temp, maxtemp);
            if (temp == maxtemp) {
                //printf("DEBUG: discrepency: %f ===> %f VS %f\n", temp, xd_subg[j], xd_disc[j]);
                printf("At point ( ");
                for (int i = 0; i < dim; i++) {
                    printf("%.4f ", x[i]);
                }
                printf(") coord %d:\n oracle = %.8f <> calc = %.8f\n", j, xd_subg[j], xd_disc[j]);
                err_cnt += 1;
            }
        }

        /* resolve lip */
        lip = sqrt(lip);
        if (lip > max_lip) {
            max_lip = lip;
        }
        tot_lip += lip;
    }

    printf("Max discrepency = %.8e\n", maxtemp);
    printf("Average / Max Lipschitz constant = %.8f / %.8f\n", tot_lip / num, max_lip);
    //printf("Test result: %d/%d erroneous values found.\n", err_cnt, num * dim);
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        perror("Expected arguments: ./gradtest <L> <eps> <num>");
        exit(EXIT_FAILURE);
    }
    const double L = stod(argv[1]);
    const double eps = stod(argv[2]);
    const int num = stoi(argv[3]);

    gradtest(L, eps, num);

    return 0;
}