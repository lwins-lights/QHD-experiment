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

using namespace std;
using namespace cnpy;

const int bar_width = 70;
void print_prog_bar(int prog) {
    const int pos = prog * bar_width / 100;
    string bar;

    bar = "[";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) bar += "=";
        else if (i == pos) bar += ">";
        else bar += " ";
    }
    cout << bar << "] " << prog << " %\r";
    cout.flush();
}

void vec_minus(const double *a, const double *b, double *c, const int size) {
    /*
        c := a - b
    */

    for (int i = 0; i < size; i++) {
        c[i] = a[i] - b[i];
    }
}

void vec_plus(const double *a, const double *b, double *c, const int size) {
    /*
        c := a + b
    */

    for (int i = 0; i < size; i++) {
        c[i] = a[i] - b[i];
    }
}

void vec_copy(const double *a, double *b, const int size, const double L) {
    /*
        b := a
        coordinates circular on [-L, L)
    */

    double temp;

    for (int i = 0; i < size; i++) {
        b[i] = a[i];

        /* map [-L, L) to [0, 1) */
        temp = (b[i] + L) / 2 / L;

        temp = temp - floor(temp);

        /* map back */
        b[i] = temp * 2 * L - L;
    }
}

void scalar_mul(double *a, const double k,  const int size) {
    /*
        a := ka
    */

    for (int i = 0; i < size; i++) {
        a[i] *= k;
    }
}

void int_to_coord(const int l, int *x, const int dim, const int len) {
    /*
        inverse of coord_to_int()
        coord_to_int(): 
            translate (x[0], x[1], ..., x[dim-1]) to x[0] + len*x[1] + len^2*x[2]+...
    */

    int temp;

    temp = l;
    for (int i = 0; i < dim; i++) {
        x[i] = temp % len;
        temp /= len;
    }
}

void add_delta(const double *x, const int i, const double d,
               double *x_new, const double L, const int size) {
    /*
        add d*\delta_i to vector x
        coordinates circular on [-L, L)
    */
    
    for (int j = 0; j < size; j++) {
        x_new[j] = x[j];
    }
    x_new[i] += d;
    if (x_new[i] >= L) {
        x_new[i] -= 2 * L;
    }
}

void grad(const int dim, const double stepsize, const double L, 
          const double *x, double *ret) {
    /*
        save \nabla V to ret[]
        V is given by get_potential()
    */

    const double V = get_potential(x);
    double x_new[dim];
    double V_new;
    
    for (int i = 0; i < dim; i++) {
        add_delta(x, i, stepsize, x_new, L, dim);
        V_new = get_potential(x_new);
        ret[i] = (V_new - V) / stepsize;
    }
}

double compute_potential(const double *pot, const int num, const int par) {
    double sorted_pot[num];
    double ret;

    #pragma omp parallel for
    for (int i = 0; i < num; i++) {
        sorted_pot[i] = pot[i];
    }

    /* parallelized sort; unnecessary if par = 1 */
    if (par > 1) {
        sort(execution::par_unseq, sorted_pot, sorted_pot + num);
    }
    ret = 0;

    #pragma omp parallel for reduction(+: ret)
    for (int i = 0; i < num; i++) {
        ret += sorted_pot[i] * (pow((double) (num - i) / num, par) - pow((double) (num - i - 1) / num, par));
    }
    return ret;
}

double prob_at_minimum(const double *pot, const int num, const double thr,
                       const int par) {
    double prob;
    //double mi;

    /* compute Pr[result < thr] */
    prob = 0;
    //mi = 100;

    #pragma omp parallel for reduction(+: prob)
    for (int i = 0; i < num; i++) {
        if (pot[i] < thr) {
            prob += (double) 1 / num;
        }
        //if (pot[i] < mi) mi = pot[i];
    }
    //cout << "probcnt=" << prob << "mi=" << mi << endl;

    /* convert to Pr[all par results < thr] */
    return 1 - pow(1 - prob, par);
}

// L, dim, tot_steps, learning_rate, par, sample_number

void subgrad(const double L, const int dim, const int tot_steps,
             const double learning_rate, const int par, const int sample_number) {
    /*
        Subgradient Method
        Update x by x - lr * f'(x), where f' is the subgradient and lr = learning_rate / sqrt(tot_steps)
    */

    const int num = sample_number;
    const double lr = learning_rate / sqrt(tot_steps);
    int prog, prog_prev;
    int v[dim];
    double x[num][dim], x_new[num][dim], temp[dim];
    double expected_pot, time_st, time_ed, thr;
    double pot[tot_steps], prob_at_min[tot_steps];
    double cur_pot[num];

    /* randomness preparation */
    default_random_engine gen;
    //normal_distribution<double> gaussian(0, noise_level * sqrt(dt) / sqrt(dim));
    uniform_real_distribution<double> uniform(0, 1);
    
    /* 
     * run subgrad with random samples in the hypercube
     */
    for (int id = 0; id < num; id++) {
        for (int j = 0; j < dim; j++) {
            x[id][j] = (uniform(gen) - 0.5) * 2 * L;
        }
    }

    /* timing */
    time_st = omp_get_wtime();

    /* main loop */
    for (int step = 0; step < tot_steps; step++) {
        expected_pot = 0;

        #pragma omp parallel for private(temp) reduction(+: expected_pot)
        for (int id = 0; id < num; id++) {
            /* evaluate potential */
            cur_pot[id] = get_potential(x[id]);

            /* load subgradient to temp */
            get_subgradient(x[id], temp);

            /* x' = x - lr * temp */
            scalar_mul(temp, lr, dim);
            vec_minus(x[id], temp, x_new[id], dim);

            /* update */
            vec_copy(x_new[id], x[id], dim, L);
        }

        /* print progress bar */
        prog = (step + 1) * 100 / tot_steps;
        if (prog != prog_prev) {
            print_prog_bar(prog);
            prog_prev = prog;
        } 

        pot[step] = compute_potential(cur_pot, num, par);
        if (step == 0) {
            thr = thr_frac * compute_potential(cur_pot, num, 1);
        }
        prob_at_min[step] = prob_at_minimum(cur_pot, num, thr, par);
    }

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("subgrad simulator runtime = %.5f s\n", time_ed - time_st);

    /* save results */
    npz_save("../result/subgrad.npz", "expected_potential", pot, {(unsigned int) tot_steps}, "w");
    npz_save("../result/subgrad.npz", "probability_at_minimum", prob_at_min, {(unsigned int) tot_steps}, "a");
    npz_save("../result/subgrad.npz", "L", &L, {1}, "a");
    npz_save("../result/subgrad.npz", "dim", &dim, {1}, "a");
    npz_save("../result/subgrad.npz", "tot_steps", &tot_steps, {1}, "a");
    npz_save("../result/subgrad.npz", "learning_rate", &learning_rate, {1}, "a");
    npz_save("../result/subgrad.npz", "par", &par, {1}, "a");
    npz_save("../result/subgrad.npz", "sample_number", &sample_number, {1}, "a");
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        perror("Expected arguments: ./subgrad <tot_steps> <learning_rate> <par> <sample_number>");
        exit(EXIT_FAILURE);
    }
    const int tot_steps = stoi(argv[1]);
    const double learning_rate = stod(argv[2]);
    const int par = stoi(argv[3]);
    const int sample_number = stoi(argv[4]);

    double L;
    int dim;
    get_potential_params(L, dim);

    printf("L=%f, dim=%d\n", L, dim);
    printf("tot_steps=%d, learning_rate=%f, par=%d, sample_number=%d\n", tot_steps, learning_rate, par, sample_number);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    subgrad(L, dim, tot_steps, learning_rate, par, sample_number);
}