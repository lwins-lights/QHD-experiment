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
        c[i] = a[i] + b[i];
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

double vec_norm(const double *a, const int size) {
    /*
        return ||a||
    */
    double sqrnorm = 0;
    for (int i = 0; i < size; i++) {
        sqrnorm += a[i] * a[i];
    }
    return sqrt(sqrnorm);
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

void lfmsgd(const double L, const int dim, const int tot_steps,
             const double noise_level, const int par, const int sample_number) {
    /*
        LFMSGD with Gaussian Noise: https://arxiv.org/pdf/2406.18287
    */

    const int num = sample_number;
    int prog, prog_prev;
    int v[dim];
    double x_init[num][dim], x[num][dim], x_new[num][dim], temp[dim], m[num][dim], m_new[num][dim], m2_sum[num], mu[num];
    double expected_pot, time_st, time_ed, thr, eta, temp_val;
    double pot[tot_steps], prob_at_min[tot_steps];
    double cur_pot[num];

    /* randomness preparation */
    default_random_engine gen;
    normal_distribution<double> gaussian(0, noise_level);   // (mean, stddev)
    uniform_real_distribution<double> uniform(0, 1);

    gen.seed(0);
    
    /* 
     * run lfmsgd with random samples in the hypercube
     */
    for (int id = 0; id < num; id++) {
        for (int j = 0; j < dim; j++) {
            x[id][j] = x_init[id][j] = (uniform(gen) - 0.5) * 2 * L;
            m[id][j] = 0;
        }
        m2_sum[id] = mu[id] = 0;
    }

    /* timing */
    time_st = omp_get_wtime();

    /* main loop */
    for (int step = 0; step < tot_steps; step++) {
        expected_pot = 0;

        #pragma omp parallel for private(temp, temp_val, eta) reduction(+: expected_pot)
        for (int id = 0; id < num; id++) {
            /* evaluate potential */
            cur_pot[id] = get_potential(x[id], L);

            /* update coefficient */
            m2_sum[id] += vec_norm(m[id], dim) * vec_norm(m[id], dim);
            vec_minus(x[id], x_init[id], temp, dim);
            temp_val = vec_norm(temp, dim);
            if (temp_val > mu[id]) {
                mu[id] = min(temp_val, 2 * L * sqrt(dim));   // cap mu by the largest distance possible in the function domain
            } 
            if (step == 0) {
                eta = lfmsgd_eps;
            } else {
                eta = mu[id] / sqrt(lfmsgd_eps + m2_sum[id]);
            }

            /* load subgradient to temp */
            get_subgradient(x[id], temp, L);

            /* add noise */
            for (int i = 0; i < dim; i++) {
                temp[i] += gaussian(gen);
            }

            /* update momentum: m_ = beta * m + (1 - beta) * f'(x) */
            memcpy(m_new[id], m[id], sizeof(m[id]));
            scalar_mul(m_new[id], lfmsgd_beta, dim);
            scalar_mul(temp, 1 - lfmsgd_beta, dim);
            vec_plus(m_new[id], temp, m_new[id], dim);

            /* update x_ */
            memcpy(temp, m_new[id], sizeof(temp));
            scalar_mul(temp, eta, dim);
            vec_minus(x[id], temp, x_new[id], dim);

            /* let x := x_ for the next rotation */
            //vec_copy(x_new[id], x[id], dim, L);
            memcpy(x[id], x_new[id], sizeof(x[id]));
            memcpy(m[id], m_new[id], sizeof(m[id]));

            /* DEBUG */
            if (id == 0 and step < 0) {
                printf("[step %d]\neta=%f; m2_sum=%f, mu=%f, m=(%f,%f), x=(%f,%f)\n", step, eta, m2_sum[id], mu[id], m[id][0], m[id][1], x[id][0], x[id][1]);
            }
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
    printf("lfmsgd simulator runtime = %.5f s\n", time_ed - time_st);

    /* save results */
    npz_save("../result/lfmsgd.npz", "expected_potential", pot, {(unsigned int) tot_steps}, "w");
    npz_save("../result/lfmsgd.npz", "probability_at_minimum", prob_at_min, {(unsigned int) tot_steps}, "a");
    npz_save("../result/lfmsgd.npz", "L", &L, {1}, "a");
    npz_save("../result/lfmsgd.npz", "dim", &dim, {1}, "a");
    npz_save("../result/lfmsgd.npz", "tot_steps", &tot_steps, {1}, "a");
    npz_save("../result/lfmsgd.npz", "noise_level", &noise_level, {1}, "a");
    npz_save("../result/lfmsgd.npz", "par", &par, {1}, "a");
    npz_save("../result/lfmsgd.npz", "sample_number", &sample_number, {1}, "a");
    npz_save("../result/lfmsgd.npz", "samples", cur_pot, {(unsigned int) sample_number}, "a");
}

int main(int argc, char **argv)
{
    if (argc != 6) {
        perror("Expected arguments: ./lfmsgd <tot_steps> <noise_level> <par> <sample_number> <L>");
        exit(EXIT_FAILURE);
    }
    const int tot_steps = stoi(argv[1]);
    const double noise_level = stod(argv[2]);
    const int par = stoi(argv[3]);
    const int sample_number = stoi(argv[4]);
    const double L = stod(argv[5]);

    int dim;
    get_potential_params(dim);

    printf("L=%f, dim=%d\n", L, dim);
    printf("tot_steps=%d, noise_level=%f, par=%d, sample_number=%d\n", tot_steps, noise_level, par, sample_number);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    lfmsgd(L, dim, tot_steps, noise_level, par, sample_number);
}