#include <iostream>
#include <algorithm>
#include <execution>
#include <complex>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "potential.hpp"
#include "config.hpp"
#include <cnpy.h>

using namespace std;
using namespace cnpy;

/* Type "comp" and "fftw_complex" will be compatible */
typedef complex<double> comp;

/* imaginary unit */
const comp iu = comp(0, 1);

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

int coord_to_int(const int *x, const int dim, const int len) {
    /*
        translate (x[0], x[1], ..., x[dim-1]) to x[0] + len*x[1] + len^2*x[2]+...
    */

    int ret, temp;

    ret = 0;
    temp = 1;
    for (int i = 0; i < dim; i++) {
        ret += temp * x[i];
        temp *= len;
    }
    return ret;
}

void int_to_coord(const int l, int *x, const int dim, const int len) {
    /*
        inverse of coord_to_int()
    */

    int temp;

    temp = l;
    for (int i = 0; i < dim; i++) {
        x[i] = temp % len;
        temp /= len;
    }
}

double sqr_norm(const int *x, const int size) {
    /*
        compute x[0]^2 + ... + x[size-1]^2
    */

    double ret;

    ret = 0;
    for (int i = 0; i < size; i++) {
        ret += (double) x[i] * x[i];
    }
    return ret;
}

double t_dep_1(const double &t) {
    return 2 / (1e-3 + t * t * t);
}

double t_dep_2(const double &t) {
    return 2 * t * t * t;
}

void hadamard_product(comp *A, comp *B, comp *C, int size) {
    /*
        entrywise product of two complex vectors
    */

    int i;

    #pragma omp parallel for private(i)
    for (i = 0; i < size; i++) {
        C[i] = A[i] * B[i];
    }
}

double frac(double x) {
    /*
        return the fractional part of x, i.e., {x}, for any x >= 0
    */
    return x - int(x);
}

void load_potential_to_array(double *V, const int len, const double L, const int dim) {
    
    int i;

    const int size = pow(len, dim);
    const double stepsize = 2 * L / len;

    int v[dim];
    double x[dim], pinned[dim], offset[dim];

    get_pinned_point(pinned, L);
    for (int i = 0; i < dim; i++) {
        offset[i] = frac((pinned[i] + L) / stepsize) * stepsize;
    }

    #pragma omp parallel for private(i, v, x)
    for (i = 0; i < size; i++) {
        int_to_coord(i, v, dim, len);
        for (int j = 0; j < dim; j++) {
            /* map [0, len) to [-L, L) */
            //x[j] = (double) v[j] * stepsize - L + offset[j];
            /* disable pinned */
            x[j] = (double) v[j] * stepsize - L;
        }
        V[i] = get_potential(x, L);
    }
}

void initialize_psi(comp *psi, const int size) {
    /*
        make psi[i] = 1/sqrt(len)
    */

    int i;
    double temp;

    temp = ((double) 1) / sqrt((double) size);

    #pragma omp parallel for private(i)
    for (i = 0; i < size; i++) {
        psi[i] = comp(temp, 0);
    }
}

void initialize_kinetic_operator(comp *op, const int dim, const int len, const double L) {
    /*
        Fourier-tranformed kinetic operator (-1/2*\nabla^2)
        discretized on the hypercube of size len^dim
    */

    int freq[len];
    int x[dim];
    int size, i;

    size = pow(len, dim);

    #pragma omp parallel for private(i)
    for (i = 0; i < len / 2; i++) {
        freq[i] = i;
        freq[len / 2 + i] = -len / 2 + i;
    }

    #pragma omp parallel for private(i, x)
    for (i = 0; i < size; i++) {
        int_to_coord(i, x, dim, len);
        for (int j = 0; j < dim; j++) {
            x[j] = freq[x[j]];
        }
        op[i] = (0.5 * M_PI * M_PI * sqr_norm(x, dim)) / (L * L);
    }
}

double expected_potential(const comp *psi, const double *V, const int size,
                          const int par) {
    /*
        compute the expected potential with wave function psi
        psi and V are discretized on {0,1,...,size-1}
    */

    struct pp_pair {
        double prob, pot;
    };

    pp_pair pp[size];
    double ret, tot_prob, tot_prob_new;

    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        pp[i].prob = (psi[i] * conj(psi[i])).real();
        pp[i].pot = V[i];
    }

    /* parallelized sort; unnecessary if par = 1 */
    if (par > 1) {
        sort(execution::par_unseq, pp, pp + size, 
             [](pp_pair a, pp_pair b) {return a.pot > b.pot;});
    }

    /* accumulate the result; this CANNOT be directly parallelized if par > 1 */
    ret = tot_prob = 0;
    if (par > 1) {
        for (int i = 0; i < size; i++) {
            tot_prob_new = tot_prob + pp[i].prob;
            ret += pp[i].pot * (pow(tot_prob_new, par) - pow(tot_prob, par));
            tot_prob = tot_prob_new;
        }
    } else {
        for (int i = 0; i < size; i++) {
            ret += pp[i].pot * pp[i].prob;
        }
    }

    return ret;
}

double prob_at_minimum(const comp *psi, const double *V, const int size,
                       const double thr, const int par) {
    /*
        compute the probability that the potential is less than thr
        proability distribution defined by the wave function psi
        psi and V are discretized on {0,1,...,size-1}
    */

    double prob;
    
    /* compute Pr[result < thr] */
    prob = 0;

    #pragma omp parallel for reduction(+: prob)
    for (int i = 0; i < size; i++) {
        if (V[i] < thr) {
            prob += (psi[i] * conj(psi[i])).real();
        }
    }

    /* convert to Pr[all par results < thr] */
    return 1 - pow(1 - prob, par);
}

void compile_distribution(const comp *psi, const int size, double *dist) {

    /* return the distribution when measuring psi */

    for (int i = 0; i < size; i++) {
        dist[i] = (psi[i] * conj(psi[i])).real();
    }
}

double get_minimum(const double *V, const int size) {
    double ret = V[0];
    for (int i = 0; i < size; i++) {
        if (V[i] < ret) {
            ret = V[i];
        }
    }
    return ret;
}

void pseudospec(const int dim, const int len, const double L, const double T, 
                const double dt, const double *V, comp *psi, const int par) {
    /*
        pseudospectral solver for time-dependent Schrodinger equation
            H(t) = t_dep_1(t)*(-1/2*\nabla^2) + t_dep_2(t)*V(x)
        defined on [-L,L]^dim discretized into len^dim hypercube nodes
        refer to Section C.2.1 in https://arxiv.org/pdf/2303.01471v1.pdf
    */

    const int size = pow(len, dim);
    const int num_steps = T / dt;

    int n[dim];
    int i, prog, prog_prev, index;
    comp kop[size], u[size], psi_new[size], temp[size];
    double time_st, time_ed, t, temp_tot, thr;
    double pot[num_steps], prob_at_min[num_steps];
    double dist[size];
    double dist_snapshot[size * n_snapshot];
    fftw_plan plan_ft, plan_ift;

    /* n for fftw later */
    for (i = 0; i < dim; i++) {
        n[i] = len;
    }

    /* load Fourier-tranformed kinetic operator */
    initialize_kinetic_operator(kop, dim, len, L);

    /* fftw setup: plan_ft for DFT and plan_ift for inverse DFT */
    
    if (fftw_init_threads() == 0)
    {
        cout << "Error with fftw_init_threads." << endl;
        exit(EXIT_FAILURE);
    }
    fftw_plan_with_nthreads(omp_get_max_threads());
    plan_ft = fftw_plan_dft(dim, n, (fftw_complex *) u, (fftw_complex *) u, FFTW_FORWARD, FFTW_MEASURE);
    plan_ift = fftw_plan_dft(dim, n, (fftw_complex *) u, (fftw_complex *) psi_new, FFTW_BACKWARD, FFTW_MEASURE);

    /* timing */
    time_st = omp_get_wtime();

    /* 
        main loop: use Trotterization to alternatively apply the potential 
        term and the kinetic term
    */
    t = 0;
    prog_prev = -1;
    for (int step = 0; step < num_steps; step++) {

        /* potential term first */

        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            temp[i] = exp(-iu * dt * t_dep_2(t) * V[i]);
        }

        hadamard_product(temp, psi, u, size);

        /* then kinetic term:  */
        fftw_execute(plan_ft);

        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            temp[i] = exp(-iu * dt * t_dep_1(t) * kop[i]);
        }

        hadamard_product(temp, u, u, size);
        fftw_execute(plan_ift);

        /* update psi and normalize */
        temp_tot = 0;

        #pragma omp parallel for private(i) reduction(+: temp_tot)
        for (i = 0; i < size; i++) {
            temp_tot += (psi_new[i] * conj(psi_new[i])).real();
        }

        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            psi[i] = psi_new[i] / sqrt(temp_tot);
        }

        /* update current time */
        t = (step + 1) * dt;

        /* print progress bar */
        prog = (step + 1) * 100 / num_steps;
        if (prog != prog_prev) {
            print_prog_bar(prog);
            prog_prev = prog;
        } 

        pot[step] = expected_potential(psi, V, size, par);
        if (step == 0) {
            thr = thr_frac * expected_potential(psi, V, size, 1);
        }
        prob_at_min[step] = prob_at_minimum(psi, V, size, thr, par);

        /* write into snapshot */
        if ((step * n_snapshot) % num_steps < n_snapshot) {
            index = (step * n_snapshot) / num_steps;
            //printf("[DEBUG] Snapshot %d\n", index);
            compile_distribution(psi, size, dist_snapshot + (index * size));
        }
    }

    compile_distribution(psi, size, dist);

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("Pseudospectral integrator runtime = %.5f s\n", time_ed - time_st);

    /* save results */
    npz_save("../result/pseudospec.npz", "expected_potential", pot, {(unsigned int) num_steps}, "w");
    npz_save("../result/pseudospec.npz", "probability_at_minimum", prob_at_min, {(unsigned int) num_steps}, "a");
    npz_save("../result/pseudospec.npz", "T", &T, {1}, "a");
    npz_save("../result/pseudospec.npz", "dt", &dt, {1}, "a");
    npz_save("../result/pseudospec.npz", "par", &par, {1}, "a");
    npz_save("../result/pseudospec.npz", "len", &len, {1}, "a");
    npz_save("../result/pseudospec.npz", "dim", &dim, {1}, "a");
    npz_save("../result/pseudospec.npz", "L", &L, {1}, "a");
    npz_save("../result/pseudospec.npz", "V", V, {(unsigned int) size}, "a");
    npz_save("../result/pseudospec.npz", "dist", dist, {(unsigned int) size}, "a");
    npz_save("../result/pseudospec.npz", "dist_snapshot", dist_snapshot, {(unsigned int)(size * n_snapshot)}, "a");
}

int main(int argc, char **argv)
{
    if (argc != 6) {
        perror("Expected arguments: ./pseudospec <len> <T> <dt> <par> <L>");
        exit(EXIT_FAILURE);
    }
    const int len = stoi(argv[1]);
    const double T = stod(argv[2]);
    const double dt = stod(argv[3]);
    const int par = stoi(argv[4]);
    const double L = stod(argv[5]);

    int dim;
    get_potential_params(dim);

    const int size = pow(len, dim);
    double V[size];
    load_potential_to_array(V, len, L, dim);
    printf("Post-Discretization Minimum: %f\n", get_minimum(V, size));

    printf("L=%f, len=%d\n", L, len);
    printf("T=%f, dt=%f\n", T, dt);

    const double stepsize = 2 * L / len;
    printf("stepsize=%f\n", stepsize);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    comp psi[size];

    initialize_psi(psi, size);
    pseudospec(dim, len, L, T, dt, V, psi, par);
    
    return 0;
}