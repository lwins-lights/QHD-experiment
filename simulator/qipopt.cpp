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

double t_dep_1(const double t) {
    //return 2 / (1e-3 + t * t * t);
    return 1;
}

double t_dep_2(const double t) {
    //return 2 * t * t * t;
    return 1;
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

void load_potential_to_array(double *V, double t, const int len, const double L, const int dim) {
    
    int i;

    const int size = pow(len, dim);

    int v[dim];
    double x[dim];

    #pragma omp parallel for private(i, v, x)
    for (i = 0; i < size; i++) {
        int_to_coord(i, v, dim, len);
        for (int j = 0; j < dim; j++) {
            /* map [0, len) to [-L, L) */
            x[j] = (double) v[j] * 2 * L / len - L;
        }
        V[i] = get_potential(x) * t_dep_2(t);
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
        #pragma omp parallel for reduction(+: ret)
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

double kinetic_energy(comp *psi, double kop_coef, comp *kop,
                       int dim, int *n, int size) {

    fftw_plan plan_ft, plan_ift;
    comp v[size], u[size];
    double tot_prob, tot_energy;

    fftw_plan_with_nthreads(omp_get_max_threads());
    plan_ft = fftw_plan_dft(dim, n, (fftw_complex *) psi, (fftw_complex *) u, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan_ft);

    tot_prob = tot_energy = 0;

    #pragma omp parallel for reduction(+: tot_prob, tot_energy)
    for (int i = 0; i < size; i++) {
        tot_prob += (u[i] * conj(u[i])).real();
        tot_energy += (u[i] * conj(u[i])).real() * kop_coef * kop[i].real();
    }

    return tot_energy / tot_prob;
}

void pseudospec(const int dim, const int len, const double L, const double T, 
                const double dt, comp *psi, const int par,
                double (*kop_coef)(double), 
                void (*pot_loader)(double*, double, int, double, int)) {
    /*
        pseudospectral solver for time-dependent Schrodinger equation
            H(t) = t_dep_1(t)*(-1/2*\nabla^2) + t_dep_2(t)*V(x)
        defined on [-L,L]^dim discretized into len^dim hypercube nodes
        refer to Section C.2.1 in https://arxiv.org/pdf/2303.01471v1.pdf
    */

    const int size = pow(len, dim);
    const int num_steps = T / dt;

    int n[dim];
    int i, prog, prog_prev;
    comp kop[size], u[size], psi_new[size], temp[size];
    double time_st, time_ed, t, temp_tot, thr;
    double pot[num_steps], kin[num_steps];
    double V[size];
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

        pot_loader(V, t, len, L, dim);

        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            temp[i] = exp(-iu * dt * V[i]);
        }

        hadamard_product(temp, psi, u, size);

        /* then kinetic term:  */
        fftw_execute(plan_ft);

        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            temp[i] = exp(-iu * dt * kop_coef(t) * kop[i]);
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
        //prob_at_min[step] = prob_at_minimum(psi, V, size, thr, par);
        kin[step] = kinetic_energy(psi, kop_coef(t), kop, dim, n, size);

        /* update current time */
        t = (step + 1) * dt;
    }

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("Pseudospectral integrator runtime = %.5f s\n", time_ed - time_st);

    /* save results */
    npz_save("../result/qipopt.npz", "expected_potential", pot, {(unsigned int) num_steps}, "w");
    npz_save("../result/qipopt.npz", "expected_kinetic", kin, {(unsigned int) num_steps}, "a");
    npz_save("../result/qipopt.npz", "T", &T, {1}, "a");
    npz_save("../result/qipopt.npz", "dt", &dt, {1}, "a");
    npz_save("../result/qipopt.npz", "par", &par, {1}, "a");
    npz_save("../result/qipopt.npz", "len", &len, {1}, "a");
    npz_save("../result/qipopt.npz", "dim", &dim, {1}, "a");
}

int main(int argc, char **argv)
{
    if (argc != 5) {
        perror("Expected arguments: ./qipopt <len> <T> <dt> <par>");
        exit(EXIT_FAILURE);
    }
    const int len = stoi(argv[1]);
    const double T = stod(argv[2]);
    const double dt = stod(argv[3]);
    const int par = stoi(argv[4]);

    double L;
    int dim;
    get_potential_params(L, dim);

    const int size = pow(len, dim);
    //double V[size];
    //load_potential_to_array(V, len, L, dim);

    printf("L=%f, len=%d\n", L, len);
    printf("T=%f, dt=%f\n", T, dt);

    const double stepsize = 2 * L / len;
    printf("stepsize=%f\n", stepsize);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    comp psi[size];

    initialize_psi(psi, size);
    pseudospec(dim, len, L, T, dt, psi, par, t_dep_1, load_potential_to_array);
    
    return 0;
}