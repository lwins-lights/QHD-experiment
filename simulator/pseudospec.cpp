#include <iostream>
#include <complex>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "potential.hpp"

using namespace std;

/* Type "comp" and "fftw_complex" will be compatible */
typedef complex<double> comp;

/* imaginary unit */
const comp iu = comp(0, 1);

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

void load_potential_to_array(double *V, const int len, const double L, const int dim) {
    
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
        V[i] = get_potential(x);
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

double expected_potential(const comp *psi, const double *V, int size) {
    /*
        compute the expected potential with wave function psi
        psi need not be normalized
        psi and V are discretized on {0,1,...,size-1}
    */

    int i;
    double pot, prob, temp;
    
    pot = prob = 0;

    #pragma omp parallel for private(i, temp) reduction(+: prob, pot)
    for (i = 0; i < size; i++) {
        temp = (psi[i] * conj(psi[i])).real();
        prob += temp;
        pot += temp * V[i];
    }

    return pot / prob;
}

void pseudospec(const int dim, const int len, const double L, const double T, 
                const double dt, const double *V, comp *psi) {
    /*
        pseudospectral solver for time-dependent Schrodinger equation
            H(t) = t_dep_1(t)*(-1/2*\nabla^2) + t_dep_2(t)*V(x)
        defined on [-L,L]^dim discretized into len^dim hypercube nodes
        refer to Section C.2.1 in https://arxiv.org/pdf/2303.01471v1.pdf
    */

    const int size = pow(len, dim);
    const int num_steps = T / dt;

    int n[dim];
    int i;
    comp kop[size], u[size], psi_new[size], temp[size];
    double time_st, time_ed, t;
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
    for (int step = 1; step <= num_steps; step++) {

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

        /* update psi; we do not strictly normalize psi for now */
        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            psi[i] = psi_new[i] / (double) size;
        }

        /* update current time */
        t = step * dt;

        /* [DEBUG] */
        double expected_value2 = expected_potential(psi, V, size);
        printf("[%d / %d] expected = %.8lf\n", step, num_steps, expected_value2);
    }

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("Pseudospectral integrator runtime = %.5f s\n", time_ed - time_st);
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        perror("Expected arguments: ./pseudospec <len> <T> <dt>");
        exit(EXIT_FAILURE);
    }
    const int len = stoi(argv[1]);
    const double T = stod(argv[2]);
    const double dt = stod(argv[3]);

    double L;
    int dim;
    get_potential_params(L, dim);

    const int size = pow(len, dim);
    double V[size];
    load_potential_to_array(V, len, L, dim);

    printf("L=%f, len=%d\n", L, len);
    printf("T=%f, dt=%f\n", T, dt);

    const double stepsize = 2 * L / len;
    printf("stepsize=%f\n", stepsize);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    comp psi[size];

    initialize_psi(psi, size);
    pseudospec(dim, len, L, T, dt, V, psi);
    
    return 0;
}