#include <iostream>
#include <algorithm>
#include <execution>
#include <complex>
#include <math.h>
#include <fftw3.h>
#include <omp.h>
#include "qipopt_F.hpp"
#include "config.hpp"
#include <cnpy.h>
#include <vector>

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

    pp_pair *pp = new pp_pair[size];
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

    delete[] pp;

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

void sample_from_psi(const comp *psi, const double rand, double *z, const double L, const int len, const int dim) {

    const int size = pow(len, dim);
    double p = 0;
    int x[dim];

    for (int i = 0; i < size; i++) {
        p += (psi[i] * conj(psi[i])).real();
        if (p > rand) {
            int_to_coord(i, x, dim, len);
            for (int j = 0; j < dim; j++)
                z[j] = x[j] * 2 * L / len - L;
            //printf("[DEBUG] i=%d rand=%.4lf", i, rand);
            return;
        }
    }
}

double kinetic_energy(comp *psi, double kop_coef, comp *kop,
                       int dim, int *n, int size) {

    return 0;

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
        //tot_energy += (u[i] * conj(u[i])).real() * kop_coef * kop[i].real();
        tot_energy += (u[i] * conj(u[i])).real() * kop[i].real();
    }

    return tot_energy / tot_prob * kop_coef;
}

double dg(const double tau) {
    /*
        the derivative of g function in the paper without normalization
    */

    if (tau <= 0 || tau >= 1) {
        return 0;
    } else {
        return exp(-1 / tau / (1 - tau));
    }
}

void initialize_mu(const int num, double *mu) {
    /*
        mu function in the paper where t \in [0,1]
        the returned vector mu[i] := mu(i/num)
    */

    mu[0] = 0;
    for (int i = 1; i < num; i++) {
        mu[i] = mu[i - 1] + dg(double(i) / num);
    }
    for (int i = 0; i < num; i++) {
        mu[i] /= mu[num - 1];
        mu[i] = 1 - (1 - mu_f) * mu[i];
    }
    //printf("%.8lf,%.8lf\n", mu[0], mu[num-1]);
}

double get_mu(double t, double *mu, int size) {
    return mu[int(size * t)];
}

/*
double get_h(double t, double *mu, int size, double h_coef) {
    return get_mu(t, mu, size) * h_coef
}
*/

double qcp_kop_coef(double t, double eta_qhd, double eta) {
    double tp; // true time rescaled
    double nominal = 0;

    if (t <= 1) {
        tp = t / eta_qhd;
        nominal = 1 / (1e-6 + tp) - 1 / (1e-6 + 1 / eta) + tp * h_coef * h_coef;
        return nominal / eta_qhd;
    } else {
        nominal = h_coef;
        return nominal / eta;
    }
}

double qcp_pot_coef(double t, double eta_qhd, double eta, double *mu, int size) {
    double tp; // true time rescaled
    double nominal = 0;
    double temp;

    if (t <= 1) {
        tp = t / eta_qhd;
        nominal = 1 * tp;
        return nominal / eta_qhd;
    } else {
        temp = get_mu(t - 1, mu, size);
        nominal = 1 / h_coef / temp / temp;
        return nominal / eta;
    }
}

void qcp_pot_loader(double *V, double *W, double t, const int len, const double L, const int dim, double eta, double *mu_data, int mu_size) {

    int size = pow(len, dim);
 
    int v[dim];
    double x[dim], F[dim];
    double mu;
    double tp = t / eta; // true time rescaled

    if (t <= 1) {
        mu = get_mu(0, mu_data, mu_size);
    } else {
        mu = get_mu(t - 1, mu_data, mu_size);
        //printf("%.4lf\n", mu);
    }
    #pragma omp parallel for private(v, x, F)
    for (int i = 0; i < size; i++) {
        int_to_coord(i, v, dim, len);
        for (int j = 0; j < dim; j++) {
            /* map [0, len) to [-L, L) */
            x[j] = (double) v[j] * 2 * L / len - L;
        }
        get_F(x, F);
        V[i] = W[i] = 0;
        for (int j = 0; j < dim; j++) {
            W[i] += F[j] / dim;
            F[j] -= mu;
            V[i] += 0.5 * F[j] * F[j];
        }
    }
}

double get_avg_comp_gap(comp *psi, double *W, int size) {

    double ret = 0;

    #pragma omp parallel for reduction(+: ret)
    for (int i = 0; i < size; i++) {
        ret += (psi[i] * conj(psi[i])).real() * W[i];
    }

    return ret;
}

void pseudospec(const int dim, const int len, const double L, const double T, 
                const double disc, comp *psi, const int par, double eta_qhd, double eta, 
                double *mu_data, int mu_size, const double rand, 
                double (*kop_coef)(double, double, double), 
                double (*pot_coef)(double, double, double, double*, int),
                void (*pot_loader)(double*, double*, double, int, double, int, double, double*, int)) {
    /*
        pseudospectral solver for time-dependent Schrodinger equation
            H(t) = t_dep_1(t)*(-1/2*\nabla^2) + t_dep_2(t)*V(x)
        defined on [-L,L]^dim discretized into len^dim hypercube nodes
        refer to Section C.2.1 in https://arxiv.org/pdf/2303.01471v1.pdf
    */

    const int size = pow(len, dim);

    int n[dim];
    int i, prog, prog_prev;
    comp *kop = new comp[size];
    comp *u = new comp[size];
    comp *psi_new = new comp[size];
    comp *temp = new comp[size];
    double time_st, time_ed, t, dt, temp_tot, thr, cur_kop_coef, cur_pot_coef;
    double *V = new double[size];
    double *W = new double[size];
    double *psi_prob = new double[size];
    double *prob = new double[101 * size];
    double z[dim];
    fftw_plan plan_ft, plan_ift;
    vector<double> time_arr, pot_arr, gap_arr;

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
    dt = dt_init;
    while (t < T) {

        /* potential term first */
        pot_loader(V, W, t, len, L, dim, eta, mu_data, mu_size);
        cur_pot_coef = pot_coef(t, eta_qhd, eta, mu_data, mu_size);

        /* compute the probability distribution */
        #pragma omp parallel for private(i)
        for (i = 0; i < size; i++) {
            psi_prob[i] = (psi[i] * conj(psi[i])).real();
        }

        /* try dt until |psi - psi_new|_1 < disc */
        dt *= dt_grow;
        while (true) {
            /* compute psi_new: first potential term */
            #pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                temp[i] = exp(-iu * V[i] * (dt * cur_pot_coef));
            }
            hadamard_product(temp, psi, u, size);

            /* then kinetic term:  */
            fftw_execute(plan_ft);
            cur_kop_coef = kop_coef(t, eta_qhd, eta);
            #pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                temp[i] = exp(-iu * kop[i] * (dt * cur_kop_coef));
            }
            hadamard_product(temp, u, u, size);
            fftw_execute(plan_ift);

            /* normalize psi_new */
            temp_tot = 0; // 2-norm of psi_new
            #pragma omp parallel for private(i) reduction(+: temp_tot)
            for (i = 0; i < size; i++) {
                temp_tot += (psi_new[i] * conj(psi_new[i])).real();
            }
            #pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                psi_new[i] /= sqrt(temp_tot);
            }

            /* compute discrepancy */
            temp_tot = 0; // discrepancy
            #pragma omp parallel for private(i) reduction(+: temp_tot)
            for (i = 0; i < size; i++) {
                temp_tot += abs((psi[i] * conj(psi[i])).real() - (psi_new[i] * conj(psi_new[i])).real());
            }

            /* shrink dt and try again if not small enough */
            if (temp_tot < disc) {
                break;
            } else {
                dt *= dt_shrink;
            }
        }

        /* update psi */
        #pragma omp parallel for private(i) 
        for (i = 0; i < size; i++) {
            psi[i] = psi_new[i];
        }

        /* print progress bar */
        prog = int(t / T * 100);
        if (prog != prog_prev) {
            print_prog_bar(prog);
            prog_prev = prog;
            /* save potential to prob */
            #pragma omp parallel for private(i)
            for (i = 0; i < size; i++) {
                prob[prog * size + i] = (psi[i] * conj(psi[i])).real();
            }
            /* debug */
            //printf("dt = %f\n", dt);
        } 

        /* store results */
        time_arr.push_back(t);
        pot_arr.push_back(expected_potential(psi, V, size, par));
        gap_arr.push_back(get_avg_comp_gap(psi, W, size));

        /* update current time */
        t += dt;
    }

    sample_from_psi(psi, rand, z, L, len, dim);

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("Pseudospectral integrator runtime = %.5f s\n", time_ed - time_st);

    /* save results */
    double *time_arr_classical = new double[time_arr.size()];
    #pragma omp parallel for private(i) 
    for (i = 0; i < time_arr.size(); i++) {
        time_arr_classical[i] = time_arr[i];
    }
    double *pot_arr_classical = new double[pot_arr.size()];
    #pragma omp parallel for private(i) 
    for (i = 0; i < pot_arr.size(); i++) {
        pot_arr_classical[i] = pot_arr[i];
    }
    double *gap_arr_classical = new double[gap_arr.size()];
    #pragma omp parallel for private(i) 
    for (i = 0; i < gap_arr.size(); i++) {
        gap_arr_classical[i] = gap_arr[i];
    }
    npz_save("../result/qipopt.npz", "timestamps", time_arr_classical, {time_arr.size()}, "w");
    npz_save("../result/qipopt.npz", "expected_potential", pot_arr_classical, {pot_arr.size()}, "a");
    npz_save("../result/qipopt.npz", "expected_duality_gap", gap_arr_classical, {gap_arr.size()}, "a");
    npz_save("../result/qipopt.npz", "probability", prob, {101, (unsigned int) size}, "a");
    npz_save("../result/qipopt.npz", "T", &T, {1}, "a");
    npz_save("../result/qipopt.npz", "L", &L, {1}, "a");
    npz_save("../result/qipopt.npz", "disc", &disc, {1}, "a");
    npz_save("../result/qipopt.npz", "par", &par, {1}, "a");
    npz_save("../result/qipopt.npz", "len", &len, {1}, "a");
    npz_save("../result/qipopt.npz", "dim", &dim, {1}, "a");
    npz_save("../result/qipopt.npz", "eta", &eta, {1}, "a");
    npz_save("../result/qipopt.npz", "sample", z, {(unsigned int) dim}, "a");
}

int main(int argc, char **argv)
{
    if (argc != 6) {
        perror("Expected arguments: ./qipopt <len> <eta_qhd> <eta> <disc> <rand>");
        exit(EXIT_FAILURE);
    }
    const int len = stoi(argv[1]);
    const double eta_qhd = stod(argv[2]);
    const double eta = stod(argv[3]);
    const double disc = stod(argv[4]);
    const double rand = stod(argv[5]);
    const int par = 1;
    const double T = 2;

    double L;
    int dim;
    get_F_params(L, dim);

    const int size = pow(len, dim);
    //double V[size];
    //load_potential_to_array(V, len, L, dim);

    printf("L=%f, len=%d\n", L, len);
    printf("T=%f, disc=%f\n", T, disc);
    printf("eta_qhd=%f, eta=%f, mu_f=%f\n", eta_qhd, eta, mu_f);

    const double stepsize = 2 * L / len;
    printf("stepsize=%f\n", stepsize);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    comp *psi = new comp[size];
    double *mu = new double[mu_size];

    initialize_mu(mu_size, mu); //printf("mu=%.8lf\n",get_mu(0.7, mu, mu_size));
    initialize_psi(psi, size);
    pseudospec(dim, len, L, T, disc, psi, par, eta_qhd, eta, mu, mu_size, rand, qcp_kop_coef, qcp_pot_coef, qcp_pot_loader);
    
    return 0;
}