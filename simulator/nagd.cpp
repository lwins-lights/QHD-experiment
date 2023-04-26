#include <iostream>
#include <complex>
#include <math.h>
#include <omp.h>
#include "potential.hpp"

using namespace std;

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

void vec_copy(const double *a, double *b, const int size) {
    /*
        b := a
    */

    for (int i = 0; i < size; i++) {
        b[i] = a[i];
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

void nagd(const int dim, const int len, const double L, const double T, 
          const double dt, const double stepsize) {
    /*
        Nesterov’s accelerated gradient descent
        refer to Section C.3.1 in https://arxiv.org/pdf/2303.01471v1.pdf
    */

    const int num = pow(len, dim);
    const int num_steps = T / dt;
    int v[dim];
    double x[num][dim], y[num][dim], x_new[num][dim], y_new[num][dim], temp[dim];
    double expected_pot, time_st, time_ed;
    
    /* run NAGD with len^dim different initial values in the hypercube */
    for (int id = 0; id < num; id++) {
        /* initialize x and y */
        int_to_coord(id, v, dim, len);
        for (int j = 0; j < dim; j++) {
            /* map [0, len) to [-L, L) */
            x[id][j] = y[id][j] = (double) v[j] * 2 * L / len - L;
        }
        //cout << x[id][0] << " " << x[id][1] << endl;
    }

    /* timing */
    time_st = omp_get_wtime();

    /* main loop */
    for (int step = 1; step <= num_steps; step++) {
        expected_pot = 0;

        for (int id = 0; id < num; id++) {
            /* x' = y - dt \nabla V(y) */
            grad(dim, stepsize, L, y[id], temp);
            //if (id == 0) cout << "temp:" << temp[0] << " " << temp[1] << endl;
            scalar_mul(temp, dt, dim);
            //if (id == 0) cout << "temp2:" << temp[0] << " " << temp[1] << endl;
            vec_minus(y[id], temp, x_new[id], dim);
            //if (id == 0) cout << "y:" << y[id][0] << " " << y[id][1] << endl;
            //if (id == 0) cout << "x_new:" << x_new[id][0] << " " << x_new[id][1] << endl;

            /* y' = x' + (k-1)/(k+2)*(x' - x) */
            vec_minus(x_new[id], x[id], temp, dim);
            //if (id == 0) cout << "temp3:" << temp[0] << " " << temp[1] << endl;
            scalar_mul(temp, (step - 1) / (step + 2), dim);
            //if (id == 0) cout << "temp4:" << temp[0] << " " << temp[1] << endl;
            vec_plus(x_new[id], temp, y_new[id], dim);
            //if (id == 0) cout << "y_new:" << y_new[id][0] << " " << y_new[id][1] << endl;

            /* update */
            vec_copy(x_new[id], x[id], dim);
            vec_copy(y_new[id], y[id], dim);

            expected_pot += get_potential(y[id]) / num;
        }
        /* [DEBUG] */
        printf("[%d / %d] expected = %.8lf\n", step, num_steps, expected_pot);
    }

    /* timing */
    time_ed = omp_get_wtime();
    printf("\n");
    printf("Pseudospectral integrator runtime = %.5f s\n", time_ed - time_st);
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        perror("Expected arguments: ./nagd <len> <T> <dt>");
        exit(EXIT_FAILURE);
    }
    const int len = stoi(argv[1]);
    const double T = stod(argv[2]);
    const double dt = stod(argv[3]);

    double L;
    int dim;
    get_potential_params(L, dim);

    printf("L=%f, len=%d\n", L, len);
    printf("T=%f, dt=%f\n", T, dt);
    const double stepsize = 2 * L / len;
    printf("stepsize=%f\n", stepsize);
    printf("Max threads: %d\n", omp_get_num_procs());
    printf("Threads: %d\n", omp_get_max_threads());

    nagd(dim, len, L, T, dt, stepsize);
}