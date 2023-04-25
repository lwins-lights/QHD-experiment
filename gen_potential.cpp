#include <iostream>

using namespace std;

const string fout_name = "pot.txt";

double pot(const double x, const double y) {
    return 4 * x * x + 4 * y * y;
}

void gen_potential(double *V, const double L, const int num_cells) {

    double x, y;

    for (int i = 0; i < num_cells; i++) {
        for (int j = 0; j < num_cells; j++) {
            x = (double) i * 2 * L / num_cells - L;
            y = (double) j * 2 * L / num_cells - L;
            V[i * num_cells + j] = pot(x, y);
        }
    }
}

void print_to_file(const double *V, const int len, const char* fn) {

    FILE *fout;

    fout = fopen(fn, "w");
    for (int i = 0; i < len; i++) {
        fprintf(fout, "%lf ", V[i]);
    }
    fclose(fout);
}

int main(int argc, char **argv) {

    //cout << argc << endl;
    if (argc != 4) {
        perror("Expected arguments: ./gen_potential <L> <num_cells> <potential_filename>");
        exit(EXIT_FAILURE);
    }
    const double L = stod(argv[1]);
    const int num_cells = stoi(argv[2]);
    const char* fn = argv[3];

    double V[num_cells * num_cells];
    gen_potential(V, L, num_cells);
    print_to_file(V, num_cells * num_cells, fn);

    return 0;
}