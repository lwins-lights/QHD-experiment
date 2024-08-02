import os, sys, time
import subprocess
import numpy as np
from termcolor import colored

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
CUR_DIR = os.getcwd()

def vmsg(s, v):
    if v == False:
        return
    print(colored(s, 'blue', 'on_black'))

def gen_F_str(n, u, M):
    s = "\n".join([
        r"#include <cmath>",
        r"using namespace std;",
        r"const double L = %f;" % (u * 0.5),
        r"const int dim = %d;" % n,
        r"void get_F_params(double &var_L, int &var_dim) {var_L = L; var_dim = dim;}",
        r"void get_F(const double *x, double *F) {",
        r"    double z[dim], s[dim];",
        r"    for (int i = 0; i < dim; i++) z[i] = x[i] + %f;" % (u * 0.5) + "\n", 
    ])
    for i in range(n):
        t = r"    s[%d] = %d" % (i, n if i == n - 1 else 0);
        for j in range(n):
            t += r" + (%f) * z[%d]" % (M[i][j], j)
        t += ";\n"
        s += t
    s += r"    for (int i = 0; i < dim; i++) F[i] = z[i] * s[i];"
    s += "\n}"
    return s;

def live_call(x, y, verbose):
    pipe = None if verbose else subprocess.PIPE
    process = subprocess.Popen(x, stdout=pipe, stderr=pipe, cwd=y)
    process.wait()
    '''
    for c in iter(lambda: process.stdout.read(1), b""):
        if verbose:
            sys.stdout.buffer.write(c)
    '''

def qcp_run(n, M, ub=2, eta_qhd=0.05, eta=0.05, res=10, disc=0.01, rand=0.5, verbose=False):

    vmsg("Dumping the LO object ...", verbose)
    with open(os.path.join(THIS_DIR, "..", "simulator", "qipopt_F.cpp"), "w") as f:
        f.write(gen_F_str(n, ub, M))

    #if verbose:
    #    time.sleep(1)
    vmsg("Compiling ...", verbose)
    live_call(['make', 'qipopt'], os.path.join(THIS_DIR, "..", "simulator"), verbose)

    #if verbose:
    #    time.sleep(1)
    vmsg("Executing ...", verbose)
    live_call(['./qipopt', str(res), str(eta_qhd), str(eta), str(disc), str(rand)], os.path.join(THIS_DIR, "..", "simulator"), verbose)

    vmsg("Retrieving results ...", verbose)
    res = np.load(os.path.join(THIS_DIR, "..", "result", "qipopt.npz"))
    z = res['sample']
    z = [z[i] + ub / 2 for i in range(len(z))]
    return z

if __name__ == '__main__':
    n = 5
    u = 2
    M = [
        [0, 0, 1, -1, 1],
        [0, 0, 0, 1, 0],
        [-1, 0, 0, 2, 0],
        [1, -1, -2, 0, 3],
        [-1, 0, 0, -3, 0]
    ]
    z = qcp_run(n, M, verbose=False, res=10)
    print(z)
