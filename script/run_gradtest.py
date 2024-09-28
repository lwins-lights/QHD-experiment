import argparse
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from seaborn import color_palette
from subprocess import run

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
simulator_path = os.path.join(root_path, "simulator")

def main(args):
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "gradtest"], cwd=simulator_path)
    run(["./gradtest", str(args.L), str(args.eps), str(args.num)], cwd=simulator_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    parser.add_argument("--L", type=float, default=1, help='specifies L for the potential function (default = 1)')
    parser.add_argument("--eps", type=float, default=1e-5, help='a small value used in the finite difference method; also used as the threshold to decide whether the finite difference matches the subgradient (default = 1e-5)')
    parser.add_argument("--num", type=int, default=10000, help='the number of random points to check (default = 10000)')
    args = parser.parse_args()
    main(args)