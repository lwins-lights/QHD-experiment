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
result_path = os.path.join(root_path, "result")
simulator_path = os.path.join(root_path, "simulator")

def main(args):
    # var init
    L_list = [float(x) for x in args.Llist.split(',')]

    # plot init
    dpi = mpl.rcParams['figure.dpi']
    f, (e_plt, p_plt) = plt.subplots(2, 1, figsize=(800/dpi, 900/dpi))
    palette = color_palette("husl", len(L_list))

    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "subgrad"], cwd=simulator_path)
    for i, L in enumerate(L_list):
        run(["./subgrad", str(args.tot), str(args.lr), str(args.par), str(args.m), str(L)], cwd=simulator_path)
        npz = np.load(os.path.join(result_path, "subgrad.npz"))
        x = np.arange(0, args.tot, 1)
        y = npz['expected_potential']
        yp = npz['probability_at_minimum']
        e_plt.plot(x, y, label=str(L), color=palette[i])
        p_plt.plot(x, yp, label=str(L), color=palette[i])

    # plot style
    e_plt.set_xlabel('step')
    e_plt.set_ylabel('expectation')
    p_plt.set_xlabel('step')
    p_plt.set_ylabel('success probability')
    e_plt.yaxis.tick_right()
    p_plt.yaxis.tick_right()
    e_plt.legend(loc="upper right", title="L", title_fontsize="8", fontsize="8")
    #p_plt.legend(loc="lower right", title="L", title_fontsize="8")
    plt_title = 'SUBGRAD on %s\n lr_coef=%f; #sample=%d;' % (args.fpath, args.lr, args.m)
    f.suptitle(plt_title)

    plt.savefig(args.output)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tot", type=int, default=10000, help='total number of updating steps')
    parser.add_argument("--lr", type=float, default=1, help='learning rate coefficient; the true learning rate will be LR/(TOT)^0.5')
    parser.add_argument("--par", type=int, default=1, help='parallelism level for calculating success probability only')
    parser.add_argument("--m", type=int, default=1000, help='number of samples for estimating average performance')
    parser.add_argument("--Llist", type=str, required=True, help='sizes of the hypercubes in which QHD runs, separated by commas')
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    parser.add_argument("--output", type=str, default=os.path.join(result_path, "qhd.png"), help='path of the output .png file')
    args = parser.parse_args()
    main(args)