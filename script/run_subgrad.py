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
    plt.rcParams["figure.figsize"] = (1600/dpi, 900/dpi)
    e_plt = plt.subplot(2, 2, 1)
    p_plt = plt.subplot(2, 2, 2)
    d_plt = plt.subplot(2, 1, 2)
    palette = color_palette("husl", len(L_list))

    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "subgrad"], cwd=simulator_path)
    for i, L in enumerate(L_list):
        if args.useqhd:
            qhd_flag = 1
        else:
            qhd_flag = 0
        run(["./subgrad", str(args.tot), str(args.lr), str(args.par), str(args.m), str(L), str(qhd_flag)], cwd=simulator_path)

        npz = np.load(os.path.join(result_path, "subgrad.npz"))
        x = np.arange(0, args.tot, 1)
        y = npz['expected_potential']
        yp = npz['probability_at_minimum']
        e_plt.plot(x, y, label=str(L), color=palette[i])
        p_plt.plot(x, yp, label=str(L), color=palette[i])

        #print(npz['expected_potential'][0])
        samples = npz['samples']
        samples.sort()
        inc = (npz['expected_potential'][0] - samples[0]) / (args.ngran - 1)
        x = []
        y = []
        cur_ind = 0
        for f_low in np.arange(samples[0], npz['expected_potential'][0] + inc / 2, inc):
            f_high = f_low + inc
            cur_prob = 0
            while cur_ind < len(samples) and samples[cur_ind] < f_high:
                cur_prob += 1 / len(samples)
                cur_ind += 1
            x.append(f_low)
            y.append(cur_prob * args.ngran)
        d_plt.plot(x, y, label=str(L), color=palette[i])

    # plot style
    e_plt.set_xlabel('time')
    e_plt.set_ylabel('expectation')
    p_plt.set_xlabel('time')
    p_plt.set_ylabel('success probability')
    d_plt.set_xlabel('returned value')
    d_plt.set_ylabel('density')
    e_plt.yaxis.tick_right()
    p_plt.yaxis.tick_right()
    d_plt.yaxis.tick_right()
    e_plt.legend(loc="upper right", title="L", title_fontsize="8", fontsize="8")
    plt_title = 'SUBGRAD on %s\n lr_coef=%f; #sample=%d\n dist_gran=%d' % (args.fpath, args.lr, args.m, args.ngran)
    plt.suptitle(plt_title)

    plt.savefig(args.output)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tot", type=int, default=10000, help='total number of updating steps')
    parser.add_argument("--lr", type=float, default=1, help='learning rate coefficient; the true learning rate will be LR/(TOT)^0.5')
    parser.add_argument("--par", type=int, default=1, help='parallelism level for calculating success probability only')
    parser.add_argument("--m", type=int, default=1000, help='number of samples for estimating average performance')
    parser.add_argument("--Llist", type=str, required=True, help='sizes of the hypercubes in which SUBGRAD runs, separated by commas')
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    parser.add_argument("--output", type=str, default=os.path.join(result_path, "subgrad.png"), help='path of the output .png file')
    parser.add_argument("--ngran", type=int, default=100, help='granularity of the final distribution graph (default = 100)')
    parser.add_argument("--useqhd", action='store_true', default=False, help='use QHD output as input')
    args = parser.parse_args()
    main(args)