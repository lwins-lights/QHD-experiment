import argparse
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pickle
import copy
from seaborn import color_palette
from subprocess import run
from datetime import datetime

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
    #f, (e_plt, p_plt, d_plt) = plt.subplots(3, 1, figsize=(800/dpi, 900/dpi))
    plt.rcParams["figure.figsize"] = (1600/dpi, 900/dpi)
    e_plt = plt.subplot(2, 2, 1)
    p_plt = plt.subplot(2, 2, 2)
    d_plt = plt.subplot(2, 1, 2)
    palette = color_palette("husl", len(L_list))

    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "pseudospec"], cwd=simulator_path)

    try:
        with open(args.assets, 'rb') as file:
            assets = pickle.load(file)
    except (FileNotFoundError, EOFError):
        assets = {}

    for i, L in enumerate(L_list):
        arglist_original = ["./pseudospec", str(args.len), str(args.T), str(args.dt), str(args.par), str(L)]
        arglist = [args.fpath] + arglist_original
        arglist_str = "; ".join(arglist)

        if arglist_str in assets and not args.force:
            print("Past result found in assets:\narglist   =  %s\ntimestamp =  %s\n" % (arglist_str, assets[arglist_str]["timestamp"]))
            npz = assets[arglist_str]["data"]
        else:
            run(arglist_original, cwd=simulator_path)
            npz = np.load(os.path.join(result_path, "pseudospec.npz"))
            assets[arglist_str] = {}
            assets[arglist_str]["timestamp"] = datetime.now()
            assets[arglist_str]["data"] = {key: copy.deepcopy(npz[key]) for key in npz.keys()} # deep copy the data
        
        x = np.arange(0, args.T, args.dt)
        y = npz['expected_potential']
        yp = npz['probability_at_minimum']
        e_plt.plot(x, y, label=str(L), color=palette[i])
        p_plt.plot(x, yp, label=str(L), color=palette[i])

        v = npz['V']
        d = npz['dist']
        vd = [[v[j], d[j]] for j in range(len(v))]
        vd.sort(key=lambda p:p[0])
        v = [p[0] for p in vd]
        d = [p[1] for p in vd]
        inc = (v[-1] - v[0]) / (args.ngran - 1)
        x = []
        y = []
        cur_ind = 0
        for f_low in np.arange(v[0], v[-1] + inc / 2, inc):
            f_high = f_low + inc
            cur_prob = 0
            while cur_ind < len(v) and v[cur_ind] < f_high:
                cur_prob += d[cur_ind]
                cur_ind += 1
            x.append(f_low)
            y.append(cur_prob * args.ngran)
        #print((len(x), len(y)))
        d_plt.plot(x, y, label=str(L), color=palette[i])

    with open(args.assets, 'wb') as file:
        pickle.dump(assets, file)

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
    #p_plt.legend(loc="lower right", title="L", title_fontsize="8")
    plt_title = 'QHD on %s\n dt=%f; gran.=%f^%d\n dist_gran=%d' % (args.fpath, args.dt, args.len, npz['dim'][0], args.ngran)
    #f.suptitle(plt_title)
    plt.suptitle(plt_title)

    plt.savefig(args.output)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--len", type=int, default=256, help='space discretization level, i.e., number of cells in one dimension')
    parser.add_argument("--T", type=float, default=10, help='total evolution time')
    parser.add_argument("--dt", type=float, default=0.001, help='time discretization level, i.e., time increment for each simulation step')
    parser.add_argument("--par", type=int, default=1, help='parallelism level for calculating success probability only')
    parser.add_argument("--Llist", type=str, required=True, help='sizes of the hypercubes in which QHD runs, separated by commas')
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    parser.add_argument("--output", type=str, default=os.path.join(result_path, "qhd.png"), help='path of the output .png file')
    parser.add_argument("--ngran", type=int, default=100, help='granularity of the final distribution graph (default = 100)')
    parser.add_argument("--assets", type=str, default=os.path.join(result_path, "qhd_assets.pkl"), help='path of the assets file')
    parser.add_argument("--force", action='store_true', default=False, help='force not loading from assets')
    args = parser.parse_args()
    main(args)