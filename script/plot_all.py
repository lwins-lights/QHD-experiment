import argparse
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from seaborn import color_palette
from subprocess import run
import pandas as pd

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
simulator_path = os.path.join(root_path, "simulator")
data_file_path = os.path.join(result_path, "best_L_automated.csv")

# Constants
SOLVER_NUM = 3                   # number of solvers to compare
TOTAL_EVOLUTION_TIME = 10        # total evolution time for QHD
DT = 0.001                       # time discretization level for QHD
T = 10000                        # total number of updating steps for lfmsgd, subgrad, and TOTAL_EVOLUTION_TIME / DT for QHD
LEN = 256                        # number of cells in one dimension
PAR = 1                          # parallelism level for calculating success probability only
NL = 1                           # intensity of noise added to the subgradient for lfmsgd
LR = 1                           # learning rate coefficient for subgrad, the true learning rate will be LR/(T)^0.5
GRAN_NUM = 100                   # granularity of the final distribution graph 
SAMPLE_NUM = 5000                # number of samples for estimating average performance

def main(args):
    # figure initialization
    func_name = os.path.splitext(os.path.basename(args.fpath))[0]
    output_path = os.path.join(result_path, func_name, f'{func_name}_integrated')
    dpi = mpl.rcParams['figure.dpi']
    plt.rcParams["figure.figsize"] = (1600/dpi, 900/dpi)
    e_plt = plt.subplot(2, 2, 1)
    p_plt = plt.subplot(2, 2, 2)
    d_plt = plt.subplot(2, 1, 2)
    palette = color_palette("husl", SOLVER_NUM)
    data = pd.read_csv(data_file_path)
    row = data[data['Problem'] == func_name]
    if row.empty:
        return f"No data found for function: {func_name}"
    subgrad_L, lfmsgd_L, qhd_L = row[['subgrad_L', 'lfmsgd_L', 'qhd_L']].values[0]
    
    # make
    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "pseudospec"], cwd=simulator_path)
    run(["make", "lfmsgd"], cwd=simulator_path)
    run(["make", "subgrad"], cwd=simulator_path)
    
    # x-axis initialization
    time_shared = np.arange(0, T, 1)         # shared time axis for suceess probability and expectation
    
    # run qhd
    run(["./pseudospec", str(LEN), str(TOTAL_EVOLUTION_TIME), str(DT), str(PAR), str(qhd_L)], cwd=simulator_path)
    npz = np.load(os.path.join(result_path, "pseudospec.npz"))
    qhd_y = npz['expected_potential']
    qhd_yp = npz['probability_at_minimum']

    v = npz['V']
    d = npz['dist']
    vd = [[v[j], d[j]] for j in range(len(v))]
    vd.sort(key=lambda p:p[0])
    v = [p[0] for p in vd]
    d = [p[1] for p in vd]
    inc = (v[-1] - v[0]) / (GRAN_NUM - 1)
    qhd_returned_value = []
    qhd_density = []
    cur_ind = 0
    for f_low in np.arange(v[0], v[-1] + inc / 2, inc):
        f_high = f_low + inc
        cur_prob = 0
        while cur_ind < len(v) and v[cur_ind] < f_high:
            cur_prob += d[cur_ind]
            cur_ind += 1
        qhd_returned_value.append(f_low)
        qhd_density.append(cur_prob * GRAN_NUM)
    
    # run lfmsgd
    run(["./lfmsgd", str(T), str(NL), str(PAR), str(SAMPLE_NUM), str(lfmsgd_L)], cwd=simulator_path)
    npz = np.load(os.path.join(result_path, "lfmsgd.npz"))
    lfmsgd_y = npz['expected_potential']
    lfmsgd_yp = npz['probability_at_minimum']
    
    samples = npz['samples']
    samples.sort()
    inc = (npz['expected_potential'][0] - samples[0]) / (GRAN_NUM - 1)
    lfmsgd_returned_value = []
    lfmsgd_density = []
    cur_ind = 0
    for f_low in np.arange(samples[0], npz['expected_potential'][0] + inc / 2, inc):
        f_high = f_low + inc
        cur_prob = 0
        while cur_ind < len(samples) and samples[cur_ind] < f_high:
            cur_prob += 1 / len(samples)
            cur_ind += 1
        lfmsgd_returned_value.append(f_low)
        lfmsgd_density.append(cur_prob * GRAN_NUM)
    
    # run subgrad
    run(["./subgrad", str(T), str(LR), str(PAR), str(SAMPLE_NUM), str(subgrad_L)], cwd=simulator_path)
    npz = np.load(os.path.join(result_path, "subgrad.npz"))
    subgrad_y = npz['expected_potential']
    subgrad_yp = npz['probability_at_minimum']
    
    samples = npz['samples']
    samples.sort()
    inc = (npz['expected_potential'][0] - samples[0]) / (GRAN_NUM - 1)
    subgrad_returned_value = []
    subgrad_density = []
    cur_ind = 0
    for f_low in np.arange(samples[0], npz['expected_potential'][0] + inc / 2, inc):
        f_high = f_low + inc
        cur_prob = 0
        while cur_ind < len(samples) and samples[cur_ind] < f_high:
            cur_prob += 1 / len(samples)
            cur_ind += 1
        subgrad_returned_value.append(f_low)
        subgrad_density.append(cur_prob * GRAN_NUM)
    
    # plot + style
    [qhd_color, lfmsgd_color, subgrad_color]= palette
    e_plt.plot(time_shared, qhd_y, label="QHD", color=qhd_color)
    e_plt.plot(time_shared, lfmsgd_y, label="LFMSGD", color=lfmsgd_color)
    e_plt.plot(time_shared, subgrad_y, label="SUBGRAD", color=subgrad_color)
    
    p_plt.plot(time_shared, qhd_yp, label="QHD", color=qhd_color)
    p_plt.plot(time_shared, lfmsgd_yp, label="LFMSGD", color=lfmsgd_color)
    p_plt.plot(time_shared, subgrad_yp, label="SUBGRAD", color=subgrad_color)
    
    d_plt.plot(qhd_returned_value, qhd_density, label="QHD", color=qhd_color)
    d_plt.plot(lfmsgd_returned_value, lfmsgd_density, label="LFMSGD", color=lfmsgd_color)
    d_plt.plot(subgrad_returned_value, subgrad_density, label="SUBGRAD", color=subgrad_color)
    
    e_plt.set_xlabel('time')
    e_plt.set_ylabel('expectation')
    p_plt.set_xlabel('time')
    p_plt.set_ylabel('success probability')
    d_plt.set_xlabel('returned value')
    d_plt.set_ylabel('density')
    e_plt.yaxis.tick_right()
    p_plt.yaxis.tick_right()
    d_plt.yaxis.tick_right()
    e_plt.legend(loc="upper right", fontsize="8")
    p_plt.legend(loc="lower right", fontsize="8")
    d_plt.legend(loc="upper center", fontsize="8")
    plt_title = f'Performance of QHD (L = {qhd_L}), LFMSGD (L = {lfmsgd_L}), SUBGRAD (L = {subgrad_L}) on {func_name[0].upper() + func_name[1:]}'
    plt.suptitle(plt_title)

    plt.savefig(output_path)
    plt.show()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    args = parser.parse_args()
    main(args)