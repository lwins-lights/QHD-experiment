import argparse
import os
from re import sub
import numpy as np
from subprocess import run
import pandas as pd

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "best_L_automated.csv")
simulator_path = os.path.join(root_path, "simulator")

# Constants
TOTAL_EVOLUTION_TIME = 10        # total evolution time for QHD
DT = 0.001                       # time discretization level for QHD
T = 10000                        # total number of updating steps for lfmsgd, subgrad, and TOTAL_EVOLUTION_TIME / DT for QHD
LEN = 256                        # number of cells in one dimension
PAR = 1                          # parallelism level for calculating success probability only
NL = 1                           # intensity of noise added to the subgradient for lfmsgd
LR = 1                           # learning rate coefficient for subgrad, the true learning rate will be LR/(T)^0.5
GRAN_NUM = 100                   # granularity of the final distribution graph 
SAMPLE_NUM = 5000                # number of samples for estimating average performance
K = 10                           # number of L to test in between L_min and L_max
ALPHA = 2                        # rate to shrink gap between L_min and L_max (gap := gap / alpha)
GAP_THRESHOLD = 1.2              # threshold to stop shrinking gap

def main(args):
    solvers = {1}    # 1: qhd, 2: subgrad, 3: lfmsgd
    qhd_L_min, qhd_L_max = 0.1, 100
    subgrad_L_min, subgrad_L_max = 1, 1000
    lfmsgd_L_min, lfmsgd_L_max = 10, 10000
    results = {"qhd_best": None, "subgrad_best": None, "lfmsgd_best": None}
    
    func_name = os.path.splitext(os.path.basename(args.fpath))[0]
    data = pd.read_csv(data_file_path)
    row = data[data['Problem'] == func_name]
    if row.empty:
        return f"No data found for function: {func_name}"
    
    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "pseudospec"], cwd=simulator_path)
    run(["make", "lfmsgd"], cwd=simulator_path)
    run(["make", "subgrad"], cwd=simulator_path)
    
    # qhd
    if 1 in solvers:
        gap = qhd_L_max / qhd_L_min
        L_best = qhd_L_min
        max_success_prob = -1
        while gap > GAP_THRESHOLD:
            L_List = np.geomspace(qhd_L_min, qhd_L_max, num=K)
            for L_temp in L_List:
                run(["./pseudospec", str(LEN), str(TOTAL_EVOLUTION_TIME), str(DT), str(PAR), str(L_temp)], cwd=simulator_path)
                npz = np.load(os.path.join(result_path, "pseudospec.npz"))
                qhd_yp = npz['probability_at_minimum']
                max_yp = qhd_yp[-1]
                if max_yp > max_success_prob:
                    L_best = L_temp
                    max_success_prob = max_yp
            gap /= ALPHA
            qhd_L_min = max(L_best / np.sqrt(gap), qhd_L_min)
            qhd_L_max = min(L_best * np.sqrt(gap), qhd_L_max)
            gap = qhd_L_max / qhd_L_min
        results["qhd_best"] = round(L_best, 4)
    
    # subgrad
    if 2 in solvers:
        gap = subgrad_L_max / subgrad_L_min
        L_best = subgrad_L_min
        max_success_prob = -1
        while gap > GAP_THRESHOLD:
            L_List = np.geomspace(subgrad_L_min, subgrad_L_max, num=K)
            for L_temp in L_List:
                run(["./subgrad", str(T), str(LR), str(PAR), str(SAMPLE_NUM), str(L_temp)], cwd=simulator_path)
                npz = np.load(os.path.join(result_path, "subgrad.npz"))
                subgrad_yp = npz['probability_at_minimum']
                max_yp = max(subgrad_yp)
                if max_yp > max_success_prob:
                    L_best = L_temp
                    max_success_prob = max_yp
            gap /= ALPHA
            subgrad_L_min = max(L_best / np.sqrt(gap), subgrad_L_min)
            subgrad_L_max = min(L_best * np.sqrt(gap), subgrad_L_max)
            gap = subgrad_L_max / subgrad_L_min
        results["subgrad_best"] = round(L_best, 4)
    
    # lfmsgd
    if 3 in solvers:
        gap = lfmsgd_L_max / lfmsgd_L_min
        L_best = lfmsgd_L_min
        max_success_prob = -1
        while gap > GAP_THRESHOLD:
            L_List = np.geomspace(lfmsgd_L_min, lfmsgd_L_max, num=K)
            for L_temp in L_List:
                run(["./lfmsgd", str(T), str(NL), str(PAR), str(SAMPLE_NUM), str(L_temp)], cwd=simulator_path)
                npz = np.load(os.path.join(result_path, "lfmsgd.npz"))
                lfmsgd_yp = npz['probability_at_minimum']
                max_yp = max(lfmsgd_yp)
                if max_yp > max_success_prob:
                    L_best = L_temp
                    max_success_prob = max_yp
            gap /= ALPHA
            lfmsgd_L_min = max(L_best / np.sqrt(gap), lfmsgd_L_min)
            lfmsgd_L_max = min(L_best * np.sqrt(gap), lfmsgd_L_max)
            gap = lfmsgd_L_max / lfmsgd_L_min
        results["lfmsgd_best"] = round(L_best, 4)
        
    row_index = data.index[data['Problem'] == func_name].tolist()[0]
    if 1 in solvers:
        data.at[row_index, 'qhd_L'] = results['qhd_best']
        print(f"qhd_best: {results['qhd_best']}")
    if 2 in solvers:
        data.at[row_index, 'subgrad_L'] = results['subgrad_best']
        print(f"subgrad_best: {results['subgrad_best']}")
    if 3 in solvers:
        data.at[row_index, 'lfmsgd_L'] = results['lfmsgd_best']
        print(f"lfmsgd_best: {results['lfmsgd_best']}")
    data.to_csv(data_file_path, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    args = parser.parse_args()
    main(args)

    