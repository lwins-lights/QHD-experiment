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

def main(args):
    solvers = {1, 2, 3}   # 1: qhd, 2: subgrad, 3: lfmsgd
    results = {"qhd_suc": None, "subgrad_suc": None, "lfmsgd_suc": None, 
               "qhd_gap": None, "subgrad_gap": None, "lfmsgd_gap": None}
    
    func_name = os.path.splitext(os.path.basename(args.fpath))[0]
    data = pd.read_csv(data_file_path)
    row = data[data['Problem'] == func_name]
    if row.empty:
        return f"No data found for function: {func_name}"
    subgrad_L, lfmsgd_L, qhd_L = row[['subgrad_L', 'lfmsgd_L', 'qhd_L']].values[0]
    
    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "pseudospec"], cwd=simulator_path)
    run(["make", "lfmsgd"], cwd=simulator_path)
    run(["make", "subgrad"], cwd=simulator_path)

    # qhd
    if 1 in solvers:
      run(["./pseudospec", str(LEN), str(TOTAL_EVOLUTION_TIME), str(DT), str(PAR), str(qhd_L)], cwd=simulator_path)
      qhd_npz = np.load(os.path.join(result_path, "pseudospec.npz"))
      qhd_y = qhd_npz['expected_potential'][-1]
      qhd_yp = qhd_npz['probability_at_minimum'][-1]
      results["qhd_gap"] = round(qhd_y, 4)
      results["qhd_suc"] = round(qhd_yp, 4)

    # subgrad
    if 2 in solvers:
      run(["./subgrad", str(T), str(LR), str(PAR), str(SAMPLE_NUM), str(subgrad_L)], cwd=simulator_path)
      subgrad_npz = np.load(os.path.join(result_path, "subgrad.npz"))
      subgrad_idx = np.argmax(subgrad_npz['probability_at_minimum'])
      subgrad_yp = subgrad_npz['probability_at_minimum'][subgrad_idx]
      subgrad_y = subgrad_npz['expected_potential'][subgrad_idx]
      results["subgrad_gap"] = round(subgrad_y, 4)
      results["subgrad_suc"] = round(subgrad_yp, 4)

    # lfmsgd
    if 3 in solvers:
      run(["./lfmsgd", str(T), str(NL), str(PAR), str(SAMPLE_NUM), str(lfmsgd_L)], cwd=simulator_path)
      lfmsgd_npz = np.load(os.path.join(result_path, "lfmsgd.npz"))
      lfmsgd_idx = np.argmax(lfmsgd_npz['probability_at_minimum'])
      lfmsgd_y = lfmsgd_npz['expected_potential'][lfmsgd_idx]
      lfmsgd_yp = lfmsgd_npz['probability_at_minimum'][lfmsgd_idx]
      results["lfmsgd_gap"] = round(lfmsgd_y, 4)
      results["lfmsgd_suc"] = round(lfmsgd_yp, 4)

    row_index = data.index[data['Problem'] == func_name].tolist()[0]
    if 1 in solvers:
        data.at[row_index, 'qhd_suc'] = results['qhd_suc'].astype(float)
        data.at[row_index, 'qhd_gap'] = results['qhd_gap'].astype(float)
        print(f"qhd_suc: {results['qhd_suc']}, qhd_gap: {results['qhd_gap']}")
    if 2 in solvers:
        data.at[row_index, 'subgrad_suc'] = results['subgrad_suc'].astype(float)
        data.at[row_index, 'subgrad_gap'] = results['subgrad_gap'].astype(float)
        print(f"subgrad_suc: {results['subgrad_suc']}, subgrad_gap: {results['subgrad_gap']}")
    if 3 in solvers:
        data.at[row_index, 'lfmsgd_suc'] = results['lfmsgd_suc'].astype(float)
        data.at[row_index, 'lfmsgd_gap'] = results['lfmsgd_gap'].astype(float)
        print(f"lfmsgd_suc: {results['lfmsgd_suc']}, lfmsgd_gap: {results['lfmsgd_gap']}")
    data.to_csv(data_file_path, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help='path of the potential function .cpp file')
    args = parser.parse_args()
    main(args)
