import argparse
import dis
import os
import pickle
from re import sub
import numpy as np
from subprocess import run
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import copy

# Path initialization
this_path = os.path.dirname(os.path.realpath(__file__))
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
output_path = os.path.join(result_path, "discretization_error.pdf")
data_file_path = os.path.join(result_path, "discretization_error.csv")
simulator_path = os.path.join(root_path, "simulator")
assets_path = os.path.join(result_path, "assets.pkl")

# Constants
TOTAL_EVOLUTION_TIME = 100      # total evolution time for QHD
DT = 0.0001                      # time discretization level for QHD
PAR = 1                          # parallelism level for calculating success probability only
GRAN_NUM = 100                   # granularity of the final distribution graph 
SAMPLE_NUM = 5000                # number of samples for estimating average performance
MAX_SUBFIGURES = 9               # Maximum subfigures for combined plot
FORCE = False                     # Force re-run if past result found in assets.pkl

def main(args):
    data = pd.read_csv(data_file_path)
    discretization_num_list = [100]     # to be changed
    discretization_err_list = []
    expectation_list = []

    x = np.arange(0, TOTAL_EVOLUTION_TIME, DT)

    func_name = os.path.splitext(os.path.basename(args.fpath))[0]
    L = None
    if func_name == "Abs":
        L = 1
    elif func_name == "Expabs":
        L = 1
    else:
        assert False, f"Unsupported function: {func_name}"

    if os.path.exists(assets_path):
        try:
            with open(assets_path, 'rb') as file:
                assets = pickle.load(file)
        except (FileNotFoundError, EOFError):
            assets = {}
    else:
        assets = {}

    run(["mkdir", "-p", result_path])
    run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
    run(["make", "pseudospec"], cwd=simulator_path)

    for discretization_num in discretization_num_list:
        arglist_original = ["./pseudospec", str(discretization_num), str(TOTAL_EVOLUTION_TIME), str(DT), str(PAR), str(L)]
        arglist = [args.fpath] + arglist_original
        arglist_str = "; ".join(arglist)

        if arglist_str in assets and not FORCE:
            print("Past result found in assets:\narglist   =  %s\ntimestamp =  %s\n" % (arglist_str, assets[arglist_str]["timestamp"]))
            npz = assets[arglist_str]["data"]
        else:
            run(arglist_original, cwd=simulator_path)
            npz = np.load(os.path.join(result_path, "pseudospec.npz"))
            assets[arglist_str] = {}
            assets[arglist_str]["timestamp"] = datetime.now()
            assets[arglist_str]["data"] = {key: copy.deepcopy(npz[key]) for key in npz.keys()} # deep copy the data

        expectation_list.append(npz['expected_potential'])
        discretization_err_list.append(npz['expected_potential'][-1])
        
        row_exists = (
            (data["Problem"] == func_name) &
            (data["discretization_num"] == discretization_num) & 
            (data["dt"] == DT)
        ).any()

        if not row_exists:
            new_row = {
                "Problem": func_name,
                "discretization_error": discretization_err_list[-1],
                "discretization_num": discretization_num,
                "L": L,
                "dt": DT
            }
            data = pd.concat([data, pd.DataFrame([new_row])], ignore_index=True)
    
    data.to_csv(data_file_path, index=False)
    with open(assets_path, 'wb') as file:
        pickle.dump(assets, file)

    # Plot subfigures with log-log scaling
    num_subfigures = min(len(discretization_num_list), MAX_SUBFIGURES)
    fig, axes = plt.subplots(1, num_subfigures, figsize=(5 * num_subfigures, 5))
    if num_subfigures == 1:
        axes = [axes]
    for idx, ax in enumerate(axes[:num_subfigures]):
        ax.plot(x, expectation_list[idx], label=f"Discretization {discretization_num_list[idx]}")
        ax.set_xscale('log')  # Set log scale for x-axis
        ax.set_yscale('log')  # Set log scale for y-axis
        ax.set_title(f"Error: {discretization_err_list[idx]:.4f}\nNum: {discretization_num_list[idx]}")
        ax.set_xlabel("Time")
        ax.set_ylabel("Expected Potential")
        ax.legend(loc="lower left")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, _: f"$10^{{{int(np.log10(val))}}}$"))

    if len(discretization_num_list) > MAX_SUBFIGURES:
        print(f"Warning: Only the first {MAX_SUBFIGURES} subfigures are plotted.")

    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()
    print(f"discretization_err_list: {discretization_err_list}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help="Path of the potential function .cpp file")
    args = parser.parse_args()
    main(args)
