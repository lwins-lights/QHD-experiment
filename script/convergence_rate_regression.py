from cProfile import label
import dis
from tracemalloc import start
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from pyparsing import line
from regex import F
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from subprocess import run
import argparse
import dis
import os
import pickle
from datetime import datetime
import copy

# Path initialization
this_path = os.path.dirname(os.path.realpath(__file__))
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
output_path = os.path.join(result_path, "convergence_rate.pdf")
data_file_path = os.path.join(result_path, "convergence_rate.csv")
simulator_path = os.path.join(root_path, "simulator")
assets_path = os.path.join(result_path, "assets.pkl")

# Constants
DISCRETIZATION_NUM = 100
DT = 0.0001
T = 100
PAR = 1
L = 1
NL = 1                           # intensity of noise added to the subgradient for lfmsgd
LR = 1                           # learning rate coefficient for subgrad, the true learning rate will be LR/(T)^0.5
SAMPLE_NUM = 5000                # number of samples for estimating average performance
FORCE = False
SOLVER = "qhd"

def main(args):
    def root_mean_squared_error(y_true, y_pred):
        return np.sqrt(mean_squared_error(y_true, y_pred))
    
    func_name = os.path.splitext(os.path.basename(args.fpath))[0]
    
    if os.path.exists(assets_path):
        try:
            with open(assets_path, 'rb') as file:
                assets = pickle.load(file)
        except (FileNotFoundError, EOFError):
            assets = {}
    else:
        assets = {}
    
    if SOLVER == "qhd":
        discretization_num = DISCRETIZATION_NUM
    else:
        discretization_num = -1
    
    # Load data and filter by PROBLEM
    data = pd.read_csv(data_file_path)
    data_row = data[(data["Problem"] == func_name) & (data["dt"] == DT) & (data["T"] == T) & (data["solver"] == SOLVER) & (data["discretization_num"] == discretization_num)]
    assert len(data_row) == 1, "Data not found or multiple data found"
    start_time = data_row["start"].values[0]
    end_time = data_row["end"].values[0]
    x = np.arange(0, T, DT)
    start_index = int(start_time // DT - 1)
    end_index = int(end_time // DT - 1)
    x_regression = x[start_index:end_index+1]
    npz = None
    
    # get objective function 
    if SOLVER == "qhd":
        run(["mkdir", "-p", result_path])
        run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
        run(["make", "pseudospec"], cwd=simulator_path)
        arglist_original = ["./pseudospec", str(DISCRETIZATION_NUM), str(T), str(DT), str(PAR), str(L)]
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
            assets[arglist_str]["data"] = {key: copy.deepcopy(npz[key]) for key in npz.keys()}  # deep copy the data
            
        with open(assets_path, 'wb') as file:
            pickle.dump(assets, file)
    elif SOLVER == "lfmsgd":
        run(["mkdir", "-p", result_path])
        run(["cp", args.fpath, os.path.join(simulator_path, "potential.cpp")])
        run(["make", "lfmsgd"], cwd=simulator_path)
        run(["./lfmsgd", str(T), str(NL), str(PAR), str(SAMPLE_NUM), str(L)], cwd=simulator_path)
        npz = np.load(os.path.join(result_path, "lfmsgd.npz"))

    # Log transformation
    y_regression = npz['expected_potential'][start_index:end_index+1]
    y = npz['expected_potential']
    x_log = np.log10(x_regression).reshape(-1, 1)
    y_log = np.log10(y_regression)

    # Regression
    model = LinearRegression()
    model.fit(x_log, y_log)
    slope = model.coef_[0]
    intercept = model.intercept_
    predict_y = np.power(10, model.predict(x_log))
    regression_error = root_mean_squared_error(y_regression, predict_y)

    # Print regression parameters and error
    print(f"Regression Slope: {slope}")
    print(f"Regression Intercept: {intercept}")
    print(f"Regression RMSE: {regression_error}")

    # Plotting
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 11  
    plt.rcParams['axes.titlesize'] = 16  
    plt.rcParams['axes.labelsize'] = 14  
    plt.rcParams['legend.fontsize'] = 10  

    _, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Plot: Objective vs. Time with bounded log-log regression 
    ax.set_xlabel("Time")
    ax.set_ylabel("Expected Potential")
    if SOLVER == "qhd":
        ax.set_title(f"{SOLVER.upper()}'s Convergence Rate at {func_name}\nwith Discretization Num = {DISCRETIZATION_NUM}")
    else:
        ax.set_title(f"{SOLVER.upper()}'s Convergence Rate at {func_name}")
    ax.plot(x, y, label="Convergence Curve", linewidth=4)
    ax.plot(x_regression, predict_y, label=f"Log-Log Regression (k = {slope:.2f})", color='red', linewidth=2)

    ax.axvline(x=start_time, color='green', linestyle='--', label="Left Boundary for Regression", linewidth=2)
    ax.axvline(x=end_time, color='orange', linestyle='--', label="Right Boundary for Regression", linewidth=2)

    # Set log scales
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Add legend
    ax.legend(loc="lower left")

    # Save and show
    plt.savefig(output_path)
    plt.show()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fpath", type=str, required=True, help="Path of the potential function .cpp file")
    args = parser.parse_args()
    main(args)
