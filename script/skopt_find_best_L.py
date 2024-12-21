from ast import List
from pickle import NONE
from util import run_subgrad, run_qhd_subgrad, run_lfmsgd, run_qhd
from skopt import gp_minimize
from skopt.space import Real
import numpy as np
from colorama import Fore, Back, Style, init
import os
from re import sub
import pandas as pd
import math
np.random.seed(0)

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "best_L_skopt.csv")

def lower_precision(value, relative_precision=1e-6):
    if value == 0:
        return 0
    # Determine the order of magnitude
    magnitude = 10 ** math.floor(math.log10(abs(value)))
    # Round to the desired relative precision
    return round(value / magnitude / relative_precision) * magnitude * relative_precision

def search_best_L(fpath: str, solvers: List, L_bounds: dict, suc_threshold: float, n_calls: int = 50):
    init()
    results = {"qhd_best_L": None, "subgrad_L": None, "lfmsgd_L": None, "minisubgrad_L": None,
               "qhd_suc": None, "subgrad_suc": None, "lfmsgd_suc": None, "pure_qhd_suc": None, 
               "qhd_gap": None, "subgrad_gap": None, "lfmsgd_gap": None, "pure_qhd_gap": None, 
               }
    #qhd_missing_subgrad_L = False
    L_list = []
    mean_List = []
    
    if 1 in solvers:  # if qhd is included

        # optimize qhd_L first

        def f(L):
            result = run_qhd(fpath, L=L[0])
            time_vs_mean = result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return lower_precision(final_mean)

        qhd_min, qhd_max = L_bounds["qhd_L_bounds"]
        search_space = [Real(qhd_min, qhd_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["qhd_L"] = round(best_L, 4)

        # run qhd again with this best L

        result = run_qhd(fpath, L=best_L)
        _, final_mean = result["time_vs_mean"][-1]
        results["pure_qhd_gap"] = round(final_mean, 4)
        results["pure_qhd_suc"] = round(sum(item["prob"] for item in result['output'] if item["value"] < suc_threshold), 4)

        # optimize minisubgrad_L next

        def g(L):
            result = run_qhd_subgrad(fpath, skip_qhd=True, subgrad_L=L[0])
            time_vs_mean = result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return lower_precision(final_mean)
        
        minisubgrad_min, minisubgrad_max = L_bounds["minisubgrad_L_bounds"]
        search_space = [Real(minisubgrad_min, minisubgrad_max, 'log-uniform')]
        optimize_result = gp_minimize(g, search_space, n_calls=n_calls, random_state=0)
        mini_best_L = optimize_result.x[0]
        results["minisubgrad_L"] = round(mini_best_L, 4)

        # run again

        result = run_qhd_subgrad(fpath, qhd_L=best_L, subgrad_L=mini_best_L)
        _, final_mean = result["time_vs_mean"][-1]
        results["qhd_gap"] = round(final_mean, 4)
        results["qhd_suc"] = round(sum(item["prob"] for item in result['output'] if item["value"] < suc_threshold), 4)
    
    if 2 in solvers:  # if subgrad is included
        subgrad_min, subgrad_max = L_bounds["subgrad_L_bounds"]
        def f(L):
            subgrad_result = run_subgrad(fpath, tot_steps=20000, L=L[0])
            time_vs_mean = subgrad_result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return final_mean
        
        search_space = [Real(subgrad_min, subgrad_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["subgrad_L"] = round(best_L, 4)
        
        # run subgrad again with this best L
        subgrad_result = run_subgrad(fpath, tot_steps=20000, L=best_L)
        _, final_mean = subgrad_result["time_vs_mean"][-1]
        results["subgrad_gap"] = round(final_mean, 4)
        results["subgrad_suc"] = round(sum(item["prob"] for item in subgrad_result['output'] if item["value"] < suc_threshold), 4)
    
    if 3 in solvers:  # if lfmsgd is included
        lfmsgd_min, lfmsgd_max = L_bounds["lfmsgd_L_bounds"]
        def f(L):
            lfmsgd_result = run_lfmsgd(fpath, tot_steps=20000, L=L[0])
            time_vs_mean = lfmsgd_result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return final_mean
        
        search_space = [Real(lfmsgd_min, lfmsgd_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["lfmsgd_L"] = round(best_L, 4)
        
        # run lfmsgd again with this best L
        lfmsgd_result = run_lfmsgd(fpath, tot_steps=20000, L=best_L)
        _, final_mean = lfmsgd_result["time_vs_mean"][-1]
        results["lfmsgd_gap"] = round(final_mean, 4)
        results["lfmsgd_suc"] = round(sum(item["prob"] for item in lfmsgd_result['output'] if item["value"] < suc_threshold), 4)

    '''
    if 4 in solvers:  # if pure qhd (no subgrad involved) is included
        qhd_min, qhd_max = L_bounds["qhd_L_bounds"]
        def f(L):
            qhd_result = run_qhd(fpath, T=20, L=L[0])
            time_vs_mean = qhd_result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return lower_precision(final_mean)
        
        search_space = [Real(qhd_min, qhd_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["pure_qhd_best_L"] = round(best_L, 4)
        
        # run qhd again with this best L
        qhd_result = run_qhd(fpath, T=20, L=best_L)
        _, final_mean = qhd_result["time_vs_mean"][-1]
        results["pure_qhd_gap"] = round(final_mean, 4)
        results["pure_qhd_suc"] = round(sum(item["prob"] for item in qhd_result['output'] if item["value"] < suc_threshold), 4)
    '''
    
    # store outputs
    func_name = os.path.splitext(os.path.basename(fpath))[0]
    data = pd.read_csv(data_file_path)
    columns_to_convert = data.columns.difference(['Problem'])  # All columns except 'Problem'
    data[columns_to_convert] = data[columns_to_convert].astype(float)  # Convert these columns to float
    row_index = data.index[data['Problem'] == func_name].tolist()[0]
    '''
    if qhd_missing_subgrad_L: 
        print(f"QHD is specified but subgrad_L is missing, skip qhd_subgrad search.")
    '''
        
    if 1 in solvers:
        for item in ["qhd_L", "minisubgrad_L", "qhd_suc", "qhd_gap", "pure_qhd_suc", "pure_qhd_gap"]:
            data.at[row_index, item] = results[item]
            print(f"{item}: {results[item]}")

    if 2 in solvers:
        data.at[row_index, 'subgrad_L'] = results["subgrad_L"]
        data.at[row_index, 'subgrad_suc'] = results['subgrad_suc']
        data.at[row_index, 'subgrad_gap'] = results['subgrad_gap']
        print(f"subgrad_L: {results['subgrad_L']}, subgrad_gap: {results['subgrad_gap']}, subgrad_suc: {results['subgrad_suc']}")
        
    if 3 in solvers:
        data.at[row_index, 'lfmsgd_L'] = results["lfmsgd_L"]
        data.at[row_index, 'lfmsgd_suc'] = results['lfmsgd_suc']
        data.at[row_index, 'lfmsgd_gap'] = results['lfmsgd_gap']
        print(f"lfmsgd_L: {results['lfmsgd_L']}, lfmsgd_gap: {results['lfmsgd_gap']}, lfmsgd_suc: {results['lfmsgd_suc']}")

    '''    
    if 4 in solvers:
        data.at[row_index, 'pure_qhd_L'] = results["pure_qhd_best_L"]
        data.at[row_index, 'pure_qhd_suc'] = results['pure_qhd_suc']
        data.at[row_index, 'pure_qhd_gap'] = results['pure_qhd_gap']
        print(f"pure_qhd_best_L: {results['pure_qhd_best_L']}, pure_qhd_gap: {results['pure_qhd_gap']}, pure_qhd_suc: {results['pure_qhd_suc']}")
    '''
        
    data.to_csv(data_file_path, index=False)
    return results

if __name__ == '__main__':
    fpath = 'func/nonsmooth/keane.cpp'
    solvers = [1, 2, 3]                                # 1: qhd, 2: subgrad, 3: lfmsgd, 4: pure qhd (no subgrad involved)
    suc_threshold = 0.05
    qhd_L_bounds = 0.1, 100
    subgrad_L_bounds = 0.1, 100
    minisubgrad_L_bounds = 0.1, 10000
    lfmsgd_L_bounds = 0.1, 100
    L_bounds = {"qhd_L_bounds": qhd_L_bounds, "subgrad_L_bounds": subgrad_L_bounds, "lfmsgd_L_bounds": lfmsgd_L_bounds, "minisubgrad_L_bounds": minisubgrad_L_bounds}
    results = search_best_L(fpath, solvers, L_bounds, suc_threshold, n_calls=50)