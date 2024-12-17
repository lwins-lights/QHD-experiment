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
np.random.seed(0)

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "best_L_skopt_0.05.csv")

def search_best_L(fpath: str, solvers: List, L_bounds: dict, suc_threshold: float, qhd_best_subgrad_L: None, n_calls: int = 50):
    init()
    results = {"qhd_best_L": None, "subgrad_best_L": None, "lfmsgd_best_L": None, "pure_qhd_best_L": None,
               "qhd_suc": None, "subgrad_suc": None, "lfmsgd_suc": None, "pure_qhd_suc": None, 
               "qhd_gap": None, "subgrad_gap": None, "lfmsgd_gap": None, "pure_qhd_gap": None, 
               "qhd_unrefined_suc": None, "qhd_unrefined_gap": None}
    qhd_missing_subgrad_L = False
    L_list = []
    mean_List = []
    
    if 1 in solvers:  # if qhd is included
        if qhd_best_subgrad_L is not None:
            qhd_min, qhd_max = L_bounds["qhd_L_bounds"]
            def f(L):
                qhd_result = run_qhd_subgrad(fpath, qhd_L=L[0], subgrad_L=qhd_best_subgrad_L)
                time_vs_mean = qhd_result["time_vs_mean"]
                _, final_mean = time_vs_mean[-1]
                return final_mean
            
            search_space = [Real(qhd_min, qhd_max, 'log-uniform')]
            optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
            best_L = optimize_result.x[0]
            results["qhd_best_L"] = round(best_L, 4)
            
            # run qhd_subgrad again with this best L
            qhd_result = run_qhd_subgrad(fpath, qhd_L=best_L, subgrad_L=qhd_best_subgrad_L)
            _, final_mean = qhd_result["time_vs_mean"][-1]
            results["qhd_gap"] = round(final_mean, 4)
            results["qhd_suc"] = round(sum(item["prob"] for item in qhd_result['output'] if item["value"] < suc_threshold), 4)
            
            # run qhd but without subgrad again with this best L
            qhd_unrefined_result = run_qhd(fpath, T=20, L=best_L)
            _, final_mean = qhd_unrefined_result["time_vs_mean"][-1]
            results["qhd_unrefined_gap"] = round(final_mean, 4)
            results["qhd_unrefined_suc"] = round(sum(item["prob"] for item in qhd_unrefined_result['output'] if item["value"] < suc_threshold), 4)
        else:
            qhd_missing_subgrad_L = True
    
    if 2 in solvers:  # if subgrad is included
        subgrad_min, subgrad_max = L_bounds["subgrad_L_bounds"]
        def f(L):
            subgrad_result = run_subgrad(fpath, tot_steps=20000, L=L[0])
            time_vs_mean = subgrad_result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            L_list.append(L[0])
            mean_List.append(final_mean)
            return final_mean
        
        search_space = [Real(subgrad_min, subgrad_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["subgrad_best_L"] = round(best_L, 4)
        
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
        results["lfmsgd_best_L"] = round(best_L, 4)
        
        # run lfmsgd again with this best L
        lfmsgd_result = run_lfmsgd(fpath, tot_steps=20000, L=best_L)
        _, final_mean = lfmsgd_result["time_vs_mean"][-1]
        results["lfmsgd_gap"] = round(final_mean, 4)
        results["lfmsgd_suc"] = round(sum(item["prob"] for item in lfmsgd_result['output'] if item["value"] < suc_threshold), 4)
    
    if 4 in solvers:  # if pure qhd (no subgrad involved) is included
        qhd_min, qhd_max = L_bounds["qhd_L_bounds"]
        def f(L):
            qhd_result = run_qhd(fpath, T=20, L=L[0])
            time_vs_mean = qhd_result["time_vs_mean"]
            _, final_mean = time_vs_mean[-1]
            return final_mean
        
        search_space = [Real(qhd_min, qhd_max, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        best_L = optimize_result.x[0]
        results["pure_qhd_best_L"] = round(best_L, 4)
        
        # run qhd again with this best L
        qhd_result = run_qhd(fpath, T=20, L=best_L)
        _, final_mean = qhd_result["time_vs_mean"][-1]
        results["pure_qhd_gap"] = round(final_mean, 4)
        results["pure_qhd_suc"] = round(sum(item["prob"] for item in qhd_result['output'] if item["value"] < suc_threshold), 4)

    
    # store outputs
    func_name = os.path.splitext(os.path.basename(fpath))[0]
    data = pd.read_csv(data_file_path)
    columns_to_convert = data.columns.difference(['Problem'])  # All columns except 'Problem'
    data[columns_to_convert] = data[columns_to_convert].astype(float)  # Convert these columns to float
    row_index = data.index[data['Problem'] == func_name].tolist()[0]
    if qhd_missing_subgrad_L: 
        print(f"QHD is specified but subgrad_best_L is missing, skip qhd_subgrad search.")
        
    if 1 in solvers and not qhd_missing_subgrad_L:
        data.at[row_index, 'qhd_L'] = results["qhd_best_L"]
        data.at[row_index, 'qhd_suc'] = results['qhd_suc']
        data.at[row_index, 'qhd_gap'] = results['qhd_gap']
        data.at[row_index, 'qhd_unrefined_suc'] = results['qhd_unrefined_suc']
        data.at[row_index, 'qhd_unrefined_gap'] = results['qhd_unrefined_gap']
        print(f"qhd_best_L: {results['qhd_best_L']}, qhd_gap: {results['qhd_gap']}, qhd_suc: {results['qhd_suc']}")
        print(f"qhd_unrefined_gap: {results['qhd_unrefined_gap']}, qhd_unrefined_suc: {results['qhd_unrefined_suc']}")

    if 2 in solvers:
        data.at[row_index, 'subgrad_L'] = results["subgrad_best_L"]
        data.at[row_index, 'subgrad_suc'] = results['subgrad_suc']
        data.at[row_index, 'subgrad_gap'] = results['subgrad_gap']
        print(f"subgrad_best_L: {results['subgrad_best_L']}, subgrad_gap: {results['subgrad_gap']}, subgrad_suc: {results['subgrad_suc']}")
        print(f"L_tested: {L_list}")
        print("---------------------------------------------------------------")
        print(f"mean_returned: {mean_List}")
        
    if 3 in solvers:
        data.at[row_index, 'lfmsgd_L'] = results["lfmsgd_best_L"]
        data.at[row_index, 'lfmsgd_suc'] = results['lfmsgd_suc']
        data.at[row_index, 'lfmsgd_gap'] = results['lfmsgd_gap']
        print(f"lfmsgd_best_L: {results['lfmsgd_best_L']}, lfmsgd_gap: {results['lfmsgd_gap']}, lfmsgd_suc: {results['lfmsgd_suc']}")
        
    if 4 in solvers:
        data.at[row_index, 'pure_qhd_L'] = results["pure_qhd_best_L"]
        data.at[row_index, 'pure_qhd_suc'] = results['pure_qhd_suc']
        data.at[row_index, 'pure_qhd_gap'] = results['pure_qhd_gap']
        print(f"pure_qhd_best_L: {results['pure_qhd_best_L']}, pure_qhd_gap: {results['pure_qhd_gap']}, pure_qhd_suc: {results['pure_qhd_suc']}")
        
    data.to_csv(data_file_path, index=False)
    return results

if __name__ == '__main__':
    fpath = 'func/nonsmooth/WF.cpp'
    solvers = [2]                                # 1: qhd, 2: subgrad, 3: lfmsgd, 4: pure qhd (no subgrad involved)
    qhd_best_subgrad_L = 12.2602                    # if qhd is included, the best L for subgrad 
    suc_threshold = 0.05
    qhd_L_bounds = 0.1, 100
    subgrad_L_bounds = 0.1, 100
    lfmsgd_L_bounds = 0.001, 10
    L_bounds = {"qhd_L_bounds": qhd_L_bounds, "subgrad_L_bounds": subgrad_L_bounds, "lfmsgd_L_bounds": lfmsgd_L_bounds}
    results = search_best_L(fpath, solvers, L_bounds, suc_threshold, qhd_best_subgrad_L=qhd_best_subgrad_L)