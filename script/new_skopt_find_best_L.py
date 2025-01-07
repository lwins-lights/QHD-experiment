from skopt import gp_minimize
from skopt.space import Real
import numpy as np
import os
import math
import csv
from util import *

# fix randomness
np.random.seed(0)

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "data.csv")

def lower_precision(value, relative_precision=1e-6):
    if value == 0:
        return 0
    # Determine the order of magnitude
    magnitude = 10 ** math.floor(math.log10(abs(value)))
    # Round to the desired relative precision
    return round(value / magnitude / relative_precision) * magnitude * relative_precision

def best_in_k(a, k):
    # a is an array of {"value":value, "prob":prob}
    prob_sum = 0
    ret = 0
    for item in reversed(a):
        ret += (pow(prob_sum + item['prob'], k) - pow(prob_sum, k)) * item['value']
        prob_sum += item['prob']
    return ret

def save_dict_to_csv_row(dict_data, csv_file_path):
    """
    Save a dictionary to a row of a CSV file, handling dynamic keys.

    Parameters:
    - dict_data (dict): The dictionary to save.
    - csv_file_path (str): Path to the CSV file.
    """
    # Check if the file already exists
    file_exists = os.path.isfile(csv_file_path)

    # If the file exists, read the existing header and data
    if file_exists:
        with open(csv_file_path, mode='r', encoding='utf-8') as csv_file:
            reader = csv.DictReader(csv_file)
            existing_fieldnames = reader.fieldnames if reader.fieldnames else []
            existing_data = list(reader)
    else:
        existing_fieldnames = []
        existing_data = []

    # Update the fieldnames to include new keys
    new_fieldnames = list(set(existing_fieldnames).union(dict_data.keys()))

    # Reorder the data to align with the new fieldnames
    updated_data = []
    for row in existing_data:
        updated_row = {key: row.get(key, '') for key in new_fieldnames}
        updated_data.append(updated_row)

    # Add the new dictionary as a row
    updated_data.append({key: dict_data.get(key, '') for key in new_fieldnames})

    # Write the updated data back to the CSV file
    with open(csv_file_path, mode='w', newline='', encoding='utf-8') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=new_fieldnames)

        # Write the new header and data
        writer.writeheader()
        writer.writerows(updated_data)

def optimize_L(
    fpath,
    optimize_qhd=True, optimize_subgrad=True, optimize_lfmsgd=True, optimize_mini=False,
    qhd_L_bounds=(0.01, 100), subgrad_L_bounds=(0.001, 1000), lfmsgd_L_bounds=(0.001, 1000), mini_L_bounds=(1, 1e6),
    qhd_resol=256, subgrad_n_samples=10000, lfmsgd_n_samples=10000, mini_n_samples=100000,
    qhd_T=10, qhd_dt=0.001, subgrad_tot_steps=10000, lfmsgd_tot_steps=10000, mini_steps=1000,
    k=1, n_calls=50,
    qhd_L=0
):

    def use_gp_minimize(f, L_bounds, L_name, gap_name, results):
        l, u = L_bounds
        search_space=[Real(l, u, 'log-uniform')]
        optimize_result = gp_minimize(f, search_space, n_calls=n_calls, random_state=0)
        results[L_name] = optimize_result.x[0]
        results[gap_name] = optimize_result.fun

    def direct_optimize(f, L, L_name, gap_name, results):
        results[L_name] = L
        results[gap_name] = f([L])

    results={"func":fpath, "k":k}

    if optimize_mini:
        def f(vars):
            L = vars[0]
            result = run_subgrad(fpath, L=L, tot_steps=mini_steps, n_sample=mini_n_samples)
            output = result['output']
            return (min(output, key=lambda x: x["value"]))["value"]
        use_gp_minimize(f, mini_L_bounds, "mini_L", "mini_gap", results)

    if optimize_qhd:
        def f(vars):
            L = vars[0]
            result = run_qhd(fpath, L=L, resol=qhd_resol, T=qhd_T, dt=qhd_dt)
            output = result['output']
            return lower_precision(best_in_k(output, k))
        if qhd_L == 0:
            use_gp_minimize(f, qhd_L_bounds, "qhd_L", "qhd_gap", results)
        else:
            direct_optimize(f, qhd_L, "qhd_L", "qhd_gap", results)

    if optimize_subgrad:
        def f(vars):
            L = vars[0]
            result = run_subgrad(fpath, L=L, tot_steps=subgrad_tot_steps, n_sample=subgrad_n_samples)
            output = result['output']
            return best_in_k(output, k)
        use_gp_minimize(f, subgrad_L_bounds, "subgrad_L", "subgrad_gap", results)

    if optimize_lfmsgd:
        def f(vars):
            L = vars[0]
            result = run_lfmsgd(fpath, L=L, tot_steps=subgrad_tot_steps, n_sample=subgrad_n_samples)
            output = result['output']
            return best_in_k(output, k)
        use_gp_minimize(f, lfmsgd_L_bounds, "lfmsgd_L", "lfmsgd_gap", results)

    return results

if __name__ == '__main__':
    np.random.seed(0)
    func_dim_list = [
        ('func/nonsmooth/Schwefel.cpp', 1, 1.47250929702301),
        ('func/nonsmooth/keane.cpp', 2, 0.251405304286405),
        ('func/nonsmooth/WF.cpp', 2, 5.59148287989168),
        ('func/nonsmooth/CrownedCross.cpp', 2, 1.20100312295498),
        ('func/nonsmooth/bukin06.cpp', 2, 7.35173280432231),
        ('func/nonsmooth/Ackley.cpp', 2, 1.80021533922279),
        ('func/nonsmooth/xinsheyang04.cpp', 2, 0.526387551268728),
        ('func/nonsmooth/CarromTable.cpp', 2, 0.815342466034414),
        ('func/nonsmooth/rana.cpp', 2, 4.25752814528056),
        ('func/nonsmooth/Damavandi.cpp', 2, 4.55998485143467),
        ('func/nonsmooth/DropWave_3d.cpp', 3, 0.316167168089332),
        ('func/nonsmooth/layeb04_3d.cpp', 3, 1.69324811265117)
    ]
    
    resol = [0, pow(int(2),18), pow(int(2),9), pow(int(2),6)]
    for fpath, dim, L in func_dim_list:
        results = optimize_L(fpath, qhd_resol=resol[dim], n_calls=100, k=30)
        results["name"] = "Regular Test v2.0"
        save_dict_to_csv_row(results, data_file_path)

    #for fpath in func_list:
    #    results = optimize_L(fpath, optimize_qhd=False, optimize_subgrad=False, optimize_lfmsgd=False, optimize_mini=True)
    #    results["name"] = "Minimum Benchmark"
    #    save_dict_to_csv_row(results, data_file_path)

    #fpath = 'func/nonsmooth/keane.cpp'
    #results = optimize_L(fpath, qhd_resol=512, subgrad_n_samples=10000, lfmsgd_n_samples=10000, n_calls=50, k=10)
    #results = optimize_L(fpath, optimize_qhd=False, optimize_subgrad=False, optimize_lfmsgd=False, optimize_mini=True)
    #results["name"] = "min benchmark"
    #print(results)
    #save_dict_to_csv_row(results, data_file_path)
