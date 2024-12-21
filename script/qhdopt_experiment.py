import numpy as np 
from sympy import symbols, sin, exp, evalf, Sum, pi, sqrt, cos, tan, log, Abs
import pandas as pd
import sys
import os
local_qhdopt_path = "/Users/zhiyuanjia/Desktop/QHDOPTResearch/QHDOPT"  
sys.path.insert(0, local_qhdopt_path)
from qhdopt.qhd import QHD

this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "real_machine_results.csv")
simulator_path = os.path.join(root_path, "simulator")

def get_total_time(res):
    if res.info["backend_time"] == 0:
        return res.info['refining_time']
    else:
        total_runtime = res.info["time_on_machine"]
        total_runtime += (res.info['refining_time'] if res.info["refine_status"] else 0)
        return total_runtime

def store_results(results, file_name):
    df = pd.DataFrame(results)
    format_sci = lambda x: f"{x:.4e}"
    
    df[['quantum_Rtime', 'quantum_min']] = \
        df[['quantum_Rtime', 'quantum_min']].map(format_sci)
    
    file_exists = os.path.isfile(file_name)
    df.to_csv(file_name, mode='a', index=False, header=not file_exists)
    
def run_qhdopt(name: str, shot: int, qhd_resolution: int, x, f, bound):
    key = "DEV-7aa39ca5c55f857048c112d91c8b819ce75b525f"
    model = QHD.SymPy(f, x, bounds=bound)
    model.dwave_setup(resolution = qhd_resolution, api_key = key, shots = shot)
        
    response = model.optimize(verbose = 1, refine = False)
    quantum_Rtime = get_total_time(response)
    quantum_min = response.coarse_minimum
    quantum_minimizer = response.coarse_minimizer
    
    result = {
        'problem': name,
        'qhd_resolution': qhd_resolution,
        'dimension': len(quantum_minimizer),
        'quantum_Rtime': quantum_Rtime,
        'quantum_min': quantum_min,
        'quantum_minimizer': quantum_minimizer
    }
    
    print(f"coarse result: {result}")
    return result

def main():
    name = "Ackley_Modified"
    r = 8
    shot = 100

    n = 10
    x = symbols(f'x0:{n+1}') 
    f = 0 
    for i in range(1, n):
        f += (-20)*exp(-0.2*Abs(x[i]) / 2)*exp(-0.2*Abs(x[i+1]) / 2) - exp(cos(2*pi*x[i]) / 2) * exp(cos(2*pi*x[i+1]) / 2)
    bound = (0, 14)
    
    result = run_qhdopt(name, shot, r, x[1:], f, bound)
    
if __name__ == '__main__':
    main()
