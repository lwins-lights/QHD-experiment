import os
import numpy as np
import filecmp
import shutil
from subprocess import run

# constant
verbose = False

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
simulator_path = os.path.join(root_path, "simulator")
potential_path = os.path.join(simulator_path, "potential.cpp")

def vprint(*params):
    if verbose:
        print(*params)

def copy_if_different(file_a, file_b):
    try:
        # Compare the files
        if not filecmp.cmp(file_a, file_b, shallow=False):  # Perform a deep comparison
            shutil.copyfile(file_a, file_b)  # Copy A to B
            vprint(f"File '{file_a}' has been copied to '{file_b}'.")
        else:
            vprint(f"Files '{file_a}' and '{file_b}' are identical. No copy needed.")
    except FileNotFoundError as e:
        vprint(f"Error: {e}")

def run_qhd(func_path, resol=256, T=10, dt=0.001, par=1, L=1):

    # compile QHD
    run(["mkdir", "-p", result_path])
    copy_if_different(func_path, potential_path)
    run(["make", "pseudospec"], cwd=simulator_path)

    # run the simulator
    arglist = ["./pseudospec", str(int(resol)), str(T), str(dt), str(int(par)), str(L)]
    run(arglist, cwd=simulator_path)

    # load the result
    npz = np.load(os.path.join(result_path, "pseudospec.npz"))
    result = {}

    # compile time-vs-mean data
    mean = npz['expected_potential']
    result['time_vs_mean'] = [(i * dt, mean[i]) for i in range(len(mean))]

    # compile distribution data
    V = npz['V']
    dist = npz['dist']
    result['output'] = [{"value": V[i], "prob": dist[i]} for i in range(len(dist))]
    result['output'].sort(key=lambda x: x["value"])

    return result

def run_subgrad(func_path, tot_steps=10000, lr_init=1, par=1, n_sample=10000, L=1):

    # compile SUBGRAD
    run(["mkdir", "-p", result_path])
    copy_if_different(func_path, potential_path)
    run(["make", "subgrad"], cwd=simulator_path)

    # run the simulator
    arglist = ["./subgrad", str(int(tot_steps)), str(lr_init), str(int(par)), str(int(n_sample)), str(L)]
    run(arglist, cwd=simulator_path)

    # load the result
    npz = np.load(os.path.join(result_path, "subgrad.npz"))
    result = {}

    # compile time-vs-mean data
    mean = npz['expected_potential']
    result['time_vs_mean'] = [(i, mean[i]) for i in range(len(mean))]

    # compile distribution data
    samples = npz['samples']
    result['output'] = [{"value": samples[i], "prob": 1.0 / len(samples)} for i in range(len(samples))]
    result['output'].sort(key=lambda x: x["value"])

    return result

def run_lfmsgd(func_path, tot_steps=10000, noise_level=1, par=1, n_sample=10000, L=1):

    # compile SUBGRAD
    run(["mkdir", "-p", result_path])
    copy_if_different(func_path, potential_path)
    run(["make", "lfmsgd"], cwd=simulator_path)

    # run the simulator
    arglist = ["./lfmsgd", str(int(tot_steps)), str(noise_level), str(int(par)), str(int(n_sample)), str(L)]
    run(arglist, cwd=simulator_path)

    # load the result
    npz = np.load(os.path.join(result_path, "lfmsgd.npz"))
    result = {}

    # compile time-vs-mean data
    mean = npz['expected_potential']
    result['time_vs_mean'] = [(i, mean[i]) for i in range(len(mean))]

    # compile distribution data
    samples = npz['samples']
    result['output'] = [{"value": samples[i], "prob": 1.0 / len(samples)} for i in range(len(samples))]
    result['output'].sort(key=lambda x: x["value"])

    return result

def run_qhd_subgrad(
    func_path, 
    resol=256, T=10, dt=0.001, par=1, qhd_L=1, 
    tot_steps=10000, lr_init=1, n_sample=10000, subgrad_L=1
):

    # run QHD first
    run_qhd(func_path, resol=resol, T=T, dt=dt, par=par, L=qhd_L)

    # compile SUBGRAD
    run(["make", "subgrad"], cwd=simulator_path)

    # run the simulator
    arglist = ["./subgrad", str(int(tot_steps)), str(lr_init), str(int(par)), str(int(n_sample)), str(subgrad_L), "1"]
    run(arglist, cwd=simulator_path)

    # load the result
    npz = np.load(os.path.join(result_path, "subgrad.npz"))
    result = {}

    # compile time-vs-mean data
    mean = npz['expected_potential']
    result['time_vs_mean'] = [(i, mean[i]) for i in range(len(mean))]

    # compile distribution data
    samples = npz['samples']
    result['output'] = [{"value": samples[i], "prob": 1.0 / len(samples)} for i in range(len(samples))]
    result['output'].sort(key=lambda x: x["value"])

    return result