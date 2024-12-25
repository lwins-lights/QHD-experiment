from util import *
from skopt import gp_minimize
from skopt.space import Real
import numpy as np
from colorama import Fore, Back, Style, init
import matplotlib.pyplot as plt

def test_keane_ternary():
    def ternary_search_min(f, a, b, epsilon=1e-3):
        """
        Perform ternary search to find the minimum of a function f(x) over the interval [a, b].
        - f: The function to maximize
        - a, b: The range to search
        - epsilon: The precision for stopping the search (smaller epsilon = more precision)
        """
        while b - a > epsilon:
            print("current range: [%.8f,%.8f]" % (a, b))

            # Calculate two mid points
            c1 = a + (b - a) / 3
            c2 = b - (b - a) / 3
            
            # Evaluate the function at these points
            f_c1 = f(c1)
            f_c2 = f(c2)
            
            # Narrow down the search range based on comparison
            if f_c1 > f_c2:
                a = c1  # Move search range to [c1, b]
            else:
                b = c2  # Move search range to [a, c2]
        
        # Return the midpoint of the final range
        return (a + b) / 2
    def f(log_L):
        result = run_qhd_subgrad('func/nonsmooth/keane.cpp', qhd_L=np.exp(log_L), subgrad_L=10)
        time_vs_mean = result["time_vs_mean"]
        time, mean = time_vs_mean[-1]
        return mean
    log_L = ternary_search_min(f, -10, 10)
    print("best L = ", np.exp(log_L))
    # Result is roughly 1/e

def test_keane_gp():
    init()
    def f(L):
        result = run_qhd_subgrad('func/nonsmooth/keane.cpp', qhd_L=L[0], subgrad_L=10)
        time_vs_mean = result["time_vs_mean"]
        time, mean = time_vs_mean[-1]
        print(Fore.GREEN + ("f(%.8lf)=%.8lf" % (L[0], mean)) + Style.RESET_ALL)
        return mean
    search_space=[Real(1e-2, 1e2, 'log-uniform')]
    result = gp_minimize(f, search_space, n_calls=50, random_state=0)
    best_lr = result.x[0]
    print(f"Optimal Learning Rate: {best_lr}")
    # 0.3843057937223915

def test_keane_and_kde():
    result = run_qhd('func/nonsmooth/keane.cpp', L=0.4)
    output = result['output']

    threshold = 0.01
    probability = sum(item["prob"] for item in result['output'] if item["value"] < threshold)

    print(f"The probability of X < {threshold} is: {probability}")

    # Extract values and probabilities
    values = np.array([item["value"] for item in result['output']])
    weights = np.array([item["prob"] for item in result['output']])

    # Create a weighted KDE
    kde = gaussian_kde(values, weights=weights)

    # Generate x-axis points
    x = np.linspace(min(values), max(values), 500)
    y = kde(x)

    # Plot the weighted KDE
    plt.plot(x, y, label="Weighted KDE", color="blue")
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.title("Kernel Density Estimation (Weighted)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.show()

def test_lfmsgd_rand():
    result = run_lfmsgd('func/nonsmooth/WF.cpp', L=6, n_sample=2, tot_steps=2)
    time_vs_mean = result['time_vs_mean']
    time, mean = time_vs_mean[-1]
    print(f'mean = {mean}')

def test():
    result = run_qhd_subgrad('func/nonsmooth/WF.cpp', qhd_L=5.8089, subgrad_L=420)
    time_vs_mean = result['time_vs_mean']
    time, mean = time_vs_mean[-1]
    suc = round(sum(item["prob"] for item in result['output'] if item["value"] < 0.05), 4)
    print(f'mean = {mean}, suc = {suc}')

def test2():
    fpath = 'func/nonsmooth/keane.cpp'
    result1 = run_qhd(fpath, dt=0.001)
    result2 = run_subgrad(fpath, tot_steps=100)
    result3 = run_lfmsgd(fpath, tot_steps=100)
    _, mean1 = result1['time_vs_mean'][0]
    _, mean2 = result2['time_vs_mean'][0]
    _, mean3 = result3['time_vs_mean'][0]
    print(f'MEAN: {mean1}; {mean2}; {mean3}')

def test_print_dist():
    fpath = 'func/nonsmooth/keane.cpp'
    #result = run_lfmsgd(fpath, L=2.1228682698967813, tot_steps=20000, n_sample=1000)
    result = run_qhd_subgrad(fpath, qhd_L=0.267101753506348, subgrad_L=6.05144432518986)
    #result = run_subgrad(fpath, L=4.81342991759568, tot_steps=20000, n_sample=1000)
    x = [item["value"] for item in result['output']]
    plt.figure(figsize=(8, 5))
    plt.plot(x, marker='o', linestyle='-', label='Linear Scale')
    plt.yscale('log')  # Set Y-axis to logarithmic scale
    plt.xlabel('Index')
    plt.ylabel('Value (log scale)')
    plt.title('1D Float Array with Logarithmic Scale')
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.show()

def test3():
    fpath = 'func/nonsmooth/WF.cpp'
    result1 = run_qhd(fpath, resol=64, L=1.7210213630447821)
    result2 = run_qhd_subgrad(fpath, skip_qhd=True, subgrad_L=1787.2900596444126)
    _, mean1 = result1['time_vs_mean'][-1]
    _, mean2 = result2['time_vs_mean'][-1]
    print(f'MEAN: {mean1}; {mean2}')

def find_min():
    fpath = 'func/nonsmooth/keane.cpp'
    def f(vars):
        L = vars[0]
        result = run_subgrad(fpath, L=L, tot_steps=1000, n_sample=100000)
        output = result['output']
        mini = (min(output, key=lambda x: x["value"]))["value"]
        print(f"L: {L}; MIN: {mini}")
        return mini
    search_space=[Real(1, 1e6, 'log-uniform')]
    optimize_result = gp_minimize(f, search_space, n_calls=50, random_state=0)
    print(f"Final L: {optimize_result.x[0]}; MIN: {optimize_result.fun}")

if __name__ == '__main__':
    #test2()
    find_min()