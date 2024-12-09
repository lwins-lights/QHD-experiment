from util import run_qhd, run_qhd_subgrad
from skopt import gp_minimize
from skopt.space import Real
import numpy as np
from colorama import Fore, Back, Style, init

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

if __name__ == '__main__':
    test_keane_gp()