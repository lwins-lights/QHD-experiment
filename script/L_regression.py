import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error


def main():
    def root_mean_squared_error(y_true, y_pred):
        return np.sqrt(mean_squared_error(y_true, y_pred))

    # Paths and data loading
    this_path = os.path.dirname(os.path.realpath(__file__))
    root_path = os.path.join(this_path, '..')
    result_path = os.path.join(root_path, "result")
    data_file_path = os.path.join(result_path, "L_lip_data.csv")
    output_path = os.path.join(result_path, "L_regression_exp_axis.pdf")

    data_avg = pd.read_csv(data_file_path).sort_values("Avg_Lipschitz")
    data_max = pd.read_csv(data_file_path).sort_values("Max_Lipschitz")
    avg_lipschitz = data_avg["Avg_Lipschitz"]
    max_lipschitz = data_max["Max_Lipschitz"]
    best_L_avg = data_avg["qhd_L"]
    best_L_max = data_max["qhd_L"]

    # Log transformation
    avg_lipschitz_log = np.log10(avg_lipschitz).values.reshape(-1, 1)
    max_lipschitz_log = np.log10(max_lipschitz).values.reshape(-1, 1)
    best_L_log_avg = np.log10(best_L_avg).values
    best_L_log_max = np.log10(best_L_max).values

    # Regression
    avg_model = LinearRegression()
    max_model = LinearRegression()
    avg_model.fit(avg_lipschitz_log, best_L_log_avg)
    max_model.fit(max_lipschitz_log, best_L_log_max)

    avg_predict_best_L = np.power(10, avg_model.predict(avg_lipschitz_log))
    max_predict_best_L = np.power(10, max_model.predict(max_lipschitz_log))

    print(f"avg_regression_error: {root_mean_squared_error(best_L_avg, avg_predict_best_L)}")
    print(f"max_regression_error: {root_mean_squared_error(best_L_max, max_predict_best_L)}")

    # Plotting
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 11  
    plt.rcParams['axes.titlesize'] = 16  
    plt.rcParams['axes.labelsize'] = 14  
    plt.rcParams['legend.fontsize'] = 10  

    _, ax = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={'wspace': 0.35})

    # Left plot: Avg Lipschitz
    ax[0].scatter(avg_lipschitz, best_L_avg, marker='o', c='r', label="Test Functions")
    ax[0].set_xlabel("Average Lipschitz")
    ax[0].set_ylabel("Best L")
    ax[0].set_title("Best L vs. Average Lipschitz Constant")
    ax[0].plot(avg_lipschitz, avg_predict_best_L, label="Log-Log Regression")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].legend(loc="lower right")

    # Right plot: Max Lipschitz
    ax[1].scatter(max_lipschitz, best_L_max, marker='o', c='r', label="Test Functions")
    ax[1].set_xlabel("Max Lipschitz")
    ax[1].set_ylabel("Best L")
    ax[1].set_title("Best L vs. Max Lipschitz Constant")
    ax[1].plot(max_lipschitz, max_predict_best_L, label="Log-Log Regression")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].legend(loc="lower right")

    # Save and show
    plt.savefig(output_path)
    plt.show()
    
if __name__ == '__main__':
    main()