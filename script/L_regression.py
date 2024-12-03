import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from sklearn.metrics import root_mean_squared_error


def main():
    # path init
    this_path = os.path.dirname(os.path.realpath(__file__)) 
    root_path = os.path.join(this_path, '..')
    result_path = os.path.join(root_path, "result")
    data_file_path = os.path.join(result_path, "L_lip_data.csv")
    output_path = os.path.join(result_path, "L_regression.png")
    
    data_avg = pd.read_csv(data_file_path).sort_values("Avg_Lipschitz")
    data_max = pd.read_csv(data_file_path).sort_values("Max_Lipschitz")
    avg_lipschitz = data_avg["Avg_Lipschitz"]
    max_lipschitz = data_max["Max_Lipschitz"]
    best_L_avg = data_avg["qhd_L"]
    best_L_max = data_max["qhd_L"]
    avg_lipschitz_log = np.log(avg_lipschitz).values.reshape(-1, 1)
    max_lipschitz_log = np.log(max_lipschitz).values.reshape(-1, 1)
    best_L_log_avg = np.log(best_L_avg).values
    best_L_log_max = np.log(best_L_max).values
    
    avg_model = LinearRegression()
    max_model = LinearRegression()
    avg_model.fit(avg_lipschitz_log, best_L_log_avg)
    max_model.fit(max_lipschitz_log, best_L_log_max)
    
    avg_predict_best_L = np.exp(avg_model.predict(avg_lipschitz_log))
    max_predict_best_L = np.exp(max_model.predict(max_lipschitz_log))
    
    print(f"avg_regression_error: {root_mean_squared_error(best_L_avg, avg_predict_best_L)}")
    print(f"max_regression_error: {root_mean_squared_error(best_L_max, max_predict_best_L)}")

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 11  
    plt.rcParams['axes.titlesize'] = 16  
    plt.rcParams['axes.labelsize'] = 14  
    plt.rcParams['legend.fontsize'] = 10  
    
    _, ax = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw={'wspace': 0.35})
    ax[0].scatter(avg_lipschitz_log, best_L_log_avg, marker='o', c='r', label="Test Functions")
    # ax[0].scatter(avg_lipschitz, best_L_avg, marker='o', c='r')
    ax[0].set_xlabel("log(average Lipschitz)")
    ax[0].set_ylabel("log(best L)")
    # ax[0].set_title("best_L vs Average Lipschitz ")
    # ax[0].plot(avg_lipschitz, avg_predict_best_L)
    ax[0].set_title("Best L vs Average Lipschitz Constant")
    ax[0].plot(avg_lipschitz_log, avg_model.predict(avg_lipschitz_log), label="Log-Log Regression")
    ax[0].legend(loc="lower right")
    
    ax[1].scatter(max_lipschitz_log, best_L_log_max, marker='o', c='r', label="Test Functions")
    # ax[1].scatter(max_lipschitz, best_L_max, marker='o', c='r')
    ax[1].set_xlabel("log(max Lipschitz)")
    ax[1].set_ylabel("log(best L)")
    # ax[1].set_title("best_L vs Max Lipschitz ")
    # ax[1].plot(max_lipschitz, max_predict_best_L)
    ax[1].set_title("Best L vs Max Lipschitz Constant")
    ax[1].plot(max_lipschitz_log, max_model.predict(max_lipschitz_log), label="Log-Log Regression")
    ax[1].legend(loc="lower right")
    
    plt.savefig(output_path)
    plt.show()
    
if __name__ == '__main__':
    main()