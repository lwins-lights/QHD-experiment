import dis
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

PROBLEM = "Abs"  # Specify the problem name to filter data
DT = 0.0001

def main():
    def root_mean_squared_error(y_true, y_pred):
        return np.sqrt(mean_squared_error(y_true, y_pred))

    # Paths and data loading
    this_path = os.path.dirname(os.path.realpath(__file__))
    root_path = os.path.join(this_path, '..')
    result_path = os.path.join(root_path, "result")
    data_file_path = os.path.join(result_path, "discretization_error.csv")
    output_path = os.path.join(result_path, "discretization_error_regression.pdf")

    # Load data and filter by PROBLEM
    data = pd.read_csv(data_file_path)
    data = data[(data["Problem"] == PROBLEM) & (data["dt"] == DT)].sort_values("discretization_num")

    # Extract filtered data
    discretization_num = data["discretization_num"].values
    discretization_err = data["discretization_error"].values

    # Log transformation
    discretization_num_log = np.log10(discretization_num).reshape(-1, 1)
    discretization_err_log = np.log10(discretization_err)

    # Regression
    model = LinearRegression()
    model.fit(discretization_num_log, discretization_err_log)
    slope = model.coef_[0]
    intercept = model.intercept_
    predict_error = np.power(10, model.predict(discretization_num_log))
    regression_error = root_mean_squared_error(discretization_err, predict_error)

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

    # Plot: Discretization Error vs Discretization Number
    ax.scatter(discretization_num, discretization_err, marker='o', c='r', label="Tested Discretization Numbers")
    ax.set_xlabel("Discretization Number")
    ax.set_ylabel("Discretization Error")
    ax.set_title("Discretization Error vs. Discretization Number")
    ax.plot(discretization_num, predict_error, label="Log-Log Regression", color='blue')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="upper right")
    
    # Save and show
    plt.savefig(output_path)
    plt.show()
    
if __name__ == '__main__':
    main()