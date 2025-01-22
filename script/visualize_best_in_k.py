from math import comb
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from seaborn import color_palette

# Path initialization
this_path = os.path.dirname(os.path.realpath(__file__))
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")
data_file_path = os.path.join(result_path, "data.csv")

# visualization setting
palette = color_palette("husl", 3)
[qhd_color, subgrad_color, lfmsgd_color] = palette

def plot_gaps_vs_k(csv_file):
    df = pd.read_csv(csv_file)
    df = df[df['name'] == 'Regular Test v2.1']
    # df = df.dropna(subset=['k', 'subgrad_gap', 'lfmsgd_gap', 'qhd_gap'])
    
    df['k'] = df['k'].astype(float)
    df['subgrad_gap'] = df['subgrad_gap'].astype(float)
    df['lfmsgd_gap'] = df['lfmsgd_gap'].astype(float)
    df['qhd_gap'] = df['qhd_gap'].astype(float)
    df['qhd_min'] = df['qhd_min'].astype(float)
    
    # Set up the plot grid
    num_instances = 12
    n_rows = 4
    n_cols = 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8.27, 11))
    # fig.suptitle('Gap vs k for NonConvex Instances', fontsize=16)
    
    for i, func in enumerate(df['func'].unique()[:num_instances]):  
        ax = axes[i // n_cols, i % n_cols]
        combined_data = df[df['func'] == func]
        theoretical_data = combined_data[combined_data['k'].isna()]
        test_data = combined_data[(combined_data['k'].notna())]
        test_data = test_data.sort_values('k')
        
        # Plot the gaps vs k in log-log scale
        ax.plot(test_data['k'], test_data['subgrad_gap'], label='Subgrad', color=subgrad_color, marker='o', markersize=4)
        ax.plot(test_data['k'], test_data['lfmsgd_gap'], label='LFMSGD', color=lfmsgd_color, marker='s', markersize=4)
        ax.plot(test_data['k'], test_data['qhd_gap'], label='QHD', color=qhd_color, marker='^', markersize=4)
        
        # Plot the theoretical qhd gap if its non-zero
        if theoretical_data.iloc[0]['qhd_min'] > 0:
            # ax.axhline(y=theoretical_data.iloc[0]['qhd_min'], color='red', linestyle='--', linewidth=1)
            handles, labels = ax.get_legend_handles_labels()
        
        # Set log scale for both axes
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        # Set title and labels
        ax.set_title(func.split('/')[-1][:-4].upper(), fontsize=10)  
        ax.set_xlabel('k', fontsize=8)
        ax.set_ylabel('Gap', fontsize=8)
    
    plt.tight_layout(rect=[0, 0, 0.95, 0.95])
    fig.legend(handles, labels, loc='upper right', fontsize=8)
    
    # Save the figure
    plt.savefig('gap_vs_k_regular_test_v2.1.pdf')
    plt.show()


if __name__ == '__main__':
    plot_gaps_vs_k(data_file_path)    