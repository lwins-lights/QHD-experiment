import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from tqdm import tqdm 
from util import *

# constant from /simulator/config.hpp
n_snapshot = 100

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")

def generate_figures(output_folder=None, snapshot_list=range(n_snapshot), n_query=None):

    print("Initializing...")

    # Data preparation for QHD
    qhd_res = np.load('./result/pseudospec.npz')
    len_qhd = qhd_res['len'][0]
    V = qhd_res['V']
    L_qhd = qhd_res['L']
    nx = ny = len_qhd
    n_snapshot = 100
    dist_qhd = qhd_res['dist_snapshot'].reshape(n_snapshot, nx, ny)
    x_qhd = np.linspace(-L_qhd, L_qhd, nx, endpoint=False)
    y_qhd = np.linspace(-L_qhd, L_qhd, ny, endpoint=False)
    x_qhd, y_qhd = np.meshgrid(x_qhd, y_qhd)

    # Data preparation for SUBGRAD
    subgrad_res = np.load('./result/subgrad.npz')
    L_subgrad = subgrad_res['L']
    sample_number = subgrad_res['sample_number'][0]
    n_snapshot = 100
    dist_subgrad = subgrad_res['dist_snapshot'].reshape(n_snapshot, sample_number, 2)
    grid_size = 256
    x_subgrad = np.linspace(-L_subgrad, L_subgrad, grid_size)
    y_subgrad = np.linspace(-L_subgrad, L_subgrad, grid_size)
    x_subgrad, y_subgrad = np.meshgrid(x_subgrad, y_subgrad)
    positions_subgrad = np.vstack([x_subgrad.ravel(), y_subgrad.ravel()])

    # Data preparation for LFMSGD
    lfmsgd_res = np.load('./result/lfmsgd.npz')
    L_lfmsgd = lfmsgd_res['L']
    sample_number_lfmsgd = lfmsgd_res['sample_number'][0]
    n_snapshot = 100
    dist_lfmsgd = lfmsgd_res['dist_snapshot'].reshape(n_snapshot, sample_number_lfmsgd, 2)
    x_lfmsgd = np.linspace(-L_lfmsgd, L_lfmsgd, grid_size)
    y_lfmsgd = np.linspace(-L_lfmsgd, L_lfmsgd, grid_size)
    x_lfmsgd, y_lfmsgd = np.meshgrid(x_lfmsgd, y_lfmsgd)
    positions_lfmsgd = np.vstack([x_lfmsgd.ravel(), y_lfmsgd.ravel()])

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fontsize=25

    print("Generating QHD figures...")
    for k in tqdm(snapshot_list, desc="QHD"):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        z_qhd = dist_qhd[k]
        ax.plot_surface(x_qhd, y_qhd, z_qhd, cmap='viridis', cstride=1, rstride=1, antialiased=False)
        ax.set_xlabel("$x_1$", fontsize=fontsize)
        ax.set_ylabel("$x_2$", fontsize=fontsize)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        plt.savefig(f"{output_folder}/QHD_{int(n_query * (k + 1) / n_snapshot)}.png")
        plt.close()

    print("Generating SUBGRAD figures...")
    for k in tqdm(snapshot_list, desc="SUBGRAD"):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        points_subgrad = dist_subgrad[k]
        kde_subgrad = gaussian_kde(points_subgrad.T)
        z_subgrad = kde_subgrad(positions_subgrad).reshape(grid_size, grid_size)
        ax.plot_surface(x_subgrad, y_subgrad, z_subgrad, cmap='viridis', cstride=1, rstride=1, antialiased=False)
        ax.set_xlabel("$x_1$", fontsize=fontsize)
        ax.set_ylabel("$x_2$", fontsize=fontsize)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        plt.savefig(f"{output_folder}/SUBGRAD_{int(n_query * (k + 1) / n_snapshot)}.png")
        plt.close()

    print("Generating LFMSGD figures...")
    for k in tqdm(snapshot_list, desc="LFMSGD"):
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='3d')
        points_lfmsgd = dist_lfmsgd[k]
        kde_lfmsgd = gaussian_kde(points_lfmsgd.T)
        z_lfmsgd = kde_lfmsgd(positions_lfmsgd).reshape(grid_size, grid_size)
        ax.plot_surface(x_lfmsgd, y_lfmsgd, z_lfmsgd, cmap='viridis', cstride=1, rstride=1, antialiased=False)
        ax.set_xlabel("$x_1$", fontsize=fontsize)
        ax.set_ylabel("$x_2$", fontsize=fontsize)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        plt.savefig(f"{output_folder}/LFMSGD_{int(n_query * (k + 1) / n_snapshot)}.png")
        plt.close()

    print(f"Figures saved to {output_folder}")

if __name__ == '__main__':

    fpath = 'func/nonsmooth/xinsheyang04.cpp'
    qhd_L = 0.526387551268728
    qhd_resol = 512
    qhd_dt = 0.001
    subgrad_L = 10.4820053914904
    lfmsgd_L = 4.57975627971285
    snapshot_list = [0, 9, 24, 99]

    for n_query in [1000, 10000]:
        run_qhd(fpath, L=qhd_L, resol=qhd_resol, T=qhd_dt*n_query, dt=qhd_dt)
        run_subgrad(fpath, L=subgrad_L, tot_steps=n_query)
        run_lfmsgd(fpath, L=lfmsgd_L, tot_steps=n_query)
        generate_figures(output_folder=os.path.join(result_path, f"visualize_dist_2d"), snapshot_list=snapshot_list, n_query=n_query)
