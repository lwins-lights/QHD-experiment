import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde

def show_side_by_side():
    # Data preparation for `show_qhd`
    qhd_res = np.load('./result/pseudospec.npz')
    len_qhd = qhd_res['len'][0]
    V = qhd_res['V']
    L_qhd = qhd_res['L']
    nx = ny = len_qhd
    k_total_qhd = 100
    dist_qhd = qhd_res['dist_snapshot'].reshape(k_total_qhd, nx, ny)
    x_qhd = np.linspace(-L_qhd, L_qhd, nx, endpoint=False)
    y_qhd = np.linspace(-L_qhd, L_qhd, ny, endpoint=False)
    x_qhd, y_qhd = np.meshgrid(x_qhd, y_qhd)

    # Data preparation for `show_sgd`
    subgrad_res = np.load('./result/subgrad.npz')
    L_sgd = subgrad_res['L']
    sample_number = subgrad_res['sample_number'][0]
    k_total_sgd = 100
    dist_sgd = subgrad_res['dist_snapshot'].reshape(k_total_sgd, sample_number, 2)
    grid_size = 100
    x_sgd = np.linspace(-L_sgd, L_sgd, grid_size)
    y_sgd = np.linspace(-L_sgd, L_sgd, grid_size)
    x_sgd, y_sgd = np.meshgrid(x_sgd, y_sgd)
    positions_sgd = np.vstack([x_sgd.ravel(), y_sgd.ravel()])

    # Set up the figure and subplots
    fig = plt.figure(figsize=(12, 6))
    ax_qhd = fig.add_subplot(121, projection='3d')
    ax_sgd = fig.add_subplot(122, projection='3d')
    plt.subplots_adjust(bottom=0.2, wspace=0.4)  # Space for sliders and between plots

    # Initial plot for `show_qhd`
    k_index_qhd = 0
    z_qhd = dist_qhd[k_index_qhd]
    surface_qhd = ax_qhd.plot_surface(x_qhd, y_qhd, z_qhd, cmap='viridis')
    ax_qhd.set_zlim(0, z_qhd.max())
    ax_qhd.set_title(f"QHD Distribution for k={k_index_qhd}")
    ax_qhd.set_xlabel("X-axis")
    ax_qhd.set_ylabel("Y-axis")
    ax_qhd.set_zlabel("Z-axis")

    # Initial plot for `show_sgd`
    k_index_sgd = 0
    points_sgd = dist_sgd[k_index_sgd]
    kde_sgd = gaussian_kde(points_sgd.T)
    z_sgd = kde_sgd(positions_sgd).reshape(grid_size, grid_size)
    surface_sgd = ax_sgd.plot_surface(x_sgd, y_sgd, z_sgd, cmap='viridis')
    ax_sgd.set_zlim(0, z_sgd.max() * 1.1)
    ax_sgd.set_title(f"SGD Distribution for k={k_index_sgd}")
    ax_sgd.set_xlabel("X-axis")
    ax_sgd.set_ylabel("Y-axis")
    ax_sgd.set_zlabel("Density")

    # Sliders
    ax_slider_qhd = plt.axes([0.15, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
    slider_qhd = Slider(ax_slider_qhd, 'QHD k', 0, k_total_qhd - 1, valinit=k_index_qhd, valstep=1)

    ax_slider_sgd = plt.axes([0.55, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
    slider_sgd = Slider(ax_slider_sgd, 'SGD k', 0, k_total_sgd - 1, valinit=k_index_sgd, valstep=1)

    # Update functions
    def update_qhd(val):
        nonlocal surface_qhd
        k_index_qhd = int(slider_qhd.val)
        z_qhd = dist_qhd[k_index_qhd]
        for collection in ax_qhd.collections:
            collection.remove()
        surface_qhd = ax_qhd.plot_surface(x_qhd, y_qhd, z_qhd, cmap='viridis')
        ax_qhd.set_zlim(0, z_qhd.max())
        ax_qhd.set_title(f"QHD Distribution for k={k_index_qhd}")
        fig.canvas.draw_idle()

    def update_sgd(val):
        nonlocal surface_sgd
        k_index_sgd = int(slider_sgd.val)
        points_sgd = dist_sgd[k_index_sgd] + np.random.normal(0, 1e-2, dist_sgd[k_index_sgd].shape)
        kde_sgd = gaussian_kde(points_sgd.T)
        z_sgd = kde_sgd(positions_sgd).reshape(grid_size, grid_size)
        for collection in ax_sgd.collections:
            collection.remove()
        surface_sgd = ax_sgd.plot_surface(x_sgd, y_sgd, z_sgd, cmap='viridis')
        ax_sgd.set_zlim(0, z_sgd.max() * 1.1)
        ax_sgd.set_title(f"SGD Distribution for k={k_index_sgd}")
        fig.canvas.draw_idle()

    slider_qhd.on_changed(update_qhd)
    slider_sgd.on_changed(update_sgd)

    plt.show()

if __name__ == '__main__':
    show_side_by_side()
