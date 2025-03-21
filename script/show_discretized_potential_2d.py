import matplotlib.pyplot as plt
import numpy as np
from util import *

# path init
this_path = os.path.dirname(os.path.realpath(__file__)) 
root_path = os.path.join(this_path, '..')
result_path = os.path.join(root_path, "result")

def generate_figure(figpath=os.path.join(result_path, "landscape.png")):

    # Load data
    qhd_res = np.load('./result/pseudospec.npz')
    grid_size = qhd_res['len'][0]  # Avoid using `len` as a variable name
    L = qhd_res['L']
    V = qhd_res['V']

    # Create coordinate grid
    x = np.linspace(-L, L, grid_size)
    y = np.linspace(-L, L, grid_size)
    surface_x, surface_y = np.meshgrid(x, y)
    surface_z = V.reshape((grid_size, grid_size))  # Ensure correct reshaping

    # Find the global minimum
    min_idx = np.unravel_index(np.argmin(surface_z), surface_z.shape)
    min_x, min_y = x[min_idx[1]], y[min_idx[0]]  # Ensure correct coordinate selection

    # Create figure with two subplots: 3D surface and 2D heatmap
    fig, axs = plt.subplots(1, 2, figsize=(11, 5))
    fig.subplots_adjust(wspace=0.2)  # Adjust spacing

    # Remove unnecessary axis ticks and labels
    axs[0].axis('off')

    # 3D Surface Plot
    ax1 = fig.add_subplot(121, projection='3d')
    surf = ax1.plot_surface(surface_x, surface_y, surface_z, cmap='jet', rstride=1, cstride=1)
    #ax1.set_title('3D Surface Plot')
    ax1.set_xlabel('$x_1$')
    ax1.set_ylabel('$x_2$')
    ax1.set_zticks([])
    #fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)  # Add color bar

    ax1.view_init(elev=15, azim=-60)

    # Reduce ticks for x_1 and x_2 axes
    ax1.set_xticks(np.linspace(x.min(), x.max(), 5))  # Reduce x_1 ticks (e.g., 5 ticks)
    ax1.set_yticks(np.linspace(y.min(), y.max(), 5))  # Reduce x_2 ticks (e.g., 5 ticks)

    # 2D Heatmap (Square aspect ratio)
    ax2 = axs[1]
    heatmap = ax2.imshow(surface_z, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', cmap='jet', aspect='equal')
    #ax2.set_title('2D Heatmap')
    ax2.set_xlabel('$x_1$')
    ax2.set_ylabel('$x_2$')

    fig.colorbar(heatmap, ax=ax2, shrink=0.5, aspect=10)

    # Mark global minimum
    ax2.scatter(min_x, min_y, color='white', marker='x', s=100, label='Global Minimum')
    legend = ax2.legend()
    legend.get_frame().set_facecolor('gray')

    plt.savefig(figpath)

if __name__ == '__main__':
    fpath = 'func/nonsmooth/xinsheyang04_nobarrier.cpp'
    run_qhd(fpath, L=10, resol=512, T=1, dt=1)
    generate_figure()
