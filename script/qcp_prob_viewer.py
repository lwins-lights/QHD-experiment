# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.widgets import Slider
from tqdm import tqdm
from itertools import product

def number_to_base(n, b):
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits

# get default dpi
dpi = mpl.rcParams['figure.dpi']

# load params from file
res = np.load('./result/qipopt.npz')
T = res['T'][0]
L = res['L'][0]
dt = res['dt'][0]
len = res['len'][0]
dim = res['dim'][0]
dt = res['dt'][0]
eta = res['eta'][0]
prob = res['probability']

# skeleton
fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.25)

def get_ly(prog, d, len, dim, prob_agg):
    size = len ** dim;
    prob = prob_agg[prog]
    ly = [0] * len
    for i in range(size):
        try:
            temp = number_to_base(i, len)[d]
        except:
            temp = 0
        ly[temp] += prob[i]
    return ly

# init: prepare all figures
ly = {}
rg = list(product(range(dim), range(101)))
for d, prog in tqdm(rg):
    ly[(d, prog)] = get_ly(prog, d, len, dim, prob)

# initial plot
lx = np.arange(-L, L, 2 * L / len)
line, = ax.plot(lx, ly[(0, 0)])
ax.set_ylim([0,1])

# slider definitions
ax_dim_slider = fig.add_axes([0.25, 0.15, 0.65, 0.03])
ax_prog_slider = fig.add_axes([0.25, 0.1, 0.65, 0.03])
dim_slider = Slider(
    ax_dim_slider, "Dim", 0, dim - 1,
    valinit=0, valstep=1
)
prog_slider = Slider(
    ax_prog_slider, "Prog", 0, 100,
    valinit=0, valstep=1
)

# slider triggers
def update(val):
    d = dim_slider.val
    prog = prog_slider.val
    line.set_ydata(ly[d, prog])
    fig.canvas.draw_idle()
dim_slider.on_changed(update)
prog_slider.on_changed(update)

# show
plt.show()