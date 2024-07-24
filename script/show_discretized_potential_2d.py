# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse

# get default dpi
dpi = mpl.rcParams['figure.dpi']

# load data
qhd_res = np.load('./result/pseudospec.npz')

len = qhd_res['len'][0]
V = qhd_res['V']
L = qhd_res['L']
num = len * len
surface_x = np.empty(shape=(len,len))
surface_y = np.empty(shape=(len,len))
surface_z = np.empty(shape=(len,len))

for i0 in range(len):
    for i1 in range(len):
        l = i0 + i1 * len
        surface_x[i0, i1] = float(i0) / len * 2 * L - L
        surface_y[i0, i1] = float(i1) / len * 2 * L - L
        surface_z[i0, i1] = V[l]

# subplots
#f, (s_plt) = plt.subplots(1, 1, figsize=(800/dpi, 900/dpi))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(surface_x, surface_y, surface_z, cmap='viridis')
 
ax.set_title('3D Surface plot')
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('obj')
plt.show()

'''
y_qhd = qhd_res['expected_potential']
y_subgrad = subgrad_res['expected_potential']
y_lfmsgd = lfmsgd_res['expected_potential']

yp_qhd = qhd_res['probability_at_minimum']
yp_subgrad = subgrad_res['probability_at_minimum']
yp_lfmsgd = lfmsgd_res['probability_at_minimum']

T = qhd_res['T'][0]
dt = qhd_res['dt'][0]
len = qhd_res['len'][0]
dim = qhd_res['dim'][0]
L = qhd_res['L'][0]
tot_steps = subgrad_res['tot_steps'][0]
learning_rate = subgrad_res['learning_rate'][0]
sample_number = subgrad_res['sample_number'][0]
noise_level = lfmsgd_res['noise_level'][0]

# init x
x_qhd = np.arange(0, T, dt)
x_subgrad = np.arange(0, T, T / tot_steps)

# subplots
f, (e_plt, p_plt) = plt.subplots(2, 1, figsize=(800/dpi, 900/dpi))
  
# plotting for expected potential
e_plt.plot(x_qhd, y_qhd, label="QHD*" + str(qhd_res['par'][0]))
e_plt.plot(x_subgrad, y_subgrad, label="SUBGRAD*" + str(subgrad_res['par'][0]))
e_plt.plot(x_subgrad, y_lfmsgd, label="LFMSGD*" + str(lfmsgd_res['par'][0]))
  
e_plt.set_xlabel('time')
e_plt.set_ylabel('expectation')

# plotting for probability at minimum
p_plt.plot(x_qhd, yp_qhd, label="QHD*" + str(qhd_res['par'][0]))
p_plt.plot(x_subgrad, yp_subgrad, label="SUBGRAD*" + str(subgrad_res['par'][0]))
p_plt.plot(x_subgrad, yp_lfmsgd, label="LFMSGD*" + str(subgrad_res['par'][0]))
  
p_plt.set_xlabel('time')
p_plt.set_ylabel('success probability')
  
# style
e_plt.yaxis.tick_right()
p_plt.yaxis.tick_right()
e_plt.legend(loc="upper right")
p_plt.legend(loc="lower right")

# show all
qhd_params = 'QHD params: dt=%f; gran.=%f^%d; L=%f' % (dt, len, dim, L)
subgrad_params = '\nSUBGRAD params: tot_steps=%d; learning_rate=%f; sample_number=%d' % (tot_steps, learning_rate, sample_number)
lfmsgd_params = '\nLFMSGD params: noise_level=%f' % (noise_level)
f.suptitle(qhd_params + subgrad_params + lfmsgd_params)
#f.suptitle('dt = ' + str(dt) + "; gran. = " + str(len) + "^" + str(dim) + ";\n tot_steps = " + str(tot_steps) + "; learning_rate = " + str(learning_rate) + "; sample_number = " + str(sample_number) + "; noise_level = " + str(noise_level))
plt.savefig(args.output)
if not args.silent:
    plt.show()
'''