# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# silent mode (no window pops up) if -s
opts, args = getopt.getopt(sys.argv[1:], "hsS", ["help", "silent", "SGDonly"])
is_silent = False
sgd_only = False
for o, a in opts:
    if o in ("-s", "--silent"):
        is_silent = True
    if o in ("-S", "--SGDonly"):
        sgd_only = True

# get default dpi
dpi = mpl.rcParams['figure.dpi']

if not sgd_only:
    qhd_res = np.load('./result/pseudospec.npz')
    nagd_res = np.load('./result/nagd.npz')
sgd_res = np.load('./result/sgd.npz')

if not sgd_only:  
    y_qhd = qhd_res['expected_potential']
    y_nagd = nagd_res['expected_potential']
y_sgd = sgd_res['expected_potential']

if not sgd_only:
    yp_qhd = qhd_res['probability_at_minimum']
    yp_nagd = nagd_res['probability_at_minimum']
yp_sgd = sgd_res['probability_at_minimum']
if not sgd_only:
    T = qhd_res['T'][0]
    dt = qhd_res['dt'][0]
    len = qhd_res['len'][0]
    dim = qhd_res['dim'][0]
    noise_level = sgd_res['noise_level'][0]
else:
    T = sgd_res['T'][0]
    dt = sgd_res['dt'][0]
    len = sgd_res['len'][0]
    dim = sgd_res['dim'][0]
    noise_level = sgd_res['noise_level'][0]

x = np.arange(0, T, dt)

# subplots
f, (e_plt, p_plt) = plt.subplots(2, 1, figsize=(800/dpi, 900/dpi))
  
# plotting for expected potential
if not sgd_only:
    e_plt.plot(x, y_qhd, label="QHD*" + str(qhd_res['par'][0]))
    e_plt.plot(x, y_nagd, label="NAGD*" + str(nagd_res['par'][0]))
e_plt.plot(x, y_sgd, label="SGD*" + str(sgd_res['par'][0]))
  
e_plt.set_xlabel('time')
e_plt.set_ylabel('expectation')

# plotting for probability at minimum
if not sgd_only:
    p_plt.plot(x, yp_qhd, label="QHD*" + str(qhd_res['par'][0]))
    p_plt.plot(x, yp_nagd, label="NAGD*" + str(nagd_res['par'][0]))
p_plt.plot(x, yp_sgd, label="SGD*" + str(sgd_res['par'][0]))
  
p_plt.set_xlabel('time')
p_plt.set_ylabel('success probability')
  
# style
e_plt.yaxis.tick_right()
p_plt.yaxis.tick_right()
e_plt.legend(loc="upper right")
p_plt.legend(loc="lower right")

# show all
f.suptitle('dt = ' + str(dt) + "; gran. = " + str(len) + "^" + str(dim) + "; noise_level = " + str(noise_level))
plt.savefig('result/result.png')
if not is_silent:
    plt.show()