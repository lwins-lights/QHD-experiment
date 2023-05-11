# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# silent mode (no window pops up) if -s
opts, args = getopt.getopt(sys.argv[1:], "hs", ["help", "silent"])
is_silent = False
for o, a in opts:
    if o in ("-s", "--silent"):
        is_silent = True

# get default dpi
dpi = mpl.rcParams['figure.dpi']

qhd_res = np.load('./result/pseudospec.npz')
nagd_res = np.load('./result/nagd.npz')
  
y_qhd = qhd_res['expected_potential']
y_nagd = nagd_res['expected_potential']
yp_qhd = qhd_res['probability_at_minimum']
yp_nagd = nagd_res['probability_at_minimum']
T = qhd_res['T'][0]
dt = qhd_res['dt'][0]
len = qhd_res['len'][0]
dim = qhd_res['dim'][0]

x = np.arange(0, T, dt)

# subplots
f, (e_plt, p_plt) = plt.subplots(2, 1, figsize=(800/dpi, 900/dpi))
  
# plotting for expected potential
e_plt.plot(x, y_qhd, label="QHD*" + str(qhd_res['par'][0]))
e_plt.plot(x, y_nagd, label="NAGD*" + str(nagd_res['par'][0]))
  
e_plt.set_xlabel('time')
e_plt.set_ylabel('expectation')
  
#e_plt.set_title('dt = ' + str(dt[0]))

e_plt.legend(loc="upper right")

# plotting for probability at minimum
p_plt.plot(x, yp_qhd, label="QHD*" + str(qhd_res['par'][0]))
p_plt.plot(x, yp_nagd, label="NAGD*" + str(nagd_res['par'][0]))
  
p_plt.set_xlabel('time')
p_plt.set_ylabel('success probability')
  
#p_plt.set_title('dt = ' + str(dt[0]))

p_plt.legend(loc="upper right")

# show all
f.suptitle('dt = ' + str(dt) + "; gran. = " + str(len) + "^" + str(dim))
plt.savefig('result/result.png')
if not is_silent:
    plt.show()