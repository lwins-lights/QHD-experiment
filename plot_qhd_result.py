# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# get default dpi
dpi = mpl.rcParams['figure.dpi']

# multiplots
fig, (plt1, plt2) = plt.subplots(2, 1, figsize=(800/dpi, 900/dpi))

# load params from file
res = np.load('./result/qipopt.npz')
T = res['T'][0]
dt = res['dt'][0]
len = res['len'][0]
dim = res['dim'][0]

#print(T)

# plot 1
x = np.arange(0, T, dt)
y = res['expected_potential']

#print(y)
  
plt1.set_xlabel('time')
plt1.set_ylabel('potential')

plt1.plot(x, y)

# plot 2
x = np.arange(0, T, dt)
y = res['expected_kinetic'] + res['expected_potential']
plt2.plot(x, y, label="Total")

y = res['expected_kinetic']
plt2.plot(x, y, label="Kinetic")

y = res['expected_potential']
plt2.plot(x, y, label="Potential")

plt2.set_xlabel('time')
plt2.set_ylabel('energy')
plt.legend()

# show all
plt.suptitle('dt = ' + str(dt) + "; gran. = " + str(len) + "^" + str(dim))
plt.savefig('result/result_qhd.png')
plt.show()