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
disc = res['disc'][0]
len = res['len'][0]
dim = res['dim'][0]

#print(T)

# plot 1
y = res['expected_potential']
x = res['timestamps']
#print(x)

#print(y)
  
plt1.set_xlabel('time')
plt1.set_ylabel('potential')
plt1.set_yscale('log')

plt1.plot(x, y)

# plot 2

y = res['expected_duality_gap']
x = res['timestamps']

plt2.set_xlabel('time')
plt2.set_ylabel('duality gap')
plt2.set_yscale('log')

plt2.plot(x, y)

# show all
plt.suptitle('disc = ' + str(disc) + "; gran. = " + str(len) + "^" + str(dim))
plt.savefig('result/result_qhd.png')
plt.show()