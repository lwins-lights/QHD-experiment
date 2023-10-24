# importing the required module
import sys, getopt
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


res = np.load('./result/qipopt.npz')
T = res['T'][0]
dt = res['dt'][0]
len = res['len'][0]
dim = res['dim'][0]

x = np.arange(0, T, dt)
y = res['expected_potential']
  
plt.xlabel('time')
plt.ylabel('expectation')

plt.plot(x, y)

# show all
plt.suptitle('dt = ' + str(dt) + "; gran. = " + str(len) + "^" + str(dim))
plt.savefig('result/result_qcp.png')
plt.show()