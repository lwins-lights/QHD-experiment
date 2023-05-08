# importing the required module
import matplotlib.pyplot as plt

import numpy as np

qhd_res = np.load('./result/pseudospec.npz')
nagd_res = np.load('./result/nagd.npz')
  
y_qhd = qhd_res['expected_potential']
y_nagd = nagd_res['expected_potential']
T = qhd_res['T']
dt = qhd_res['dt']

x = np.arange(0, T, dt)
  
# plotting the points 
plt.plot(x, y_qhd, label="QHD*" + str(qhd_res['par'][0]))
plt.plot(x, y_nagd, label="NAGD*" + str(nagd_res['par'][0]))
  
plt.xlabel('time')
plt.ylabel('expectation')
  
plt.title('dt = ' + str(dt[0]))

plt.legend(loc="upper right")
  
plt.show()