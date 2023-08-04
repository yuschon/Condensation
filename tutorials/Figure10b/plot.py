import numpy as np
import math
import matplotlib.pyplot as plt
file='LiquidAccumulation.dat'
a0=np.loadtxt(file)
xData0=a0[:,0]
yData0=a0[:,2]
yData1=np.sqrt(2*xData0*0.8/958.8/4200/(0.5+2.26e6/4200/20))

plt.xlabel('Time/s', size=14)
plt.ylabel('Interface position/m', size=14)
plt.plot(xData0, yData1, color='b', linestyle='--', label='Analytical solution')
plt.plot(xData0, yData0, color='r', linestyle='-', label='mcor=30')

plt.legend(loc='lower right')
plt.savefig('Stefan.png')
plt.show()

