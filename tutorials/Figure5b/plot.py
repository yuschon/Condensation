import numpy as np
import math
import matplotlib.pyplot as plt
file='LiquidAccumulation.dat'

a0=np.loadtxt(file)

xData0=a0[:,0]
yData0=a0[:,2]
yData1=2*0.1355*np.sqrt(0.681*xData0/1000/4200)
yData2=np.sqrt(2*xData0*0.681/1000/4200/(0.5+2.4e6/4200/20))


plt.xlabel('Time/s', size=14)
plt.ylabel('Interface position/m', size=14)
plt.plot(xData0, yData1, color='b', linestyle='--', label='Analytical solution')
plt.plot(xData0, yData0, color='r', linestyle='-', label='Sun model')

plt.legend(loc='lower right')
plt.savefig('plot.png')
plt.show()

