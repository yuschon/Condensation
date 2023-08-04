import numpy as np
import math
import matplotlib.pyplot as plt
file1='Bubble_Condensation.dat'
a1=np.loadtxt(file1)

xData=a1[:,0]

yData0=a1[:,1]
yData1=a1[:,4]
yData2=a1[:,5]

plt.xlabel('x-axis', size=14)
plt.ylabel('y-axis', size=14)
plt.plot(xData, yData0, color='black', linestyle='-', label='Total')
plt.plot(xData, yData1, color='r', linestyle='-', label='Bottom')
plt.plot(xData, yData2, color='blue', linestyle='-', label='Top')

plt.savefig('plot.png')
plt.show()

