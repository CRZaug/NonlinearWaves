import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS

sd ='Aug2Data'
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

labels = ['Tutuila','Palmyra', 'Honolulu', 'Yakutat']
plt.figure(figsize=(6,5))
j = 0
u=0
files = listdirNH(sd+'/Raw')
print(files)
for f in files:
    datavals = np.transpose(np.loadtxt(f).view(float))
    N = len(datavals[1])
    # Make sure there are an even number of data points
    if N % 2 != 0:
        x = datavals[0][1:]
        y = datavals[1][1:]
        ly = 10**(y/10)
        N = len(x)
        L = x[-1]-x[0]
        
        m = max(y)
        i = np.where(y == m)
        print(i,m,x[i])
        
        i = i[0]
        
        plt.plot(x,ly,label = labels[u])
        plt.plot(np.mean(x[i]),np.mean(ly[i]),'k.', markersize = 10)
        u=u+1
    else:
        x = datavals[0]
        y = datavals[1]
        ly = 10**(y/10)
        N = len(x)
        L = x[-1]-x[0]
        
        
        m = max(y)
        i = np.where(y == m)
        print(i,m,x[i])
        i = i[0]
        
        plt.plot(x,ly,label = labels[u])
        plt.plot(np.mean(x[i]),np.mean(ly[i]),'k.',markersize = 10)
        u=u+1
        
plt.title('Aug 13.7 Swell Data')
plt.legend()
plt.ylabel('Energy Density (cm'+'$^2$'+'/mHz)')
plt.xlabel('Frequency (mHz)')
#plt.show()
plt.savefig('Aug2TalkFig.png',dpi=480)
j=j+1