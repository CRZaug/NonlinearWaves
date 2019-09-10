"""
This code makes the data even in length, but that's not really necessary. See RescaleFiles1.py
"""
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS

subdirs = ['Aug1Data','Aug2Data','JulyData']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

j = 0
for sd in subdirs:
    files = listdirNH(sd+'/Raw')
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
            
            np.savetxt(sd+'/'+'Rescaled/'+f[-10:],np.transpose(np.append([x],[ly],axis=0)))
            
            m = max(y)
            i = np.where(y == m)
            print(i,m,x[i])
            
            
            plt.plot(x,ly)
        else:
            x = datavals[0]
            y = datavals[1]
            ly = 10**(y/10)
            N = len(x)
            L = x[-1]-x[0]
            
            np.savetxt(sd+'/'+'Rescaled/'+f[-10:],np.transpose(np.append([x],[ly],axis=0)))
            
            m = max(y)
            i = np.where(y == m)
            print(i,m,x[i])
            
            plt.plot(x,ly)
    plt.title(sd)
    plt.show()
    j=j+1