"""
I realize there is no reason to make the data even in length. The only reason to make them even is for FFTs,
but when I define B I can choose an even number of points.

"""
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS

subdirs = ['Aug3Data','Aug4Data']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

j = 0
for sd in subdirs:
    files = listdirNH(sd+'/Raw')
    fig, ax = plt.subplots()
    for f in files:
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        # Make sure there are an even number of data points
        x = datavals[0][1:]
        y = datavals[1][1:]
        ly = 10**(y/10)
        N = len(x)
        L = x[-1]-x[0]
        
        for v in range(len(ly)):
            print(ly[v])
            if ly[v]<0:
                
                ly[v]==0
        
        
        np.savetxt(sd+'/'+'Rescaled/'+f[-10:-4]+'.out',np.transpose(np.append([x],[ly],axis=0)))
        
        m = max(y)
        i = np.where(y == m)

        
        
        ax.plot(x,ly)
    
    
    ax.set_title(sd)
    ax.set_xlabel('Frequency (mHz)')
    ax.set_ylabel('Energy Density (cm'+r'$^2$'+'/mHz)')
    #plt.show()
    plt.savefig(sd+'/Processed Data Figs/'+sd+'.png',dpi = 500)
    j=j+1