import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
import random as rand

# seed random number generator


subdirs = ['Aug1Data','Aug2Data','JulyData']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

j = 0
for sd in subdirs:
    files = listdirNH(sd+'/Rescaled')
    n = 0
    for f in files:
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        ly = datavals[1] # Magnitudes
        L = 3600*3
        x = x*0.001 # Gets me to milliHertz (but the 2 pi version...)
        
        scalevals = np.sqrt(ly/2)
        randvals = np.zeros(len(ly))
        randpn = np.zeros(len(ly))
        pn = [-1,1]
        for k in range(len(ly)):
            rand.seed(k)
            randpn[k] = rand.choice(pn)
            rand.seed(k+1)
            value = rand.random()*rand.choice(pn)
            randvals[k] = value
        ascale = randvals*scalevals

        bscale = 1j*np.sqrt(ly-ascale**2)*randpn

        fakevals = ascale+bscale
        
        # Need to make a new k vector that actually has the frequencies we want
        # Find the average spacing between the freqs
        av = 0
        for ind in range(1,len(x)):
            diff = x[ind]-x[ind-1]
            av = av+diff
        aav = av/len(x)
        
        
        xappend = np.linspace(0,x[0],int(x[0]/aav))
        halfx = np.append(xappend,x)
        fullx = np.append(halfx,-np.flip(halfx[1:-1],0))*1/L
        
        yappend1 = np.zeros(len(xappend),dtype=complex)
        yappend2 = np.zeros(len(halfx[1:-1]),dtype=complex)
        fully = np.concatenate((yappend1,fakevals,yappend2))
    
        
        if n == 0:
            m = max(fully)
            i = np.where(fully == m)
            if len(i[0]) > 1:
                newi = i[0][len(i[0])//2]
                carriermode = np.array(newi)
                carrierloc = fullx[carriermode]
            else:
                newi = i[0][0]
                carriermode = np.array(newi)
                carrierloc = fullx[carriermode]
            print(f)
            print(carriermode)
            print(carrierloc)
        n = n+1
        
        # Rearrange data so the carrier wave is at mode zero
        
        # First, find the carrier mode
        #loc = np.where(fullx == carrierloc)
        loc = np.where(np.logical_and(fullx>carrierloc*0.999, fullx<carrierloc*1.009))
        print(carrierloc,loc,fullx[loc])
        
        if len(loc[0])>1:
            loc = loc[0][0]
        else:
            loc = loc[0][0]
        #print(carrierloc,loc,fullx[loc],fully[loc])
        
        
        
        
        yhata = fully[0:loc] 
        yhatb = fully[loc:] 
        yhatnew = np.append(yhatb,yhata) # Putting the carrier frequency in the 0 location 
        
        #plt.plot(fullx,np.abs(yhatnew),'.',markersize = 5)
        
        
        
        
        
        
        # Define new t data (and a new number of data points)
        NNN = 1024
        tnew = np.linspace(0,10800-10800/NNN,NNN) # ADD TO ICTEMPORALDATA
        
        # Find B
        B=0 # ADD TO ICTEMPORALDATA
        for v in range(N):
            newvalue = yhatnew[v]*np.exp(2j*np.pi*fullx[v]*tnew) # Sum a new fourier series
            B = B+newvalue
        
        plt.title(f)
        #plt.plot(tnew,np.real(B))
        
        plt.plot(fullx,np.abs(yhatnew),'.',label = 'act')
        plt.plot(1/L*NLS.kvec(NNN)*0.001,1/NNN*np.abs(fft(B)),'.',label = 'B')
        plt.legend()
        
        
        plt.show()

    
    # plt.title(sd)
    # plt.show()
          