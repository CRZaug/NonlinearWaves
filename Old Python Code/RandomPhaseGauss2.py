"""
This is some test code to play around with generating surfaces from just Fourier magnitudes.

Generates random PHASE by using angles, tries to copy what happens in the Feb Data (albeit poorly)

"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
import random as rand

subdirs = ['Aug1Data','Aug2Data','JulyData']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

# Read in the data
j = 0
for sd in subdirs:
    files = listdirNH(sd+'/Rescaled')
    n = 0
    for f in files:
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        sly = datavals[1] # Magnitudes
        
        #ly = np.sqrt(sly*x)*0.01 # Get Amplitudes in units of CM
        ly  = sly
        
        L = 3600*3
        
        k = (0.001*x)//(2*np.pi/L) # Gets me to milliHertz THIS IS THE K VECTOR
        
        # Generate angles
        aang = []
        for h in ly:
            while len(aang) < N:
                value = rand.gauss(0, 2)
                if -0.002 < value < 0.002:
                    aang.append(value)
        
        ascale = np.zeros(N)
        bscale = np.zeros(N,dtype=complex)
        for p in range(N):
            if p % 2 ==0:
                ascale[p]= ly[p]*np.sin(aang[p])
                bscale[p] = 1j*ly[p]*np.cos(aang[p])
            else:
                ascale[p]= ly[p]*np.cos(aang[p])
                bscale[p] = 1j*ly[p]*np.sin(aang[p])
        
        fakevals = ascale+bscale
        
        
    
        
        # FIND THE LOCATION OF THE CARRIER WAVE (defined by the first file)
        if n == 0:
            m = max(ly)
            i = np.where(ly == m)
            if len(i[0]) > 1:
                newi = i[0][len(i[0])//2]
                carriermode = np.array(newi)
                carrierloc = x[carriermode]
            else:
                newi = i[0][0]
                carriermode = np.array(newi)
                carrierloc = x[carriermode]
        n = n+1
        
        # Rearrange data so the carrier wave is at mode zero
        
        # First, find the carrier mode in ANY file, not just the first one
        loc = np.where(np.logical_and(x>carrierloc*0.999, x<carrierloc*1.009))
     
        # Be a little more restrictive
        if len(loc[0])>1:
            loc = loc[0][0]
        else:
            loc = loc[0][0]
        
        # REDEFINE THE K VECTOR SO THAT THE CARRIER MODE IS AT 0
        knew = k-k[loc]
        xnew = x-x[loc]
        
    
        
        # RESTRICT THE FREQUENCY RANGE
        # fakevals[:loc-10]=0+0j
        # fakevals[loc+10:]=0+0j
        
        # REMOVE DUPLICATE VALUES
        ndk = np.array(())
        ndy = np.array(())
        for fi in range(len(knew)):
            num = knew[fi]
            yum = fakevals[fi]
            if num not in ndk:
                ndk = np.append(ndk,num)
                ndy = np.append(ndy,yum)
        
        
        # plt.plot(xnew, 1/N*np.abs(fakevals),'.')
        # plt.show()

        # GET B
        # Define new t data (and a new number of data points)
        NNN = 1024
        tnew = np.linspace(-L/2,L/2-L/NNN,NNN) # ADD TO ICTEMPORALDATA
        # Find B
        B=np.zeros(NNN,dtype=complex) # ADD TO ICTEMPORALDATA
        for v in range(len(ndy)):
            newvalue = 1/N*ndy[v]*np.exp(2j*np.pi*ndk[v]/L*tnew) # Sum a new fourier series
            B = B+newvalue
        
        # Save Results
        np.savetxt(sd + '/Processed/'+f[-10:], np.transpose(np.append(tnew,B,axis=0)).view(float))
        
        # Plot results
        fig, ax = plt.subplots(2,1,figsize = (8,6))
        fig.suptitle(sd+' '+f[-10:-4],fontsize = 20)
        ax[0].set_title('Fourier Space')
        ax[0].plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.',markersize = 5, label = 'B')
        ax[0].plot(ndk,1/N*np.abs(ndy),'.',markersize = 3,label = 'actual data')
        ax[0].legend()
        ax[0].set_xlabel('K value (1)')
        ax[0].set_ylabel('Fourier Amplitude (m)')
        
        ax[1].set_title('Temporal Space')
        ax[1].plot(tnew,np.real(B))
        ax[1].set_xlabel('time (s)')
        ax[1].set_ylabel('Amplitude (m)')
        
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.88)
        plt.show()
        

          