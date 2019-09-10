"""
This is some test code to play around with generating surfaces from just Fourier magnitudes.
1. Read in rescaled data (scaled as per your MATLAB code)
2. Get the frequencies (these are equivalent to k, but I call them x just to make things confusing)
3. Generate random complex values with the correct magnitude
4. Find the location of the carrier wave
5. Shift over the k data so that the 0th mode lines up with the carrier with the max amp.
6. Restrict the frequency range to +/- 10 modes (for now)
7. Get B using the shifted k values and the restricted amplitudes
8. Plot results
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
        
        ly = np.sqrt(sly*x)*0.01
        
        
        L = 3600*3
        k = (0.001*x)//(2*np.pi/L) # Gets me to milliHertz THIS IS THE K VECTOR
        
        
        scalevals = np.sqrt(ly/2)
        randvals = np.zeros(len(ly))
        randpn = np.zeros(len(ly))
        pn = [-1,1]
        for l in range(len(ly)):
            rand.seed(l) # Seed random number generator
            randpn[l] = rand.choice(pn)
            rand.seed(l+1)
            value = rand.random()*rand.choice(pn) # Get a random value with a random sign
            randvals[l] = value
        ascale = randvals*scalevals # Scale the value

        bscale = 1j*np.sqrt(ly-ascale**2)*randpn # Get the complex part

        fakevals = ascale+bscale # THESE ARE NOW THE FOURIER AMPLITUDES
        
        
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
            print(f)
            print(carriermode)
            print(carrierloc)
        n = n+1
        
        # Rearrange data so the carrier wave is at mode zero
        
        # First, find the carrier mode in ANY file, not just the first one
        loc = np.where(np.logical_and(x>carrierloc*0.999, x<carrierloc*1.009))
     
        # Be a little more restrictive
        if len(loc[0])>1:
            loc = loc[0][0]
        else:
            loc = loc[0][0]
        print(carrierloc,loc,x[loc],fakevals[loc])
        
        # REDEFINE THE K VECTOR SO THAT THE CARRIER MODE IS AT 0
        knew = k-k[loc]
        xnew = x-x[loc]

        # REMOVE DUPLICATE VALUES
        ndk = np.array(())
        ndy = np.array(())
        for f in range(len(knew)):
            num = knew[f]
            yum = fakevals[f]
            if num not in ndk:
                ndk = np.append(ndk,num)
                ndy = np.append(ndy,yum)

        plt.plot(ndk,np.abs(ndy))
        plt.show()
        

        # RESTRICT THE FREQUENCY RANGE
        fakevals[:loc-10]=0+0j
        fakevals[loc+10:]=0+0j
        

        # GET B
        # Define new t data (and a new number of data points)
        NNN = 1024
        tnew = np.linspace(-L/2,L/2-L/NNN,NNN) # ADD TO ICTEMPORALDATA
        # Find B
        B=np.zeros(NNN,dtype=complex) # ADD TO ICTEMPORALDATA
        for v in range(N):
            newvalue = 1/N*fakevals[v]*np.exp(0.001*xnew[v]*1j*tnew) # Sum a new fourier series
            B = B+newvalue
        
        
        # Plot results
        # fig, ax = plt.subplots(2,1,figsize = (8,6))
        # fig.suptitle(sd+' '+f[-10:-4],fontsize = 20)
        # ax[0].set_title('Fourier Space')
        # ax[0].plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.',label = 'B')
        # ax[0].plot(xnew,1/N*np.abs(fakevals),'.',label = 'actual data')
        # ax[0].legend()
        # ax[0].set_xlabel('K value (1)')
        # ax[0].set_ylabel('Fourier Amplitude (m)')
        # 
        # ax[1].set_title('Temporal Space')
        # ax[1].plot(tnew,np.real(B))
        # ax[1].set_xlabel('time (s)')
        # ax[1].set_ylabel('Amplitude (m)')
        # 
        # 
        # fig.tight_layout()
        # fig.subplots_adjust(top=0.88)
        # plt.show()

          