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
        
        ly = np.sqrt(sly*x)*0.01 # Get Amplitudes in units of CM
        
        L = 3600*3
        x = x*0.001*L # Gets me to milliHertz (but the 2 pi version...) THIS IS THE K VECTOR
        
        scalevals = np.sqrt(ly/2)
        randvals = np.zeros(len(ly))
        randpn = np.zeros(len(ly))
        pn = [-1,1]
        
        ascale = []
        randpn = np.zeros(N)
        ki = 0
        for h in scalevals:
            randpn[ki] = rand.choice(pn)
            while len(ascale) < N:
                value = rand.gauss(0.00001*h, 0.1)
                if -h < value < h:
                    ascale.append(value)
            ki=ki+1
                    
        ascale = np.array(ascale)
          
        bscale = 1j*np.sqrt(ly-ascale**2)*randpn # Get the complex part
        
        fakevals = ascale+bscale 
        
        
        # for k in range(len(ly)):
        #     rand.seed(k) # Seed random number generator
        #     randpn[k] = rand.choice(pn)
        #     rand.seed(k+1)
        #     value = rand.random()*rand.choice(pn) # Get a random value with a random sign
        #     randvals[k] = value
        # ascale = randvals*scalevals # Scale the value
        # 
        # bscale = 1j*np.sqrt(ly-ascale**2)*randpn # Get the complex part
        # 
        # fakevals = ascale+bscale # THESE ARE NOW THE FOURIER AMPLITUDES
        
        
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
        xnew = x-x[loc]

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
            newvalue = 1/N*fakevals[v]*np.exp(2j*np.pi*np.ceil(xnew[v])*tnew/L) # Sum a new fourier series
            B = B+newvalue
        
        
        # Plot results
        # plt.title(f) 
        # plt.plot(xnew,1/N*np.abs(fakevals),'.',label = 'actual data')
        # plt.plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.',label = 'B')
        # plt.legend()
        
        plt.title(f)
        plt.plot(tnew,np.abs(B))
        plt.xlabel('time (s)')
        plt.ylabel('Amplitude (m)')
        plt.show()

          