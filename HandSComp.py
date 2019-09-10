"""

~~~ IMPORT EXPERIMENTAL DATA, PROCESS, AND NONDIMENSIONALIZE ~~~

This code reads in the rescaled Snodgrass data and compares parameters
to known parameters found in the Henderson and Segur paper.

1. Get distances
2. Read in the gauge data for each event (get frequencies and Fourier magnitudes)
3. Adjust the y axis units
4. Get the k vector using integer division and clean up
5. Get the carrier wave location (requires some restricting)
6. Factor out carrier wave
7. Get the energies at each gauge
8. Get nondimensionalization constants

"""
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
import random as rand
from scipy import interpolate

### STEP 1: Get distance information
distv = np.array([0.0,2400000.0,4200000.0,8700000.0]) # Distances between gauges in METERS



### STEP 2: Read in information at each gauge for each event
subdirs = ['Aug1Data','Aug2Data','JulyData']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))


# Read in the data
j = 0

for sd in subdirs:
    files = listdirNH(sd+'/Rescaled')
    
    # Initialize some values
    n = 0
    pi =0
    fig1,ax1 = plt.subplots(4,1)
    plt.suptitle(sd)
    
    # Get files
    Deltavals = []
    for f in files:
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        sly = datavals[1] # Magnitudes
        
        #ly = np.sqrt(sly*x)*0.01 #MULTIPLY VERSION (get the amplitude in meters)
        
        mns = []
        for w in range(N-1):
            mns.append(np.abs(x[w+1]-x[w]))
        #mns.append(np.mean(mns))
        
        
        
### STEP 3: Adjust the y axis units
        ly = np.sqrt(sly*np.mean(mns))*0.01 # INTEGRATED VERSION 

### STEP 4: Get the k vector using integer division and clean up
        L = 3600*3 # The period
        k = (x*0.001)//(2*np.pi/L) # Convert to mHz, then divide by 2pi/L to get the k vector

        # REMOVE DUPLICATE VALUES
        ndk = np.array(())
        
        for fi in range(len(k)):
            num = k[fi]
            if num not in ndk:
                ndk = np.append(ndk,num)
        
        lll =[]
        for h in ndk:
            l1=np.where(k==h)[0]
            lll.append(l1)
        ndy = np.array(())
        for ar in lll:
            val = np.mean(ly[ar])

            ndy=np.append(ndy,val)
                    
### STEP 5: Get the location of the carrier wave (defined by the first gauge)
        if n == 0:
            m = max(ndy)
            i = np.where(ndy == m)
            if len(i[0]) > 1:
                newi = i[0][len(i[0])//2]

                carriermode = np.array(newi)
                carrierloc = ndk[carriermode]
            else:
                newi = i[0][0]
                carriermode = np.array(newi)
                carrierloc = ndk[carriermode]

        
        # First, find the carrier mode in ANY file, not just the first one
        loc = np.where(np.logical_and(ndk>carrierloc*0.99, ndk<carrierloc*1.001))
        #loc = np.where(np.logical_and(ndk>carrierloc-1, ndk<carrierloc+1))
        
        # Be a little more restrictive
        if len(loc[0])>1:
            loc = loc[0][0]
        else:
            loc = loc[0][0]
            
            
            
### STEP 6: Redefine the k vector so that the carrier mode is at 0 (factor it out)
        knew = ndk-ndk[loc]
        xnew = x-x[loc]
        
### STEP 7: Get the "Energy" integrals
        fnc = interpolate.interp1d(x, sly,kind ='cubic')
        longx = np.linspace(x[0],x[-1],1000)
        newy = fnc(longx)
        A0 = np.sqrt(2*NLS.trapezoid(newy,(x[-1]-x[0])))*0.01
        
        
        
        figg,axx = plt.subplots()
        axx.plot(x,sly,'.',markersize=7)
        axx.plot(longx,newy)
        plt.show()


        M000 = NLS.trapezoid(newy[np.where(np.logical_and(longx>41.2,longx<75.6))],(74.6-41.2))
     
        Deltavals.append(M000)
        
### STEP 8: Get nondimensionalization constants
        g = 9.81 #(m/s^2)
        if n==0:
            w0 = (2*np.pi)**2/L*ndk[loc] # Get the value from the integer
            k0 = w0**2/g # The carrier wavenumber
            m = max(ndy)
            epsilon = 2*m*k0 # The nondimensionalization constant epsilon
            heps = A0*k0
            
            
            
            
            
            print(f,'Special Values')
            print('2A0',A0)
            print('Maximum value',m)
            print('Carrier frequency',w0)
            print('Wavenumber',k0)
            print('MY epsilon',epsilon)
            print('HENDERSON EPSILON', heps)
            print('period',L)

            
        n = n+1
    
    M0 = Deltavals[0]
    
    MX = Deltavals/M0
    energyy = np.log(MX)
    
    # Get the fit and define a new y vector
    
    A = np.vstack([distv, np.ones(len(distv))]).T
    m, b = np.linalg.lstsq(A, energyy,rcond=-1)[0] # m is delta
    
    hdeltab = -m
    
    hdeltas = hdeltab/(2*heps**2*k0)
    
    xplot = np.linspace(distv[0],distv[-1],100)
    newy = m*xplot+b
    print('HENDERSON BIG Delta ',hdeltab, 'b ', b)
    print('HENDERSON LITTLE delta', hdeltas)
    print()
