"""

~~~ IMPORT EXPERIMENTAL DATA, PROCESS, AND NONDIMENSIONALIZE ~~~

This code reads in the rescaled Snodgrass data, finds important values, creates a
new surface B, and nondimensionalizes the data for simulation input.

1. Get the distances between gauges (from the Snodgrass gauge distances)
2. Read in each file for each storm event in a big loop
3. Interpolate the data and restrict it down based on a given frequency to get new data
4. Get the k vector
5. Get real and imaginary phase (using random phase assumption)
6. Get Hermitian amplitudes to examine the plain surface (optional)
7. Find the location of the carrier wave
8. Find the constants for nondimensionalization
9. Factor out the carrier wave
10. Sum to get back to temporal space and save
11. Nondimensionalize and save

"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
import random as rand
from scipy import interpolate

# Choose whether or not to show/save any plots (either "no" or "go")
doIplot = 'no'
doIplot_chi = 'go'



### STEP 1: Get distance information
distv = np.array([0.0,2400000.0,4200000.0,8700000.0]) # Distances between gauges in METERS



### STEP 2: Read in information at each gauge for each event
subdirs = ['Aug1Data','Aug2Data','JulyData']
# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

j = 0
for sd in subdirs:
    files = listdirNH(sd+'/Rescaled')
    
    # Initialize some values
    n = 0
    pi =0
    gaugenumber = 0
    
    if doIplot_chi == 'go':
        fig1,ax1 = plt.subplots(4)
    
    # Get files
    for f in files:
        print()
        
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        sly = datavals[1] # Magnitudes
        
        
        
### STEP 3: Interpolate the data and get new y values
        
        # Interpolation function   
        fnc = interpolate.interp1d(x, sly,kind ='cubic')
        sly2 = fnc(x)
        slymax = max(sly2)
        
        xmaxloc = np.where(slymax==sly2)[0][0]
    
        aa = x[0]
        bb = x[-1]
    
        # Choose the spacing of the grid based on the first maximum value
        if gaugenumber == 0:
        
            #newn = 25 # An arbitrary number of points that provides a range of x values
            newn = 50
            deltaf=x[xmaxloc]*2/(newn+1) # The new spacing between points
            
            # Create the new grid of x points based on the spacing deltaf and the location of the max
            xn = np.zeros(newn)
            xn[newn//2]=x[xmaxloc]
            for p in range(newn//2-1):
                xn[newn//2+p] = x[xmaxloc]+deltaf*(p)
                xn[newn//2-p] = x[xmaxloc]-deltaf*(p)

        # Restrict the x range to the x range given by the actual Snodgrass data
        xfinal = np.array([])
        for ixn in xn:
            if np.logical_and(bb>ixn,aa<ixn):
                xfinal=np.append(xfinal,ixn)
        
        
        # Get the new values
        ly = fnc(xfinal)*deltaf
        
        # Adjust the units on the x axis to get to mHz
        xfinal = xfinal*0.001
        
        gaugenumber +=1      
        
        
        
### STEP 4: Get the k vector using integer division and clean up
        lenly = len(ly)
        L = 3600*5 # The chosen period
        k = (xfinal)//(2*np.pi/L) # Then divide by 2pi/L (using integer division) to get the k vector
        
        
        
### STEP 5: Generate random phase and real and imaginary parts
        randvals = np.random.rand(lenly)*2*np.pi
        ascale = np.cos(randvals)*np.sqrt(ly)*0.01 # Rescale y axis
        bscale = np.sin(randvals)*np.sqrt(ly)*0.01 # Rescale y axix
        
        # Add the real and imaginary parts together
        fakevals = (ascale+1j*bscale)



### STEP 6: Remove duplicate values and generate a 2-sided Hermitian spectrum for testing
        ndk = np.arange(k[0],k[-1]+1)
        ndy = np.zeros(len(ndk),dtype=complex)
        
        ndk= np.append(ndk,ndk[-1]+1)
        ndy = np.append(ndy,0)
        
        pip = 0
        for ele in range(len(ndk)):
            if ndk[ele] in k:
                ndy[ele] = fakevals[pip]
                pip+=1

        extrazeros = np.zeros(int(ndk[0]))
        positivey = np.append(extrazeros,ndy)
        negativey = np.conj(np.append(np.flip(ndy[:-1],axis=0),extrazeros[1:]))
        
        # New y values
    
        ynew2=np.append(positivey,negativey)
        #print('ynew',max(np.abs(ynew2)),len(ynew2))
        ynew2 = ynew2*len(ynew2)
        extrak = np.arange(0,ndk[0])
        
        # New x axis values
        k1 = np.append(np.append(np.append(extrak,ndk),-np.flip(ndk[:-1],axis=0)),-np.flip(extrak,axis=0))[:-1]
    
        # Optional plotting
        if doIplot=='go':
            
            bbb = np.real(ifft(ynew2))
            timex = np.linspace(0,L,len(bbb))
            
            f1,a=plt.subplots(3)
            
            a[0].plot(k1,np.real(ynew2),'.')
            a[1].plot(k1,np.imag(ynew2),'.')
            a[2].plot(k1,ynew2*np.conj(ynew2),'.',markersize=5)
            a[2].plot(k,ly*len(ynew2)/2*0.01**2*len(ynew2)/2,'.')
            
            f1.suptitle('Period: '+ str(L) + ' s')
            f1.subplots_adjust(top=0.88)
            
            
            g,b=plt.subplots(3)
            b[0].plot(timex,np.real(bbb))
            b[1].plot(timex,np.imag(bbb))
            b[2].plot(timex,np.abs(bbb))
            b[2].set_ylabel('Wave Height (m)')
            b[1].set_ylabel('Wave Height (m)')
            b[0].set_ylabel('Wave Height (m)')
            b[2].set_xlabel('Time (s)')
            b[1].set_xlabel('Time (s)')
            b[0].set_xlabel('Time (s)')
            g.suptitle(f)
            g.subplots_adjust(top=0.88) 
        # plt.show()
        
        
### STEP 7: Locate the carrier wave
        #halfy = ynew2[:len(ynew2)//2]
        carrierfreqs = k1*(2*np.pi/L)
        if n==0:
            
            carriermax = max(np.abs(ynew2))
            carrierloc = np.where(np.logical_and(np.abs(ynew2)>=carriermax, np.abs(ynew2)<=carriermax))[0][0]
            CarFreq = k1[carrierloc]
        
        loc = np.where(k1==CarFreq) # The location of the carrier wave
        
        
### STEP 8: Get nondimensionalization constants    
        g = 9.81 #(m/s^2)
        if n==0:
            w0 = carrierfreqs[loc] # Get the value from the integer
            k0 = w0**2/g # The carrier wavenumber
            m = max(np.abs(ynew2)/len(ynew2))
            epsilon = 2*m*k0 # The nondimensionalization constant epsilon
        
            print(f,'Special Values')
            print('delta f', deltaf)
            print('period',L)
            print('Maximum value',m)
            print('Wavenumber',float(k0))
            print('Carrier frequency',float(w0))
            print('epsilon',float(epsilon))
            
            file = open(sd+'/SpecialVals.txt','w') 
 
            file.write('delta f, '+str(float(deltaf))+'\n') 
            file.write('period, '+str(float(L))+'\n') 
            file.write('maximum value, '+str(float(m))+'\n') 
            file.write('wavenumber k0, '+str(float(k0))+'\n')
            file.write('carrier frequency w0, '+ str(float(w0))+'\n')
            file.write('epsilon, '+str(float(epsilon))+'\n')
 
            file.close() 
           
           
           
            print()
            
            n = n+1
  
        
        
### STEP 9: Factor out the carrier wave

        # Shorten the old spectrums
        yhat1 = ynew2[:len(ynew2)//2] 
        k_new_1 = k1[:len(ynew2)//2]
        newcarrierfreqs = carrierfreqs[:len(ynew2)//2]
        k_new_2 = k_new_1 - k_new_1[loc]

    
    
### STEP 10: Sum to get back to temporal space and save
    
        # Define new t data (and a new number of data points)
        NNN = 1024
        tnew = np.linspace(0,L-L/NNN,NNN) # ADD TO ICTEMPORALDATA
        yN = len(yhat1)
        
        # Find B, the new surface modulating envelope
        B=np.zeros(NNN) 
        for v in range(len(yhat1)):
            newvalue = yhat1[v]/yN*np.exp(1j*(2*np.pi/L)*k_new_2[v]*tnew) # Sum a new fourier series
            B = B+newvalue
        
        # Optional plotting
        if doIplot=='go':
            bfig,bax = plt.subplots(2)
            bax[0].plot(tnew,np.real(B))
            bax[1].plot(tnew,np.imag(B))
            bfig.suptitle(str(f))
            bfig.subplots_adjust(top=0.88)
        
        
            figg,axx = plt.subplots()
            axx.set_title(f)
            axx.plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.',markersize=10)
            axx.plot(k_new_2,1/len(yhat1)*np.abs(yhat1),'.',markersize = 5) ###################
            plt.show()
        #plt.close('all')
        
        # Save the x and y values for t and B
        #np.savetxt(sd + '/Processed/'+f[-10:], np.transpose(np.append(tnew,B,axis=0)).view(float))
  
    

### STEP 11: Nondimensionalize values, plot, and save results
        ndB = k0/epsilon*B # new "y"
        xi = w0*epsilon*tnew # new "x"
        
        # Fix up xi so that it is periodic and has the correct number of points to match B
        xi_L = (xi[-1]-xi[0])+(xi[1]-xi[0])
        xi_new = np.linspace(0,xi_L-xi_L/NNN,NNN,dtype=complex)
        
        chi = epsilon**2*k0*distv # new "t"
    
        # Optional plotting
        if doIplot_chi=='go':
            ax1[pi].set_title(r'$\chi$'+' = '+str(chi[pi]))
            ax1[pi].plot(xi_new,np.real(ndB))
            ax1[pi].plot(xi_new,np.imag(ndB))
            
            pi+=1
        
        # Save the nondimensional x and y axis, chi and non dim B
        #np.savetxt(sd + '/NonDim Data/ND'+f[-10:],np.append([xi_new],[ndB],axis=0).view(float))
        
    # Save just the nondimensonal time vector
   # np.savetxt(sd + '/chi.txt', chi.view(float))
    
    # Optional plotting of all the dimensionless surfaces together
    n=0
    if doIplot_chi =='go':
        fig1.suptitle(sd+', L = '+str(L))
        fig1.tight_layout()
        fig1.subplots_adjust(top=0.88)
        #plt.show()
       # plt.savefig(sd + '/NonDim Figs/allgauges.png',dpi=500)