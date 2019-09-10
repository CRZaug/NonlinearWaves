"""

~~~ IMPORT EXPERIMENTAL DATA, PROCESS, AND NONDIMENSIONALIZE ~~~

This code reads in the rescaled Snodgrass data, finds important values, creates a
new surface B, and nondimensionalizes the data for simulation input.

1. Get the distances between gauges
2. Read in the gauge data for each event (get frequencies and Fourier magnitudes)
3. Adjust the y axis units
4. Get the k vector using integer division and clean up
5. Get the carrier wave location (requires some restricting)
6. Factor out carrier wave
7. Get nondimensionalization constants
8. Determine the phases of the data from the Feb Data
9. Restrict frequency range (optional)
10. Determine the new surface, B, and save it
11. Nondimensionalize, plot, and save the nondimensional values
"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
import random as rand
############## GET PHASE DATA

phasedata = np.transpose(np.loadtxt('Henderson Phase Data/phasefebx050.out').view(float))
nphase = len(phasedata[0])
kphase = phasedata[0]
pphase = phasedata[1]
#print(phasedata)
knegs = kphase[nphase//2+1:]
kpos = kphase[:nphase//2]
pnegs = pphase[nphase//2+1:]
ppos = pphase[:nphase//2]

kp = np.append(knegs,kpos)
pp = np.append(pnegs,ppos)

phasezero = np.where(kp==0)[0]


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
        ly = np.sqrt(sly*np.mean(mns))*0.01 # INTEGRATED VERSION (get the amplitude in meters) (This is a lot messier as a result)



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
        
        
        
### STEP 7: Get nondimensionalization constants    
        g = 9.81 #(m/s^2)
        if n==0:
            w0 = (2*np.pi)**2/L*ndk[loc] # Get the value from the integer
            k0 = w0**2/g # The carrier wavenumber
            m = max(ndy)
            epsilon = 2*m*k0 # The nondimensionalization constant epsilon
            print(f,'Special Values')
            print('Maximum value',m)
            print('Carrier frequency',w0)
            print('Wavenumber',k0)
            print('epsilon',epsilon)
            print('period',L)
            print()
            
        n = n+1

          
### STEP 8: Determine the phases of the data
        # GENERATE "RANDOM" ANGLES
        lndy = len(ndy)
        zloc = int(np.where(knew==0)[0])
        #print(zloc)
        if lndy<nphase:
            #print('The Snodgrass is smaller')
            fakevals = np.zeros(lndy,dtype=complex)

            for z in range(zloc):
                
                phase1 = phasezero-z
                targetloc1 = zloc - z
                ang1 = pp[phase1]
                #print('part1',knew[targetloc1],kp[phase1],pp[phase1])
                ascale = ndy[targetloc1]*np.sin(ang1)
                bscale = 1j*ndy[targetloc1]*np.cos(ang1)
                fakevals[targetloc1]=ascale+bscale

                
            for z in range(1,lndy-zloc):    
                phase2 = phasezero+z
                targetloc2 = zloc+z
                ang2 = pp[phase2]
                #print('part2',knew[targetloc2],kp[targetloc2],pp[phase2])
                ascale = ndy[targetloc2]*np.sin(ang2)
                bscale = 1j*ndy[targetloc2]*np.cos(ang2)
                fakevals[targetloc2]=ascale+bscale
    
                
        else:
            #print('The Snodgrass is larger')
            aa = zloc-50
            for z in range(zloc-50,zloc):
                
                phase1 = phasezero-z
                targetloc1 = zloc - z
                ang = pp[phase1]
                
                ascale = ndy[targetloc1]*np.sin(ang1)
                bscale = 1j*ndy[targetloc1]*np.cos(ang1)
                fakevals[aa]=ascale+bscale
                aa=aa+1
            for z in range(zloc+1,zloc+50):    
                phase2 = phasezero+z
                targetloc2 = zloc+z
                ang2 = pp[phase2]
                
                ascale = ndy[targetloc2]*np.sin(ang2)
                bscale = 1j*ndy[targetloc2]*np.cos(ang2)
                fakevals[aa]=ascale+bscale
                aa=aa+1
            # need to add the vales before or after the feb data values expire
            aa=0
            for v in range(zloc-50):
                ang = rand.random()*rand.choice(pn)
                scale = ang*2*np.pi
                ascale = ndy[v]*np.cos(scale)
                bscale = ndy[v]*np.sin(scale)*1j
                fakevals[aa] = ascale+bscale
                aa=aa+1
            aa=zloc+50
            for v in range(zloc+50,lndy):
                ang = rand.random()*rand.choice(pn)
                scale = ang*2*np.pi
                ascale = ndy[v]*np.cos(scale)
                bscale = ndy[v]*np.sin(scale)*1j
                fakevals[aa] = ascale+bscale
                aa=aa+1
            
            print()  

        
        ####################
        knew = ndk
        fakevals=ndy

### STEP 9 (optional): Restrict the frequency range
        #fakevals[:loc-10]=0+0j
        #fakevals[loc+10:]=0+0j
        
        
        
### STEP 10: Get B, the new surface
        NNN = 1024 # New number of grid points
        tnew = np.linspace(0,L-L/NNN,NNN) # new temporal x axis
        
        # Find B
        B=np.zeros(NNN,dtype=complex)
        for v in range(len(fakevals)):
            newvalue = 1/N*fakevals[v]*np.exp(2j*np.pi*knew[v]/L*tnew)
            B= B+newvalue
        
        #print(max(np.abs(1/NNN*fft(B)))-max(1/N*np.abs(fakevals)))
        # plt.plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.')
        # plt.plot(knew,1/N*np.abs(fakevals),'.')
        # plt.show()
        
        # SAVE
        np.savetxt(sd + '/Processed/'+f[-10:], np.transpose(np.append(tnew,B,axis=0)).view(float))
        
        
     
### STEP 11: Nondimensionalize values, plot, and save results
        ndB = k0/epsilon*B # new "y"
        xi = w0*epsilon*tnew # new "x"
        
        # Fix up xi so that it is periodic and has the correct number of points to match B
        xi_L = (xi[-1]-xi[0])+(xi[1]-xi[0])
        xi_new = np.linspace(0,xi_L-xi_L/NNN,NNN,dtype=complex)
        
        chi = epsilon**2*k0*distv # new "t"
    
        #Plot results
        ax1[pi].set_title(r'$\chi$'+' = '+str(chi[pi]))
        ax1[pi].plot(xi_new,np.real(ndB))
        ax1[pi].plot(xi_new,np.imag(ndB))
        pi+=1
        
        # Save Results
        np.savetxt(sd + '/NonDim Data/ND'+f[-10:],np.append([xi_new],[ndB],axis=0).view(float))
        
    
    np.savetxt(sd + '/chi.txt', chi.view(float))
    fig1.tight_layout()
    fig1.subplots_adjust(top=0.88)
    #plt.show()
    plt.savefig(sd + '/NonDim Figs/allgauges.png',dpi=500)