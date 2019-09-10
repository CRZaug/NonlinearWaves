"""

~~~ IMPORT EXPERIMENTAL DATA, PROCESS, AND NONDIMENSIONALIZE ~~~

This code reads in the rescaled Snodgrass data, finds important values, creates a
new surface B, and nondimensionalizes the data for simulation input.

1. Get the distances between gauges
2. Read in the gauge data for each event (get frequencies and Fourier magnitudes)
3. Interpolate data and adjust the y axis units
4. Get the k vector using integer division and clean up
5. Get the carrier wave location (requires some restricting)
6. Factor out carrier wave
7. Get nondimensionalization constants
8. Determine the phases of the data
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
    
    # fig1,ax1 = plt.subplots(4,1)
    # plt.suptitle(sd)
    
    # Get files
    for f in files:
        print()
        print(f)
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        sly = datavals[1] # Magnitudes
        
        slymax = max(sly)
        
        print(max(sly))
    
        # Interpolation function   
        fnc = interpolate.interp1d(x, sly,kind ='cubic')

        
        
    
        deltaf = 4 #mHz

        aa = x[0]
        bb = x[-1]
    
        newn = (bb-aa)/deltaf
        xn = np.arange(int(newn)+1)*deltaf+aa
        
        #sly= fnc(x)
        ly = fnc(xn)
       
 
    ### STEP 3: Adjust the y axis units
        
        mmax = max(ly)
        maxi = np.where(ly==mmax)[0]
        
        xn = xn*0.001
        
    
### STEP 5: Get the location of the carrier wave (defined by the first gauge)
        if n == 0:
            m = max(ly)
            newi = np.where(ly == m)[0]
            
            carriermode = np.array(newi)
            carrierloc = xn[carriermode]
            

        
        # First, find the carrier mode in ANY file, not just the first one
        loc = np.where(np.logical_and(xn>carrierloc*0.9, xn<carrierloc*1.1))

        # Be a little more restrictive
        #print('loc',loc)
        
        if len(loc[0])>1:
            loc = loc[0][0]
        else:
            loc = loc[0][0]
            
        #print('new',loc)
        
        ####################################################### DO THE STUFF IN OCEAN OPTICS WEBBOOK
        
        lenly = len(ly)
        # print(lenly)
        # if lenly % 2 == 0:
        #     ly = ly
        # else:
        #     print('fix it')
        #     ly = ly[1:]
        #     xn = xn[1:]
        #     lenly=len(ly)
        # print(lenly)
        # print()
        

        
### STEP 4: Get the k vector using integer division and clean up
        L = 3*3600 # The period
        k = (xn)//(2*np.pi/L) # Convert to mHz, then divide by 2pi/L to get the k vector
        
        randvals = np.random.rand(lenly)*2*np.pi
        ascale = np.cos(randvals)*np.sqrt(ly)*0.01
        bscale = np.sin(randvals)*np.sqrt(ly)*0.01
        
        fakevals = ascale+1j*bscale
        
        
        
        # REMOVE DUPLICATE VALUES
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
        negativey = np.conj(np.append(np.flip(ndy[:-1],axis=0),extrazeros[1:]))  ### issue of meeting in the middle +/-
        ynew2=np.append(positivey,negativey)
        
        
        
        #ynew2 = 1/2*ynew2
        
        
        #input(fakevals*np.conj(fakevals)-ly*0.01**2)
        print('maxes')
        print(max(np.abs(ynew2)),max(np.sqrt(ly)*0.01))
      

        
        
        extrak = np.arange(0,ndk[0])
        k1 = np.append(np.append(np.append(extrak,ndk),-np.flip(ndk[:-1],axis=0)),-np.flip(extrak,axis=0))[:-1]
    
        # print(ndk,ndy)
        # print()
        
        for r in range(len(k1)):
            print(k1[r],ynew2[r])
            
        print(len(k1),len(ynew2))
      
        
        bbb = ifft(ynew2)
        timex = np.linspace(0,100,len(bbb))
    
  
        f,a=plt.subplots(3)
        
        a[0].plot(k1,np.real(ynew2),'.')
        a[1].plot(k1,np.imag(ynew2),'.')
        a[2].plot(k1,np.abs(ynew2),'.',markersize=5)
        a[2].plot(k,np.sqrt(ly)*0.01,'.')
        
        
        g,b=plt.subplots(3)
        b[0].plot(timex,np.real(bbb))
        b[1].plot(timex,np.imag(bbb))
        b[2].plot(timex,np.abs(bbb))
        plt.show()
        


        #input(ndy)
                
        
        # ndk = np.array(())
        #ndy = np.array(())
        # for fi in range(len(k)):
        #     num = k[fi]
        #     if num not in ndk:
        #         ndk = np.append(ndk,num)
        # 
        # lll =[]
        # for h in ndk:
        #     l1=np.where(k==h)[0]
        #     lll.append(l1)
        # 
        # 
        # for ar in lll:
        #     val = np.mean(fakevals[ar])
        # 
        #     ndy=np.append(ndy,val)
            

        
        # plt.plot(ndk,ndy)
        # plt.plot(k,fakevals)
        # plt.show()
        # 
        
            
### STEP 6: Redefine the k vector so that the carrier mode is at 0 (factor it out)
        knew = ndk-ndk[loc]
        xnew = x-x[loc]
        
#         
#         
# ### STEP 7: Get nondimensionalization constants    
#         g = 9.81 #(m/s^2)
#         if n==0:
#             w0 = (2*np.pi)**2/L*ndk[loc] # Get the value from the integer
#             k0 = w0**2/g # The carrier wavenumber
#             m = max(ndy)
#             epsilon = 2*m*k0 # The nondimensionalization constant epsilon
#             print(f,'Special Values')
#             print('Maximum value',m)
#             print('Carrier frequency',w0)
#             print('Wavenumber',k0)
#             print('epsilon',epsilon)
#             print('period',L)
#             print()
#             
#         n = n+1
#             
#           
#         
#         
#         
#         
#         knew = ndk-ndk[loc]
#         xnew = x-x[loc]
#         
# 
# ### STEP 9 (optional): Restrict the frequency range
#         #fakevals[:loc-10]=0+0j
#         #fakevals[loc+10:]=0+0j
#         
#         
#         
# ### STEP 10: Get B, the new surface
#         NNN = 4096 # New number of grid points
#         tnew = np.linspace(0,L-L/NNN,NNN) # new temporal x axis
#         
#         # Find B
#         fakevals = np.append(np.flip(fakevals,axis=0),fakevals)
#         ndk = np.append(-np.flip(ndk,axis=0),ndk)
#         B=np.zeros(NNN,dtype=complex)
#         for v in range(len(fakevals)):
#             #newvalue = 1/N*fakevals[v]*np.exp(2j*np.pi*knew[v]/L*tnew)
#             newvalue = 1/N*fakevals[v]*np.exp(2j*np.pi*ndk[v]/L*tnew)
#             B= B+newvalue
#         
#         
#         # test = np.append(np.flip(fakevals,axis=0),fakevals)
# 
#         ggg,hhh=plt.subplots(2,1)
#         hhh[0].plot(tnew,np.real(B))
#         hhh[1].plot(tnew,np.imag(B))
#         plt.show()
#         
#         
#         
#         
#         
#         # print(max(np.abs(1/NNN*fft(B)))-max(1/N*np.abs(fakevals)))
#         #figg,axx = plt.subplots()
#         # axx.set_title(f)
#         #axx.plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.')
#         # axx.plot(knew[np.where(np.abs(fakevals)==max(np.abs(fakevals)))],max(1/N*np.abs(fakevals)),'*',markersize=10)
#         #axx.plot(knew,1/N*np.abs(fakevals),'.')
# 
#         # axx.plot(knew,1/N*ndy,'^',markersize='0.25')
#         # 
#         #plt.show()
#         
#         # SAVE
#         np.savetxt(sd + '/Processed/'+f[-10:], np.transpose(np.append(tnew,B,axis=0)).view(float))
#         
#         
#      
# ### STEP 11: Nondimensionalize values, plot, and save results
#         ndB = k0/epsilon*B # new "y"
#         xi = w0*epsilon*tnew # new "x"
#         
#         # Fix up xi so that it is periodic and has the correct number of points to match B
#         xi_L = (xi[-1]-xi[0])+(xi[1]-xi[0])
#         xi_new = np.linspace(0,xi_L-xi_L/NNN,NNN,dtype=complex)
#         
#         chi = epsilon**2*k0*distv # new "t"
#     
#         #Plot results
#         ax1[pi].set_title(r'$\chi$'+' = '+str(chi[pi]))
#         ax1[pi].plot(xi_new,np.real(ndB))
#         ax1[pi].plot(xi_new,np.imag(ndB))
#         pi+=1
#         
#         # Save Results
#         np.savetxt(sd + '/NonDim Data/ND'+f[-10:],np.append([xi_new],[ndB],axis=0).view(float))
#         
#     
#     np.savetxt(sd + '/chi.txt', chi.view(float))
#     fig1.tight_layout()
#     fig1.subplots_adjust(top=0.88)
#     #plt.show()
#     plt.savefig(sd + '/NonDim Figs/allgauges.png',dpi=500)