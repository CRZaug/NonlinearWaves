"""

~~~ FIND CONSERVED QUANTITIES AND SIDEBANDS OF THE DIMENSIONLESS INITIAL DATA ~~~

This code reads in the dimensionless initial data, finds the values M, P, omega_m, omega_p, and the sidebands.

1. Read in the dimensionless data
2. Finds the locations of the sidebands by looking at gauge 50 (the first gauge)
3. Finds M, P, omega_m, omega_p, and sideband values at each gauge (each chi value)
4. Saves the data
5. Plots the values against chi
6. Fits the decay of M to find delta (also get the error and plot the results)


"""


import NLS
import numpy as np
from numpy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
import glob
import os

# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))


### STEP 1: Read in the xi, B, and chi data

whichset  = 'JulyData/'


# Choose the name of the file the data will be pulled from
dir = whichset+'NonDim Data'
dirnames = listdirNH(dir)
dirlength = len(os.listdir(dir))


tvector = np.loadtxt(whichset+'chi.txt').view(float) #times
values = {} # The x and y values at each location
for name  in dirnames:

    vdatavals = np.loadtxt(name).view(complex)# xi and B

    title = name[-10:-4]
    print(title)
    N = len(vdatavals[0])

    # Save values to dictionaries
    values[title]=np.append([vdatavals[0]],[vdatavals[1]],axis=0) 



### STEP 2: Find the sideband values and the carrier wave location

# Find the sideband values to track through the rest of the program
x=values['gauge2'][0]
y=values['gauge2'][1]

# Perform an FFT of the y values
yhat =fft(y) # Fourier amplitudes
yhat1 = 1/N*np.abs(yhat) # Normalized fourier amplitudes


# Define some constants/constant vectors
L = (x[-1]-x[0])+(x[1]-x[0]) # The period
k1=NLS.kvec(N)


k=1/L*k1

#sv = np.array([-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12])

# Find max Fourier amplitude and location
m = max(yhat1) # The max amplitude (m)
i = np.where(yhat1 == m)
print(i)
i = i[0][0]

#print(i)
carrier = k[i]

# plt.plot(k,yhat1,'.')
# plt.plot(carrier,m,'.',markersize = 5)
# plt.show()


#Shift over the sideband vector so that the carrier wave is in the middle (carrier is 0 so this doesn't really matter)
sv = np.array([-3,-2,-1,0,1,2,3]) # The sideband vector


plt.plot(k1,yhat1,'.')
plt.plot(k1[sv],yhat1[sv],'.',markersize=10)
plt.show()


sv=sv+i

#print(i,sv,k1[i])


#print('sb',k[sv])
#print(m,carrier,L)



### STEP 3: Find P, M, omega_m, omega_p, and the sidebands at each gauge location (each xi)

# Preallocate space for the data values to be stored
Pvals = np.zeros(len(tvector),dtype = float)
Mvals = np.zeros(len(tvector),dtype = float)
PMvals = np.zeros(len(tvector),dtype = float)
wpeak = np.zeros(len(tvector),dtype = float)
sideband7 = np.zeros((len(tvector),7),dtype = float)
j= 0
for n in values:
    
    x=values[n][0]
    y=values[n][1]
    
    
    # Define some constants/constant vectors
    L = (x[-1]-x[0])+(x[1]-x[0])
    # Perform an FFT of the y values
    yhat =fft(y) # Fourier amplitudes
    yhat1 = 1/N*np.abs(yhat) # Normalized fourier amplitudes
    
    # Find max Fourier amplitude and location
    m = max(yhat1) # The max amplitude (m)
    i = np.where(yhat1 == m)
    i = i[0][-1]

    
    carrier = k[i] # Find the location of the peak frequency
    
    
    P = NLS.CQ_P1(y,L,k1) # Find P
    M = NLS.CQ_M(y,L) # Find M
    PM = P/M # Find omega_m
    

    sidebands = yhat1[sv] # Find fourier amplitudes at sidebands
    #plt.plot(k,yhat1,'.')
    #plt.plot(k[sv],sidebands,'.',markersize = 10)
    #plt.show()
    
 
    Pvals[j] = np.real(P)
    Mvals[j] = np.real(M)
    PMvals[j] = np.real(PM)
    wpeak[j] = np.real(carrier)
    sideband7[j]=np.real(sidebands)
    
    j=j+1
    
print('Pvals',Pvals)
print('Mvals',Mvals)
    
    
### STEP 4: Save the Data

dir = whichset+'Data CQs/NonDim CQ Values/'
savedata = np.append([tvector],[Pvals],axis=0)
np.savetxt(dir+'P'+'.txt', np.transpose(savedata).view(float))

savedata = np.append([tvector],[Mvals],axis=0)
np.savetxt(dir+'M'+'.txt', np.transpose(savedata).view(float))

savedata = np.append([tvector],[PMvals],axis=0)
np.savetxt(dir+'PM'+'.txt', np.transpose(savedata).view(float))

savedata = np.append([tvector],[wpeak],axis=0)
np.savetxt(dir+'wp'+'.txt', np.transpose(savedata).view(float))

savedata = np.insert(sideband7, 0, tvector, axis=1)
np.savetxt(dir+'sb'+'.txt', savedata.view(float))



### STEP 5: Plot the values at each time chi

plotem = 1
if plotem == 1:

    # Plotting vectors
    fig1, ax1 = plt.subplots(2,2,figsize = (10,6.5))
    fig2, ax2 = plt.subplots(7,sharex=True,figsize = (10,6.5))
    fig1.suptitle('Quantities of Interest',fontsize=16)
    fig2.suptitle('Select Fourier Amplitudes',fontsize=16)
    
    plotter1 = [Pvals,Mvals,PMvals,wpeak]
    titles1 = ['CQ P', 'CQ M', r'$\omega_m$', r'$\omega_p$']
    
    titles2 = [r'$a_{-3}$',r'$a_{-2}$',r'$a_{-1}$',r'$a_0$',r'$a_1$',r'$a_2$',r'$a_3$']
    
    # CQ figure
    ax1 = ax1.flatten()
    for i in range(len(plotter1)):
        ax1[i].plot(tvector,plotter1[i],'.',markersize =15)
        ax1[i].set_title(titles1[i])
        ax1[i].set_xlabel(r'$\chi$')
    vp = sideband7[:,3]
    
    # Sideband figure
    ax2[3]=plt.subplot2grid((4,2),(0,0),colspan=2,fig=fig2)
    ax2[3].plot(tvector,vp,'.',markersize =10)
    ax2[3].set_title(titles2[3])
    for i in range(0,3):
        vp = sideband7[:,i]
        print(vp)
        ax2[i]=plt.subplot2grid((4,2),(3-i,0),fig=fig2)
        ax2[i].plot(tvector,vp,'.', markersize =10)
        ax2[i].set_title(titles2[i])
        i2 = i+4
        vp = sideband7[:,i2]
        print(vp)
        ax2[i2]=plt.subplot2grid((4,2),(i+1,1),fig=fig2)
        ax2[i2].plot(tvector,vp,'.', markersize =10)
        ax2[i2].set_title(titles2[i2])
    
    # Fix configurations
    fig1.tight_layout()
    fig1.subplots_adjust(top=0.88)
    fig2.tight_layout()
    fig2.subplots_adjust(top=0.88)
    plt.show()



# STEP 6: Get a fit of M

fitx = -2*tvector # The x vals of the fit   
fity = np.log(Mvals) # The y vals of the fit

# Get the fit and define a new y vector
A = np.vstack([fitx, np.ones(len(fitx))]).T
m, b = np.linalg.lstsq(A, fity,rcond=-1)[0] # m is delta
fittedy = m*fitx+b
newy = np.exp(fittedy)
print(tvector)
print('delta ',m, 'b ', b)

# Get the error of the fit
error = np.sum(np.abs(newy-Mvals)**2)/len(newy)
print('error ', error)
 
# Plot results of the fit
fig3,ax3 = plt.subplots(2,1)
fig3.suptitle('Fitting M for Delta')
plt.text(.7, .7,'y='+str(m)[1:7]+'x'+str(b)[1:6])
ax3[0].plot(fitx,fity,'.',label = 'Linearized values')
ax3[0].plot(fitx,fittedy, label ='Linear fit')
ax3[0].set_xlabel('-2t')
ax3[0].set_ylabel('ln(M)')
ax3[0].legend(loc='upper left')

ax3[1].plot(tvector,Mvals,'.', label = 'Actual values')
ax3[1].plot(tvector,newy,label = 'Fit curve')
ax3[1].set_xlabel('t')
ax3[1].set_ylabel('M')
ax3[1].legend(loc='upper right')
fig3.tight_layout()
fig3.subplots_adjust(top=0.88)
plt.show()