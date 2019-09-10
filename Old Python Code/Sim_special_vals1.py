"""

~~~ FIND CONSERVED QUANTITIES AND SIDEBANDS OF THE SIMULATION RESULTS' ~~~

This code reads in the dimensionless initial data, finds the values M, P, omega_m, omega_p, and the sidebands.

1. Read in the dimensionless simulation results
2. Finds the locations of the sidebands by looking at gauge 50 (the first gauge)
3. Finds M, P, omega_m, omega_p, and sideband values at each gauge (each chi value)
4. Plots the values against chi
5. Saves the data


"""



import NLS
import numpy as np
from numpy.fft import fft, ifft, fftfreq
import matplotlib.pyplot as plt
import os
import glob

# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))



### STEP 1: Load in simulation data

# Load time values
whichset = 'JulyData/'
tvector = np.loadtxt(whichset+'Simulations/SimTime.txt').view(float)

# Choose the name of the file the data will be pulled from
masterdir = whichset+'Simulations/'
dir = ['dNLS Sim', 'Dysthe Sim','NLS Sim','vDysthe Sim']

NLSd = {} 
dNLSd = {}
Dysthed = {}
vDysthed = {}
Dictionaries = [dNLSd,Dysthed,NLSd,vDysthed] # Alphabetized

# Read in the intital data
IC = np.loadtxt(whichset+'NonDim Data/NDgauge2.out').view(complex)
x = IC[0]
y = IC[1]

h = 0
for m in dir:
    dict = Dictionaries[h]
    dirnames = listdirNH(masterdir+m)
    dirlength = len(dirnames)
    kk = 0
    for name in dirnames:
        if os.path.isfile(name) == True:
            title = tvector[kk]
            vdatavals = np.loadtxt(name).view(complex)
            N = len(vdatavals)
            dict[title]=np.append([x],[vdatavals],axis=0)
            kk=kk+1
    h=h+1



### STEP 2: Find the sideband values and the carrier wave location

# Perform an FFT of the y values
yhat =fft(y) # Fourier amplitudes
yhat1 = 1/N*np.abs(yhat) # Normalized fourier amplitudes

# Define some constants/constant vectors
L = x[-1]-x[0]+(x[1]-x[0])
k=NLS.kvec(N)
#sv = np.array([-3,-2,-1,0,1,2,3]) # The sideband vector
sv = [] # the sideband vector
for j in range(len(yhat1)):
    if yhat1[j]>0.00000000001:
        sv.append(j)

lll = len(sv)

# Find max Fourier amplitude and location
mt = max(yhat1) # The max amplitude (m)
i = np.where(yhat1 == mt)[0][0] ################################
carrier = 1/L*k[i]

### STEP 3: Find P, M, omega_m, omega_p, and the sidebands at each gauge location (each xi)

NLSCQ = {}
dNLSCQ = {}
DystheCQ = {}
vDystheCQ = {}

CQDict = [dNLSCQ, DystheCQ, NLSCQ, vDystheCQ]
keys = ['P', 'M', 'PM', 'wp', 'sb']
dname = ['dNLS CQ','Dysthe CQ','NLS CQ','vDysthe CQ']

cid = 0
for dict in Dictionaries:
    Pvals = np.zeros(len(tvector))
    Mvals = np.zeros(len(tvector))
    PMvals = np.zeros(len(tvector))
    wpeak = np.zeros(len(tvector))
    sideband7 = np.zeros((len(tvector),lll))
    j=0
    CQs = CQDict[cid]
    for n in dict:
        
        x=dict[n][0]
        y=dict[n][1]
        # Perform an FFT of the y values
        yhat =fft(y) # Fourier amplitudes
        yhat1 = 1/N*np.abs(yhat) # Normalized fourier amplitudes
        
        # Find max Fourier amplitude and location
        m = max(yhat1) # The max amplitude (m)
        i = np.where(yhat1 == m)
        
        carrier = 1/L*k[i] # omega_p
        P = NLS.CQ_P1(y,L,k) #P
        M = NLS.CQ_M(y,L) # M 
        PM = P/M # omega_m
        sidebands = yhat1[sv] # sidebands 
        
        Pvals[j] = np.real(P)
        Mvals[j] = np.real(M)
        PMvals[j] = np.real(PM)
        wpeak[j] = np.real(carrier[0])
        sideband7[j]=sidebands
        j=j+1
    
        valuevect = [Pvals,Mvals,PMvals,wpeak,sideband7]
    
    for val in range(len(keys)):
        CQs[keys[val]] = valuevect[val]

    # Get the number on each sideband for later labeling purposes
    svlabp = []
    iop  = 0
    for gg in range(len(sv)):
        if sv[gg] <N//2:
            svlabp.append(iop)
            iop +=1
    
    svlabn = []
    iop  = 1
    for gg in range(len(sv)):
        if np.flip(sv,axis =0)[gg] >N//2:
            svlabn.append(-iop)
            iop +=1
    svlab = np.append(svlabp,np.flip(svlabn,axis=0),axis=0) 
    
### STEP 4: Plot the results

    plotem = 1
    if plotem == 1:
        tn = dname[cid]
        
        # Plotting vectors
        fig1, ax1 = plt.subplots(2,2,figsize = (10,6.5))
        fig1.suptitle(tn+' Quantities of Interest',fontsize=16)

        
        plotter1 = [Pvals,Mvals,PMvals,wpeak]
        titles1 = ['CQ P', 'CQ M', r'$\omega_m$', r'$\omega_p$']
                
        ax1 = ax1.flatten()
        for i in range(len(plotter1)):
            ax1[i].plot(tvector,plotter1[i])
            ax1[i].set_title(titles1[i])
            ax1[i].set_xlabel(r'$\chi$')
        fig1.tight_layout()
        fig1.subplots_adjust(top=0.88)
        plt.savefig(whichset+'Simulations/CQ Figs/CQ '+tn+'.png')
        
        fig2, ax2 = plt.subplots(lll,sharex=True,figsize = (7,1.625*lll))
        fig2.suptitle(tn+' Select Fourier Amplitudes',fontsize=16)
        for po in range(lll):
            vp = sideband7[:,po]
            ax2[po].plot(tvector,vp)
            ax2[po].set_ylabel('a'+ r'$_{'+str(svlab[po])+'}$')
        
        
        fig2.tight_layout()
        fig2.subplots_adjust(top=0.97)
        plt.savefig(whichset+'Simulations/CQ Figs/Sidebands '+tn+'.png')


    cid = cid +1



### STEP 5: Save the Data

dval = 0
for dict in CQDict:
    dctnm = dname[dval]
    m=0
    for cons in dict:
        filename = whichset+'Simulations/Conserved Quantities/'+dctnm+'/'+cons+'.txt'

        consval =  dict[cons]

        if consval.shape==tvector.shape:
            savedata = np.append([tvector],[consval],axis=0)
            np.savetxt(filename, np.transpose(savedata).view(float))
        else:
            savedata = np.insert(consval, 0, tvector, axis=1)
            np.savetxt(filename, savedata.view(float))
        m=m+1
        
    dval = dval +1