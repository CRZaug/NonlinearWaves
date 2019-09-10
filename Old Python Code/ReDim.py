"""

~~~ REDIMENSIONALIZE THE SIMULATION AND DATA FILES ~~~

Re dimensionalizing and plotting the CQ Data (sims AND experimental)
will be done in a few steps

1. Read in the 3 experimental data files
2. Read in the 3*4 = 12 simulation data files
3. Dimensionalize P, M, wp, wm, and sb for both the data and the simulations
4. Use the old plotting program to plot everything together in 2 figures.
5. Save the results to data files.

"""


import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft

# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))


### STEP 1: READ IN THE EXPERIMENTAL DATA FILES

# Define the dictionaries
P = {}
M = {}
PM = {}
sb = {}
wp = {}
MasterDict = [M,P,PM,sb,wp]

# Start reading in the data
whichset = 'JulyData/'
dir = whichset+'Data CQs/NonDim CQ Values'
files = listdirNH(dir)
key1 = 'Data'
l = 0
for k in files:
    Dict = MasterDict[l]
    data = np.transpose(np.loadtxt(k).view(float))
    Dict[key1]=data
    l=l+1



### STEP 2: READ IN THE SIMULATION DATA
dir = whichset+'Simulations/Conserved Quantities'
key2 = ['dNLS CQ', 'Dysthe CQ', 'NLS CQ', 'vDysthe CQ']
dk = 0
for subdir in key2:
    files = listdirNH(dir+'/'+subdir)
    l = 0
    for k in files:
        Dict = MasterDict[l]
        data = np.transpose(np.loadtxt(k).view(float))
        Dict[key2[dk]]=data
        l=l+1
    dk = dk+1

### STEP 3: DIMENSIONALIZE THE DATA

# Define dimensionalization constants
g = 9.81
epsilon = 7.84869668208853e-05
w0 =   0.05864306286700947
k0 = w0**2/g


# Dim P
dim_P = {}
for key in P:
    ent = P[key]
    xi = ent[0]
    p = ent[1]
    x = xi/(epsilon**2*k0)
    dim_p = epsilon**3*w0/k0**2*p
    dim_P[key] = np.append([x],[dim_p],axis = 0)
    
# Dim M
dim_M = {}
for key in M:
    ent = M[key]
    xi = ent[0]
    m = ent[1]
    x = xi/(epsilon**2*k0)
    dim_m = (epsilon/k0)**2*m
    dim_M[key] = np.append([x],[dim_m],axis = 0)

# Dim PM
dim_PM = {}
for key in PM:
    ent = PM[key]
    Pent = dim_P[key]
    Ment = dim_M[key]
    m = Ment[1]
    p = Pent[1]
    xi = ent[0]
    x = xi/(epsilon**2*k0)
    dim_pm = (p/m)*1000 # Gives mHz
    dim_PM[key] = np.append([x],[dim_pm],axis = 0)

# Dim sb
dim_sb = {}
for key in sb:
    ent = sb[key]
    xi = ent[0]
    x = xi/(epsilon**2*k0)
    sbv = np.zeros((len(ent)-1,len(xi)))
    for j in range(1,len(ent)):
        sideband = ent[j]
        dim_sideband = (epsilon/k0)*sideband
        sbv[j-1]=dim_sideband
    dim_sb[key] = np.vstack([x,sbv])
    

# Dim wp
dim_wp = {}
for key in wp:
    ent = wp[key]
    xi = ent[0]
    peak = ent[1]
    x = xi/(epsilon**2*k0)
    dim_peak = peak+w0*1000 # Gives mHz
    dim_wp[key] = np.append([x],[dim_peak],axis = 0)



### STEP 4: PLOT THE RESULTS

# Initialize for plotting
plotter1 = [dim_M,dim_P,dim_PM,dim_wp]
key2[:0] = [key1]
titles1 = ['CQ M', 'CQ P', r'$\omega_m$', r'$\omega_p$']
titles2 = np.loadtxt(os.getcwd()+'/'+whichset+'sidebandnums.txt').view(float)
y1 = ['M (m'+r'$^2$'+')','P (m'+r'$^2$'+'/s)',r'$\omega_m$'+' (mHz)',r'$\omega_p$'+' (mHz)']
#y2 = [r'$|a_{-3}|$'+' (m)',r'$|a_{-2}|$'+' (m)',r'$|a_{-1}|$'+' (m)',r'$|a_0|$'+' (m)',r'$|a_1|$'+' (m)',r'$|a_2|$'+' (m)',r'$|a_3|$'+' (m)']
disp = ['.k',':m','-.g','--r','-c']
sizes = [13,1,1,1,1]

# Begin plotting
fig1, ax1 = plt.subplots(4,1,figsize = (11,6.5))
fig1.suptitle('Quantities of Interest',fontsize=16)

dispind = 0
for key in key2:
    ax1 = ax1.flatten()
    for i in range(len(plotter1)):
        dict = plotter1[i]
        VALUES = dict[key]
        x = VALUES[0]
        y = VALUES[1]
        ax1[i].plot(x,y,disp[dispind],markersize = sizes[dispind])
        ax1[i].set_title(titles1[i])
        ax1[i].set_ylabel(y1[i])
        ax1[i].set_xlabel('Location (m)')
        ax1[i].ticklabel_format(style='sci',scilimits=(-1,1),axis='both')
    dispind += 1
    ax1[0].legend(key2,bbox_to_anchor=(1, 1))
    fig1.tight_layout()
    fig1.subplots_adjust(top=0.88)
plt.savefig(whichset+'Final Figures/CQResultFig.png',dpi=500)

fig2, ax2 = plt.subplots(len(titles2),sharex=True,figsize = (7,1.625*len(titles2)))
fig2.suptitle('Select Fourier Amplitudes',fontsize=16)
dispind = 0
for key in key2: 
    sbvals = dim_sb[key]
    x = sbvals[0,:]
    sideband7=np.delete(sbvals, 0, 0)
    for po in range(len(titles2)):
        ax2[po].plot(x,sideband7[po],disp[dispind],markersize = sizes[dispind])
        ax2[po].set_ylabel('a'+ r'$_{'+str(int(titles2[po]))+'}$')
    fig2.tight_layout()
    fig2.subplots_adjust(top=0.97)
    dispind += 1
plt.savefig(whichset+'Final Figures/FAResultFig.png',dpi=500)


### STEP 5: SAVE THE RESULTS

# Save P, M, wp, wm
md =[dim_P,dim_M,dim_PM,dim_wp]
val = ['dimP','dimM','dimPM','dimwp','dimsb']
#key2 = ['dNLS CQ', 'Dysthe CQ', 'NLS CQ', 'vDysthe CQ']
o = 0
for cqval in md:
    # Save the Data Values
    for ky in cqval:
        if ky == 'Data':
            np.savetxt(whichset+'Data CQs/Dim CQ Values/'+val[o]+'.txt',np.transpose(cqval[ky]).view(float))
        else:
            np.savetxt(whichset+'Simulations/Dimensional Results/'+str(ky)[:-3]+' dimCQ/'+val[o]+'.txt',np.transpose(cqval[ky]).view(float))
    o=o+1

# Save sidebands
for ky in dim_sb:
    if ky == 'Data':
        np.savetxt(whichset+'Data CQs/Dim CQ Values/'+val[-1]+'.txt',np.transpose(dim_sb[ky]).view(float))
    else:
        np.savetxt(whichset+'Simulations/Dimensional Results/'+str(ky)[:-3]+' dimCQ/'+val[-1]+'.txt',np.transpose(dim_sb[ky]).view(float))
