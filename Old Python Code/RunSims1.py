"""

~~~ RUNS SIMULATIONS OF THE DIMENSIONLESS PROCESSED INITIAL DATA ~~~

This code reads in the dimensionless initial data and runs simulations of NLS, Dysthe, dNLS, and vDysthe.

1. Get the dimensionless chi values, define a new time vector, and save it
2. Get the dimensionless xi and B values
3. Run simulations and save data to a grid using already found parameters
4. Save data to text files


"""
import os
import numpy as np
import matplotlib.pyplot as plt
import NLS
import dNLS
import Dysthe as dy
import vDysthe as vdy
from numpy.fft import fft, ifft



### STEP 1: Read in dimensionless data to create time values

# Choose the name of the file the data will be pulled from
whichset = 'Aug2Data/'
dir = whichset+'NonDim Data'


Num_O_Times = 100
times  = np.loadtxt(whichset+'chi.txt').view(float)
simtimes = np.linspace(times[0],times[-1],Num_O_Times) # The time vector: the times the code will step out to


# Save time vector
np.savetxt(whichset+'Simulations/SimTime.txt', simtimes.view(float))



### STEP 2: Get the xi and B initial data

# The initial data
IC = np.loadtxt(dir+'/NDgauge2.out').view(complex)
xspace = IC[0] # xi
#u0 = IC[1] # B
u0 = max(IC[1])*np.sin(xspace)
# plt.plot(xspace,np.real(u0),xspace,np.imag(u0))
# plt.show()


### STEP 3: Run simulations and save data to grid

# Enter parameters already in ProcessNonDim and Data_special_vals1
epsilon =  0.006482153402815046
Del = 0.07007178564680976



# Collect values needed to perform operator splitting
L = xspace[-1]-xspace[0]+(xspace[2]-xspace[1]) # Since x[-1] isn't L, add on deltax to get to L
gridnum = len(xspace)
k,expconsts = dy.kvec(gridnum,L)
starttime = simtimes[0]
stoptime = simtimes[-1]

rk4steps = 1

# operator splitting parameters


# Run simulations 
r_NLS = np.zeros((len(simtimes),gridnum),dtype = complex)
r_dNLS = np.zeros((len(simtimes),gridnum),dtype = complex)
r_Dysthe = np.zeros((len(simtimes),gridnum),dtype = complex)
r_vDysthe = np.zeros((len(simtimes),gridnum),dtype = complex)

r_NLS[0,:] = u0
r_dNLS[0,:] = u0
r_Dysthe[0,:] = u0
r_vDysthe[0,:] = u0


for t in range(1,Num_O_Times):
    steps = t
    stoptime = simtimes[t]
    deltat = (stoptime-starttime)/steps
    print(steps,stoptime,deltat)
    sim_NLS = NLS.sixth(u0,deltat,steps,expconsts)
    sim_dNLS = dNLS.sixth(u0,deltat,steps,expconsts,Del)
    sim_Dysthe = dy.sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon)
    sim_vDysthe = vdy.sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del)
    
    r_NLS[t,:] = sim_NLS
    r_dNLS[t,:] = sim_dNLS
    r_Dysthe[t,:] = sim_Dysthe
    r_vDysthe[t,:] = sim_vDysthe
    


    

### STEP 3: Save the data to text files (make sure the director is correct)

p=0
for data in range(1,Num_O_Times):
    NLSd = r_NLS[data,:]
    dNLSd = r_dNLS[data,:]
    Dysthed = r_Dysthe[data,:]
    vDysthed = r_vDysthe[data,:]
    # np.savetxt(whichset+'Simulations/NLS Sim/SimNLS'+str(simtimes[data])+'.txt', NLSd.view(float))
    # np.savetxt(whichset+'Simulations/dNLS Sim/SimdNLS'+str(simtimes[data])+'.txt', dNLSd.view(float))
    # np.savetxt(whichset+'Simulations/Dysthe Sim/SimDysthe'+str(simtimes[data])+'.txt', Dysthed.view(float))
    # np.savetxt(whichset+'Simulations/vDysthe Sim/SimvDysthe'+str(simtimes[data])+'.txt', vDysthed.view(float))
    np.savetxt(whichset+'Simulations/NLS Sim/SimNLS'+str(p)+'.txt', NLSd.view(float))
    np.savetxt(whichset+'Simulations/dNLS Sim/SimdNLS'+str(p)+'.txt', dNLSd.view(float))
    np.savetxt(whichset+'Simulations/Dysthe Sim/SimDysthe'+str(p)+'.txt', Dysthed.view(float))
    np.savetxt(whichset+'Simulations/vDysthe Sim/SimvDysthe'+str(p)+'.txt', vDysthed.view(float))
    p=p+1