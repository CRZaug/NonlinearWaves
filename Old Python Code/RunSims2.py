import os
import numpy as np
import matplotlib.pyplot as plt
import NLS
import dNLS
import Dysthe as dy
import vDysthe as vdy
from numpy.fft import fft, ifft

# Define master directory

whichset = 'JulyData/'

# Read in x and y data
IC = np.loadtxt(whichset+'NonDim Data/NDgauge2.out').view(complex)
xspace = IC[0]
u0 = IC[1]

# Read in time data
times = np.loadtxt(whichset+'chi.txt').view(float)
starttime = times[0]
stoptime = times[-1]
num_o_times = 300
simtimes = np.linspace(starttime,stoptime,num_o_times)

# Save time data
np.savetxt(whichset+'/Simulations/SimTime.txt',simtimes.view(float))

# Set operator splitting parameters
L = xspace[-1]-xspace[0]+(xspace[1]-xspace[0])
gridnum = len(xspace)
k, expconsts = vdy.kvec(gridnum,L)
epsilon = 7.84869668208853e-05
Del =  58884.92711968474
per = 3/32
rk4steps = 1

# Perform operator splitting
r_NLS = np.zeros((num_o_times,gridnum),dtype=complex)
r_dNLS = np.zeros((num_o_times,gridnum),dtype=complex)
r_Dysthe = np.zeros((num_o_times,gridnum),dtype=complex)
r_vDysthe = np.zeros((num_o_times,gridnum),dtype=complex)

r_NLS[0,:]=u0
r_dNLS[0,:]=u0
r_Dysthe[0,:]=u0
r_vDysthe[0,:]=u0

for t in range(1,num_o_times):
    steps = t*1
    endtime = simtimes[t]
    deltat = (endtime-starttime)/steps
    #print(steps, endtime, deltat)
    
    sim_NLS = NLS.sixth(u0,deltat,steps,expconsts,per)
    sim_dNLS = dNLS.sixth(u0,deltat,steps,expconsts,Del,per)
    sim_Dysthe = dy.sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon,per)
    sim_vDysthe = vdy.sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del,per)
    
    if np.isnan(sim_vDysthe[0]):
        os.system('say "Error"')
        input()
    
    r_NLS[t,:]=sim_NLS
    r_dNLS[t,:]=sim_dNLS
    r_Dysthe[t,:]=sim_Dysthe
    r_vDysthe[t,:]=sim_vDysthe

# Save data
for s in range(num_o_times):
    if s < 10:
        ss = '00'+str(s)
    elif 9<s<100:
        ss='0'+str(s)
    else:
        ss = str(s)
    np.savetxt(whichset+'Simulations/NLS Sim/simNLS'+ss+'.txt',r_NLS[s].view(float))
    np.savetxt(whichset+'Simulations/dNLS Sim/simdNLS'+ss+'.txt',r_dNLS[s].view(float))
    np.savetxt(whichset+'Simulations/Dysthe Sim/simDysthe'+ss+'.txt',r_Dysthe[s].view(float))
    np.savetxt(whichset+'Simulations/vDysthe Sim/simvDysthe'+ss+'.txt',r_vDysthe[s].view(float))