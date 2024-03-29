"""
This code performs recurring functions necessary for NLS analysis.

It has linear and nonlinear subroutines.
It has 1st, 2nd, and 6th order operator splitting.
It has a function to create intial data (sech and cn) and the period for cn.
It has an error analysis function.
It has a plot function which plots real and imaginary parts of a profile as well as the magnitude.

When using this module, it is useful to define the following in the script:
- gridnum (number of grid points in x)
- kappa (a constant for the initial profile, both sech and cn, try 1.0/3.0)
- k0 (a constnat for cn initial profile, 0<k0<1)
- L (the period. Define L for sech, and use cn_L for Cn)
- starttime (float)
- endtime (float)



"""

import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

#______________________ LINEAR _______________________________
def linear(initial_condition,deltat,expconsts,Del,per): 
    uhat0 = fft(initial_condition)
    uhatnew = uhat0*np.exp(deltat*(1j*expconsts**2-Del))
    
    N = len(initial_condition)
    A = int(per*N/2)
    uhatnew[N//2-A:N//2] = np.zeros(A)

    return ifft(uhatnew)


#______________________NONLINEAR_________________________________
def nonlinear(initial_condition,deltat):
    return initial_condition*np.exp(4j*np.conj(initial_condition)*initial_condition*deltat)


#__________ FIRST ORDER OPERATOR SPLITTING _______________________
def first(initial_condition,deltat,steps,expconsts,Del,per):
    
    unew=initial_condition
    for j in range(0,steps):
        unew=linear(unew,deltat,expconsts,Del,per)
        unew=nonlinear(unew,deltat)
    return unew


#______________SECOND ORDER OPERATOR SPLITTING__________________
def second(initial_condition,deltat,steps,expconsts,Del,per):
    # Preparation before the operator splitting
    sdt = deltat*0.5 # Creates a new timestep size
    # The first part
    unew = linear(initial_condition,sdt,expconsts,Del,per) # linear part
    unew = nonlinear(unew,deltat) # nonlinear part, which outputs unew to be used in the next loop
    
    for i in range(0,steps-1): # The second part; the bulk of the loop (includes every timestep except first and last)
        unew = linear(unew, deltat,expconsts,Del,per) # Uses unew from previous loop 
        unew = nonlinear(unew,deltat) # Outputs a new unew
    ufinal = linear(unew,sdt,expconsts,Del,per) # This is the final result

    return ufinal


#________________SIXTH ORDER OPERATOR SPLITTING___________________
def sixth(initial_condition,deltat,steps,expconsts,Del,per): # Take in the initial space profile and timevector 
    # Preparation for the operator splitting
        
    w3 = 0.7845136104775600
    w2 = 0.2355732133593570
    w1 = -1.177679984178870
    w0 = 1.3151863206839060
    
    lt = deltat*0.5
    nlt = deltat
    
    
    # Begin operator splitting
    unew= linear(initial_condition,lt*w3,expconsts,Del,per) # linear part, the first step
    
    for i in range(0,steps-1): # The second step, loop through everything except the first and last timestep
        unew = nonlinear(unew, nlt*w3)
        unew = linear(unew, lt*(w3+w2),expconsts,Del,per)
        unew = nonlinear(unew, nlt*w2)
        unew = linear(unew, lt*(w1+w2),expconsts,Del,per)
        unew = nonlinear(unew,nlt*w1)
        unew = linear(unew, lt*(w1+w0),expconsts,Del,per)
        unew = nonlinear(unew, nlt*w0)
        unew = linear(unew, lt*(w0+w1),expconsts,Del,per)
        unew = nonlinear(unew, nlt*w1)
        unew = linear(unew, lt*(w1+w2),expconsts,Del,per)
        unew = nonlinear(unew, nlt*w2)
        unew = linear(unew, lt*(w2+w3),expconsts,Del,per)
        unew = nonlinear(unew, nlt*w3)
        unew = linear(unew, nlt*w3,expconsts,Del,per) # Note different timestep
        
    # The final step: perform the operatons one last time and output the result    
    unew = nonlinear(unew, nlt*w3)
    unew = linear(unew, lt*(w3+w2),expconsts,Del,per)
    unew = nonlinear(unew, nlt*w2)
    unew = linear(unew, lt*(w1+w2),expconsts,Del,per)
    unew = nonlinear(unew,nlt*w1)
    unew = linear(unew, lt*(w1+w0),expconsts,Del,per)
    unew = nonlinear(unew, nlt*w0)
    unew = linear(unew, lt*(w0+w1),expconsts,Del,per)
    unew = nonlinear(unew, nlt*w1)
    unew = linear(unew, lt*(w1+w2),expconsts,Del,per)
    unew = nonlinear(unew, nlt*w2)
    unew = linear(unew, lt*(w2+w3),expconsts,Del,per)
    unew = nonlinear(unew, nlt*w3)
    ufinal = linear(unew, lt*w3,expconsts,Del,per) # This is the final result
    return ufinal


#________________________K vector_____________________________
def kvec(gridnum):
    split = int(gridnum/2)
    k = np.zeros(gridnum)
    k[0:split+1] = np.arange(0,split+1,1)
    k[split+1:gridnum] = np.arange(-split+1,0,1)
    return k


#___________________INITIAL DATA _____________________________
def dNLS_id(gridnum,u0,k0,zeta,L):
    xspace = np.linspace(-L/2.0, L/2.0-L/gridnum, gridnum)
    initialprofile = u0*np.exp(1j*k0*xspace+1j*zeta)
    k = kvec(gridnum)
    expconsts=1j*(2.0*np.pi*k/L)**2
    return (xspace,initialprofile,expconsts)


#___________________ SOLUTION _____________________________
def dNLSsoln(xspace,u0,k0,zeta,Del,t,L):
    return u0*np.exp(1j*k0*xspace+1j*(-t*k0**2+(2*u0**2/Del)*(1-np.exp(-2*t*Del)))-t*Del+1j*zeta)


#______________________OPERATOR: Mth DERIVATIVE________________
def der_fft(u,L,k,n=1):
    # n is the number of derivatives to take
    consts = 2*np.pi*1j*k/L
    return ifft(consts**n*fft(u))


#____________________OPERATOR: INTEGRATE__________________________
def trapezoid(u,L):
    return (L/len(u))*np.sum(u)


#_______________________CONSERVED QUANTITY: M______________________
def CQ_M(u,L):
    return 1/L*trapezoid(u*np.conj(u),L) # u should be defined along one period


#______________________CONSERVED QUANTITY: P1________________________
def CQ_P1(u,L,k): # more efficient
    u_x = der_fft(u,L,k,1)
    ccu_x = np.conj(u_x)
    return -1/(L)*trapezoid(np.imag(u*ccu_x),L)


#______________________CONSERVED QUANTITY: P2________________________
def CQ_P2(u,L,k):
    u_x = der_fft(u,L,k,1) # Derivative of u
    ccu_x = np.conj(u_x) # Complex conjugate of the derivative of u
    ccu = np.conj(u) # Complex conjugate of u
    return 1j/(2*L)*trapezoid((u*ccu_x)-(u_x*ccu),L)


#_________________________MAGNITUDE ERROR___________________________
def mag_er(an_ufinal, ufinal):
    return np.sum(np.abs(an_ufinal*np.conj(an_ufinal)-(ufinal*np.conj(ufinal)))**2)/len(ufinal)


#__________________________REAL/IMAG ERROR__________________________
def ri_er(an_ufinal,ufinal): #At any given time t
    N = len(ufinal)
    ruf = np.real(ufinal)
    iuf = np.imag(ufinal)
    ranuf = np.real(an_ufinal)
    ianuf = np.imag(an_ufinal)
    realerror = np.sum(np.abs(ranuf-ruf)**2)/N
    imagerror = np.sum(np.abs(ianuf-iuf)**2)/N
    return (realerror,imagerror)


#_______________________SUBPLOTS_________________________________
def rim_plt(xspace,initialprofile,ufinal,showplot=1):
    # real, imaginary, magnitude plot
    
    final_mag = ufinal*np.conj(ufinal)
    
    fig, (ax1, ax2) = plt.subplots(2)

    ax1.plot(xspace,np.imag(ufinal),'r',label = 'Imaginary')
    ax1.plot(xspace,np.real(ufinal),'b',label = 'real')
    ax1.set(title='Real and Imaginary Parts')
    ax1.legend(loc="upper right")
    
    ax2.plot(xspace,initialprofile**2,label = 'Initial Profile')
    ax2.plot(xspace, final_mag,'.',label = 'Final Profile')
    ax2.set(title='Initial and Final Magnitudes')
    ax2.legend(loc="upper right")
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
    if showplot == 1:
        plt.show()



#_______________________RUN SIMULATIONS_________________________________
def rundNLS(PARAMS,simtimes,u0,expconsts,Del,per):
    # PARAMS = [whichset,num_o_times,starttime,endtime,rk4steps,gridnum]
    # whichset is the directory to put into
    # num_o_times is an integer telling how many data points will be produced
    # starttime is ususally 0, the beginning of the simulation
    # endtime is the last time to step out to
    # rk4steps is an integer telling how many steps of rk4 to perform in the nonlinear part (usually 1)
    # gridnum is the x spatial resolution
    
    whichset = PARAMS[0]
    num_o_times = PARAMS[1]
    starttime = PARAMS[2]
    rk4steps = PARAMS[3]
    gridnum = PARAMS[4]
    
    # Initialize the matrix
    r_dNLS = np.zeros((num_o_times,gridnum),dtype=complex)
    r_dNLS[0,:]=u0
    
    # During iteration, deltat (the size of the step remains the same)
    # The initial profile u0 is always the input
    # The final time to step out increases in each loop, thus the number of steps increases in each iteration
    
    
    # y = np.abs(fft(u0))
    # N=gridnum
    # plt.plot(np.linspace(0,gridnum,gridnum),np.append(y[N//2:],y[:N//2]),'.')
    # plt.show()
    
    # Use True to run the nonlinear sim
    nln = True
    if nln == True:
        print('nonlinear')
        for t in range(1,num_o_times):
            steps = t*100
            endtime = simtimes[t]
            deltat = (endtime-starttime)/steps
            #print(steps, endtime, deltat,flush=True)
            
            sim_dNLS = sixth(u0,deltat,steps,expconsts,Del,per)
            
            
            # y = np.abs(fft(sim_dNLS))
            # plt.plot(np.linspace(0,gridnum,gridnum),np.append(y[N//2:],y[:N//2]),'.')
            # plt.show()
            # 
            r_dNLS[t,:]=sim_dNLS
    else:
        print('linear')
        #### THIS RUNS THE LINEAR PARTS OF dNLS ONLY
        uhat0 = fft(u0)    
        for t in range(1,num_o_times):
            deltat = simtimes[t]
            uhatnew = uhat0*np.exp(deltat*(1j*expconsts**2-Del))
            r_dNLS[t,:]=ifft(uhatnew)
    
    
    
    
    # Save data
    for s in range(num_o_times):
        if s < 10:
            ss = '00'+str(s)
        elif 9<s<100:
            ss='0'+str(s)
        else:
            ss = str(s)
        np.savetxt(whichset+'Simulations/dNLS Sim/simdNLS'+ss+'.txt',r_dNLS[s].view(float))
