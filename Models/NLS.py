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
from scipy import special
import matplotlib.pyplot as plt

#______________________ LINEAR _______________________________
def linear(initial_condition,deltat,expconsts): 
    uhat0=fft(initial_condition)
    return ifft(uhat0*np.exp(1j*expconsts**2*deltat))

#______________________NONLINEAR_________________________________
def nonlinear(initial_condition,deltat):
    return initial_condition*np.exp(4j*np.conj(initial_condition)*initial_condition*deltat)
    
#__________ FIRST ORDER OPERATOR SPLITTING _______________________
def first(initial_condition,deltat,steps,expconsts):
    
    unew=initial_condition
    for j in range(0,steps):
        unew=linear(unew,deltat,expconsts)
        unew=nonlinear(unew,deltat)
    return unew

#______________SECOND ORDER OPERATOR SPLITTING__________________
def second(initial_condition,deltat,steps,expconsts):
    # Preparation before the operator splitting
    sdt = deltat*0.5 # Creates a new timestep size
    # The first part
    unew = linear(initial_condition,sdt,expconsts) # linear part
    unew = nonlinear(unew,deltat) # nonlinear part, which outputs unew to be used in the next loop
    
    for i in range(0,steps-1): # The second part; the bulk of the loop (includes every timestep except first and last)
        unew = linear(unew, deltat,expconsts) # Uses unew from previous loop 
        unew = nonlinear(unew,deltat) # Outputs a new unew
    ufinal = linear(unew,sdt,expconsts) # This is the final result

    return ufinal

#________________SIXTH ORDER OPERATOR SPLITTING___________________

def sixth(initial_condition,deltat,steps,expconsts): # Take in the initial space profile and timevector 
    # Preparation for the operator splitting
    w3 = 0.7845136104775600
    w2 = 0.2355732133593570
    w1 = -1.177679984178870
    w0 = 1.3151863206839060
    
    lt = deltat*0.5
    nlt = deltat
    
    # Begin operator splitting
    unew= linear(initial_condition,lt*w3,expconsts) # linear part, the first step
    
    for i in range(0,steps-1): # The second step, loop through everything except the first and last timestep
        unew = nonlinear(unew, nlt*w3)
        unew = linear(unew, lt*(w3+w2),expconsts)
        unew = nonlinear(unew, nlt*w2)
        unew = linear(unew, lt*(w1+w2),expconsts)
        unew = nonlinear(unew,nlt*w1)
        unew = linear(unew, lt*(w1+w0),expconsts)
        unew = nonlinear(unew, nlt*w0)
        unew = linear(unew, lt*(w0+w1),expconsts)
        unew = nonlinear(unew, nlt*w1)
        unew = linear(unew, lt*(w1+w2),expconsts)
        unew = nonlinear(unew, nlt*w2)
        unew = linear(unew, lt*(w2+w3),expconsts)
        unew = nonlinear(unew, nlt*w3)
        unew = linear(unew, nlt*w3,expconsts) # Note different timestep
        
    # The final step: perform the operatons one last time and output the result    
    unew = nonlinear(unew, nlt*w3)
    unew = linear(unew, lt*(w3+w2),expconsts)
    unew = nonlinear(unew, nlt*w2)
    unew = linear(unew, lt*(w1+w2),expconsts)
    unew = nonlinear(unew,nlt*w1)
    unew = linear(unew, lt*(w1+w0),expconsts)
    unew = nonlinear(unew, nlt*w0)
    unew = linear(unew, lt*(w0+w1),expconsts)
    unew = nonlinear(unew, nlt*w1)
    unew = linear(unew, lt*(w1+w2),expconsts)
    unew = nonlinear(unew, nlt*w2)
    unew = linear(unew, lt*(w2+w3),expconsts)
    unew = nonlinear(unew, nlt*w3)
    ufinal = linear(unew, lt*w3,expconsts) # This is the final result
    return ufinal

#________________________K vector_____________________________
def kvec(gridnum):
    split = int(gridnum/2)
    k = np.zeros(gridnum)
    k[0:split+1] = np.arange(0,split+1,1)
    k[split+1:gridnum] = np.arange(-split+1,0,1)

    return k

#___________________INITIAL DATA (SECH)_______________________
def sech_id(gridnum,kappa,L):
    
    xspace=np.linspace(-L/2.0, L/2.0-L/gridnum, gridnum) # The x values at time t=0
    initialprofile=1/np.sqrt(2)*kappa*1/np.cosh(kappa*xspace) # The y values at time t=0
    
    k = kvec(gridnum)
    expconsts=-1j*(2.0*np.pi/L*k)**2
    
    return (xspace,initialprofile,expconsts)

#___________________INITIAL DATA (CN)_______________________

def cn_L(kappa,k0):
    L=4.0*special.ellipk(k0**2)/kappa # period, whole interval
    return L

def cn_id(gridnum,kappa,k0,L):
    
    xspace=np.linspace(-L/2.0, L/2.0-L/gridnum, gridnum)
    _,cn,_,asd=special.ellipj(kappa*(xspace),k0**2)
    initialprofile=k0*kappa/np.sqrt(2.0)*cn
    
    k = kvec(gridnum)
    expconsts=-1j*(2.0*np.pi/L*k)**2

    return (xspace,initialprofile,expconsts)

#______________________NLS SOLN: SECH ___________________________
def sechsoln(xspace,t,kappa):
    return 1/np.sqrt(2)*kappa*1/np.cosh(kappa*xspace)*np.exp(1j*kappa**2*t)

#_______________________NLS SOLN: CN____________________________
def cnsoln(xspace,t,kappa,k0):
    _,cn,_,asd=special.ellipj(kappa*(xspace),k0**2)
    return (kappa*k0*cn/np.sqrt(2))*np.exp(1j*(kappa**2)*((2*k0**2)-1)*t)


#______________________OPERATOR: Mth DERIVATIVE________________
def der_fft(u,L,k,n):
    # n is the number of derivatives to take
    consts = (2*np.pi*1j*k/L)
    return ifft((consts**n)*fft(u))

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
    return np.sum(np.abs(an_ufinal*np.conj(an_ufinal)-(ufinal*np.conj(ufinal))))

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
