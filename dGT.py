import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

#____________________K VECTOR AND CONSTANTS__________________________
def kvec(gridnum,L):
    split = int(gridnum/2)
    k = np.zeros(gridnum,dtype = complex)
    k[0:split+1] = np.arange(0,split+1,1)
    k[split+1:gridnum] = np.arange(-split+1,0,1)
    expconsts=(1j*2.0*np.pi*k)/L
    return (k, expconsts)


#____________________HILBERT TRANSFORM AND DERIVATIVE__________________________
def der_hil(u0,expconsts,k,n):
    #n is the number of derivatives to take; n=1 is a first derivative
    return ifft(-1j*np.sign(k)*expconsts**n*fft(u0))


#____________________OPERATOR: Nth DERIVATIVE__________________________
def der_fft(u,expconsts,n):
    # n is the number of derivatives to take
    return ifft((expconsts**n)*fft(u))


#________________________VECTOR RK4__________________________________
def rk4(u0, dt, rk4steps,k,expconsts,epsilon):
    # deltat is the final timestep, esentially.
    t = 0
    for i in range(rk4steps):
        km1 = drhs(u0, t,k,expconsts,epsilon)
        km2 = drhs(u0+(dt/2)*km1, t+dt/2,k,expconsts,epsilon)
        km3 = drhs(u0+(dt/2)*km2, t+dt/2,k,expconsts,epsilon)
        km4 = drhs(u0+(dt)*km3, t+dt,k,expconsts,epsilon)
        u0 = u0 +(dt/6)*(km1+2*km2+2*km3+km4)
        t= t+dt

    return u0


#__________________________RK4 RHS_________________________________
def drhs(u0, t,k,expconsts,epsilon):
    u0p = der_fft(u0,expconsts,1) # First derivative
    U0 = np.conj(u0) # Complex conjugate 
    U0p = np.conj(u0p) # Complex conjugate of the derivative
    modu2 = U0*u0 # Modulus squared
    return 4j*modu2*u0+epsilon*(32*modu2*u0p-8j*u0*der_hil(modu2,expconsts,k,1))
    

#____________________________NONLINEAR_____________________________
def nonlinear(u0,deltat,rk4steps,k,expconsts,epsilon):
    dt = deltat/rk4steps
    return rk4(u0, dt, rk4steps,k,expconsts,epsilon)


#_____________________________LINEAR________________________________
def linear(u0,deltat,linexpconsts,per):
    uhat0=fft(u0)
    uhatnew = uhat0*np.exp(linexpconsts*deltat)
    
    N = len(u0)
    A=int(per*N/2)
    uhatnew[N//2-A:N//2] = np.zeros(A)
    return ifft(uhatnew)


#__________________FIRST ORDER OPERATOR SPLITTING____________________
def first(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del,per):
    
    linexpconsts=expconsts**2*(1j+10*epsilon**2*Del)+5*epsilon*Del*1j*expconsts-Del
    unew = u0
    
    for j in range(steps):
        unew=linear(unew,deltat,linexpconsts,per)
        unew=nonlinear(unew,deltat,rk4steps,k,expconsts,epsilon)
        
    return unew


#__________________SECOND ORDER OPERATOR SPLITTING____________________
def second(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del,per):
    
    sdt = deltat*0.5
    linexpconsts=expconsts**2*(1j+10*epsilon**2*Del)+5*epsilon*Del*1j*expconsts-Del
    
    unew = linear(u0,sdt,linexpconsts,per) # linear part
    unew = nonlinear(unew,deltat,rk4steps,k,expconsts,epsilon) # nonlinear part, which outputs unew to be used in the next loop
    
    for i in range(0,steps-1): # The second part; the bulk of the loop (includes every timestep except first and last)
        unew = linear(unew,deltat,linexpconsts,per) # Uses unew from previous loop 
        unew = nonlinear(unew,deltat,rk4steps,k,expconsts,epsilon) # Outputs a new unew

    ufinal = linear(unew,sdt,linexpconsts,per) # This is the final result
 
    return ufinal


#__________________SIXTH ORDER OPERATOR SPLITTING____________________
def sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del,per):
    
    # Preparation for the operator splitting
    w3 = 0.7845136104775600
    w2 = 0.2355732133593570
    w1 = -1.177679984178870
    w0 = 1.3151863206839060
    
    lt = deltat*0.5
    nlt = deltat
    
    linexpconsts=expconsts**2*(1j+10*epsilon**2*Del)+5*epsilon*Del*1j*expconsts-Del
    
    # Begin operator splitting
    unew= linear(u0,lt*w3,linexpconsts,per)
    
    for i in range(0,steps-1): # The second step, loop through everything except the first and last timestep
        unew = nonlinear(unew, nlt*w3,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w3+w2),linexpconsts,per)
        unew = nonlinear(unew, nlt*w2,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w1+w2),linexpconsts,per)
        unew = nonlinear(unew,nlt*w1,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w1+w0),linexpconsts,per)
        unew = nonlinear(unew, nlt*w0,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w0+w1),linexpconsts,per)
        unew = nonlinear(unew, nlt*w1,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w1+w2),linexpconsts,per)
        unew = nonlinear(unew, nlt*w2,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, lt*(w2+w3),linexpconsts,per)
        unew = nonlinear(unew, nlt*w3,rk4steps,k,expconsts,epsilon)
        unew = linear(unew, nlt*w3,linexpconsts,per) # Note different timestep
        
    # The final step: perform the operatons one last time and output the result    
    unew = nonlinear(unew, nlt*w3,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w3+w2),linexpconsts,per)
    unew = nonlinear(unew, nlt*w2,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w1+w2),linexpconsts,per)
    unew = nonlinear(unew,nlt*w1,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w1+w0),linexpconsts,per)
    unew = nonlinear(unew, nlt*w0,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w0+w1),linexpconsts,per)
    unew = nonlinear(unew, nlt*w1,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w1+w2),linexpconsts,per)
    unew = nonlinear(unew, nlt*w2,rk4steps,k,expconsts,epsilon)
    unew = linear(unew, lt*(w2+w3),linexpconsts,per)
    unew = nonlinear(unew, nlt*w3,rk4steps,k,expconsts,epsilon)
    ufinal = linear(unew, lt*w3,linexpconsts,per) # This is the final result

    return ufinal

#________________________EXACT SOLUTION__________________________
def vDysthesoln(xspace,b0,l,rho,epsilon,Del,t):
    
    wr = -Del*(1+5*epsilon*l)*t
    wi = -t*l**2+(2*b0**2*(1+6*epsilon*l))/(Del*(1+5*epsilon*l))*(1-np.exp(-2*Del*(1+5*epsilon*l)*t))
    
    return b0*np.exp(1j*l*xspace+wr+1j*wi+1j*rho)

#________________________INITIAL DATA__________________________
def vDysthe_id(gridnum,L,b0,l,rho,epsilon,Del):
    
    xspace=np.linspace(-L/2.0, L/2.0-L/gridnum, gridnum,dtype=complex)
    initialprofile = vDysthesoln(xspace,b0,l,rho,epsilon,Del,0)
    
    return (xspace,initialprofile)

#____________________OPERATOR: INTEGRATE__________________________
def trapezoid(u,L):
    return (L/len(u))*np.sum(u)


#_______________________CONSERVED QUANTITY: M______________________
def CQ_M(u,L):
    return 1/L*trapezoid(u*np.conj(u),L) # u should be defined along one period


#______________________CONSERVED QUANTITY: P1________________________
def CQ_P1(u,L,expconsts): # more efficient
    u_x = der_fft(u,expconsts,1)
    ccu_x = np.conj(u_x)
    return -1/(L)*trapezoid(np.imag(u*ccu_x),L)


#______________________CONSERVED QUANTITY: P2________________________
def CQ_P2(u,L,expconsts):
    u_x = der_fft(u,expconsts,1) # Derivative of u
    ccu_x = np.conj(u_x) # Complex conjugate of the derivative of u
    ccu = np.conj(u) # Complex conjugate of u
    return 1j/(2*L)*trapezoid((u*ccu_x)-(u_x*ccu),L)


#__________________MINI CONSERVED QUANTITY: Q________________________
def CQ_Q(u,L,expconsts):
    u_x = der_fft(u,expconsts,1) # Derivative of u
    ccu_x = np.conj(u_x) # Complex conjugate of the derivative of u
    return 1/(L)*trapezoid(u_x*ccu_x,L)


#__________________MINI CONSERVED QUANTITY: R________________________
def CQ_R(u,L,expconsts):
    u_xx = der_fft(u,expconsts,2) # Derivative of u
    ccu = np.conj(u) # Complex conjugate of u
    return 1/(L)*trapezoid(np.imag(u*ccu**2*u_xx),L)


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




#_______________________RUN SIMULATIONS_________________________________

def rundGT(PARAMS,simtimes,u0,k,expconsts,epsilon,Del,per):
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
    
    # Initialize matrix
    r_dGT = np.zeros((num_o_times,gridnum),dtype=complex)
    r_dGT[0,:]=u0
    
    # During iteration, deltat (the size of the step remains the same)
    # The initial profile u0 is always the input
    # The final time to step out increases in each loop, thus the number of steps increases in each iteration
    
    # Use True to run the nonlinear sim
    nln = True
    if nln == True:
        print('nonlinear')
        for t in range(1,num_o_times):
            steps = t*10
            endtime = simtimes[t]
            deltat = (endtime-starttime)/steps
        
            sim_dGT = sixth(u0,deltat,steps,rk4steps,k,expconsts,epsilon,Del,per)
            r_dGT[t,:]=sim_dGT
    else:
        print('linear')
        #### THIS RUNS THE LINEAR PARTS OF dGT ONLY
        uhat0 = fft(u0)
        linexpconsts=expconsts**2*(1j+10*epsilon**2*Del)+5*epsilon*Del*1j*expconsts-Del
        for t in range(1,num_o_times):
            deltat = simtimes[t]
            uhatnew = uhat0*np.exp(linexpconsts*deltat)
            r_dGT[t,:]=ifft(uhatnew)
    
    # Save data
    for s in range(num_o_times):
        if s < 10:
            ss = '00'+str(s)
        elif 9<s<100:
            ss='0'+str(s)
        else:
            ss = str(s)
        np.savetxt(whichset+'Simulations/dGT Sim/simdGT'+ss+'.txt',r_dGT[s].view(float))