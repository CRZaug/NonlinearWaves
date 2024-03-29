import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import random as rand
from scipy import interpolate
import NLS
import dNLS
import Dysthe as dy
import vDysthe as vdy
import dGT
import IslasSchober as IS
import csv
from sys import exit


def createstruct(refdir, masterdir, subdirs):
    # This code looks at a sample directory and recreates that structure to store
    # the results of each new run. The code assumes this structure in order to work!
    # Refdir is the reference directory structure, stored on my computer for EACH data set
    # Masterdir is the directory under which all the other work will be stored.
    # Data_sets currently include Aug1Data, Aug2Data, and JulyData. Will be
    # given by a list of indices (I think)
    
    
    
    for sd in subdirs:
                    
        if os.path.isdir(sd+refdir)==False:
            print('This is not a valid reference directory. Shutting down.')
            exit()
        
        inputpath = sd+refdir
        outputpath = sd+masterdir 
        
        for dirpath, dirnames, filenames in os.walk(inputpath):
            structure = os.path.join(outputpath, dirpath[len(inputpath):])
            if not os.path.isdir(structure):
                os.mkdir(structure) 
            else:
                print("Folder already exits!")



def choosesims(SIMULATIONS):
    # Given a vector of "y" and "n" values (Ex: ['y','n','n','n','n','n'], this function creates a
    # vector representing an index that pulls the correct simulations in other parts of the script
    # SIMULATIONS controls which simulations to be considered. ORDER MATTERS.
        # Order: NLS, dNLS, Dysthe, vDysthe, dGT, IS (Just the order in which I implemented the models)
    
    
    # Translate the vector into values corresponding to each simulation. ORDER MATTERS
    y_NLS = SIMULATIONS[0]
    y_dNLS = SIMULATIONS[1]
    y_Dysthe = SIMULATIONS[2]
    y_vDysthe = SIMULATIONS[3]
    y_dGT = SIMULATIONS[4]
    y_IS = SIMULATIONS[5]

    # Create an array that says which simulations we should draw from    
    whichsims = np.array([],dtype = int)

    if y_dGT == 'y':
        whichsims = np.append(whichsims,0)
    
    if y_dNLS =='y':
        whichsims = np.append(whichsims,1)
    
    if y_Dysthe =='y':
        whichsims = np.append(whichsims,2)
        
    if y_IS == 'y':
        whichsims = np.append(whichsims,3)
        
    if y_NLS =='y':
        whichsims = np.append(whichsims,4)

    if y_vDysthe == 'y':
        whichsims = np.append(whichsims,5)
    
    return whichsims



def readsvals(whichset):
    # This function reads and returns values (as a tuple) from the Special Values text file associated with each sim
    # Specifically, it returns w_0, epsilon, and delta
    # whichset is the directory to pull from
    
    
    
    # Open the text file that contains the desired information and read it
    with open(whichset+'SpecialVals.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
    
        vals=[row for idx, row in enumerate(csv_reader) if idx in (4,5,6)]
    
        # Find omega_0, epsilon, and delta. 
        w0 = vals[0][1]
        epsilon = vals[1][1]
        Del = vals[2][1]
    
    # Convert the results to floats so that they are actually intrepreted as numbers
    return (float(w0),float(epsilon),float(Del))



def processnondim(masterdir, deltaf, L, subdirs, doIplot = 'no',doIplot_chi = 'go'):
    # This procecess the data, finds special values, and returns nondimensionalized versions of it
    # Masterdir is the place all of the new information will be sent to
    # Deltaf is the width of the discretized spectrum, should be 1
    # L is the period, should be three hours
    # Subdirs is a vector of the correct subdirectors (out of Aug1, Aug2, July) to pull from
        # Should be of the form subdirs = ['Aug1Data/']
    # doIplot shows the data processing process, such as generating the real surface of the wave
    # doIplot_chi plots all the dimensionless surfaces together
    
    
    
    ### STEP 1: Get distance information
    #distv = np.array([0.0,2400000.0,4200000.0,8700000.0]) # Distances from gauge 1 to gauge n in METERS
    distv = np.array([0.0,2400000.0,4200000.0,8700000.0])
    #distv = np.array([0.0,3300000.0, 5700000.0, 7500000.0]) #USE ONLY FOR AUG4
    
    # Define something that will list directories that are not hidden
    def listdirNH(path):
        return sorted(glob.glob(os.path.join(path, '*')))
    
    j = 0
    for sd in subdirs:
        
        # Find the files to study
        whichset = sd+masterdir
        files = listdirNH(sd+'Rescaled')
        
        
        # Initialize some values
        n = 0
        pi =0
        gaugenumber = 0
        
        if doIplot_chi == 'go':
            fig1,ax1 = plt.subplots(4)
        
        # Get files
        for f in files:
            datavals = np.transpose(np.loadtxt(f).view(float))
            N = len(datavals[1])
            x = datavals[0] # Frequencies in mHz
            sly = datavals[1] # Magnitudes in cm^2/mHz
            
            
    ### STEP 3: Interpolate the data and get new y values
            
            # Interpolation function   
            fnc = interpolate.interp1d(x, sly,kind ='cubic')
            sly2 = fnc(x)
            slymax = max(sly2)
            
            
            xmaxloc = np.where(slymax==sly2)[0][0]
        
            # The first and last values of the Snodgrass data
            aa = x[0]
            bb = x[-1]
        
            # Choose the spacing of the grid based on the first maximum value
            if gaugenumber == 0:
            
                newn = x[xmaxloc]*2/deltaf-1
                newn = int(round(newn,0)) # The new spacing between points
                
                # Create the new grid of x points based on the spacing deltaf and the location of the max
                xn = np.zeros(newn)
                xn[newn//2]=x[xmaxloc]
                for p in range(newn//2-1):
                    xn[newn//2+p] = x[xmaxloc]+deltaf*(p)
                    xn[newn//2-p] = x[xmaxloc]-deltaf*(p)
    
            # Restrict the x range to the x range given by the actual Snodgrass data
            xfinal = np.array([])
            for ixn in xn:
                if np.logical_and(bb>ixn,aa<ixn):
                    xfinal=np.append(xfinal,ixn)
        
            # Get the new values
            ly = fnc(xfinal)*deltaf

            # Get rid of potential negative values
            negind = np.where(ly<0)
            ly[negind]=0
            
            # Find the significant wave height, if desired
            #wvht = 4*np.sqrt(sum(ly))
           
            # Adjust the units on the x axis to get to get from mHz = 0.001Hz to Hz
            xfinal = xfinal*0.001 # now in Hz
            
            gaugenumber +=1
            
    ## STEP 5: Generate random phase and real and imaginary parts
            lenly = len(ly)
            np.random.seed(1) # This ensures that the random phase will be the same 2 runs in a row
            randvals = np.random.rand(lenly)*2*np.pi
            ascale = np.cos(randvals)*np.sqrt(ly)*0.01 # Rescale y axis to m
            bscale = np.sin(randvals)*np.sqrt(ly)*0.01 # Rescale y axix to m
            
            # Add the real and imaginary parts together
            fakevals = (ascale+1j*bscale) # Amplitudes in m
            
            
    ### STEP 7: Locate the carrier wave
            #halfy = ynew2[:len(ynew2)//2]
            if n==0:
                
                carriermax = max(np.abs(fakevals))
                carrierloc = np.where(np.abs(fakevals)==carriermax)[0][0] # The location of the carrier wave
                CarFreq = xfinal[carrierloc] # The freq. at which the carrier wave is located


            
    ### STEP 8: Get nondimensionalization constants    
            g = 9.81 #(m/s^2)
            if n==0:
                w0 = CarFreq*2*np.pi # Get the value from the integer (angular frequency)
                k0 = w0**2/g # The carrier wavenumber
                #m = max(np.abs(ynew2)/len(ynew2)) # a0 in the original FD paper
                m = max(np.abs(fakevals)) # a0 in the original FD paper
                epsilon = 2*m*k0 # The nondimensionalization constant epsilon
                #epsilon = 0.0022
            
                print(f,'Special Values')
                print('delta f', deltaf)
                print('period',L)
                print('Maximum value',m)
                print('Wavenumber',float(k0))
                print('Carrier frequency',float(w0))
                print('epsilon',float(epsilon))
                
                file = open(whichset+'/SpecialVals.txt','w') 
     
                file.write('delta f, '+str(float(deltaf))+'\n') 
                file.write('period, '+str(float(L))+'\n') 
                file.write('maximum value, '+str(float(m))+'\n') 
                file.write('wavenumber k0, '+str(float(k0))+'\n')
                file.write('carrier frequency w0, '+ str(float(w0))+'\n')
                file.write('epsilon, '+str(float(epsilon))+'\n')
     
                file.close() 
               
                
                n = n+1
      
            
                    
                        
    ### STEP 9: Factor out the carrier wave
    
            # Shorten the old spectrums
            #yhat1 = ynew2[:len(ynew2)//2] 
            #k_new_1 = k1[:len(ynew2)//2]
            xfinal_n = xfinal - xfinal[carrierloc]


            
            # PUT INTO k SLOTS, just to be safe. Should be unecessary
            k =np.round(xfinal_n*L,0)
            
            print(xfinal)
            print(xfinal_n*L)
            print(k)
            
            
    ### STEP 6: Make sure there are no gaps 
            ndk = np.arange(k[0],k[-1]+1)
            ndy = np.zeros(len(ndk),dtype=complex)
            
            ndk= np.append(ndk,ndk[-1]+1)
            ndy = np.append(ndy,0)
            
            # The result: ndy will be the spectrum with zeros when necessary
            pip = 0
            for ele in range(len(ndk)):
                if ndk[ele] in k:
                    ndy[ele] = fakevals[pip]
                    pip+=1
            
            print(ndk)
            # hig,hx = plt.subplots(1)
            # hx.plot(ndk,np.abs(ndy),'.')
            #plt.show()
            

    
        
    ### STEP 10: Sum to get back to temporal space and save
        
            # Define new t data (and a new number of data points)
            NNN = 1024
            tnew = np.linspace(0,L-L/NNN,NNN) # ADD TO ICTEMPORALDATA
            
            
            # Find B, the new surface modulating envelope
            B=np.zeros(NNN) 
            for v in range(len(ndy)):
                newvalue = ndy[v]*np.exp(1j*(2*np.pi/L)*ndk[v]*tnew) # Sum a new fourier series
                #newvalue = yhat1[v]*np.exp(1j*(2*np.pi/L)*k_new_2[v]*tnew) # Sum a new fourier series
                B = B+newvalue
            print(B)  
            # Optional plotting
            if doIplot=='go':
                bfig,bax = plt.subplots(2)
                bax[0].plot(tnew,np.real(B))
                bax[1].plot(tnew,np.imag(B))
                bax[0].set_title('Real(B), m')
                bax[1].set_title('Imag(B), m')
                bfig.suptitle(str(f))
                bfig.subplots_adjust(top=0.88)
            
            
                figg,axx = plt.subplots()
                axx.set_title(f+' checking preservation of norm')
                axx.set_ylabel('FFT of factored B')
                axx.plot(NLS.kvec(NNN),1/NNN*np.abs(fft(B)),'.',markersize=10)
                axx.plot(ndk,np.abs(ndy),'*',markersize = 5) ###################
                
                
                #plt.show()
            #plt.close('all')
            
            # Save the x and y values for t and B
            np.savetxt(whichset + '/Processed/'+f[-10:], np.append([tnew],[B],axis=0).view(float))
      
        
    
    ### STEP 11: Nondimensionalize values, plot, and save results
            ndB = k0/epsilon*B # new "y"
            xi = w0*epsilon*tnew # new "x"
            
            
            
            # for j in 1/NNN*fft(ndB):
            #     print(np.abs(j),j)
            # input()
            
            # Fix up xi so that it is periodic and has the correct number of points to match B
            xi_L = (xi[-1]-xi[0])+(xi[1]-xi[0])
            xi_new = np.linspace(0,xi_L-xi_L/NNN,NNN,dtype=complex)
            
            
            chi = epsilon**2*k0*distv # new "t"
        
            # Optional plotting
            if doIplot_chi=='go':
                ax1[pi].set_title(r'$\chi$'+' = '+str(chi[pi]))
                ax1[pi].plot(xi_new,np.real(ndB))
                ax1[pi].plot(xi_new,np.imag(ndB))
                
                pi+=1
            
            # Save the nondimensional x and y axis, chi and non dim B
            np.savetxt(whichset + '/NonDim Data/ND'+f[-10:],np.append([xi_new],[ndB],axis=0).view(float))
            
        # Save just the nondimensonal time vector
        np.savetxt(whichset + '/chi.txt', chi.view(float))
        
        # Optional plotting of all the dimensionless surfaces together
        n=0
        if doIplot_chi =='go':
            fig1.suptitle('Dimensionless '+sd+', L = '+str(L))
            fig1.tight_layout()
            fig1.subplots_adjust(top=0.88)
            plt.show()
            #plt.savefig(whichset + '/NonDim Figs/allgauges.png',dpi=500)
    plt.close('all')


def dataspecialvals(masterdir, subdirs, showplots='no',plotem = 'no'):
    # This script takes in the processed data and gets out the conserved quantities
    # Masterdir is the place all of the new information will be sent to
    # Subdirs is a vector of the correct subdirectors (out of Aug1, Aug2, July) to pull from
        # Should be of the form subdirs = ['Aug1Data/']
    # showplots shows which sidebands are getting counted, shows the CQs, and shows the fit of m
    # plotem controls whether the CQs even get plotted at all 
    
    
    # Define something that will list directories that are not hidden
    def listdirNH(path):
        return sorted(glob.glob(os.path.join(path, '*')))
    
    
    ### STEP 1: Read in the xi, B, and chi data

    
    for set in subdirs:
        
        whichset = set + masterdir
        
        
        # Choose the name of the file the data will be pulled from
        dir = whichset+'NonDim Data'
        dirnames = listdirNH(dir)
        dirlength = len(os.listdir(dir))
        
        tvector = np.loadtxt(whichset+'/chi.txt').view(float) #times
        values = {} # The x and y values at each location
        for name  in dirnames:
            vdatavals = np.loadtxt(name).view(complex)# xi and B
        
            title = name[-10:-4]
            N = len(vdatavals[0])
            
            # Save values to dictionaries
            values[title]=np.append([vdatavals[0]],[vdatavals[1]],axis=0) 
        
        
        ### STEP 2: Find the sideband values and the carrier wave location
        
        # Find the sideband values to track through the rest of the program
        x=values['gauge2'][0]
        y=values['gauge2'][1]
        
        # Perform an FFT of the y values
        yhat =fft(y) # Fourier amplitudes # used to be yhat
        yhat1 = 1/N*np.abs(yhat) # Normalized fourier amplitudes
        
        # Define some constants/constant vectors
        L = (x[-1]-x[0])+(x[1]-x[0]) # The period
        k1=NLS.kvec(N)
        
        
        k=1/L*k1
        
        
        # Find max Fourier amplitude and location
        m = max(yhat1) # The max amplitude (m)
        i = np.where(yhat1 == m)
        
        i = i[0][0]
        
        carrier = k[i]
        
        # plt.plot(k,yhat1,'.')
        # plt.plot(carrier,m,'.',markersize = 5)
        # plt.show()
        
        
        #Shift over the sideband vector so that the carrier wave is in the middle (carrier is 0 so this doesn't really matter)
        
        sv = [] # the sideband vector
        for j in range(len(yhat1)):
            if yhat1[j]>0.00000000001:
                sv.append(j)
       
        for gauge in values:
            
            tempyhat = np.abs(fft(values[gauge][1]))
            # rfnc = interpolate.interp1d(k1, tempyhat,kind ='cubic') ### Need to make the curve make sense before it can be interpolated
            # yhat3 = rfnc()
            
            # Check which values are getting counted as sidebands:
            plt.plot(k1,tempyhat,'.')
            plt.plot(k1[sv],tempyhat[sv],'.',markersize=10)
            if showplots =='go':
                plt.show()
        
        
        
        ### STEP 3: Find P, M, omega_m, omega_p, and the sidebands at each gauge location (each xi)
        
        # Preallocate space for the data values to be stored
        Pvals = np.zeros(len(tvector),dtype = float)
        Mvals = np.zeros(len(tvector),dtype = float)
        PMvals = np.zeros(len(tvector),dtype = float)
        wpeak = np.zeros(len(tvector),dtype = float)
        sideband7 = np.zeros((len(tvector),len(sv)),dtype = float)
        
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
            #print(max(y*np.conj(y)),L)
            PM = P/M # Find omega_m
            
        
            sidebands = yhat1[sv] # Find fourier amplitudes at sidebands
            # plt.plot(k,yhat1,'.')
            # plt.plot(k[sv],sidebands,'.',markersize = 10)
            # plt.show()
            
            
        
         
            Pvals[j] = np.real(P)
            Mvals[j] = np.real(M)
            PMvals[j] = np.real(PM)
            wpeak[j] = np.real(carrier)
            sideband7[j]=np.real(sidebands)
            
            j=j+1
        
        # Get the number on each sidebands for labeling purposes
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
        np.savetxt(whichset+'/sidebandnums'+'.txt', svlab.view(int))
            
        ### STEP 4: Save the Data
        
        dir = whichset+'/Data CQs/NonDim CQ Values/'
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
        
        if plotem == 'go':
        
            # Plotting vectors
            fig1, ax1 = plt.subplots(2,2,figsize = (10,6.5))
            fig2, ax2 = plt.subplots(len(sv),sharex=True,figsize = (10,6.5))
            fig1.suptitle('Quantities of Interest',fontsize=16)
            fig2.suptitle('Select Fourier Amplitudes',fontsize=16)
            
            plotter1 = [Pvals,Mvals,PMvals,wpeak]
            titles1 = ['CQ P', 'CQ M', r'$\omega_m$', r'$\omega_p$']
            
            # CQ figure
            ax1 = ax1.flatten()
            for i in range(len(plotter1)):
                ax1[i].plot(tvector,plotter1[i],'.',markersize =15)
                ax1[i].set_title(titles1[i])
                ax1[i].set_xlabel(r'$\chi$')
                
        
            # Sideband figure
            for k in range(len(sv)):
                vp = sideband7[:,k]
                ax2[k].plot(tvector,vp,'.',markersize=10)
                ax2[k].set_ylabel('a'+ r'$_{'+str(svlab[k])+'}$')
            
            # Fix configurations
            fig1.tight_layout()
            fig1.subplots_adjust(top=0.88)
            fig2.tight_layout()
            fig2.subplots_adjust(top=0.88)
            if showplots =='go':
                plt.show()
            
        
        
        
        # STEP 6: Get a fit of M
        
        fitx = -2*tvector # The x vals of the fit   
        fity = np.log(Mvals) # The y vals of the fit
        
        # Get the fit and define a new y vector
        A = np.vstack([fitx, np.ones(len(fitx))]).T
        m, b = np.linalg.lstsq(A, fity,rcond=-1)[0] # m is delta
        fittedy = m*fitx+b
        newy = np.exp(fittedy)
        
        
        file = open(whichset+'SpecialVals.txt','r')
        last = file.readlines(-1)
        line =str(last[-1][:-3])
        if line[0:5] != 'delta':
            file.close()
            file = open(whichset+'SpecialVals.txt','a')
            file.write('delta, '+str(float(m))+'\n')
        file.close() 
        
        
        # Get the error of the fit
        error = np.sum(np.abs(newy-Mvals)**2)/len(newy)
    
    
        # Plot results of the fit
        fig3,ax3 = plt.subplots(2,1)
        fig3.suptitle('August 2: Fitting ' r'$\mathcal{M}$'+' to find '+r'$\delta$')
        plt.text(.7, .7,'y='+str(m)[1:7]+'x'+str(b)[1:6])
        ax3[1].plot(fitx,fity,'.',label = 'Linearized values')
        ax3[1].plot(fitx,fittedy, label ='Linear fit')
        ax3[1].set_xlabel('-2'+r'$\chi$')
        ax3[1].set_ylabel('ln('+r'$\mathcal{M}$'+')')
        ax3[1].legend(loc='upper left')
        ax3[1].text(-8e-05,-0.58,r'$\delta$'+' = '+str(float(m))[:9])
        
        ax3[0].plot(tvector,Mvals,'.', label = 'Actual values')
        ax3[0].plot(tvector,newy,label = 'Fit curve')
        ax3[0].set_xlabel(r'$\chi$')
        ax3[0].set_ylabel(r'$\mathcal{M}$')
        ax3[0].legend(loc='upper right')
        fig3.tight_layout()
        fig3.subplots_adjust(top=0.88)      
        if showplots =='go':
            #plt.savefig('fitM.png',dpi=500)
            plt.show()
    plt.close('all')


def runsims(SIMULATIONS, masterdir,subdirs, num_o_times, bet, per):
    # SIMULATIONS is the vector of simulations to run
    # Masterdir is the directory everything will be placed into
    # Subdirs control what set of data (Aug1, Aug2, July) are simulated
    # num_o_times gives the number of points between the first and last gauge to output
    # bet is the parameter for the Islas Schober Eqn
    # per controls the 3/16ths rule; ususally set to 0
    
    y_NLS = SIMULATIONS[0]
    y_dNLS = SIMULATIONS[1]
    y_Dysthe = SIMULATIONS[2]
    y_vDysthe = SIMULATIONS[3]
    y_dGT = SIMULATIONS[4]
    y_IS = SIMULATIONS[5]

    
    
    # Define master directory
    for set in subdirs:
        
        whichset = set+masterdir
        
        w0,epsilon,Del = readsvals(whichset)
        
        
        # Read in x and y data
        IC = np.loadtxt(whichset+'NonDim Data/NDgauge2.out').view(complex)
        xspace = IC[0]
        u0 = IC[1]
        
        # Read in time data
        times = np.loadtxt(whichset+'chi.txt').view(float)
        starttime = times[0]
        stoptime = times[-1]
        simtimes = np.linspace(starttime,stoptime,num_o_times)
        
        simtimes = np.sort(np.append(simtimes[1:-1],times)) # Include the ACTUAL data positions in the simtimes
        num_o_times = len(simtimes) # Update num_o_times to reflect the new addition
    
        
        # Save time data
        np.savetxt(whichset+'/Simulations/SimTime.txt',simtimes.view(float))
        
        # Set operator splitting parameters
        L = xspace[-1]-xspace[0]+(xspace[1]-xspace[0])
        gridnum = len(xspace)
        k, expconsts = vdy.kvec(gridnum,L)
        endtime = simtimes[-1]
        rk4steps = 1
        
        PARAMS = [whichset,num_o_times,starttime,rk4steps,gridnum]
        
        if y_NLS =='y':
            print('NLS')
            NLS.runNLS(PARAMS,simtimes,u0,expconsts,per)
            
        if y_dNLS =='y':
            print('dNLS')
            dNLS.rundNLS(PARAMS,simtimes,u0,expconsts,Del,per)
            
        if y_Dysthe =='y':
            print('Dysthe')
            dy.runDysthe(PARAMS,simtimes,u0,k,expconsts,epsilon,per)
            
        if y_vDysthe == 'y':
            print('vDysthe')
            vdy.runvDysthe(PARAMS,simtimes,u0,k,expconsts,epsilon,Del,per)
        
        if y_dGT == 'y':
            print('dGT')
            dGT.rundGT(PARAMS,simtimes,u0,k,expconsts,epsilon,Del,per)
            
        if y_IS == 'y':
            print('IS')
            IS.runIS(PARAMS,simtimes,u0,k,expconsts,bet,epsilon,Del,per)
            
        

def simspecialvals(SIMULATIONS,masterdir, subdirs):
    # SIMULATIONS are the simulations to be considered
    # Masterdir is the directory all new results will be put in
    # Subdirs controls which set(s) of data will be considered (Aug1, Aug2, July)
    
    
    # Define something that will list directories that are not hidden
    def listdirNH(path):
        return sorted(glob.glob(os.path.join(path, '*')))
    
    ### STEP 0: Determine what simulations should be dealt with
    whichsims = choosesims(SIMULATIONS)
    
    ### STEP 1: Load in simulation data
    
    # Load time values
 
    for set in subdirs:
        whichset = set + masterdir
        tvector = np.loadtxt(whichset+'Simulations/SimTime.txt').view(float)
    
        # Choose the name of the file the data will be pulled from
        masterdir1 = whichset+'Simulations/'
        dir = np.array(['dGT Sim','dNLS Sim', 'Dysthe Sim','IS Sim', 'NLS Sim','vDysthe Sim'])[whichsims]
        
        NLSd = {} 
        dNLSd = {}
        Dysthed = {}
        ISd = {}
        vDysthed = {}
        dGTd = {}
        Dictionaries = [dGTd,dNLSd,Dysthed,ISd, NLSd,vDysthed] # Alphabetized
        Dictionaries = np.array(Dictionaries)[whichsims]
        
        # Read in the intital data
        IC = np.loadtxt(whichset+'NonDim Data/NDgauge2.out').view(complex)
        x = IC[0]
        y = IC[1]
        
        h = 0
        for m in dir:
            dict = Dictionaries[h]
            dirnames = listdirNH(masterdir1+m)
            dirlength = len(dirnames)
            kk = 0
            for name in dirnames:
                if os.path.isfile(name) == True:
                    title = tvector[kk]
                    vdatavals = np.loadtxt(name).view(complex)
                    N = len(vdatavals)
                    print(N)
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
        ISCQ  = {}
        vDystheCQ = {}
        dGTCQ = {}
        
        CQDict = np.array([dGTCQ, dNLSCQ, DystheCQ, ISCQ, NLSCQ, vDystheCQ])[whichsims]
        keys = ['P', 'M', 'PM', 'wp', 'sb']
        dname = np.array(['dGT CQ','dNLS CQ','Dysthe CQ', 'IS CQ','NLS CQ','vDysthe CQ'])[whichsims]
        
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
                yhat =fft(y) # Fourier amplitudes-
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
            plotem = 'go'
            if plotem == 'go':
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
            plt.close('all')
        
        
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
        


def redim(SIMULATIONS,masterdir,subdirs):
    # SIMULATIONS are the simulations to be considered
    # Masterdir is the directory all new results will be put in
    # Subdirs controls which set(s) of data will be considered (Aug1, Aug2, July)
    
    # Define something that will list directories that are not hidden
    def listdirNH(path):
        return sorted(glob.glob(os.path.join(path, '*')))
    
    
     ### STEP 0: Determine what simulations should be dealt with
    whichsims = choosesims(SIMULATIONS)
    
    ### STEP 1: READ IN THE EXPERIMENTAL DATA FILES
    
    for set in subdirs:
        
        whichset = set + masterdir    
        
        # Define the dictionaries
        P = {}
        M = {}
        PM = {}
        sb = {}
        wp = {}
        MasterDict = [M,P,PM,sb,wp]
        
        
        
        # Start reading in the data
    
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
        key2 = np.array(['dGT','dNLS', 'Dysthe', 'IS', 'NLS', 'vDysthe'])[whichsims]
        dk = 0
        for subdir in key2:
            files = listdirNH(dir+'/'+subdir+' CQ')
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
        
        w0,epsilon,Del = readsvals(whichset)
        
        k0 = w0**2/g
        
        
        # Dim P
        dim_P = {}
        for key in P:
            ent = P[key]
            xi = ent[0]
            p = ent[1]
            x = xi/(epsilon**2*k0) #gives m
            dim_p = epsilon**3*w0/k0**2*p #gives m^2/s
            dim_P[key] = np.append([x],[dim_p],axis = 0)
            
        # Dim M
        dim_M = {}
        for key in M:
            ent = M[key]
            xi = ent[0]
            m = ent[1]
            x = xi/(epsilon**2*k0) #gives m
            dim_m = (epsilon/k0)**2*m #gives m^2
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
            x = xi/(epsilon**2*k0) #gives m
            dim_pm = (p/m)+w0/(2*np.pi)*1000 # Gives mHz
            dim_PM[key] = np.append([x],[dim_pm],axis = 0)
        
        # Dim sb
        dim_sb = {}
        for key in sb:
            ent = sb[key]
            xi = ent[0]
            x = xi/(epsilon**2*k0) #gives m
            sbv = np.zeros((len(ent)-1,len(xi)))
            for j in range(1,len(ent)):
                sideband = ent[j]
                dim_sideband = (epsilon/k0)*sideband #gives m
                sbv[j-1]=dim_sideband
            dim_sb[key] = np.vstack([x,sbv])

        # Dim wp
        dim_wp = {}
        for key in wp:
            ent = wp[key]
            xi = ent[0]
            peak = ent[1]
            x = xi/(epsilon**2*k0) #gives m
            dim_peak = peak+w0/(2*np.pi)*1000 # Gives mHz
            dim_wp[key] = np.append([x],[dim_peak],axis = 0)
        
        
        
        ### STEP 4: PLOT THE RESULTS
        
        # Initialize for plotting
        plotter1 = [dim_M,dim_P,dim_PM,dim_wp]
        #key2[:0] = [key1]
        key2 = np.append(key1,key2)
        titles1 = ['Quantity '+r'$\mathcal{M}$', 'Quantity '+r'$\mathcal{P}$', 'Quantity '+r'$\omega_m$', 'Quantity '+r'$\omega_p$']
        titles2 = np.loadtxt(os.getcwd()+'/'+whichset+'sidebandnums.txt').view(float)
        y1 = [r'$\mathcal{M}$'+' (m'+r'$^2$'+')',r'$\mathcal{P}$'+' (m'+r'$^2$'+'/s)',r'$\omega_m$'+' (mHz)',r'$\omega_p$'+' (mHz)']
        #y2 = [r'$|a_{-3}|$'+' (m)',r'$|a_{-2}|$'+' (m)',r'$|a_{-1}|$'+' (m)',r'$|a_0|$'+' (m)',r'$|a_1|$'+' (m)',r'$|a_2|$'+' (m)',r'$|a_3|$'+' (m)']
        # https://matplotlib.org/devdocs/gallery/lines_bars_and_markers/linestyles.html
        disp = [0, (0, (1, 1)), (0, (5, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, ()), (0, (3, 1, 1, 1)),(0, (1,5))]
        colors = np.array(['k','#BF8EDE','#EfA0A0','#84E3BE','#EFD7B0','#443E9D','#D65050'])[np.append(0,whichsims+1)]
        sizes = np.array([13,1.5,1.5,1.5,1.5,1.5,1.5])[np.append(0,whichsims+1)]
        
        # Begin plotting
        fig1, ax1 = plt.subplots(4,1,sharex=True,figsize = (11,6.5))
        fig1.suptitle('Quantities of Interest',fontsize=16)
        
        dispind = 0
        for key in key2:
            ax1 = ax1.flatten()
            for i in range(len(plotter1)):
                dict = plotter1[i]
                VALUES = dict[key]
                x = VALUES[0]
                y = VALUES[1]
                if dispind ==0:
                    l0=ax1[i].plot(x,y,'.',color = colors[dispind],markersize = sizes[dispind])
                else:
                    ax1[i].plot(x,y,linestyle = disp[dispind],color = colors[dispind],linewidth = sizes[dispind])
                ax1[i].set_title(titles1[i])
                ax1[i].set_ylabel(y1[i])
                #ax1[i].ticklabel_format(style='sci',scilimits=(-1,1),axis='both')
            dispind += 1
            fig1.legend([l0],     # The line objects
            borderaxespad=0.1,    # Small spacing around legend box
            bbox_to_anchor=(.98,.85),   # Position of legend
            labels=key2,   # The labels for each line
            )
            ax1[-1].set_xlabel('Location (m)')
            fig1.tight_layout()
            fig1.subplots_adjust(top=0.88)
            fig1.subplots_adjust(right=.87)
        plt.savefig(whichset+'Final Figures/CQResultFig.png',dpi=500)
        
        fig2, ax2 = plt.subplots(len(titles2),sharex=True,figsize = (7,1.625*len(titles2)))
        fig2.suptitle('Select Fourier Amplitudes',fontsize=16)
        dispind = 0
        for key in key2: 
            sbvals = dim_sb[key]
            x = sbvals[0,:]
            sideband7=np.delete(sbvals, 0, 0)
            for po in range(len(titles2)):
                if dispind ==0:
                    # print(sideband7)
                    # print(colors[dispind])
                    # print(sizes[dispind])
                    # print()
                    ax2[po].plot(x, sideband7[po], '.', color = colors[dispind], markersize = sizes[dispind])
                else:
                    ax2[po].plot(x,sideband7[po],linestyle = disp[dispind],color = colors[dispind],linewidth = sizes[dispind])
                ax2[po].set_ylabel('a'+ r'$_{'+str(int(titles2[po]))+'}$'+' (m)')

            fig2.tight_layout()
            fig2.subplots_adjust(top=0.88)
            dispind += 1
            ax2[0].legend(key2,bbox_to_anchor=(1, 1))
            ax2[-1].set_xlabel('Location (m)')
        plt.savefig(whichset+'Final Figures/FAResultFig.png',dpi=500)
        
        plt.close(fig2)
        plt.close(fig1)
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
                    np.savetxt(whichset+'Simulations/Dimensional Results/'+str(ky)+' dimCQ/'+val[o]+'.txt',np.transpose(cqval[ky]).view(float))
            o=o+1
        
        # Save sidebands
        for ky in dim_sb:
            if ky == 'Data':
                np.savetxt(whichset+'Data CQs/Dim CQ Values/'+val[-1]+'.txt',np.transpose(dim_sb[ky]).view(float))
            else:
                np.savetxt(whichset+'Simulations/Dimensional Results/'+str(ky)+' dimCQ/'+val[-1]+'.txt',np.transpose(dim_sb[ky]).view(float))

def sberror(SIMULATIONS,masterdir,subdirs):
    # This code computes the errors for the sidebands, comparing simulations to data
    # Note: It does not count error for when there is no sideband information at a certain gauge
    # SIMULATIONS are the simulations to be considered
    # Masterdir is the directory all new results will be put in
    # Subdirs controls which set(s) of data will be considered (Aug1, Aug2, July)
    
    #Determine what simulations should be dealt with
    whichsims = choosesims(SIMULATIONS)

    
    
    # Since our simtime and actual time vectors don't exactly match, we will use this to find the closest values.
    ######## STILL NEED TO STOP PENALIZING THE FUNCTION IF THE SIDEBAND DROPS TO 0??!!
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    # Begin finding errors for each set of data

    for set in subdirs:
        whichset = set + masterdir
        # Load in the data sideband values
        datacqvalues = np.loadtxt(whichset+'Data CQs/Dim CQ Values/dimsb.txt')
        
        # Load in the data and simulations time vectors
        tvector = np.loadtxt(whichset+'Simulations/SimTime.txt').view(float)
        chi = np.loadtxt(whichset+'chi.txt').view(float)
        
        # Choose the name of the file the data will be pulled from
        masterdir1 = whichset+'Simulations/Dimensional Results/'
        dir = np.array(['dGT dimCQ/', 'dNLS dimCQ/', 'Dysthe dimCQ/','IS dimCQ/', 'NLS dimCQ/','vDysthe dimCQ/'])[whichsims]
        
        # Set up some values for file writing
        titles = np.array(['dGT', 'dNLS', 'Dysthe', 'IS', 'NLS', 'vDysthe'])[whichsims]
        errorvect = np.array([])
        
        # Open a file to write the error in
        file = open(whichset+'Simulations/Errors/SidebandError.txt','w+') 
        k = 0
        for d in dir:
       
            # Read in the intital data
            sb = np.loadtxt(masterdir1+d+'dimsb.txt').view(float)

            # Find the time values that won't exactly match
            dif1 = chi[1]
            dif2 = chi[2]
            
            # Find next closest time values
            v1 = find_nearest(tvector,dif1)
            i1 = np.where(tvector==v1)[0][0]
            v2 = find_nearest(tvector,dif2)
            i2 = np.where(tvector==v2)[0][0]
            
            # Now sum error over each of these times.
            usable_times = [0,i1,i2,-1]
            simulationtotalerror = 0
            j=0
            for t in usable_times:
                sim_val = sb[t][1:] # First entry is the distance in the ocean; don't want it
                data_val = datacqvalues[j][1:] # First entry is the distance in the ocean; don't want it
                
                # Get rid of zero values
                zeroindex = np.where(data_val<10**-15)[0] # Where values less than 10^-10 are located; these are effectively 0
                data_val = np.delete(data_val, zeroindex)
                sim_val = np.delete(sim_val, zeroindex)
                error = 1/(len(data_val)-1)*np.sum(np.abs(sim_val-data_val)**2) # Error at each gauge
                simulationtotalerror += error
                j=j+1
                
            # The final total error at all gauges, summed    
            FinalError = np.sum(simulationtotalerror)
            
            # Add this value to a vector of all the errors
            errorvect=np.append(errorvect,FinalError)
            
            # Write the error to a file
            file.write(titles[k]+' Sideband Error, '+str(float(FinalError))+'\n')
            k+=1
         
        file.close()
        
        # Find the minimum error and the simulation that accomplished that
        m = min(errorvect)
        h = np.where(errorvect==m)[0][0]
        
        # Write the minimum error to a file
        file = open(whichset+'Simulations/Errors/MinErrors.txt','w+')
        file.write('Sideband ' +', '+ titles[h]+', '+str(float(m))+'\n')
        file.close()


def cqerror(SIMULATIONS,masterdir,subdirs):
    # This code computes the errors for conserved quantities
        # M
        # P
        # w_m
        # w_p
    # Using the standard error norm. It compares data to simulations
    # SIMULATIONS are the simulations to be considered
    # Masterdir is the directory all new results will be put in
    # Subdirs controls which set(s) of data will be considered (Aug1, Aug2, July)
    
    
    #Determine what simulations should be dealt with
    whichsims = choosesims(SIMULATIONS)
    
    # Since our simtime and actual time vectors don't exactly match, we will use this to find the closest values
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]

    # Begin finding errors for each set of data

    for set in subdirs:
        whichset = set + masterdir
        # Load in the data and simulations time vectors
        tvector = np.loadtxt(whichset+'Simulations/SimTime.txt').view(float)
        chi = np.loadtxt(whichset+'chi.txt').view(float)
        # Find the time values that won't exactly match
        dif1 = chi[1]
        dif2 = chi[2]
    
        # Find next closest time values
        v1 = find_nearest(tvector,dif1)
        i1 = np.where(tvector==v1)[0][0]
        v2 = find_nearest(tvector,dif2)
        i2 = np.where(tvector==v2)[0][0]

        usable_times = [0,i1,i2,-1]

            
        cqnamevector = ['M','P','PM','wp']
        g=0
        for type in ['dimM.txt','dimP.txt','dimPM.txt','dimwp.txt']:
            datacqvalues = np.loadtxt(whichset+'Data CQs/Dim CQ Values/'+type).view(float)[:,1]
 
            # Choose the name of the file the data will be pulled from
            masterdir1 = whichset+'Simulations/Dimensional Results/'
            dir = np.array(['dGT dimCQ/','dNLS dimCQ/', 'Dysthe dimCQ/','IS dimCQ/', 'NLS dimCQ/','vDysthe dimCQ/'])[whichsims]
            
            # Set up some values for file writing
            titles = np.array(['dGT', 'dNLS', 'Dysthe', 'IS','NLS', 'vDysthe'])[whichsims]
            errorvect = np.array([])
            
            # Open a file to write the error in
            file = open(whichset+'Simulations/Errors/'+cqnamevector[g]+' Error.txt','w+') 
            k = 0
            for d in dir:
                # Read in the intital data
                sim_val = np.loadtxt(masterdir1+d+type).view(float)[:,1]

                simulationtotalerror = 1/(len(usable_times)-1)*np.abs(np.sum(sim_val[usable_times]-datacqvalues))
                # Add this value to a vector of all the errors
                errorvect=np.append(errorvect,simulationtotalerror)
                
                # Write the error to a file
                file.write(titles[k]+' '+cqnamevector[g]+' Error, '+str(float(simulationtotalerror))+'\n')
                k+=1
            
            file.close()
            
            # Find the minimum error and the simulation that accomplished that
            m = min(errorvect)
            h = np.where(errorvect==m)[0][0]
            
            # Write the minimum error to a file
            cqnamevector[g]
            file = open(whichset+'Simulations/Errors/MinErrors.txt','a')
            file.write(cqnamevector[g]+', ' + titles[h]+', ' +str(float(m))+ '\n')
            file.close()
            g+=1
            