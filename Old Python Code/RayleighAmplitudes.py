"""

This code produces Rayleigh distributed amplitudes for the Snodgrass data.
We are likely not going to use random amplitudes for the simulations because
1. Swell isn't really random
2. How would we compare random amplitudes at gauge 1 vs gauge 2 across sims?

"""
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy import interpolate

newn = 130
L = 3600*3

# p = np.zeros(1000)
# for i in range (0,1000):
#     p[i]=np.random.rayleigh(10)
# plt.hist(p,bins = 100)
# plt.show()


### STEP 1: Get distance information
distv = np.array([0.0,2400000.0,4200000.0,8700000.0]) # Distances between gauges in METERS



### STEP 2: Read in information at each gauge for each event
subdirs = ['Aug1Data/','Aug2Data/','JulyData/']


# Define something that will list directories that are not hidden
def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))

j = 0
for sd in subdirs:
    

    files = listdirNH(sd+'Rescaled')
    
    # Initialize some values
    n = 0
    pi =0
    gaugenumber = 0
    

    
    # Get files
    for f in files:
        datavals = np.transpose(np.loadtxt(f).view(float))
        N = len(datavals[1])
        x = datavals[0] # Frequencies
        sly = datavals[1] # Magnitudes
        
        
        
### STEP 3: Interpolate the data and get new y values
        
        # Interpolation function   
        fnc = interpolate.interp1d(x, sly,kind ='cubic')
        sly2 = fnc(x)
        slymax = max(sly2)
        
        xmaxloc = np.where(slymax==sly2)[0][0]
    
        aa = x[0]
        bb = x[-1]
    
        # Choose the spacing of the grid based on the first maximum value
        if gaugenumber == 0:
        
        
            deltaf = 1
            
            newn = x[xmaxloc]*2/deltaf-1
            newn = int(round(newn,0))
            print(newn)
            
            #deltaf=x[xmaxloc]*2/(newn+1) # The new spacing between points
            
            #print(deltaf)
            
            # deltaf = 1 mc/s or 1 mHz from the Snodgrass paper.
            #deltaf = 1
            
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
        
        # Find the significant wave height
        
        #wvht = 4*np.sqrt(sum(ly))
       
        
        # Adjust the units on the x axis to get to mHz
        xfinal = xfinal*0.001
        
        gaugenumber +=1      
        
        
        
### STEP 4: Get the k vector using integer division and clean up
        lenly = len(ly)
        k = np.round((xfinal)/(2*np.pi/L),0) # Then divide by 2pi/L (using integer division) to get the k vector
        
        
### STEP 5: Generate random phase and amplitude
        # Random phase
        randvals = np.random.rand(lenly)*2*np.pi
        
        # Random amplitude
        lypre = ly
        for i in range(100):
            ly = np.random.rayleigh(ly*np.sqrt(2/np.pi)) # Rayleigh distribute
            
            # Construct the complex result
            ascale = np.cos(randvals)*np.sqrt(ly)*0.01 # Rescale y axis
            bscale = np.sin(randvals)*np.sqrt(ly)*0.01 # Rescale y axix
            
            # Add the real and imaginary parts together
            fakevals = (ascale+1j*bscale)
            
            plt.plot(k,fakevals*np.conj(fakevals),'.')
        plt.plot(k,(np.sqrt(lypre)*0.01)**2)
        #plt.plot((xfinal)/(2*np.pi/L),fakevals*np.conj(fakevals),'.')
    # plt.plot(k,np.zeros(len(k)),'.')
        plt.show()
        
        


