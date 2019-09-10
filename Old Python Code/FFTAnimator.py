
import numpy as np
from numpy.fft import fft, ifft, fftfreq
from scipy import special
import matplotlib.pyplot as plt
from matplotlib import animation
import glob
import os
import NLS


def listdirNH(path):
    return glob.glob(os.path.join(path, '*'))



filename = 'vDysthe'
Writer = animation.writers['ffmpeg']
writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=-1)

#Set framerate. 100 = 10 frames per second
fr = 1

# Takes an array where the first row is the x data and the following rows are
# y data at different time intervals.
whichset = 'JulyData/'

NNN = 1024

filelist = []
for fname in listdirNH(whichset+'Simulations/'+filename+' Sim'):
    if os.path.isfile(fname) == True:
        # skip directories
        filelist.append(fname)

datavals = np.zeros((len(filelist),NNN),dtype=complex)
g = 0
for file in filelist:
    vals = np.loadtxt(file).view(complex)
    datavals[g]=np.abs(1/NNN*fft(vals))
    g=g+1


# Pull out x and y data
L=18000
#xspace = np.linspace(0,0.04373794708462155-0.04373794708462155/1024,1024)
xspace = 1/L*NLS.kvec(NNN)


#Find the number of timesteps
tvs = len(datavals)
print(tvs)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
#ax = plt.axes(xlim=(0, L), ylim=(-max(datavals[1])*1.5, max(datavals[1])*1.5))
ax = plt.axes(xlim=(-max(xspace), max(xspace)), ylim=(0, max(datavals[1])*1.5))
line1, = ax.plot([], [], lw=2)



# initialization function: plot the background of each frame
def init():
    line1.set_data([], [])
    line1.set_marker('.')
    ax.set(title ='FFT')
    return line1,

# animation function.  This is called sequentially
def animate(i):
    yi = datavals[i,:]
    line1.set_data(xspace, yi)
    line1.set_marker('.')
 
    return line1, 

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, interval = fr,
                               frames=tvs, blit=True)

anim.save(whichset+'Simulations/Animations/FFT'+filename+'.mp4', writer=writer)

#plt.show()