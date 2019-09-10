
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
    datavals[g]=vals
    g=g+1

NNN = 1024
# Pull out x and y data
L=0.019739236375746372
xspace = np.linspace(0,L-L/NNN,NNN)



#Find the number of timesteps
tvs = len(datavals)
print(tvs)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
#ax = plt.axes(xlim=(0, L), ylim=(-max(datavals[1])*1.5, max(datavals[1])*1.5))
ax = plt.axes(xlim=(0, L), ylim=(-max(datavals[1])*5, 5*max(datavals[1])))
line1, = ax.plot([], [], lw=2,label = 'imaginary')
line2, = ax.plot([],[], lw=2,label = 'real')
ax.legend()

# initialization function: plot the background of each frame
def init():
    line1.set_data([], [])
    line2.set_data([],[])
    ax.set(title ='Numerical: Real (o) and Imaginary (b) Parts')
    ax.legend()
    return line1, line2,

# animation function.  This is called sequentially
def animate(i):
    yr = np.real(datavals[i,:])
    yi = np.imag(datavals[i,:])
    line1.set_data(xspace, yi)
    line2.set_data(xspace,yr)
    return line1, line2,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init, interval = fr,
                               frames=tvs, blit=True)

anim.save(whichset+'Simulations/Animations/'+filename+'.mp4', writer=writer)

#plt.show()