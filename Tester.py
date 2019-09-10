import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS

# lenly = 10
# for i in range(10):
#     np.random.seed(1)
#     randvals = np.random.rand(lenly)*2*np.pi
#     print(randvals,flush=True)
# 
xspace = np.linspace(0,8,1024)
u = 0.48994682264670886*np.exp(1j*xspace)+0.48994682264670886*np.exp(-1j*xspace)
v = u

up = u+1j*v
fup = 1/1024*fft(up)
#print(fup)

#plt.plot(xspace,np.real(up))
#plt.plot(xspace,np.imag(up))
#plt.plot(NLS.kvec(1024),np.angle(fft(up)),'.')
plt.plot(NLS.kvec(1024),np.angle(fup),'.')
#plt.show()



whichset = 'Aug1Data/dGT_Test_3/Processed/'
filename = 'gauge2.out'
data = np.loadtxt(whichset+filename).view(complex)


t = data[0]
u = data[1]

k = NLS.kvec(len(t))


# plt.plot(t,np.real(u))
# plt.plot(t,np.imag(u))
#plt.show()

fu = 1/len(t)*np.abs(fft(u))

plt.plot(k,1/len(t)*np.abs(fft(u)),'.')


w = np.where(fu == max(fu))[0][0]
k0 = k[w]
fu0 = fu[w]
plt.plot(k0,np.abs(fu0),'*',markersize = 6)
k1 = k[w-2]
fu1 = fu[w-2]

u1 = fup[w-2]
print(u1)
print(k1)
plt.plot(k1,np.abs(fu1),'^',markersize = 3)
plt.plot(k1,np.angle(u1),'^',markersize = 6)
plt.show()
print(np.angle(u1))
# 
# profiledata = np.transpose(np.vstack((np.array([np.real(t)],dtype=float),
#                                     np.array([np.real(u)],dtype = float),
#                                     np.array([np.imag(u)],dtype = float))))
# np.savetxt(whichset+'profile.txt',profiledata)
# 
# fftdata = np.transpose(np.vstack((np.array([np.real(k)],dtype=float),
#                                     np.array([np.real(fu)],dtype = float),
#                                     np.array([np.imag(fu)],dtype = float))))
# np.savetxt(whichset+'fftdata.txt',fftdata)