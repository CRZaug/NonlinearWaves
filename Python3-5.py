import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft
import NLS
from scipy.fftpack import dct

# xspace = np.linspace(-50)
# uhat = 0.13 * np.exp(-(xspace-(63/1000))**2/100)
# 
# NNN = 1024
# 
# mu, sigma = 0, 0.00001 # mean and standard deviation
# amp = np.random.normal(mu, sigma, NNN)
# 
# mu, sigma = 63, 0.00001
# 
# freq = np.random.normal(mu, sigma, NNN)
# phase = np.random.rand(NNN)*2*np.pi
# 
# L = 10800
# 
# 
# 
# x = np.linspace(0,L,1024)
# amp = 0.1*np.cos(-6*x/L)
# y = amp*np.cos(2*np.pi*x/L+phase)
# 
# 
# plt.plot(x,y)
# plt.show()
# 
# plt.plot(x, 1/NNN*np.abs(fft(y)))
# plt.show()

L = 10800
N = 1024
nu = 1/L
x = np.linspace(0,L,N)
def A(j,N):
    return 0.1*np.exp(-3*j/N)
def phi(j):
    np.random.seed(1)
    return np.random.rand(1)*2*np.pi

def z(x):
    z=0
    for j in range(N//2):
        z =+ A(j,N)*np.cos(2*np.pi*j*x*nu+phi(j))
    return z

k = NLS.kvec(N)

z = z(x)
plt.plot(x,z)
plt.show()

#plt.plot(k,1/N*np.abs(fft(z)),'.')
print(max(1/N*np.abs(fft(z))))
#plt.show()
########################

#z = z-np.mean(z)




# def rho(z,x,tvect):
#     rhovect = np.zeros(len(tvect))
#     print(len(tvect))
#     for i in range(len(tvect)):
#         print(i)
#         t = tvect[i]
#         rhovect[i]=np.mean(z(x)*z(x+t))
#         print(rhovect[i],flush=True)
#     return rhovect
# 
# tvect = np.linspace(-1000,1000,1024)
# 
# rh= rho(z,x,tvect)
# 
# plt.plot(tvect,rh)
# plt.show()
# 
# plt.plot(k, dct(rh))
# plt.show()
