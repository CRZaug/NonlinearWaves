"""
THIS IS REALLY TERRIBLE CODE. ReadISErrors.py is much better and accomplishes the same goals.

Goal: Read in errors created during the optimization process for IS

"""

import csv
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import glob
import os
import itertools

def readsvals(whichset):
    with open(whichset+'SpecialVals.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        deltaf = next(csv_reader)[1]
    return deltaf

def listdirNH(path):
    return sorted(glob.glob(os.path.join(path, '*')))


tests = np.array([])
xvector = np.array([])
#for i in range(1,101):
files = listdirNH('JulyDataIS')

for f in files:
    #tests = np.append(tests,'L_18000_Test_'+str(i+6))
    tests = np.append(tests,f)
    xvector = np.append(xvector,float(f[31:]))

# lists = sorted(itertools.izip(*[xvector, tests]))
# xvector, new_y = list(itertools.izip(*lists))    

errorvectors = np.zeros((1,199))

k = 0
for test in tests:
    
    modelchoice = []
    with open(test+'/Simulations/Errors/SidebandError.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            errorvectors[0,k]=row[1]
            print(row[1])
    k+=1
# plt.title(dir+' beta vs IS error')
# plt.plot(np.linspace(-250,-50,100),errorvectors[i,:],'.')
# plt.show()
#i+=1

set1 = np.vstack((xvector,np.array(errorvectors[0,:99])))

np.savetxt('JulyData/IS Optimization Errors/run_50_to_250.txt',np.transpose(set1))
    
    
    
tests2 = np.array([])  
# 
dirs = ['July']
# #,'Aug2Data','JulyData']
# tests = np.array([])
for iy in range(1,101):
#for i in np.linspace(-250,-50,100):
    #tests = np.append(tests,'L_18000_Test_'+str(i+6))
    tests2 = np.append(tests2,'IS_optimization_bet_'+str(iy))

for dir in dirs:
    print('-'*20)
    print(dir)
    for test in tests2:
        deltaf=readsvals(dir+'/'+test+'/')
        
        modelchoice = []
        with open(dir+'/'+test+'/Simulations/Errors/SidebandError.txt') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                errorvectors[0,k]=row[1]
                print(row[1])
        k+=1
    plt.title(dir+' beta vs IS error')
    plt.plot(np.linspace(-50,50,100),errorvectors[0,99:],'k.')
    plt.plot(xvector,errorvectors[0,:99],'k.')
    #plt.show()
    
set2 = np.vstack((np.linspace(-50,50,100),np.array(errorvectors[0,99:])))
print(set2)

np.savetxt('JulyData/IS Optimization Errors/run_-50_to_50.txt',np.transpose(set2))