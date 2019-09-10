"""

Use this code!

"""

import csv
import numpy as np
from collections import Counter
import operator
import glob
import os
import matplotlib.pyplot as plt

def readsvals(whichset):
    with open(whichset+'SpecialVals.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        deltaf = next(csv_reader)[1]
    return deltaf

def listdirNH(path):
    return sorted(glob.glob(os.path.join(path, '*')))

dir = 'JulyData'

dGT = np.array([])
dNLS = np.array([])
Dysthe = np.array([])
IS = np.array([])
NLS = np.array([])
vDysthe = np.array([])


files = listdirNH(dir)

i=0
for test in files: #Test is a set of simulation runs
    testerror = np.array([])
    
    if test.find('Ensemble_') != -1:
        i+=1
        with open(test+'/Simulations/Errors/SidebandError.txt') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                testerror=np.append(testerror,float(row[1]))

        dGT = np.append(dGT, testerror[0])
        dNLS = np.append(dNLS, testerror[1])
        #Dysthe = np.append(Dysthe, testerror[2])
        #IS = np.append(IS, testerror[3])
        NLS = np.append(NLS, testerror[2])
        vDysthe = np.append(vDysthe, testerror[3])
        print(i)
        if i>99:
            break

#sims = [dGT,dNLS,Dysthe,IS,NLS,vDysthe]
sims = [dGT,dNLS,NLS,vDysthe]
for i in sims:
    m = np.mean(i)
    s = np.std(i)
    print(m,s)



