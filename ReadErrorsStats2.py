import csv
import numpy as np
from collections import Counter
import operator

def readsvals(whichset):
    with open(whichset+'SpecialVals.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        deltaf = next(csv_reader)[1]
    return deltaf

dirs = ['Aug1Data','Aug2Data','JulyData']
tests = np.array([])
for i in range(100):
    #tests = np.append(tests,'L_18000_Test_'+str(i+6))
    tests = np.append(tests,'L_18000_smalldfcheck_'+str(i))


names = ['dGT','dNLS','Dysthe','NLS','vDysthe']





set = np.zeros((5,5))
for dir in dirs:
    print('-'*20)
    print(dir)
    
    for k in range(5):
        sberror = np.array([]) # The 5 model errors
        with open(dir+'/'+tests[k]+'/Simulations/Errors/SidebandError.txt') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                #print(row)
                sberror=np.append(sberror, float(row[1]))
         
        set[k,:] = sberror

    set =  np.matrix(set)
    
    error = []
    stdev = []
    for j in range(5):
    
        print(names[j],np.average(set[:,j]), '+/-', np.std(set[:,j])/10)
        error.append(np.average(set[:,j]))
        stdev.append(np.std(set[:,j]))
 
    index1, value = max(enumerate(error), key=operator.itemgetter(1))
    index2, value = min(enumerate(error), key=operator.itemgetter(1))
    print('max:', names[index1])
    print('min:', names[index2])
        
    # dGT=set[0,:]
    # dNLS = set[1,:]
    # Dysthe = set[2,:]
    # NLS = set[3,:]
    # vDysthe = set[4,:]
    # 
    # 
    # errorsets = [dGT, dNLS,Dysthe,NLS,vDysthe]
    # 
    # print(vDysthe)
    # averagemat = []
    # h=0
    # for set in errorsets:
    #  
    #     average = np.mean(set)
    # 
    #     averagemat.append(average)
    #     h+=1
    # averagemat = np.transpose(averagemat)
    # print(averagemat)
    # print()
    # for r in range(len(averagemat)):
    #     index1, value = max(enumerate(averagemat[r]), key=operator.itemgetter(1))
    #     index2, value = min(enumerate(averagemat[r]), key=operator.itemgetter(1))
    #     print('delta f value:', readsvals(dir+'/L_18000_Groupso7_'+str(r+1)+'/'))
    #     print('max:', names[index1])
    #     print('min', names[index2])
    #     print()