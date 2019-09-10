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
# for i in range(100):
#     #tests = np.append(tests,'L_18000_Test_'+str(i+6))
#     tests = np.append(tests,'L_18000_StatRun_'+str(i))
tests = 'L_3hour_'+str(0)

names = ['dGT','dNLS','Dysthe','NLS','vDysthe']






for dir in dirs:
    print('-'*20)
    print(dir)
    
    smallerror = []


    sberror = np.array([]) # The 5 model errors
    with open(dir+'/'+tests+'/Simulations/Errors/SidebandError.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            #print(row)
            sberror=np.append(sberror, float(row[1]))
        index2, value = min(enumerate(sberror), key=operator.itemgetter(1))
        smallerror.append(names[index2])

    Counter1 = Counter(smallerror) 
    most_occur = Counter1.most_common(4)
    print(most_occur)