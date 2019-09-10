import csv
import numpy as np
from collections import Counter 

def readsvals(whichset):
    with open(whichset+'SpecialVals.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        deltaf = next(csv_reader)[1]
    return deltaf

dirs = ['Aug1Data','Aug2Data','JulyData']
tests = np.array([])
for i in range(1,101):
    #tests = np.append(tests,'L_18000_Test_'+str(i+6))
    tests = np.append(tests,'IS_optimization_bet_'+str(i))


for dir in dirs:
    print('-'*20)
    print(dir)
    for test in tests:
        print(test)
        deltaf=readsvals(dir+'/'+test+'/')
        print(deltaf)
        
        modelchoice = []
        with open(dir+'/'+test+'/Simulations/Errors/MinErrors.txt') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                print(row)
                modelchoice.append(row[1])
        Counter1 = Counter(modelchoice) 
  
        # most_common() produces k frequently encountered 
        # input values and their respective counts. 
        most_occur = Counter1.most_common(4)
        print(most_occur)
                
        print()
        print()
        print()
