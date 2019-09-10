"""

This code reads in old and new information from IS simulations for the purposes of
beta optimization.

It saves the "new" information to a file, which can then be accessed in the future
as "old" information.

THIS CODE IS MEANT TO BE USED WITH THE SCRIPT "RUNSCRIPTIS"

"""
import csv
import numpy as np
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

pullfrom = 'JulyData/' # This contains new error information from recent sims
subdir = 'JulyData/' # This contains old error information from past sims

#Access existing data
xvector = np.array([])
errors = np.array([])
eerr = listdirNH(subdir+'IS Optimization Errors')
i = 0
for e in eerr:
    print(e,i)
    data = np.transpose(np.loadtxt(e))
    for d in range(len(data[0])):
        #if data[0][d]>-8700:
        if data[0][d]>-1000000000:
            xvector = np.append(xvector,data[0][d])
            errors = np.append(errors,data[1][d])
    i+=1

# Plot old data

z= np.polyfit(xvector,errors,2)
print(z)
p = np.poly1d(z)
px=np.linspace(0,6e9,2000)
py = p(px)

mn = np.where(py==min(py))
mnx = px[mn]
print(mnx)
plt.title(subdir[:-1])
plt.plot(xvector,errors,'k.',markersize = 10)
plt.plot(px,py,'.',markersize = 5)
plt.plot(mnx,min(py),'*',markersize = 7)

# 
plt.show()



# ##Get new information
# pullfrom = 'JulyData/'
# 
# subdir = pullfrom
# newerror = np.array([])
# newxvector = np.array([])
# files = listdirNH(pullfrom)
# 
# for test in files: # Test is a set of simulation runs
#     print(test)
#     if test.find('IS_optimization_bet') != -1:
#         if test[-11:].find('_') == -1: # This may need to be ==
#             print('test',test[-11:])
#             with open(test+'/Simulations/Errors/SidebandError.txt') as csv_file:
#                 csv_reader = csv.reader(csv_file, delimiter=',')
#                 for row in csv_reader:
#                     newerror=np.append(newerror,float(row[1]))
#             newxvector = np.append(newxvector,float(test[-11:]))
#             
# #newxvector = np.linspace(0,35000,100)
# print(newerror)
# #Plot the new error information
# plt.plot(newxvector,newerror,'k.')
# plt.show()
# 
# # Save the new error information to the overall data directory
# newset = np.vstack((newxvector,newerror))
# mn = min(newxvector)
# mx = max(newxvector)
# print(pullfrom,mn,mx)
# print()
# np.savetxt(subdir+'IS Optimization Errors/run_'+str(mn)+'_to_'+str(mx)+'.txt',np.transpose(newset))
# print('it saved')
