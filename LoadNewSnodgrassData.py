"""
This loads new Snodgrass data and creates figures of the data.
"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import csv

def listdirNH(path):
    return sorted(glob.glob(os.path.join(path, '*')))

subdir = 'Aug4Data/'

# THIS ASSUMES THAT THE GAUGES ARE CORRECTLY ORDERED BY gauge1, gauge2...
eerr = listdirNH(subdir+'Raw/')

# labels for the plot legend
hand =  ['New Zealand', 'Tutuila', 'Palmyra','Honolulu','FLIP','Yakutat']
i = 0
for e in eerr:
    
    
    ### USE THIS TO LOAD EXISTING TEXT FILES
    data = np.transpose(np.loadtxt(e))
    x=data[0]
    y=data[1]
    
    
    ### USE THIS TO CONVERT CSV FILES TO TXT FILES
    # x = np.array([])
    # y = np.array([])
    # with open(e) as csv_file:
    #     csv_reader = csv.reader(csv_file, delimiter=',')
    #     for row in csv_reader:
    #         x = np.append(x,float(row[0]))
    #         y = np.append(y,float(row[1]))
    
    # Save the text, if desired
    np.savetxt(subdir+'Raw/'+e[-10:-4]+'.txt',np.transpose(np.vstack((x,y))))
    
    # Plot and save the figures 
    plt.plot(x,y,label = hand[i])
    i+=1
    
plt.legend()
plt.title(subdir[:-1])
plt.savefig(subdir+'Processed Data Figs/'+subdir[:-1]+' Log Plot.png',dpi=500)
#plt.show()