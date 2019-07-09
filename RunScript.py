#import Master as ms
import MasterNEW_IS as ms
import os
from sys import exit
 
############################
# IMPORTANT::::: DEL IS SET TO 0 FOR DNLS AND VDYSTHE!!!!!! 6TH ORDER OS

ref = 'Reference Directory Structure 6 Sims/'



# The period in seconds. (3 h period)
period = 10800

# Delta f (SHOULD BE 1)
deltaf = 1

# The data sets to run. Must be enclosed in brackets (AKA a list)
#subdirs = ['Aug1Data/','Aug2Data/','JulyData/']
subdirs = ['Aug2Data/']

# The simulations to run. Put 'y' if it is to be run:
SIMULATIONS = ['y', #NLS
               'y', #dNLS
               'y', #Dysthe
               'y', #vDysthe
               'y', #dGT
               'y', #IS
                    ]


# Define the name of the directory
masterdir = '190709_Test/'
#masterdir = 'L_'+str(period)+'_New/'

print('-'*20)
print('Creating file structure', flush=True)
print('. . .')
ms.createstruct(ref,masterdir,subdirs)

print('Nondimensionalizing', flush=True)
print('. . . ',flush=True)
ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')

print('Getting data special values',flush=True)
print('. . . ',flush=True)
ms.dataspecialvals(masterdir,subdirs,'no')

print('Running simulations',flush=True)
print('. . . ',flush=True)
ms.runsims(SIMULATIONS,masterdir,subdirs,100,1,0)

print('Getting simulation special values',flush=True)
print('. . . ',flush=True)
ms.simspecialvals(masterdir,subdirs)

print('Redimensionalizing',flush=True)
print('. . . ',flush=True)
ms.redim(masterdir,subdirs)

print('Finding error',flush=True)
print('. . . ',flush=True)
ms.sberror(masterdir,subdirs)
ms.cqerror(masterdir,subdirs)

print('PROCESSING COMPLETE')
print()
print('*'*20)
print('ALL PROCESSING COMPLETE')
# os.system('say "Processing Complete"')
