#import Master as ms
import MasterSept as ms
import os
from sys import exit
import numpy as np
 

ref = 'Reference Directory Structure 6 Sims/'



# Delta f (SHOULD BE 1)
deltaf = 1

# The period in seconds. (3 h period)
#period = 10800
period = 1/(0.001*deltaf)

# The data sets to run. Must be enclosed in brackets (AKA a list)
#subdirs = ['Aug1Data/','Aug2Data/','JulyData/']
subdirs = ['JulyData/']

# The simulations to run. Put 'y' if it is to be run:
SIMULATIONS = ['y', #NLS
               'y', #dNLS
               'n', #Dysthe
               'y', #vDysthe
               'y', #dGT
               'n', #IS
                    ]
i = 1
#for bet in np.linspace(-8300,8300,100):
for bet in np.array([0]):

    masterdir = 'New_Test_' + str(i)+ '/'
    #bet = -8263.91319566
    
    # print('-'*20)
    # print('Creating file structure', flush=True)
    # print('. . .')
    # ms.createstruct(ref,masterdir,subdirs)
    # 
    # print('Nondimensionalizing', flush=True)
    # print('. . . ',flush=True)
    # ms.processnondim(masterdir,deltaf,period,subdirs,'go','go')
    # 
    # print('Getting data special values',flush=True)
    # print('. . . ',flush=True)
    # ms.dataspecialvals(masterdir,subdirs,'no','no')
      
    # print('Running simulations',flush=True)
    # print('. . . ',flush=True)
    # ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0) #5/16 for aug2 worked
    # 
    print('Getting simulation special values',flush=True)
    print('. . . ',flush=True)
    ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
    
    print('Redimensionalizing',flush=True)
    print('. . . ',flush=True)
    ms.redim(SIMULATIONS,masterdir,subdirs)
    
    print('Finding error',flush=True)
    print('. . . ',flush=True)
    ms.sberror(SIMULATIONS,masterdir,subdirs)
    ms.cqerror(SIMULATIONS,masterdir,subdirs)
            
    print('PROCESSING COMPLETE')
    i+=1
print()
print('*'*20)
print('ALL PROCESSING COMPLETE')
    ##os.system('say "Processing Complete"'
    
    
#     
#     
#     
#     
#     
#     
#     
# subdirs = ['Aug2Data/']
# 
# # The simulations to run. Put 'y' if it is to be run:
# SIMULATIONS = ['n', #NLS
#                'n', #dNLS
#                'n', #Dysthe
#                'n', #vDysthe
#                'n', #dGT
#                'y', #IS
#                     ]
#     
# 
# i = 0
# for bet in np.linspace(-4100,4100,100):
#     masterdir = 'Ensemble_' + str(i)+ '/'
#     #bet = -4091.34567284
#     
#     print('-'*20)
#     print('Creating file structure', flush=True)
#     print('. . .')
#     ms.createstruct(ref,masterdir,subdirs)
#     
#     print('Nondimensionalizing', flush=True)
#     print('. . . ',flush=True)
#     ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')
#     
#     print('Getting data special values',flush=True)
#     print('. . . ',flush=True)
#     ms.dataspecialvals(masterdir,subdirs,'no','no')
#     
#     print('Running simulations',flush=True)
#     print('. . . ',flush=True)
#     ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0)
#     
#     print('Getting simulation special values',flush=True)
#     print('. . . ',flush=True)
#     ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
#     
#     print('Redimensionalizing',flush=True)
#     print('. . . ',flush=True)
#     ms.redim(SIMULATIONS,masterdir,subdirs)
#     
#     print('Finding error',flush=True)
#     print('. . . ',flush=True)
#     ms.sberror(SIMULATIONS,masterdir,subdirs)
#     ms.cqerror(SIMULATIONS,masterdir,subdirs)
#             
#     print('PROCESSING COMPLETE')
#     i+=1
# print()
# print('*'*20)
# print('ALL PROCESSING COMPLETE')
# 
# 




# 
# subdirs = ['Aug1Data/','Aug2Data/','JulyData/'] #this used to be only july
# 
# # The simulations to run. Put 'y' if it is to be run:
# SIMULATIONS = ['y', #NLS
#                'y', #dNLS
#                'n', #Dysthe
#                'y', #vDysthe
#                'y', #dGT
#                'n', #IS
#                     ]
#     
# 
# for i in range(0,100): 
#     masterdir = 'Linearized_Ensemble_' + str(i)+ '/'
#     #bet = 34891.04455223 # large amp beta
#     bet = 2.47023512e+09 # small amp beta
#     
#     print('-'*20)
#     print('Creating file structure', flush=True)
#     print('. . .')
#     ms.createstruct(ref,masterdir,subdirs)
#     
#     print('Nondimensionalizing', flush=True)
#     print('. . . ',flush=True)
#     ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')
#     
#     print('Getting data special values',flush=True)
#     print('. . . ',flush=True)
#     ms.dataspecialvals(masterdir,subdirs,'no','no')
#     
#     print('Running simulations',flush=True)
#     print('. . . ',flush=True)
#     ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0)
#     
#     print('Getting simulation special values',flush=True)
#     print('. . . ',flush=True)
#     ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
#     
#     print('Redimensionalizing',flush=True)
#     print('. . . ',flush=True)
#     ms.redim(SIMULATIONS,masterdir,subdirs)
#     
#     print('Finding error',flush=True)
#     print('. . . ',flush=True)
#     ms.sberror(SIMULATIONS,masterdir,subdirs)
#     ms.cqerror(SIMULATIONS,masterdir,subdirs)
#             
#     print('PROCESSING COMPLETE')
# 
#     print()
#     print('*'*20)
#     print('ALL PROCESSING COMPLETE')
    