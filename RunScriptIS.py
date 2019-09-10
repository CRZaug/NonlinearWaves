#import Master as ms
import MasterNEW_IS as ms
import os
from sys import exit
import numpy as np

ref = 'Reference Directory Structure 6 Sims/'

# The period in seconds. (3 h period)
period = 10800

# Delta f (SHOULD BE 1)
deltaf = 1 #mc/s

# The simulations to run. Put 'y' if it is to be run:
SIMULATIONS = ['n', #NLS
               'n', #dNLS
               'n', #Dysthe
               'n', #vDysthe
               'n', #dGT
               'y', #IS
                     ]
# SIMULATIONS = ['y', #NLS
#                'y', #dNLS
#                'y', #Dysthe
#                'y', #vDysthe
#                'y', #dGT
#                'y', #IS
#                      ]

# RUN THREE SEPARATE SIMULATION SETS, ONE FOR EACH SUB DIRECTORY GIVEN PRELIMINARY RESULTS


###################################
###################################
###################################
#           July
###################################
###################################
###################################

# #The data sets to run. Must be enclosed in brackets (AKA a list)
# subdirs = ['JulyData/']
# 
# # Define the name of the directory (should end in a /)
# #masterdir = 'IS_optimization/'
# #masterdir = 'L_'+str(period)+'_New/'
# 
# for bet in np.linspace(0,6e9,50):
# #for bet in np.linspace(-50,50,100):
# #for bet in np.array([1]):   
#     print(bet,flush = True)
#     masterdir = 'IS_optimization_bet_'+str(bet)[:11]+'/'
# 
#     
#     
#    # print('-'*20)
#    # print('Creating file structure', flush=True)
#    # print('. . .')
#     ms.createstruct(ref,masterdir,subdirs)
#     
#    # print('Nondimensionalizing', flush=True)
#    # print('. . . ',flush=True)
#     ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')
# 
#    # print('Getting data special values',flush=True)
#    # print('. . . ',flush=True)
#     ms.dataspecialvals(masterdir,subdirs,'no','go')
#     
#    # print('Running simulations',flush=True)
#    # print('. . . ',flush=True)
#     ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0)
#     
#    # print('Getting simulation special values',flush=True)
#    # print('. . . ',flush=True)
#     ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
#     
#    # print('Redimensionalizing',flush=True)
#    # print('. . . ',flush=True)
#     ms.redim(SIMULATIONS,masterdir,subdirs)
#     
#    # print('Finding error',flush=True)
#    # print('. . . ',flush=True)
#     ms.sberror(SIMULATIONS,masterdir,subdirs)
#     ms.cqerror(SIMULATIONS,masterdir,subdirs)
# # print('PROCESSING COMPLETE')
# # print()
# # # print('*'*20)
# # print('ALL PROCESSING COMPLETE')
# # os.system('say "Processing Complete"')


# 
# ###################################
# ###################################
# ###################################
# #           AUG 1
# ###################################
# ###################################
# ###################################
# 
# #The data sets to run. Must be enclosed in brackets (AKA a list)
# subdirs = ['Aug1Data/']
# 
# # Define the name of the directory (should end in a /)
# #masterdir = 'IS_optimization/'
# #masterdir = 'L_'+str(period)+'_New/'
# 
# 
# enddd = 1000000
# for bet in np.linspace(0,enddd,100):
#     print(bet)
#     masterdir = 'new_IS_optimization_bet_'+str(bet)[:7]+'/'
# 
#     
#    # print('-'*20)
#    # print('Creating file structure', flush=True)
#    # print('. . .')
#     ms.createstruct(ref,masterdir,subdirs)
#     
#    # print('Nondimensionalizing', flush=True)
#    # print('. . . ',flush=True)
#     ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')
# 
#    # print('Getting data special values',flush=True)
#    # print('. . . ',flush=True)
#     ms.dataspecialvals(masterdir,subdirs,'no','go')
#     
#    # print('Running simulations',flush=True)
#    # print('. . . ',flush=True)
#     ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0)
#     
#    # print('Getting simulation special values',flush=True)
#    # print('. . . ',flush=True)
#     ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
#     
#    # print('Redimensionalizing',flush=True)
#    # print('. . . ',flush=True)
#     ms.redim(SIMULATIONS,masterdir,subdirs)
#     
#    # print('Finding error',flush=True)
#    # print('. . . ',flush=True)
#     ms.sberror(SIMULATIONS,masterdir,subdirs)
#     ms.cqerror(SIMULATIONS,masterdir,subdirs)
#     input('stopping for now...')
#     print('PROCESSING COMPLETE')    
#     print()
# print('*'*20)
# print('ALL PROCESSING COMPLETE')
# os.system('say "Processing Complete"')

###################################
###################################
###################################
#           AUG 2
###################################
###################################
###################################

#The data sets to run. Must be enclosed in brackets (AKA a list)
subdirs = ['Aug2Data/']

#Define the name of the directory (should end in a /)
masterdir = 'TEST/'
masterdir = 'L_'+str(period)+'_New/'

#for bet in np.linspace(-15000,2000,200):
#for bet in np.linspace(-50,50,100):
for bet in np.array([1]):   

    #masterdir = 'new_IS_optimization_bet_'+str(bet)[:7]+'/'
    masterdir = 'TEST'+str(bet)[:7]+'/'


    
   # print('-'*20)
   # print('Creating file structure', flush=True)
   # print('. . .')
    ms.createstruct(ref,masterdir,subdirs)
    
   # print('Nondimensionalizing', flush=True)
   # print('. . . ',flush=True)
    ms.processnondim(masterdir,deltaf,period,subdirs,'no','no')

   # print('Getting data special values',flush=True)
   # print('. . . ',flush=True)
    ms.dataspecialvals(masterdir,subdirs,'no','go')
    
   # print('Running simulations',flush=True)
   # print('. . . ',flush=True)
    ms.runsims(SIMULATIONS,masterdir,subdirs,100,bet,0)
    
   # print('Getting simulation special values',flush=True)
   # print('. . . ',flush=True)
    ms.simspecialvals(SIMULATIONS,masterdir,subdirs)
    
   # print('Redimensionalizing',flush=True)
   # print('. . . ',flush=True)
    ms.redim(SIMULATIONS,masterdir,subdirs)
    
   # print('Finding error',flush=True)
   # print('. . . ',flush=True)
    ms.sberror(SIMULATIONS,masterdir,subdirs)
    ms.cqerror(SIMULATIONS,masterdir,subdirs)
print('PROCESSING COMPLETE')
# print()
# # # print('*'*20)
# # # print('ALL PROCESSING COMPLETE')
# # # os.system('say "Processing Complete"')


