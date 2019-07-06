#import Master as ms
import MasterNEW as ms
import os
from sys import exit
 
############################
# IMPORTANT::::: DEL IS SET TO 0 FOR DNLS AND VDYSTHE!!!!!! 6TH ORDER OS

ref = 'Reference Directory Structure New/'



# The period in seconds. (3 h period)
period = 10800

# Constant df value
newn = 130

# Define the name of the directory
masterdir = 'L_18000_runforfigs_'+str(0)+'/'
#masterdir = 'L_'+str(period)+'_New/'
# 
# print('-'*20)
# print('Creating file structure', flush=True)
# print('. . .')
# ms.createstruct(ref,masterdir)
# 
print('Nondimensionalizing', flush=True)
print('. . . ',flush=True)
ms.processnondim(masterdir,newn,period,'no','no')
# 
# print('Getting data special values',flush=True)
# print('. . . ',flush=True)
# ms.dataspecialvals(masterdir,'yes')

# print('Running simulations',flush=True)
# print('. . . ',flush=True)
# ms.runsims(masterdir,100,0)
# 
# print('Getting simulation special values',flush=True)
# print('. . . ',flush=True)
# ms.simspecialvals(masterdir)

# print('Redimensionalizing',flush=True)
# print('. . . ',flush=True)
# ms.redim(masterdir)
# 
# print('Finding error',flush=True)
# print('. . . ',flush=True)
# ms.sberror(masterdir)
# ms.cqerror(masterdir)
# 
# print('PROCESSING COMPLETE')
# print()
# print('*'*20)
# print('ALL PROCESSING COMPLETE')
# #os.system('say "Processing Complete"')
