import csv
from PY_INPUT import *
from datetime import datetime
from scipy.spatial import distance
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

##########################################################################
######################## RUN FIRST SIMULATION ############################
##########################################################################

### RUN SIMULATION UNTIL FIRST FAILURE

now = datetime.now().time() # time object
print("First simulation:", now)

if geometry=="triangular":
    jobname='T'+str(int(reldensity*1000))+'_S'+str(nx)+'_M'+str(mesh)
if geometry=="kagome":
    jobname="K"+str(int(reldensity*1000))+"_S"+str(nx)+"_M"+str(mesh)

stepNR=1
with open('Manager.txt', 'w') as doc:
    doc.write(str(stepNR)+','+str(cracksize)+','+str(0)+','+str(0)+'\n')
#
#os.system('abaqus job='+jobname+'_step1.inp interactive cpus=3 double')

for i in range(50):
    now = datetime.now().time() # time object
    print("Sumulation "+str(i+1)+":", now)
    ##########################################################################
    ########################### INTERPRET RESULTS ############################
    ##########################################################################
    
    
    os.system('python PY_CRACK_UPDATE.py')
    
    
    ## python ~/code/hw01/script.py
    
    ##########################################################################
    ########################### CREATE NEW INPUT ############################
    ##########################################################################
    
    
    os.system('python PY_GEN_new_step.py')
    
    ##########################################################################
    ####################### RESTART WITH NEW STEP ############################
    ##########################################################################
    
    stepNR=stepNR+1
#    os.system('abaqus int cpus=3 job='+jobname+'_step'+str(stepNR)+' oldjob='+jobname+'_step'+str(stepNR-1))
    #abaqus interactive cpus=$NSLOTS mp_mode=mpi job=$INPUT.$JOB_ID double input=$INPUT \
    #       scratch=$ABAQUS_PARALLELSCRATCH $ABAQUS_ARGS
    #abaqus job=s6mod-step2 oldjob=s6mod cpus=$RESCALE_CORES_PER_SLOT mp_mode=mpi interactive
    
    os.system('abaqus restartjoin originalodb='+jobname+'_step1.odb restartodb='+jobname+'_step'+str(stepNR)+'.odb')