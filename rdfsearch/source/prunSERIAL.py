import numpy as np
import os
if False:
    nfiles=100
    for i in range(nfiles):
        astring=str(i).zfill(4)
        bstring='POSCAR_'+astring
        os.system('cp POSCAR_0001 '+bstring)
if True:
    arr = os.listdir()
    inputs=[]
    i=0
    for i0 in range(len(arr)):
        astring=arr[i]
#       if astring[-4:] =='.txt' :
#           print(arr[i])
#       if astring[-3:] =='.py' :
#           print(arr[i])
        if astring[:7] =='POSCAR_' :
#           print(arr[i])
            inputs.append(astring)
        i=i+1
    inputs=sorted(inputs)
    nfiles=len(inputs)
#   print(inputs)
if True:
    ndir=20
    for i in range(ndir):
        astring=str(i).zfill(4)
        print(astring)
        bstring='mkdir '+astring
        os.system(bstring)
if True:
    for i in range(ndir):
        astring=str(i).zfill(4)
        bstring='cp SERIAL.pbs ./'+astring
        os.system(bstring)
        bstring='cp INCAR_rlx  ./'+astring
        os.system(bstring)
        bstring='cp INCAR_rlxall  ./'+astring
        os.system(bstring)
        bstring='cp INCAR_bs  ./'+astring
        os.system(bstring)
        bstring='cp POTCAR  ./'+astring
        os.system(bstring)
if True:
    dirlist=[ False for i in range(ndir)]
    i=0
    for k in range(nfiles):
        i0=int(np.mod(i,ndir))
        os.system('cp '+inputs[k]+' ./'+str(i0).zfill(4))
        dirlist[i0]=True
        i=i+1
if True:
    for i in range(ndir):
        astring=str(i).zfill(4)
        if dirlist[i] :
            bstring='cd ./'+astring+' ; ls * ; cd ../'
#           bstring='cd ./'+astring+' ; qsub SERIAL.pbs ; cd ../'
            os.system(bstring)
