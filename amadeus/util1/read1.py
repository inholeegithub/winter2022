# Written by In-Ho Lee, KRISS, February 12, 2020.
import numpy as np
def datafile():
    afile=open('fort.1',"r")
    lines=afile.readlines()
    w=[]
    for line in lines:
        for i in range(len(line.split())):
            w.append(line.split()[i])
    afile.close()
    k=0
    ndeg=int(w[k])
    k=k+1
    npop=int(w[k])
    k=k+1
    npop1=int(w[k])
    k=k+1
    posi1_r=np.zeros((npop1,ndeg))
    for i in range(npop1):
        for j in range(ndeg):
            posi1_r[i,j]=float(w[k])
            k=k+1
    posi1_energy=np.zeros((npop1,))
    for i in range(npop1):
        posi1_energy[i]=float(w[k])
        k=k+1
    posi_r=np.zeros((npop,ndeg))
    for i in range(npop):
        for j in range(ndeg):
            posi_r[i,j]=float(w[k])
            k=k+1
    posi_energy=np.zeros((npop,))
    for i in range(npop):
        posi_energy[i]=float(w[k])
        k=k+1
    posi_best_r=np.zeros((ndeg,))
    for i in range(ndeg):
        posi_best_r[i]=float(w[k])
        k=k+1
    best_energy=float(w[k])
    k=k+1
    davg_r=float(w[k])
#   print(ndeg,npop,npop1)
#   print(posi1_r)
#   print(posi1_energy)
#   print(posi_r)
#   print(posi_energy)
#   print(posi_best_r)
#   print(best_energy)
#   print(davg_r)
    return ndeg,npop,posi_r,posi_energy
