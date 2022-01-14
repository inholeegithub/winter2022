"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/mds0.py
To visualize conformations in two-dimensional space, multidimensional scaling (MDS), was performed. 
In this way, relative descriptor-distances are maintained, 
and an accurate two-dimensional representation of the many-dimensional conformations is obtained. 
Thus, distances between points in the plot quantify the similarity between the various conformations.
The closer, the more similar.
"""
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

kount=-1
isam=-1
nsamples=-1
with open('outputreadfort1', 'r') as afile:
     for line in afile:
         if len(line.split()) == 4 and 'original' == line.split()[0] and 'data' == line.split()[1] :
            if isam == -1 and kount == -1:
               nsamples=int(line.split()[2]) 
               m=int(line.split()[3]) 
               print(nsamples,m)
               Y=np.zeros((nsamples,m+1)) 
               continue
         if len(line.split()) > 0 and 'VBLE.' == line.split()[0] and 'PROJ-1' == line.split()[1] :
            kount=0
            isam=0
            continue
         if len(line.split()) > 0 and '------' == line.split()[0] :
            continue
         if len(line.split()) > 0 and  kount ==  0 and isam >= 0 :
                isam=isam+1
#               print(line.split()[0])
                for j in range(7):
                    Y[isam-1,j]=float(line.split()[j+1])
         if isam ==  nsamples and isam > 0:
            kount=-1
            break

# similarity map of the conformations, based on bond-length distribution
txtset=[]
for k in range(nsamples):
    txtset.append(k+1)

color=iter(cm.rainbow(np.linspace(0,1,nsamples)))
fig, ax = plt.subplots(figsize=(10,10))
plt.margins(0.1, 0.1)
for k in range(nsamples):
    c=next(color)
    s=(nsamples-k)
    ax.plot(Y[k,0],Y[k,1],'o',lw=1, color=c, markersize=s, alpha=0.3)
for k, txt in enumerate(txtset):
#   ax.annotate(txt,(Y[k,0],Y[k,1]))
    fs=12
    if k+1 > 10:
       fs=10
    if k+1 > 20:
       fs=8
    if k+1 > 30:
       fs=6
    if k+1 > 40:
       fs=4
    plt.text(Y[k,0],Y[k,1], "{}".format(txt) ,fontsize=fs)
ax.set_aspect('equal', 'box')
ax.grid()
fig.tight_layout()
fig.savefig('mdstest.pdf')
plt.show()

