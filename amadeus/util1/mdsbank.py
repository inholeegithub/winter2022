"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/mdsbank.py
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

def mds(d, ndimensions = 2):
    """
    Author: Jeremy M. Stober
    function: mds.py
    Description: Multidimensional Scaling
    Multidimensional Scaling - Given a matrix of interpoint distances,
    find a set of low dimensional points that have similar interpoint
    distances.
    """
    (n,n) = d.shape
    E = (-0.5 * d**2)
# Use mat to get column and row means to act as column and row means.
    Er = np.mat(np.mean(E,1))
    Es = np.mat(np.mean(E,0))
# From Principles of Multivariate Analysis: A User's Perspective (page 107).
    F = np.array(E - np.transpose(Er) - Es + np.mean(E))
    [U, S, V] = svd(F)
    Y = U * np.sqrt(S)
    return (Y[:,0:ndimensions], S)

def norm(vec):
    return np.sqrt(sum(vec**2))

def disthist(point,natoms,rcut):
    a11=point[3*natoms+0]
    a12=point[3*natoms+1]
    a13=point[3*natoms+2]
    a21=point[3*natoms+3]
    a22=point[3*natoms+4]
    a23=point[3*natoms+5]
    a31=point[3*natoms+6]
    a32=point[3*natoms+7]
    a33=point[3*natoms+8]
#   print(a11,a12,a13)
#   print(a21,a22,a23)
#   print(a31,a32,a33)
    k=int(rcut/0.1)+1
    dhist=np.zeros(k)
    for i in range(natoms-1):
        d1i=point[i*3+0]
        d2i=point[i*3+1]
        d3i=point[i*3+2]
        for j in range(i+1,natoms):
            d1=d1i-point[j*3+0]
            d2=d2i-point[j*3+1]
            d3=d3i-point[j*3+2]
            d1=d1-np.rint(d1)
            d2=d2-np.rint(d2)
            d3=d3-np.rint(d3)
            x=d1*a11+d2*a21+d3*a31
            y=d1*a12+d2*a22+d3*a32
            z=d1*a13+d2*a23+d3*a33
            r=np.sqrt(x*x+y*y+z*z)
            if r <= rcut :
               k=int(r/0.1)
               dhist[k]=dhist[k]+1.
    return dhist

n1=3
n2=3
n3=3
rcut=6.

kount1=0
with open('csa.in', 'r') as afile:
     for line in afile:
         kount1=kount1+1
         if kount1 == 1:
            nspecies=int(line.split()[0])
         if kount1 == 3:
             natoms=0
             for k in range(nspecies):
                 natoms=natoms+int(line.split()[k])
             natoms=natoms*(n1*n2*n3)
         if len(line.split()) ==  4:
             if line.split()[2] == 'npop,' and line.split()[3] == 'npop1':
                nsamples=int(line.split()[0])
         if len(line.split()) ==  5:
             if line.split()[2] == 'npop' and line.split()[4] == 'npop1':
                nsamples=int(line.split()[0])
print(nspecies,natoms,nsamples)
ninputd=3*natoms+9
points=np.zeros((nsamples,ninputd))
for isam in range(nsamples):
    filename='QOSCAR_'+str(isam+1).zfill(4)
    kount1=0
    with open(filename, 'r') as afile:
         for line in afile:
             kount1=kount1+1
             if kount1 == 2:
                 if len(line.split()) == 1 :
                    scale0=float(line.split()[0])
             if kount1 == 3:
                a11=float(line.split()[0])
                a12=float(line.split()[1])
                a13=float(line.split()[2])
                a11=a11*scale0
                a12=a12*scale0
                a13=a13*scale0
             if kount1 == 4:
                a21=float(line.split()[0])
                a22=float(line.split()[1])
                a23=float(line.split()[2])
                a21=a21*scale0
                a22=a22*scale0
                a23=a23*scale0
             if kount1 == 5:
                a31=float(line.split()[0])
                a32=float(line.split()[1])
                a33=float(line.split()[2])
                a31=a31*scale0
                a32=a32*scale0
                a33=a33*scale0
             if kount1 == 7:
                 nspecies=len(line.split())
                 natoms=0
                 for k in range(nspecies):
                     natoms=natoms+int(line.split()[k])
             if kount1 == 8:
                direct=line.split()[0]
                direct=direct.lower()
                direct=direct[0]
                if direct == 's':
                   kount1=kount1-1
             if kount1 == 8:
                direct=line.split()[0]
                direct=direct.lower()
                direct=direct[0]
                data=[]
             if kount1 > 8 and len(line.split()) > 0 :
                if direct == 'd':
                   d1=float(line.split()[0])
                   d2=float(line.split()[1])
                   d3=float(line.split()[2])
                   data.append(d1)
                   data.append(d2)
                   data.append(d3)
    data.append(a11)
    data.append(a12)
    data.append(a13)
    data.append(a21)
    data.append(a22)
    data.append(a23)
    data.append(a31)
    data.append(a32)
    data.append(a33)
    data=np.array(data)
    points[isam,:]=data[:]

mdistance = np.zeros((nsamples,nsamples))
valmax=-1.
for (i, pointi) in enumerate(points):
    veci=disthist(pointi,natoms,rcut)
    for (j, pointj) in enumerate(points):
        if i > j :
           vecj=disthist(pointj,natoms,rcut)
           mdistance[i,j] = norm(veci - vecj)
           mdistance[j,i] = mdistance[i,j] 
           if valmax <  mdistance[i,j] :
              valmax =  mdistance[i,j] 
mdistance=mdistance/valmax
Y, eigs = mds(mdistance)

# similarity map of the conformations, based on bond-length distribution
txtset=[]
for k in range(nsamples):
    txtset.append(k+1)
color=iter(cm.rainbow(np.linspace(0,1,nsamples)))
fig, ax = plt.subplots(figsize=(10,10))
plt.margins(0.1, 0.1)
#ax.plot(Y[:,0],Y[:,1],'D',lw=1, color='red', markersize=5, markerfacecoloralt='lawngreen')
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





