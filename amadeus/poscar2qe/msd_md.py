import matplotlib
import matplotlib.pyplot as plt
import numpy as np
def cart2drt(b,uec):
    vec=np.matmul(b,uec)
    for i in range(3):
        vec[i]=vec[i]-round(vec[i])
        if vec[i] < 0.:
           vec[i]=vec[i]+1.
    return vec
def cart2drtunit(a1,a2,a3,uec):
    aa1=a1
    aa2=a2
    aa3=a3
    dvd=aa1[0]*aa2[1]*aa3[2]-aa1[1]*aa2[0]*aa3[2]-aa1[0]*aa2[2]*aa3[1]
    dvd=dvd+aa1[2]*aa2[0]*aa3[1]+aa1[1]*aa2[2]*aa3[0]-aa1[2]*aa2[1]*aa3[0]
    b=np.zeros((3,3))
    b[0,0]=-aa2[2]*aa3[1]+aa2[1]*aa3[2]
    b[1,0]= aa1[2]*aa3[1]-aa1[1]*aa3[2]
    b[2,0]=-aa1[2]*aa2[1]+aa1[1]*aa2[2]
    b[0,1]= aa2[2]*aa3[0]-aa2[0]*aa3[2]
    b[1,1]=-aa1[2]*aa3[0]+aa1[0]*aa3[2]
    b[2,1]= aa1[2]*aa2[0]-aa1[0]*aa2[2]
    b[0,2]=-aa2[1]*aa3[0]+aa2[0]*aa3[1]
    b[1,2]= aa1[1]*aa3[0]-aa1[0]*aa3[1]
    b[2,2]=-aa1[1]*aa2[0]+aa1[0]*aa2[1]
    b[:,:]=b[:,:]/dvd
    return cart2drt(b,uec)
afile=open('POSCAR','r')
line=afile.readline()
line=afile.readline()
ascale0=float(line.split()[0])
a1=np.zeros(3)
a2=np.zeros(3)
a3=np.zeros(3)
line=afile.readline()
a1[0]=float(line.split()[0])
a1[1]=float(line.split()[1])
a1[2]=float(line.split()[2])
line=afile.readline()
a2[0]=float(line.split()[0])
a2[1]=float(line.split()[1])
a2[2]=float(line.split()[2])
line=afile.readline()
a3[0]=float(line.split()[0])
a3[1]=float(line.split()[1])
a3[2]=float(line.split()[2])
afile.close()
print(a1)
print(a2)
print(a3)
afile=open('OUTCAR','r')
kount1=0
kount=0
natot=0
energy0=0.
lfirst=True
vec=np.zeros(3)
wec=np.zeros(3)
t=[]
y=[]
z=[]
for line in afile:
    if len(line.split()) == 12 and line.split()[8] == 'ions'  and line.split()[9] == 'NIONS'  :
       natot=int(line.split()[11])
       print(natot)
       drt=np.zeros((natot,3))
    if len(line.split()) == 6 and line.split()[0] == 'direct' and line.split()[1] == 'lattice' :
       kount1=kount1+1
    if len(line.split()) == 8 and line.split()[2] == 'entropy' and line.split()[6] == '='  :
       energy0=float(line.split()[7])
    if len(line.split()) == 3 and line.split()[0] == 'POSITION' and line.split()[1] == 'TOTAL-FORCE' :
       kount=kount+1
#      if kount > 2:
#         break
       if True:
          line=afile.readline()
          for i in range(natot):
              line=afile.readline()
              drt[i,0]=float(line.split()[0])          
              drt[i,1]=float(line.split()[1])          
              drt[i,2]=float(line.split()[2])          
              drt[i,:]=cart2drtunit(a1,a2,a3,drt[i,:])
              drt[i,0]=drt[i,0]-round(drt[i,0])
              drt[i,1]=drt[i,1]-round(drt[i,1])
              drt[i,2]=drt[i,2]-round(drt[i,2])
              if drt[i,0] < 0.:
                 drt[i,0]=drt[i,0]+1.
              if drt[i,1] < 0.:
                 drt[i,1]=drt[i,1]+1.
              if drt[i,2] < 0.:
                 drt[i,2]=drt[i,2]+1.
          if kount == 1:
             drt0=np.zeros((natot,3))
             for i in range(natot):
                 drt0[i,:]=drt[i,:]
       if True:
          rsq=0.
          for i in range(natot):
              vec[0]=drt[i,0]-drt0[i,0] 
              vec[1]=drt[i,1]-drt0[i,1] 
              vec[2]=drt[i,2]-drt0[i,2] 
              wec[0]=vec[0]*a1[0]+vec[1]*a2[0]+vec[2]*a3[0]           
              wec[1]=vec[0]*a1[1]+vec[1]*a2[1]+vec[2]*a3[1]           
              wec[2]=vec[0]*a1[2]+vec[1]*a2[2]+vec[2]*a3[2]           
              rsq=rsq+wec[0]**2+wec[1]**2+wec[2]**2
          rsq=rsq/natot
          z.append(rsq)
          t.append(kount*1.e-3)
          y.append(energy0)
afile.close()
print(natot)
print(kount1)
print(kount)
#print(t)
#print(y)
#print(z)

fig, ax = plt.subplots()
if True:
   ax.set(xlabel='Time (ps)', ylabel='MSD '+r'$(\mathrm{\AA})}$', title='')
   ax.plot(t, z)
if False:
   ax.set(xlabel='Time (fs)', ylabel='Energy (eV)', title='')
   ax.plot(t, y)
#ax.grid()
#fig.savefig("test.png")
plt.show()
