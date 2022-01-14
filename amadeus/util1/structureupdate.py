"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/structureupdate.py
"""
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext
import numpy as np

oldiid=[]
newiid=[]
oldenergy=[]
newenergy=[]
kount=0
with open("csa.out", 'r') as afile:
     for line in afile:
         if(len(line.split()) == 6 and line.split()[1] == 'type,'):
             kount=kount+1
             if(line.split()[0] == 'old'):
                   oldiid.append(kount)
                   oldenergy.append(float(line.split()[3]))
             if(line.split()[0] == 'new'):
                   newiid.append(kount)
                   newenergy.append(float(line.split()[3]))

oldiid=np.array(oldiid)
newiid=np.array(newiid)
oldenergy=np.array(oldenergy)
newenergy=np.array(newenergy)
#print kount
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw={'facecolor': "#ebf5ff"}, figsize=(7,6))
ax.plot(oldiid, oldenergy, 'o', ms=5, lw=2, alpha=0.7, mfc='orange')
ax.plot(newiid, newenergy, '*', ms=5, lw=2, alpha=0.7, mfc='blue')
ax.grid()
tmp=np.amin(oldiid)
if tmp > np.amin(newiid) :
   tmp = np.amin(newiid) 
tmq=np.amax(oldiid)
if tmq < np.amax(newiid) :
   tmq=np.amax(newiid)
ax.set_xlim(tmp,tmq)
tmp=np.amin(oldenergy)
if tmp > np.amin(newenergy) :
   tmp = np.amin(newenergy) 
tmq=np.amax(newenergy)
if tmq < np.amax(newenergy) :
   tmq = np.amax(newenergy) 
tmp=tmp-1.
if tmq > 1.e18:
    tmq=tmp+np.abs(tmp)*0.1
ax.set_ylim(tmp,tmq)
ax.set_xlabel('structure update', fontsize=20)
ax.set_ylabel('energy (eV)', fontsize=20)
fig.savefig("strup.pdf")
plt.show()
