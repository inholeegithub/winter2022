"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/idenergy.py
"""
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext
import numpy as np

iid=[]
volume=[]
energy=[]
kount=0
with open("csa.out", 'r') as afile:
     for line in afile:
         if(len(line.split()) == 11 and line.split()[3] == 'outcar'):
#            iid.append(int(line.split()[0]))
             iid.append(kount)
             volume.append(float(line.split()[1]))
             energy.append(float(line.split()[2]))
             kount=kount+1

volume=np.array(volume)
energy=np.array(energy)
iid=np.array(iid)
#print kount
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw={'facecolor': "#ebf5ff"}, figsize=(7,6))
#ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='orange')
#ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='yellow')
#ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='cyan')
#ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='lawngreen')
#ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='forestgreen')
ax.plot(iid, energy, 'o', ms=4, lw=2, alpha=0.8, mfc='dodgerblue')
ax.grid()
# modify the energy limits below.
tmp=np.amin(energy)-1.
tmq=np.amax(energy)+0.
ax.set_ylim(tmp,tmq)
tmp=np.amin(iid)
tmq=np.amax(iid)
ax.set_xlim(tmp,tmq)
ax.set_xlabel('structure', fontsize=20)
ax.set_ylabel('energy (eV)', fontsize=20)
fig.savefig("id_energy.pdf")
plt.show()
