"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/energyvolume.py
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
             iid.append(int(line.split()[0]))
             volume.append(float(line.split()[1]))
             energy.append(float(line.split()[2]))
             kount=kount+1

iid=np.array(iid)
volume=np.array(volume)
energy=np.array(energy)
#print kount
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw={'facecolor': "#ebf5ff"}, figsize=(7,6))
#ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='orange')
#ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='lawngreen')
#ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='lightgreen')
#ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='dodgerblue')
#ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='coral')
ax.plot(volume, energy, 'o', ms=5, lw=2, alpha=0.7, mfc='chartreuse')
ax.grid()
tmp=np.amin(volume)
tmq=np.amax(volume)
ax.set_xlim(tmp,tmq)
tmp=np.amin(energy)-1.
tmq=np.amax(energy)
ax.set_ylim(tmp,tmq)
ax.set_xlabel('volume (\AA$^3$)', fontsize=20)
ax.set_ylabel('energy (eV)', fontsize=20)
fig.savefig("energyvol.pdf")
plt.show()
