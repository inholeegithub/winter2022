# Wrritten by In-Ho Lee, KRISS, February 20, 2020.
def get_nedos():
    afile=open('INCAR',"r")
    for line in afile:
        if len(line.split()) >= 2 :
          if line.split()[0].lower() == 'nedos':
             npt=int(line.split()[2])
    afile.close()
    print(npt)
    return npt

def get_fermi():
    afile=open('OUTCAR',"r")
    for line in afile:
        if len(line.split()) == 8 :
          if line.split()[0] == 'E-fermi':
             fermi=float(line.split()[2])
    afile.close()
    print(fermi)
    return fermi

def get_nstn():
    afile=open('KPOINTS',"r")
    jline=0
    for line in afile:
        jline=jline+1
        if jline == 2:
           nstn=int(line.split()[0])
           break
    afile.close()
    print(nstn)
    return nstn

import numpy as np
def get_nenknbandi():
    fermi=get_fermi()
    afile=open('EIGENVAL',"r")
    jline=0
    for line in afile:
        jline=jline+1
        if jline == 6:
           ne=int(line.split()[0])
           nk=int(line.split()[1])
           nbandi=int(line.split()[2])
           ei=np.zeros((nk,nbandi))
           w=np.zeros((nk,))
           ik=0
        if jline > 6:
           lvl=-1
           if len(line.split()) == 0:
              continue
           if len(line.split()) == 4:
              w[ik]=float(line.split()[3])
              continue
           if len(line.split()) == 3:
              lvl=int(line.split()[0])
              if lvl >= 1 or lvl <= nbandi:
                 ei[ik,lvl-1]=float(line.split()[1])
           if lvl == nbandi:
               ik=ik+1
    afile.close()
    for ik in range(nk):
        for lvl in range(nbandi):
            ei[ik,lvl]=ei[ik,lvl]-fermi
    fermi=0.
    print(ne,nk,nbandi)
    return ne,nk,nbandi,ei

ne,nk,nbandi,ei=get_nenknbandi()
nstn=get_nstn()
xgrd=[]
xxr=0.
for ik in range(nk):
    xxr=xxr+0.01
    xgrd.append(xxr)
xgrd=np.array(xgrd)
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator


plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_xlabel('',fontsize=22)
ax.set_ylabel('Energy (eV)',fontsize=22)
#majorLocator= MultipleLocator(10)
majorLocator= MultipleLocator(1)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%d')
minorFormatter= FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator= MultipleLocator(0.50)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%4.2f')
minorFormatter= FormatStrFormatter('%4.2f')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
ax.tick_params(which='major', length=2, color='black')
ax.tick_params(which='minor', length=4, color='brown')
plt.xticks([], [])
#plt.grid(True)
#ax.set_facecolor("ivory")
#ax.set_facecolor("cornsilk")
#ax.set_facecolor("gainsboro")
ax.set_facecolor("lightcyan")
#ax.set_facecolor("beige")
plt.xlim(0, 6 )
plt.ylim(-1, 1)
#plt.text(-3, 8.0, r"T'-WSe$_2$", {'color' : 'blue', 'fontsize' : 16} )
for ib in range(nbandi):
    plt.plot(xgrd,ei[:,ib], '-.', lw=1, alpha=0.9, color='blue')
plt.tight_layout()
plt.savefig('band.eps',dpi=1200)
plt.show()
#plt.close()
