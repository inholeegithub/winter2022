from os import listdir
from os.path import isfile, join
from pathlib import Path
from scipy.stats import pearsonr
import numpy as np
import random
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from datetime import datetime
import time

def database():
    afile=open('OUTCAR',"r")
    for line in afile:
        if len(line.split()) == 8 :
          if line.split()[0] == 'E-fermi':
             fermi=float(line.split()[2])
    afile.close()
    print(fermi)
    xgrd=[]
    data=[]
    kount=0
    koun1=0
    afile=open('DOSCAR',"r")
    for line in afile:
        kount=kount+1
        if kount == 6 :
           npt=int(line.split()[2])
           print(npt)
        if kount >  6 :
           koun1=koun1+1
           xgrd.append(float(line.split()[0]))
           data.append(float(line.split()[1]))
           print(float(line.split()[0]),float(line.split()[1]))
           if koun1 == npt:
              break
    afile.close()
#   sys.exit()
    xgrd=np.array(xgrd)
    data=np.array(data)
    xgrd[:]=xgrd[:]-fermi
    data[:]=data[:]/2.
    print('data from materials project (red)')
    return xgrd,data

def database4():
    afile=open('../../b4/soc/OUTCAR',"r")
    for line in afile:
        if len(line.split()) == 8 :
          if line.split()[0] == 'E-fermi':
             fermi=float(line.split()[2])
    afile.close()
    print(fermi)
    xgrd=[]
    data=[]
    kount=0
    koun1=0
    afile=open('../../b4/soc/DOSCAR',"r")
    for line in afile:
        kount=kount+1
        if kount == 6 :
           npt=int(line.split()[2])
           print(npt)
        if kount >  6 :
           koun1=koun1+1
           xgrd.append(float(line.split()[0]))
           data.append(float(line.split()[1]))
           print(float(line.split()[0]),float(line.split()[1]))
           if koun1 == npt:
              break
    afile.close()
#   sys.exit()
    xgrd=np.array(xgrd)
    data=np.array(data)
    xgrd[:]=xgrd[:]-fermi
    data[:]=data[:]/2.
    print('data from materials project (red)')
    return xgrd,data

xgrd,xdata=database()
xgrd4,xdata4=database4()

plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_xlabel('Energy',fontsize=22)
ax.set_ylabel('Density of States'+' (1/eV/[f.u.])',fontsize=22)
#majorLocator= MultipleLocator(10)
majorLocator= MultipleLocator(1)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%d')
minorFormatter= FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator= MultipleLocator(1.00)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%4.2f')
minorFormatter= FormatStrFormatter('%4.2f')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
ax.tick_params(which='major', length=2, color='black')
ax.tick_params(which='minor', length=4, color='brown')
plt.grid(True)
#ax.set_facecolor("ivory")
#ax.set_facecolor("cornsilk")
#ax.set_facecolor("gainsboro")
ax.set_facecolor("lightcyan")
#ax.set_facecolor("beige")
plt.xlim(-3, 3 )
plt.ylim(-1, 15)
plt.text(-3, 8.0, r"T'-WSe$_2$", {'color' : 'blue', 'fontsize' : 16} )
plt.plot(xgrd,xdata, '-.', lw=2, alpha=0.9, color='brown')
plt.plot(xgrd4,xdata4, '-.', lw=2, alpha=0.9, color='gray')
plt.tight_layout()
plt.savefig('dos46.eps',dpi=1200)
plt.show()
#plt.close()
