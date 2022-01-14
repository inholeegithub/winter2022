# Wrritten by In-Ho Lee, KRISS, February 20, 2020.
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

def get_nedos():
    afile=open('INCAR_bs_dos',"r")
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

def database():
    fermi=get_fermi()
    npt=get_nedos()
    xgrd=[]
    data=[]
    kount=0
    koun1=0
    afile=open('DOSCAR',"r")
    for line in afile:
        kount=kount+1
#       if kount == 6 :
#          npt=int(line.split()[2])
#          print(npt)
        if kount >  6 :
           koun1=koun1+1
           xgrd.append(float(line.split()[0]))
           data.append(float(line.split()[1]))
#          print(float(line.split()[0]),float(line.split()[1]))
           if koun1 == npt:
              break
    afile.close()
#   sys.exit()
    xgrd=np.array(xgrd)
    data=np.array(data)
    xgrd[:]=xgrd[:]-fermi
    data[:]=data[:]/2.
    print('data from materials project (red)')
    for i in range(npt):
        print(xgrd[i],data[i])
    return xgrd,data

xgrd,xdata=database()

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
plt.plot(xgrd,xdata, '-.', lw=2, alpha=0.9, color='blue')
plt.tight_layout()
plt.savefig('dos.eps',dpi=1200)
plt.show()
#plt.close()
