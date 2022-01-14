from os import listdir,system,getpid
from os.path import isfile, join
from pathlib import Path
from scipy.stats import pearsonr
import numpy as np
import random
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
from datetime import datetime
aa=[]
bb=[]
cc=[]
dd=[]
afile=open("POSCAR_tmq_fin","r")
ii=0
for line in afile:
    if len(line.split()) == 4:
       continue
    if len(line.split()) == 1:
       ii=ii+1
       aa.append(ii)
       bb.append(float(line.split()[0]))
afile.close()
afile=open("POSCAR_tmp_fin","r")
ii=0
for line in afile:
    if len(line.split()) == 4:
       continue
    if len(line.split()) == 1:
       ii=ii+1
       cc.append(ii)
       dd.append(float(line.split()[0]))
afile.close()
aa=np.array(aa)
bb=np.array(bb)
cc=np.array(cc)
dd=np.array(dd)
plt.figure(figsize=(6,6))
plt.rcParams['axes.linewidth'] = 1.4 #set the value globally
ax=plt.axes()
ax.set_xlabel('Index',fontsize=22)
ax.set_ylabel('Finger print (arb. unit)',fontsize=22)
majorLocator= MultipleLocator(200)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%d')
minorFormatter= FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator= MultipleLocator(0.02)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%4.2f')
minorFormatter= FormatStrFormatter('%4.2f')
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)
ax.tick_params(which='major', length=2, color='black')
ax.tick_params(which='minor', length=4, color='brown')
ax.set_facecolor("cornsilk")
ax.tick_params(which='major', length=2, color='black')
ax.tick_params(which='minor', length=4, color='brown')
ax.set_facecolor("cornsilk")
corr, _ = pearsonr(bb,dd)
ctsty0=corr
plt.text(200., 0.06, "Pearson's distance = "+f'{ctsty0: .3f}', {'color' : 'black', 'fontsize' : 16} )
plt.xlim(1., 800)
plt.ylim(0., 0.08)
plt.plot(aa,bb, '--', lw=1, alpha=0.8, color='green')
plt.plot(cc,dd, ':', lw=1, alpha=0.6, color='blue')
plt.tight_layout()
plt.savefig('finger_diff.eps')
plt.show()
plt.close()
