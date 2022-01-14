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
colors = []
afile=open("fort.11","r")
ii=0
for line in afile:
    if len(line.split()) == 4:
       continue
    if len(line.split()) == 2:
       ii=ii+1
       bb.append(float(line.split()[1]))
       aa.append(int(line.split()[0]))
       colors.append(float(line.split()[0])/230.)
afile.close()
aa=np.array(aa)
bb=np.array(bb)
colors=np.array(colors)
plt.figure(figsize=(6,6))
plt.rcParams['axes.linewidth'] = 1.4 #set the value globally
ax=plt.axes()
ax.set_xlabel('Index',fontsize=22)
ax.set_ylabel('Finger print (arb. unit)',fontsize=22)
majorLocator= MultipleLocator(30)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%d')
minorFormatter= FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator= MultipleLocator(0.10)
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
#corr, _ = pearsonr(bb,dd)
#ctsty0=corr
#plt.text(200., 0.06, "Pearson's distance = "+f'{ctsty0: .3f}', {'color' : 'black', 'fontsize' : 16} )
N=len(aa)
areas = (30 * np.random.rand(N))*0.0 +2.  # 0 to 15 point radii
#plt.scatter(aa, bb, s, c="g", alpha=0.5, marker=r'$\clubsuit$', label="Luck")
plt.scatter(aa, bb, s=areas, c=colors, cmap='hsv', alpha=0.8)
plt.ylabel("Correlation")
plt.xlabel("Space group")
plt.xlim(1., 231)
plt.ylim(0.6, 1.0)
#plt.plot(aa,bb, '--', lw=1, alpha=0.8, color='green')
plt.tight_layout()
plt.savefig('spg_rrfort11.eps')
plt.show()
plt.close()
