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

rr=np.zeros((200,))
aa=np.zeros((200,9))
colors = []
afile=open("fort.7","r")
ii=-1
for line in afile:
    if len(line.split()) >= 3:
       colors.append(ii)
       ii=ii+1
       jj=0
       continue
    if len(line.split()) == 1:
       continue
    if len(line.split()) == 2:
       rr[jj]=(float(line.split()[0]))
       aa[jj,ii]=(float(line.split()[1]))
       jj=jj+1
afile.close()
aa=np.array(aa)
colors=np.array(colors)
plt.figure(figsize=(8,8))
plt.rcParams['axes.linewidth'] = 1.4 #set the value globally
plt.rcParams['text.usetex'] = True
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
ax=plt.axes()
ax.set_xlabel('Distance'+' (\u00c5)',fontsize=22)
ax.set_ylabel('Fingerprint (arb. unit)',fontsize=22)
majorLocator= MultipleLocator(1)
minorLocator= AutoMinorLocator()
majorFormatter= FormatStrFormatter('%d')
minorFormatter= FormatStrFormatter('%d')
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
majorLocator= MultipleLocator(0.005)
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
#plt.text(10., 0.95, "Sn"+r'$_4$'+"Se"+r'$_4$',fontsize=26)
plt.text(2.5, 0.09, "Sn"+r'$_4$'+"Se"+r'$_4$',fontsize=26)
#plt.text(200., 0.06, "Pearson's distance = "+f'{ctsty0: .3f}', {'color' : 'black', 'fontsize' : 16} )
#N=9
#area = (30 * np.random.rand(N))*0.0 +1.  # 0 to 15 point radii
#plt.scatter(aa, bb, s, c="g", alpha=0.5, marker=r'$\clubsuit$', label="Luck")
#plt.scatter(aa, bb, s=area, c=colors, alpha=0.5)
plt.xlim(1., 6.)
plt.ylim(-0.005, 0.02)
if False:
   plt.plot(rr,aa[:,0], ':' , lw=1, alpha=0.8, color='blue')
   plt.plot(rr,aa[:,1], '-.', lw=1, alpha=0.8, color='green')
   plt.plot(rr,aa[:,2], '--', lw=1, alpha=0.8, color='red')
   plt.plot(rr,aa[:,3], '-',  lw=1, alpha=0.8, color='cyan')
   plt.plot(rr,aa[:,4], ':' , lw=1, alpha=0.8, color='magenta')
   plt.plot(rr,aa[:,5], '-.', lw=1, alpha=0.8, color='olive')
   plt.plot(rr,aa[:,6], '--', lw=1, alpha=0.8, color='black')
   plt.plot(rr,aa[:,7], '-' , lw=1, alpha=0.8, color='violet')
   plt.plot(rr,aa[:,8], ':' , lw=1, alpha=0.8, color='orange')
if True:
   plt.rc('legend',fontsize=18)
   plt.plot(rr,aa[:,0], ':' , lw=1, alpha=0.8, color='blue',    label=r'$G^1(r)$')
   plt.plot(rr,aa[:,1], '-.', lw=1, alpha=0.8, color='green',   label=r'$G^2(r)$')
   plt.plot(rr,aa[:,2], '--', lw=1, alpha=0.8, color='red',     label=r'$G^3(r)$')
   plt.plot(rr,aa[:,3], '-',  lw=1, alpha=0.8, color='cyan',    label=r'$G^4(r)$')
   plt.plot(rr,aa[:,4], ':' , lw=1, alpha=0.8, color='magenta', label=r'$G^5(r)$')
   plt.plot(rr,aa[:,5], '-.', lw=1, alpha=0.8, color='olive',   label=r'$G^6(r)$')
   plt.plot(rr,aa[:,6], '--', lw=1, alpha=0.8, color='black',   label=r'$G^7(r)$')
   plt.plot(rr,aa[:,7], '-' , lw=1, alpha=0.8, color='violet',  label=r'$G^8(r)$')
   plt.plot(rr,aa[:,8], ':' , lw=1, alpha=0.8, color='orange',  label=r'$G^9(r)$')
   plt.legend(loc='upper left', borderaxespad=0.13)
plt.tight_layout()
plt.savefig('gifort7.eps')
plt.show()
plt.close()
