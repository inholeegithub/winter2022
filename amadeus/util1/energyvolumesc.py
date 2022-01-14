"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/energyvolumesc.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import matplotlib.mathtext as mathtext

# Fixing random state for reproducibility
np.random.seed(19680801)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# the random data
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
x=np.asarray(volume)
y=np.asarray(energy)

nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]

# start with a rectangular Figure
plt.figure(1, figsize=(9, 8), facecolor=  "#ebf5ff", edgecolor=  "#ebf5ff")

axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx)
axHisty = plt.axes(rect_histy)

# no labels
axHistx.xaxis.set_major_formatter(nullfmt)
axHisty.yaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.scatter(x, y, color='blueviolet', s=0.9,  alpha=0.9)
axScatter.grid(color='silver', linestyle='-', linewidth=1)
axScatter.set_xlabel('volume (\AA$^3$)', fontsize=20)
axScatter.set_ylabel('energy (eV)', fontsize=20)

# now determine nice limits by hand:
binwidth = 0.25
binwidth = 1.00

# manipulate the numbers below: 550., 1000., -175., -150., 300, 200
tmp=np.amin(x)
tmq=np.amax(x)
axScatter.set_xlim((tmp,tmq))
tmp=np.amin(y)
tmq=np.amax(y)
axScatter.set_ylim((tmp,tmq))
axHistx.hist(x, bins=300, color='blueviolet')
axHisty.hist(y, bins=200, color='blueviolet', orientation='horizontal')
axHistx.grid(color='silver', linestyle='-', linewidth=1)
axHisty.grid(color='silver', linestyle='-', linewidth=1)

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())

plt.tight_layout()
plt.savefig("sc.pdf")
plt.show()
