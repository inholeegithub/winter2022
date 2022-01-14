"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/energyhistogram.py
"""
import matplotlib.mathtext as mathtext
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path


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
#fig, ax = plt.subplots()
fig, ax = plt.subplots(subplot_kw={'facecolor': "#ebf5ff"}, figsize=(6,6))
#print kount
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

data=energy
# histogram our data with numpy

n, bins = np.histogram(data, 100)

# get the corners of the rectangles for the histogram
left = np.array(bins[:-1])
right = np.array(bins[1:])
bottom = np.zeros(len(left))
top = bottom + n

# we need a (numrects x numsides x 2) numpy array for the path helper
# function to build a compound path
XY = np.array([[left, left, right, right], [bottom, top, top, bottom]]).T

# get the Path object
barpath = path.Path.make_compound_path_from_polys(XY)

# make a patch out of it
patch = patches.PathPatch(barpath)
ax.add_patch(patch)

# update the view limits
ax.set_xlim(left[0], right[-1])
#ax.set_xlim(-175,-140)
ax.set_ylim(bottom.min(), top.max())

ax.set_xlabel('energy (eV)', fontsize=20)
fig.savefig("energyhist.pdf")
plt.show()
