"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
Pie chart showing portions related to discard/new/old types of searches.
$python ~/csa_vasp/util1/status.py
"""
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext
import numpy as np

iid=[]
volume=[]
energy=[]
kount1=0
with open("csa.out", 'r') as afile:
     for line in afile:
         if(len(line.split()) == 11 and line.split()[3] == 'outcar'):
#            iid.append(int(line.split()[0]))
             iid.append(kount1)
             volume.append(float(line.split()[1]))
             energy.append(float(line.split()[2]))
             kount1=kount1+1
volume=np.array(volume)
energy=np.array(energy)

oldiid=[]
newiid=[]
oldenergy=[]
newenergy=[]
kount2=0
kount21=0
kount22=0
with open("csa.out", 'r') as afile:
     for line in afile:
         if(len(line.split()) == 6 and line.split()[1] == 'type,'):
             kount2=kount2+1
             if(line.split()[0] == 'old'):
                   oldiid.append(kount2)
                   kount21=kount21+1
                   oldenergy.append(float(line.split()[3]))
             if(line.split()[0] == 'new'):
                   newiid.append(kount2)
                   kount22=kount22+1
                   newenergy.append(float(line.split()[3]))
oldiid=np.array(oldiid)
newiid=np.array(newiid)
oldenergy=np.array(oldenergy)
newenergy=np.array(newenergy)
print(kount1,kount2)
print(kount21,kount22)

kount21=round(kount21/float(kount1)*100.)
kount22=round(kount22/float(kount1)*100.)
kount11=round((kount1-kount2)/float(kount1)*100.)
print(kount21,kount22,kount11)

# Pie chart, where the slices will be ordered and plotted counter-clockwise:
labels = 'Discard', 'Old', 'New'
sizes = [kount11, kount21, kount22]
explode = (0.11, 0, 0)  # only "explode" the 2nd slice (i.e. 'Old')

fig1, ax1 = plt.subplots()
ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.2f%%', shadow=True, startangle=90, radius=0.5, textprops={'size': 'large'} )
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
fig1.savefig("percent.pdf")
plt.show()
