"""
A part of AMADEUS protocol
Written by In-Ho Lee, KRISS, July 14, 2018.
$python ~/csa_vasp/util1/abcalbega.py
"""
import matplotlib.pyplot as plt
import matplotlib.mathtext as mathtext
import numpy as np

iid=[]
alpha=[]
beta=[]
gamma=[]
aa=[]
bb=[]
cc=[]
kount1=0
with open("csa.out", 'r') as afile:
     for line in afile:
         if(len(line.split()) == 11 and line.split()[3] == 'outcar'):
#            iid.append(int(line.split()[0]))
             iid.append(kount1)
             alpha.append(float(line.split()[7]))
             beta.append(float(line.split()[8]))
             gamma.append(float(line.split()[9]))
             aa.append(float(line.split()[4]))
             bb.append(float(line.split()[5]))
             cc.append(float(line.split()[6]))
             kount1=kount1+1
alpha=np.array(alpha)
beta=np.array(beta)
gamma=np.array(gamma)
aa=np.array(aa)
bb=np.array(bb)
cc=np.array(cc)

fig, axs = plt.subplots(2, 3, sharey=True, tight_layout=True)
n_bins=20
# We can set the number of bins with the `bins` kwarg
axs[0,0].hist(aa, bins=n_bins)
axs[0,0].set_xlabel(r'${\it a}$ $ ({\rm  \AA})$')
axs[0,1].hist(bb, bins=n_bins)
axs[0,1].set_xlabel(r'${\it b}$ $ ({\rm  \AA})$')
axs[0,2].hist(cc, bins=n_bins)
axs[0,2].set_xlabel(r'${\it c}$ $ ({\rm  \AA})$')
axs[1,0].hist(alpha, bins=n_bins)
axs[1,0].set_xlabel(r'$\alpha$ (deg.)')
axs[1,1].hist(beta, bins=n_bins)
axs[1,1].set_xlabel(r'$\beta$ (deg.)')
axs[1,2].hist(gamma, bins=n_bins)
axs[1,2].set_xlabel(r'$\gamma$ (deg.)')

fig.savefig("abcdist.pdf")
plt.show()
