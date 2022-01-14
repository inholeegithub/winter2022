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

def get_corr():
   nrref=200
   ninv=9
   aa=np.zeros((nrref,ninv))
   afile=open("fort.7","r")
   j=-1
   for line in afile:
       if len(line.split()) == 3:
          ii=0
          j=j+1
          continue
       if len(line.split()) == 1 and line.split() == '&' :
          continue
       if len(line.split()) == 2:
          aa[ii,j]=float(line.split()[1])
          ii=ii+1
   afile.close()
   bb=np.zeros((nrref,))
   cc=np.zeros((nrref,))
   for i in range(ninv):
       for j in range(ninv):
           if i > j :
              bb[:]=aa[:,i]        
              cc[:]=aa[:,j]        
              corr, _ = pearsonr(bb,cc)
              print(i,j,corr)

get_corr()
