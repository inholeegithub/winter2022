# Written by In-Ho Lee, KRISS, September 25, 2020.
from os import system,getpid
import time
import math

def get_rp():
    afile=open("outdiff","r")
    for line in afile:
        if len(line.split()) == 4:
           if line.split()[2] == 'r,'  :
              rp=float(line.split()[0])
    afile.close()
    if math.isnan(rp):
       rp=0.
    return rp

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

alist=[32, 50, 54, 55, 100, 106, 117, 125, 127]
alist=[i for i in range(1,231)]
print(alist)
bfile=open('job_similarcrystal'+'_kill-9_'+str(getpid()),"w")
bfile.close()
# target crystal structure:  POSCAR_tmq
new=0
system('rm POSCAR_tmq_fin')
system('rm fort.7_tmq')
system('rm POSCAR_tmp_fin')
system('rm outdiff')
system("nice /home/ihlee/rdfsearch/rcry.x")
print('We will continue to use this file as input. rdfsearch.in')
#for i in range(230*1000*len(alist)):
for i in range(100000):
    system("nice /home/ihlee/rdfsearch/rcry.x >/dev/null")
    time.sleep(0.300)
    system("nice /home/ihlee/rdfsearch/fingerdiff.x > outdiff")
    time.sleep(0.300)
    rp=get_rp()
    if rp > 0.80:
       print('{:12.8f} a similar thing'.format(rp))
       jline=0
       afile=open('POSCAR_tmp',"r")
       for line in afile:
           jline=jline+1
           if jline == 1:
              ispgr=int(line.split()[0])
              if ispgr in alist:
                 new=new+1
                 system('cp POSCAR_tmp '+'POSCAR_'+str(ispgr)+'_'+str(new))
       afile.close()
    else:
       print('{:12.8f} a thing of no resemblance'.format(rp))
print(new)
system('rm POSCAR_tmq_fin')
system('rm POSCAR_tmp_fin')
system('rm fort.7_tmq')
system('rm outdiff')
