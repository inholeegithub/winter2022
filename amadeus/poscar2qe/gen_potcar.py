import os
aa=[]
afile=open("POSCAR","r")
ii=0
for line in afile:
    ii=ii+1
    if ii == 6 :
       nspecies=len(line.split())
       for j in range(nspecies):
           aa.append(line.split()[j])
afile.close()
print(aa)
for j in range(nspecies):
    if j == 0:
       cmd='cat /TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND/'+aa[j].strip()+'/POTCAR   > ./POTCAR'
       cmd='cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/'+aa[j].strip()+'/POTCAR   > ./POTCAR'
    else :
       cmd='cat /TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND/'+aa[j].strip()+'/POTCAR  >> ./POTCAR'
       cmd='cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/'+aa[j].strip()+'/POTCAR  >> ./POTCAR'
    print(cmd)
    returned_value = os.system(cmd)
