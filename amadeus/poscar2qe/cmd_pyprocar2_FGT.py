import os
import numpy as np
import pyprocar
if False:
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
   for j in range(nspecies):
       if j == 0:
#         cmd='cat /TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND/'+aa[j].strip()+'/POTCAR   > ./POTCAR'
          cmd='cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/'+aa[j].strip()+'/POTCAR   > ./POTCAR'
       else :
#         cmd='cat /TGM/Apps/VASP/POTCAR/2.POTPAW.PBE.54.RECOMMEND/'+aa[j].strip()+'/POTCAR  >> ./POTCAR'
          cmd='cat /TGM/Apps/VASP/POTCAR/1.POTPAW.LDA.54.RECOMMEND/'+aa[j].strip()+'/POTCAR  >> ./POTCAR'
       print(cmd)
       returned_value = os.system(cmd)
if False:
   if os.path.isfile('KPOINTS'):
      os.system('rm KPOINTS')
   os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_rlx  INCAR')
   os.system('sleep 0.3')
   os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
   os.system('cp CONTCAR  POSCAR')
   for kter in range(3):
       os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_rlxall  INCAR')
       os.system('sleep 0.3')
       os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
       os.system('cp CONTCAR  POSCAR')
       os.system('cp CONTCAR  POSCAR_')
if False:
   if os.path.isfile('KPOINTS'):
      os.system('rm KPOINTS')
   os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_mechanical_stability    INCAR')
   os.system('sleep 0.3')
   os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
   os.system('sleep 0.3')
   os.system('python /home/ihlee/csa_vasp/poscar2qe/MechElastic_V2.py -d=3D > output_mech')
if False:
   if os.path.isfile('KPOINTS'):
      os.system('rm KPOINTS')
   os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_rlxall0  INCAR')
   os.system('sleep 0.3')
   os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
   os.system('sleep 0.3')
   pyprocar.kpath('POSCAR','KPOINTS',40,True, 'hpkot',1.e-7,1.e-5,-1.0, np.eye(3))
   os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_bands  INCAR')
   os.system('sleep 0.3')
   os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
   os.system('sleep 0.3')
if False:
#  pyprocar.repair('PROCAR','PROCAR_repaired')
   pyprocar.filter('PROCAR','PROCAR_up',spin=[0])
   pyprocar.filter('PROCAR','PROCAR_dn',spin=[1])
   pyprocar.bandsplot('PROCAR_up',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[0] , savefig='sbands_up.eps')
   pyprocar.bandsplot('PROCAR_dn',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[0] , savefig='sbands_dn.eps')
   pyprocar.bandsplot('PROCAR_up',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[1,2,3] , savefig='pbands_up.eps')
   pyprocar.bandsplot('PROCAR_dn',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[1,2,3] , savefig='pbands_dn.eps')
   pyprocar.bandsplot('PROCAR_up',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[4,5,6,7,8] , savefig='dbands_up.eps')
   pyprocar.bandsplot('PROCAR_dn',outcar='OUTCAR',elimit=[-2,2],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[4,5,6,7,8] , savefig='dbands_dn.eps')
#  pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],mode='plain',color='blue',kpointsfile='KPOINTS' )
#  pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[0] )
#  pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[1,2,3] , savefig='bands.eps')
#  pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', orbitals=[4,5,6,7,8] )
if False:
   pyprocar.repair('PROCAR','PROCAR_repaired')
   pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', atoms=[0,1] )
   pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', atoms=[2,3] )
   pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', atoms=[4,5,6,7,8,9,10,11,12,13,14,15] )
if False:
   pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', atoms=[0],orbitals=[0] )
   pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-4,3],kpointsfile='KPOINTS',cmap='jet',mode='parametric', atoms=[1],orbitals=[1,2,3] )
if False:
   if os.path.isfile('KPOINTS'):
      os.system('rm KPOINTS')
   os.system('cp /home/ihlee/csa_vasp/poscar2qe/INCAR_dos_tetrahedron  INCAR')
   os.system('sleep 0.3')
   os.system('mpirun -np 8 /TGM/Apps/VASP/bin/5.4.4/O3/NORMAL/vasp.5.4.4.pl2.O3.NORMAL.ncl.x > stdout.log')
   os.system('sleep 0.3')
if False:
   pyprocar.repair('PROCAR','PROCAR_repaired')
   pyprocar.bandgap(procar="PROCAR_repaired", outcar="OUTCAR", code="vasp")
   pyprocar.dosplot(filename='vasprun.xml', mode='stack',colors=['green','orange','cyan'],items=dict(Fe=[4,5,6,7,8],Ge=[0,1,2,3],Te=[0,1,2,3]), elimit=[-2, 2], orientation='vertical', plot_total=True, savefig='dos.eps')
#  pyprocar.dosplot(filename='vasprun.xml', mode='plain', elimit=[-2, 2], spins=[0], orientation='vertical' )
#  pyprocar.dosplot(filename='vasprun.xml', mode='plain', elimit=[-2, 2], orientation='vertical' )
#  pyprocar.dosplot(filename='vasprun.xml', mode='parametric', elimit=[-2, 2], orientation='vertical',orbitals=[0],atoms=[4,5,6,7,8,9,10,11,12,13,14,15] )
#  pyprocar.dosplot(filename='vasprun.xml', mode='parametric', elimit=[-2, 2], orientation='vertical',orbitals=[0], plot_total=True)
#  pyprocar.dosplot(filename='vasprun.xml', mode='parametric', elimit=[-2, 2], orientation='vertical',orbitals=[0], plot_total=True)
#  pyprocar.dosplot(filename='vasprun.xml', mode='stack',colors=['green','orange','cyan'],items=dict(C=[0,1,2,3],S=[0,1,2,3],H=[0]), elimit=[-4, 3], orientation='vertical',orbitals=[0], plot_total=True, savefig='dos.eps')
#  pyprocar.dosplot(filename='vasprun.xml', mode='parametric', elimit=[-2, 2], orientation='vertical',atoms=[4,5,6,7,8,9,10,11,12,13,14,15] )
#  pyprocar.dosplot(filename='vasprun.xml', mode='parametric', elimit=[-2, 2], orientation='vertical',orbitals=[1,2,3] )
if True:
   pyprocar.repair('PROCAR','PROCAR_repaired')
   pyprocar.fermi3D(procar='PROCAR_repaired',outcar='OUTCAR',mode='plain')
if False:
   pyprocar.repair('PROCAR','PROCAR_repaired')
   pyprocar.filter('PROCAR','PROCAR_up',spin=[0])
   pyprocar.filter('PROCAR','PROCAR_dn',spin=[1])
   pyprocar.fermi3D(procar='PROCAR_up',outcar='OUTCAR',mode='plain' )
