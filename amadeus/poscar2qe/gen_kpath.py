import os
import numpy as np
import pyprocar
if os.path.isfile('KPOINTS'):
   os.system('rm KPOINTS')
if False:
   pyprocar.kpath('POSCAR','KPOINTS',40,True, 'hpkot',1.e-7,1.e-5,-1.0, np.eye(3))
if True:
   pyprocar.kpath('POSCAR','KPOINTS',40,True, 'hpkot',1.e-5,1.e-3,-1.0, np.eye(3))
