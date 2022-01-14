import os
import numpy as np
import pyprocar
if True:
#   pyprocar.dosplot(filename='vasprun.xml', mode='stack',colors=['green','orange','cyan'],items=dict(C=[0,1,2,3],S=[0,1,2,3],H=[0]), elimit=[-4, 3], orientation='vertical',orbitals=[0], plot_total=True)
   pyprocar.dosplot(filename='vasprun.xml', mode='plain', elimit=[-4, 3], spins=[0], orientation='vertical', savefig='dos.eps' )
