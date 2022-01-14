import pyprocar
pyprocar.repair('PROCAR','PROCAR_repaired')
pyprocar.bandsplot('PROCAR_repaired',outcar='OUTCAR',elimit=[-3,3],mode='plain',color='blue',kpointsfile='KPOINTS', savefig='bands.eps' )
