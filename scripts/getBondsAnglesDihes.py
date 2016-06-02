################################################
from MDAnalysis import *
import numpy as np
PSF   = '../../../data/cg_average-4-All-D.psf'
#PDB   = '../trajectories/ref_structures/stator_ps2-CA.pdb'
DCD   = '../../../trajectories/cg-2.dcd'
#PDB2  = '../trajectories/ref_structures/stator_ps3-CA.pdb'

# OUTOUT FILE
writeBonds   = open('bonds-2.txt',     'w')
writeAngles  = open('angles-2.txt',    'w')
writeDihes   = open('dihedrals-2.txt', 'w')

# Calculate the difference vector.
u = Universe(PSF,DCD)

# The Universe contains all the information of the system 
# including bonds, angles and dihedrals.
# MDAnalysis.core.AtomGroup.Universe()._psf --> the topology

bondlist     = u._psf['_bonds']
anglelist    = u._psf['_angles']
dihedrallist = u._psf['_dihe'] 

for ts in u.trajectory[:]:

   dr = []; an = []; dihe = []
   for bond in bondlist:
      ag = u.atoms[list(bond)]
      dr.append(ag.bond()) 
   
   for angle in anglelist:
      ag = u.atoms[list(angle)]
      an.append(ag.angle())   

   for dihedral in dihedrallist:
      ag = u.atoms[list(dihedral)]
      dihe.append(ag.dihedral())

   dr_strg = '' 
   for ii in dr: dr_strg += '{} '.format(ii)  

   an_strg = '' 
   for ii in an: an_strg += '{} '.format(ii)  
 
   dihe_strg = '' 
   for ii in dihe: dihe_strg += '{} '.format(ii)  
 
   writeBonds.write("{} {}\n".format(ts.frame, dr_strg))
   writeAngles.write("{} {}\n".format(ts.frame, an_strg))
   writeDihes.write("{} {}\n".format(ts.frame, dihe_strg))
   
exit()
