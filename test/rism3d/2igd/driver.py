#!/usr/bin/env python3

import sander
from parmed.amber.readparm import AmberParm, Rst7

# Initialize the topology object with coordinates
parm = AmberParm("2igd.parm7")
rst = Rst7.open("2igd.rst7")

# Set up input options 
inp = sander.pme_input()
inp.irism = 1

rism_inp = sander.RismInputOptions()
rism_inp.solvcut = 8.1
print("setting solvcut to %8.3f" % rism_inp.solvcut)
rism_inp.grdspc = 0.8
print("setting grdspc to %8.3f" % rism_inp.grdspc)

sander.setup(parm, rst.coordinates, rst.box, inp, rism_inp)

# Compute the energies and forces
ene, frc = sander.energy_forces()

# Do whatever you want with the energies and forces
print( "Energy    = %10.5f" % ene.tot )
print( "  vdw     = %10.5f" % ene.vdw )
print( "  elec    = %10.5f" % ene.elec )
print( "  vdw_14  = %10.5f" % ene.vdw_14 )
print( "  elec_14 = %10.5f" % ene.elec_14 )
print( "  bond    = %10.5f" % ene.bond )
print( "  angle   = %10.5f" % ene.angle )
print( "  dihed   = %10.5f" % ene.dihedral )
print( "  rism    = %10.5f" % ene.rism )

# Free up our memory
sander.cleanup()
