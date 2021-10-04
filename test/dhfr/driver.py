#!/usr/bin/env python3

import sander
from parmed.amber.readparm import AmberParm, Rst7

# Initialize the topology object with coordinates
parm = AmberParm("prmtop_an")
rst = Rst7.open("md12.x")

# Set up input options 
inp = sander.pme_input()
inp.irism = 1
print("inp.rism is %d\n" % inp.irism)

qm_inp = sander.QmInputOptions()


rism_inp = sander.RismInputOptions()
rism_inp.solvcut = 10.0

sander.setup(parm, rst.coordinates, rst.box, inp, None, rism_inp)
# Compute the energies and forces
ene, frc = sander.energy_forces()

# Do whatever you want with the energies and forces
print( "Energy = %10.5f\n" % ene.tot )

# Free up our memory
sander.cleanup()
