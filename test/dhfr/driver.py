#!/usr/bin/env python3

import sander
from parmed.amber.readparm import AmberParm, Rst7

# Initialize the topology object with coordinates
parm = AmberParm("prmtop_an")
rst = Rst7.open("md12.x")

# Set up input options 
inp = sander.pme_input()
sander.setup(parm, rst.coordinates, rst.box, inp)

# Compute the energies and forces
ene, frc = sander.energy_forces()

# Do whatever you want with the energies and forces
print( "Energy = %10.5f\n" % ene.tot )

# Free up our memory
sander.cleanup()
