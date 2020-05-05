#!/usr/bin/env python
from __future__ import division

import os
import sys
try:
    import sander
    import parmed
    from parmed.amber.readparm import Rst7
    from parmed.utils.six.moves import range
    import numpy as np
except ImportError:
    print('Could not import sander. Skipping test')
    sys.exit(0)

# For duck-typing test with set_positions
class Vec3(tuple):

    def __new__(cls, a, b, c):
        return tuple.__new__(cls, (a, b, c))

def compare(computed, regression, desc):
    """
    Compares expected and computed values. Returns True if the comparison
    failed, False otherwise
    """
    if abs(computed - regression) > 2e-4:
        print("%s failed: Expected %15.4f got %15.4f" %
              (desc, regression, computed))
        return True
    return False

def compare_forces(computed, regression):
    """ Compares two forces for the same components """
    if abs(computed - regression) > 3e-4:
        print('Force comparison failed: %s vs. %s (%s)' % (computed, regression,
            abs(computed-regression)))
        return True
    return False

def print_result(failed):
    if failed:
        print("Possible FAILURE")
    else:
        print("PASSED")
    print("="*62)

# Run the first test
print('Testing GB sander interface')
options = sander.gas_input(7)
options.cut = 9999.0
options.saltcon = 0.2
options.gbsa = 1
rst = Rst7.open("../gb7_trx/trxox.2.4ns.x")
sander.setup("../gb7_trx/prmtop_an", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 631.8993, 'Bond')
failed = compare(e.angle, 898.2543, 'Angle') or failed
failed = compare(e.dihedral, 566.4453, 'Dihedral') or failed
failed = compare(e.vdw_14, 348.8246, '1-4 vdW') or failed
failed = compare(e.elec_14, 5980.5047, '1-4 Elec') or failed
failed = compare(e.vdw, -768.3629, 'van der Waals') or failed
failed = compare(e.elec, -7874.4913, 'Electrostatic') or failed
failed = compare(e.gb, -1943.0838, 'EGB') or failed
failed = compare(e.surf, 33.8338, 'SASA (GBSA)') or failed

print_result(failed)
assert len(f) == 3 * sander.natom()
sander.cleanup()

failed2 = failed
failed = False

# Run the second test
print('Testing PME sander interface')
options = sander.pme_input()
options.cut = 8.0
rst = Rst7.open('../4096wat/eq1.x')
sander.setup("../4096wat/prmtop", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0, "Bond") or failed
failed = compare(e.angle, 0.0, "Angle") or failed
failed = compare(e.dihedral, 0.0, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0, "1-4 Elec") or failed
failed = compare(e.vdw, 6028.9517, "van der Waals") or failed
failed = compare(e.elec, -45371.5995, "Electrostatic") or failed

print_result(failed)

assert len(f) == 3 * sander.natom()

failed2 = failed2 or failed
failed = False

# Run the third test
print('Checking for consistent forces')
fi = open("../4096wat/mdfrc_cmp.save", 'r')
try:
    line = fi.readline()
    while line and not line.strip().startswith('forces ='):
        line = fi.readline()

    for i in range(sander.natom()):
        i3 = i * 3
        frcs = [float(x) for x in fi.readline().split()]
        failed = failed or compare_forces(f[i3  ], frcs[0])
        failed = failed or compare_forces(f[i3+1], frcs[1])
        failed = failed or compare_forces(f[i3+2], frcs[2])
finally:
    fi.close()

print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run the 4th test
print('Testing the QM/MM non-periodic interface')
options = sander.gas_input(1)
options.cut = 99.0
options.ifqnt = 1
qmmm_options = sander.qm_input()
qmmm_options.iqmatoms[:3] = [8, 9, 10]
qmmm_options.qm_theory = "PM3"
qmmm_options.qmcharge = 0
qmmm_options.qmgb = 2
qmmm_options.adjust_q = 0

rst = Rst7.open("../qmmm2/lysine_PM3_qmgb2/lysine.crd")
sander.setup("../qmmm2/lysine_PM3_qmgb2/prmtop",
             rst.coordinates, rst.box, options, qmmm_options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0016, "Bond") or failed
failed = compare(e.angle, 0.3736, "Angle") or failed
failed = compare(e.dihedral, 0.0026, "Dihedral") or failed
failed = compare(e.vdw_14, 3.7051, "1-4 vdW") or failed
failed = compare(e.elec_14, 65.9137, "1-4 Elec") or failed
failed = compare(e.vdw, 0.1908, "van der Waals") or failed
failed = compare(e.elec, -4.1241, "Electrostatic") or failed
failed = compare(e.gb, -80.1406, "EGB") or failed
failed = compare(e.scf, -11.9100, "QM Escf") or failed


print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run the 5th test
print('Testing the QM/MM periodic interface (PM3-PDDG)')
options = sander.pme_input()
options.cut = 8.0
options.ifqnt = 1
options.jfastw = 4

qmmm_options = sander.QmInputOptions()
qmmm_options.qm_theory = "PDDG-PM3"
qmmm_options.qmmask = ":1-2"
qmmm_options.qmcharge = 0
qmmm_options.scfconv = 1e-10
qmmm_tight_p_conv = 1
qmmm_options.qmmm_int = 5

rst = Rst7.open("../qmmm2/MechEm_nma-spcfwbox/inpcrd")
sander.setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
             rst.coordinates, rst.box, options, qmmm_options)
e, f = sander.energy_forces()

failed = compare(e.bond, 605.7349, "Bond") or failed
failed = compare(e.angle, 331.7679, "Angle") or failed
failed = compare(e.dihedral, 0.0000, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0000, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0000, "1-4 Elec") or failed
failed = compare(e.vdw, 1281.8450, "van der Waals") or failed
failed = compare(e.elec, -7409.7167, "Electrostatic") or failed
failed = compare(e.scf, -37.1277, "QM Escf") or failed

print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run 6th test
print('Testing the QM/MM periodic interface (DFTB)')
if not os.path.exists('../../dat/slko/C-C.skf'):
    print('Could not find the SLKO files. Skipping this test.')
    print("="*62)
else:
    qmmm_options.qm_theory = "DFTB"
    qmmm_options.dftb_3rd_order = "NONE"
    sander.setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
                 rst.coordinates, rst.box, options, qmmm_options)
    e, f = sander.energy_forces()
    failed = compare(e.bond, 605.7349, "Bond") or failed
    failed = compare(e.angle, 331.7679, "Angle") or failed
    failed = compare(e.dihedral, 0.0000, "Dihedral") or failed
    failed = compare(e.vdw_14, 0.0000, "1-4 vdW") or failed
    failed = compare(e.elec_14, 0.0000, "1-4 Elec") or failed
    failed = compare(e.vdw, 1281.8450, "van der Waals") or failed
    failed = compare(e.elec, -7409.7167, "Electrostatic") or failed
    failed = compare(e.scf, -1209.0254, "QM Escf") or failed
    
    print_result(failed)
    sander.cleanup()

# Run the 7th test
print('Testing the broader API functionality')

failed2 = failed or failed2
failed = False

options = sander.gas_input(7)
options.cut = 9999.0
options.saltcon = 0.2
options.gbsa = 1
rst = Rst7.open("../gb7_trx/trxox.2.4ns.x")
sander.setup("../gb7_trx/prmtop_an", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 631.8993, 'Bond')
failed = compare(e.angle, 898.2543, 'Angle') or failed
failed = compare(e.dihedral, 566.4453, 'Dihedral') or failed
failed = compare(e.vdw_14, 348.8246, '1-4 vdW') or failed
failed = compare(e.elec_14, 5980.5047, '1-4 Elec') or failed
failed = compare(e.vdw, -768.3629, 'van der Waals') or failed
failed = compare(e.elec, -7874.4913, 'Electrostatic') or failed
failed = compare(e.gb, -1943.0838, 'EGB') or failed
failed = compare(e.surf, 33.8338, 'SASA (GBSA)') or failed

inpcrd = Rst7.open('../gb7_trx/trxox.2.4pns.x')
sander.set_positions(inpcrd.coordinates)
e, f = sander.energy_forces()

failed = compare(e.bond, 342.0589, "Bond") or failed
failed = compare(e.angle, 877.9924, "Angle") or failed
failed = compare(e.dihedral, 580.6551, "Dihedral") or failed
failed = compare(e.vdw_14, 357.5908, "1-4 vdW") or failed
failed = compare(e.elec_14, 5973.9713, "1-4 Elec") or failed
failed = compare(e.vdw, -770.7973, "van der Waals") or failed
failed = compare(e.elec, -7883.3799, "Electrostatic") or failed
failed = compare(e.gb, -1938.5736, "EGB") or failed
failed = compare(e.surf, 33.8944, "SASA (GBSA)") or failed

sander.cleanup()

options = sander.pme_input()

rst = Rst7.open('../4096wat/eq1.x')
sander.setup('../4096wat/prmtop', rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0, "Bond") or failed
failed = compare(e.angle, 0.0, "Angle") or failed
failed = compare(e.dihedral, 0.0, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0, "1-4 Elec") or failed
failed = compare(e.vdw, 6028.9517, "van der Waals") or failed
failed = compare(e.elec, -45371.5995, "Electrostatic") or failed

sander.set_box(51, 51, 51, 90, 90, 90)

e, f = sander.energy_forces()

failed = compare(e.bond, 0., "Bond") or failed
failed = compare(e.angle, 0., "Angle") or failed
failed = compare(e.dihedral, 0., "Dihedral") or failed
failed = compare(e.vdw_14, 0., "1-4 vdW") or failed
failed = compare(e.elec_14, 0., "1-4 Elec") or failed
failed = compare(e.vdw, 12567.4552, "van der Waals") or failed
failed = compare(e.elec, -42333.7294, "Electrostatic") or failed

box = sander.get_box()

failed = compare(box[0], 51.0, "a-vector length") or failed
failed = compare(box[1], 51.0, "b-vector length") or failed
failed = compare(box[2], 51.0, "c-vector length") or failed
failed = compare(box[3], 90.0, "alpha angle") or failed
failed = compare(box[4], 90.0, "beta angle") or failed
failed = compare(box[5], 90.0, "gamma angle") or failed

print_result(failed)

failed = False
print('Checking get_positions call')
coords = [Vec3(*xyz) for xyz in rst.coordinates.reshape((sander.natom(), 3))]
sander.set_positions(coords)

failed = not np.allclose(rst.coordinates.flatten(),
                     sander.get_positions(as_numpy=True)) or failed

print_result(failed)
sander.cleanup()

failed = False
print('Checking context manager functionality')
with sander.setup('../4096wat/prmtop', rst.coordinates, rst.box, options) as context:
    # Check that context.positions is a descriptor shortcut for get_positions
    failed = not np.allclose(rst.coordinates.flatten(),
                             np.array(context.positions).flatten()) or failed
    # Same for context.natom
    failed = context.natom != sander.natom() or failed
    # Check the context.positions setter
    context.positions = np.zeros((context.natom, 3))
    failed = not np.allclose(sander.get_positions(as_numpy=True),
            np.zeros(context.natom*3)) or failed
    # Check the unit cell
    failed = context.box != sander.get_box() or failed
    # Check the unit cell setter
    context.box = (51, 51, 51, 90, 90, 90)
    box = sander.get_box()
    failed = compare(box[0], 51.0, "a-vector length") or failed
    failed = compare(box[1], 51.0, "b-vector length") or failed
    failed = compare(box[2], 51.0, "c-vector length") or failed
    failed = compare(box[3], 90.0, "alpha angle") or failed
    failed = compare(box[4], 90.0, "beta angle") or failed
    failed = compare(box[5], 90.0, "gamma angle") or failed
    # Make sure sander.is_setup() is True
    failed = not bool(context) or failed
    failed = not sander.is_setup() or failed

# Make sure that closing the context manager cleaned everything up
failed = bool(context) or failed
failed = sander.is_setup() or failed

print_result(failed)

# Test belly
options = sander.pme_input()
options.cut = 8.0
options.ibelly = 1
options.bellymask = ':WAT'
rst = Rst7.open('../ubiquitin/inpcrd')

print('Testing belly = True')
with sander.setup("../ubiquitin/prmtop", rst.coordinates, rst.box, options):
  e, f = sander.energy_forces()
  failed = compare(e.tot, -53832.2117, "Potential energy") or failed
  print_result(failed)

print('Testing belly = False')
options.ibelly = 0
options.bellymask = ''
with sander.setup("../ubiquitin/prmtop", rst.coordinates, rst.box, options):
    e, f = sander.energy_forces()
    failed = compare(e.tot, -47567.4371, "Potential energy") or failed
    print_result(failed)

# Test restraintmask
print('Testing restraintmask')
options = sander.pme_input()
options.cut = 8.0
options.ntr = 1
options.refc = '../ubiquitin/inpcrd'
options.restraintmask = '!:WAT'
options.restraint_wt = 1.

rst = parmed.load_file('../ubiquitin/inpcrd')
with sander.setup("../ubiquitin/prmtop", rst.coordinates, rst.box, options):
    e, f = sander.energy_forces()
    failed = compare(e.tot, -47502.5121, "Potential energy") or failed
    print_result(failed)
