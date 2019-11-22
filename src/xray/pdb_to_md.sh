#!/bin/sh

#  Outline of a general routine to go from a PDB entry to an MD
#  simulation using X-ray structure factor restraints.

$pdbid='xxxx'

if false; then
#----------------------------------------------------------------------------
#  1.  get the structure-factor files, starting from ${pdbid}-sf.cif
#----------------------------------------------------------------------------

#  run a short (even zero-cycle) phenix refinement to get
#  ${pdbid}_data.mtz

phenix.mtz.dump -c -f s ${pdbid}_data.mtz | tr ',' '\t' > $pdbid.fmtz

#  edit $pdbid.fmtz: should be six fields: h k l Fobs SIGF-obs R-free-factor
#  remove the column headers, first line should contain the number of
#    reflections  (will try to simplify this soon).

fi

if false; then
#----------------------------------------------------------------------------
#  2.  run phenix.AmberPrep to get the unit cell Amber setup
#----------------------------------------------------------------------------

#  (not many details here...see the phenix instructions

fi

if false; then
#----------------------------------------------------------------------------
#  3.  (Optional): add waters to the simulation
#                  this could be done before step 2, above....
#----------------------------------------------------------------------------

UnitCell -p 4phenix_$pdbid.pdb -o ${pdbid}_uc.pdb

AddToBox -c ${pdbid}_uc.pdb -a spce.pdb -na 25000 -RP 2.8 -RW 2.8 \
   -G 0.1 -V 1 -o ${pdbid}_ucw.pdb

cat <<EOF > leap.in
set default nocenter on
set default reorder_residues off
source leaprc.protein.ff14SB
source leaprc.water.spce

x = loadPdb ${pdbid}_ucw.pdb
setBox x vdw 10.0

saveAmberParm x ${pdbid}_ucw.parm7 ${pdbid}_ucw.rst7a
quit
EOF

tleap -f leap.in

ChBox -c ${pdbid}_ucw.rst7a -o ${pdbid}_ucw.rst7 \
    -X 138.916 -Y 138.916 -Z 63.907 -al 90. -bt 90. -gm 120.

/bin/rm leap.in ${pdbid}_ucw.rst7a
fi

if false; then
#----------------------------------------------------------------------------
#  4.  prepare the prmtop revisions that are needed
#----------------------------------------------------------------------------

$AMBERHOME/bin/add_pdb -i ${pdbid}_ucw.parm7 -p ${pdbid}_ucw.pdb -o foo1.parm7
$AMBERHOME/bin/add_xray -i foo1.parm7 -o foo2.parm7
/bin/mv foo2.parm7 ${pdbid}_ucw.parm7
/bin/rm -f foo1.parm7

fi
