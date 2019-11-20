
This directory tree contains "msander", a "modern" version of parts of
sander.

This project is generally licensed under the GNU General Public License,
version 3 (GPL v3).  Some components use different, but compatible, open
source licenses.  See the LICENSE file for more information.

# Design goals:

* This project is a fork of the sander code in AmberTools.  It tries to
(greatly) simplify the code base, choosing the best and most useful parts of
the code, and to serve as a test bed for how modern Fortran coding techniques
can be used.  Key application areas are expected to be in structure
refinements using NMR, cryoEM or Xray diffraction information.

* Some pieces are missing from the sander program in AmberTools:

  * Things that should be easy to re-introduce later: emil, sebomd, pbsa, APBS

  * Things are are problably gone for good, but which don't represent the best
current program practice: Path-integral methods, thermostats that don't follow
the "middle" scheme, Berendsen barostat

  * Things that might be useful, but really complicate the code: evb
potentials, some parts of adaptive QM/MM, nudged elastic band

  * Non-periodic 3D-RISM has been removed for now, in an attempt to get the
simplest possible RISM code, perhaps as a basis for future GPU work.

* Key pieces of code that are still there, and being emphasized:

  * Periodic and non-periodic simulations, with all of Amber's GB models

  * QM/MM, including hooks to external codes

  * NMR, cryoEM and Xray restraints (including quite a bit of new code)

  * Thermodynamic integration and non-equilibrium sampling methods

  * Replica exchange capabilities, and constant pH and redox potential
simulations

# Building the code:

*Conda build:

   conda build recipe  [N.B: not yet tested!]

*Non-conda build:  (MacOSX, Linux)

   ./configure --help   #  then choose the options you want
   make install
   make test

