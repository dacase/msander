# UPDATE: December 2024:
This repository has now been merged into
`https://github.com/Amber-MD/AmberClassic.`  Futher development will take
place there.  I will leave this repository here for now, but at some point
it will go away, to avoid confusion.

# Overview

This directory tree contains `msander`, a "modern" version of parts
of sander.  Also included are various X-ray and cryoEM-related code
and utilities, as well as versions of a number of the “classic”
(and most-used) parts of AmberTools: `tleap, antechamber, sqm, metatwist`
and `paramfit`. With these tools, many systems can be set up for simulation in
msander.

The documentation and authorship credits are in the *doc/msander.pdf* file.

# Warning

This is a work in progress, and may not always be in a stable state
(although that is my goal for the main branch).  I may not be able to
respond to requests for support, but please create a github issue if you
have comments or suggestions.

This code is probably most useful to those who are already familiar with
AmberTools, and access to either AmberTools (https://ambermd.org), or at
least to `cpptraj` (https://github.com/Amber-MD) will be needed by most
users.

# Design goals

* This project began as a fork of the `sander` code in `AmberTools`.  
It tries to (greatly) simplify the code base, choosing the best and 
most useful parts of the code, and to serve as a test bed for how 
modern Fortran coding techniques can be used.  Key application areas 
are expected to be in structure refinements using NMR, cryoEM or 
Xray diffraction information.  This version has a fair amount of OpenMP
support, especially for Xray and 3D-RISM calculations.  Parts of the Xray
code uses GPU acceleration.

* One additional goal of this collection is to make compiling and
installation as simple as possible. There is a pretty simple configure
script, and minimal dependencies on external packages.

# Key differences in functionality versus sander

* Some pieces are missing from the sander program in AmberTools:

  * Things that should be easy to re-introduce later: emil, sebomd, pbsa, APBS

  * Things are are problably gone for good, but which don't represent the best
current practice: Path-integral methods, thermostats that don't follow
the "middle" scheme, Berendsen barostat

  * Things that might be useful, but really complicate the code: evb
potentials, QM/MM, nudged elastic band, constant pH
and constant redox potential simulations.  The API interface has also been
removed.

  * Non-periodic 3D-RISM has been removed for now, in an attempt to get the
simplest possible RISM code, perhaps as a basis for future GPU work.

(If you need some of these deleted pieces, use *sander* from AmberTools
instead.)

* Key pieces of code that are still there, and being emphasized:

  * Periodic and non-periodic simulations, with all of Amber's GB models

  * 3D-RISM in periodic boundary conditions

  * NMR, cryoEM and Xray restraints (including quite a bit of new code; Xray
    restraints include NVIDIA GPU-enabled capabilities)

  * Thermodynamic integration and non-equilibrium sampling methods,
    including adaptively biassed sampling and self-guided Langevin dynamics

  * Sampling and minimization using the lmod and xmin approaches; these
    can now be used in conjunction with SHAKE and SETTLE.

  * Replica exchange capabilities, except for constant pH and redox potential
    simulations

# Execution speed

Force field evaluation is still slow compared to many other codes.  This
project thus focusses on systems where other parts of the simulation (such
as QM, RISM, or Xray/Cryoem restraints) are the bottleneck, so that force
field speed is not the limiting component.

# Building the code

* MacOSX, Linux, probably WSL:
```
   ./configure --help   #  then choose the options you want
   make install
   make test
```

# License

This project is licensed under the GNU General Public License, 
version 2, or (at your option) any later version.   Some components use 
different, but compatible, open source licenses.  See the LICENSE file 
for more information.

