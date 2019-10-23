
This directory tree contains "amber_phenix", a subset of AmberTools.  It
is somewhat experimental, but is intended to be a simple way for users to
get access to the Amber capabilities needed for phenix.

This project is generally licensed under the Lesser GNU General Public License,
version 3 (LGPL v3).  Some components use different, but compatible, open
source licenses.  See the LICENSE file for more information.

Note: As of phenix release 1.16, this package is a dependency of the
phenix build, and so is automatically included in the phenix release.
Hence, this package is for developers, not (regular) users.

** Overview:

*Conda build:

   conda build ./recipe_ap

*Non-conda build:  (MacOSX, Linux)

   export AMBERHOME=`pwd`
   ./configure_ambermini --no-netcdf --no-rism
   cd AmberTools/src && make -f Makefile.ambermini install

Requires C and Fortran compilers, bison, flex, and python development headers.
On ubuntu (and lubuntu, etc.) the following should work (many of these 
may already be installed, depending on which distribution you have):

     sudo apt-get install gfortran bison flex python-dev

On MacOSX, you need to install XCode and a fortran compiler.  Instructions 
are here:

     http://ambermd.org/Installation.php

