#
#************************************************************************
#                             AMBER                                    **
#                                                                      **
#                Copyright (c) 1986, 1991, 1995                        **
#             Regents of the University of California                  **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#

.SUFFIXES : .c .F .F90 .o
SHELL=/bin/sh
NSRC= 	nucgen.F90 ngfil.F90
NOBJ=	nucgen.o ngfil.o
ASRC=   ambpdb.F90 forFIRSTtype.F90  forFIRSTbond.F90 \
        forFIRSThbond.F90 forFIRSTteth.F90 \
        forFIRSToutp.F90  forFIRSTsugar.F90
AOBJ=   ambpdb.o forFIRSTtype.o forFIRSTbond.o forFIRSThbond.o \
                 forFIRSTteth.o forFIRSToutp.o forFIRSTsugar.o
include ../config.h

# trick nab to find header files in $(AMBER_PREFIX)/include
# At install time, AMBERHOME points to $(AMBER_SOURCE), which 
# does not have include folder (if using --prefix)
export AMBERHOME=$(AMBER_PREFIX)

PROGS= elsize$(SFX) molsurf$(SFX) resp$(SFX) \
       tinker_to_amber$(SFX) new_crd_to_dyn$(SFX) new_to_old_crd$(SFX) \
       fix_new_inpcrd_vel$(SFX)

PERLPROGS= process_mdout.perl process_minout.perl mdout2pymbar.pl

PYTHONPROGS= cpinutil.py ceinutil.py cpeinutil.py mdout_analyzer.py finddgref.py \
             fitpkaeo.py fixremdcouts.py genremdinputs.py softcore_setup.py charmmlipid2amber.py

SHELLPROGS= sgldinfo.sh sgldwt.sh

all:		$(PROGS) python

nukit$(SFX):	nukit.o
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o nukit nukit.o $(LDFLAGS) $(AMBERLDFLAGS)

nucgen$(SFX):	nucgen_modules.o $(NOBJ) syslib
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o nucgen$(SFX) \
		nucgen_modules.o $(NOBJ) ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

ambpdb$(SFX): 	$(AOBJ) syslib
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o ambpdb$(SFX) \
		$(AOBJ) ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

add_pdb$(SFX): add_pdb.o syslib
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o add_pdb$(SFX) \
		add_pdb.o ../lib/sys2.a $(LDFLAGS) $(AMBERLDFLAGS)

calcpka$(SFX):
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o calcpka$(SFX) calcpka.F90

elsize$(SFX): elsize.o
	@echo "[ETC]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o elsize$(SFX) elsize.o $(LM)

molsurf$(SFX):  molsurf.nab
	@echo "[ETC]  NAB $@"
	$(VB)$(BINDIR)/nab$(SFX) -o molsurf$(SFX) molsurf.nab

resp$(SFX):	lapack.o resp.o
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) $(LDFLAGS) $(AMBERLDFLAGS) \
		lapack.o resp.o -o resp$(SFX)

tinker_to_amber$(SFX): tinker_to_amber.o cspline.o
	@echo "[ETC]  FC $@"
	$(VB)$(FC) -o tinker_to_amber$(SFX) tinker_to_amber.o cspline.o

new_crd_to_dyn$(SFX): new_crd_to_dyn.o nxtsec
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(LDFLAGS) -o new_crd_to_dyn$(SFX) new_crd_to_dyn.o ../lib/nxtsec.o

new_to_old_crd$(SFX): new_to_old_crd.o nxtsec
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(LDFLAGS) -o new_to_old_crd$(SFX) new_to_old_crd.o ../lib/nxtsec.o

fix_new_inpcrd_vel$(SFX): fix_new_inpcrd_vel.o nxtsec
	@echo "[ETC]  FC $@"
	$(VB)$(FC) $(LDFLAGS) -o fix_new_inpcrd_vel$(SFX) fix_new_inpcrd_vel.o ../lib/nxtsec.o

nxtsec:
	cd ../lib; $(MAKE) nxtsec.o

syslib::
	cd ../lib && $(MAKE) sys2.a

install:	all
	/bin/mv $(PROGS) $(BINDIR)
	/bin/cp $(PERLPROGS) $(BINDIR)
	/bin/cp $(SHELLPROGS) $(BINDIR)

parallel:
	$(FC) $(FPPFLAGS) -o numprocs numprocs.F90
	/bin/cp numprocs ../../test
	/bin/cp numprocs ../../../test

python:
	$(VB)(if [ "$(SKIP_PYTHON)" = "no" ]; then \
		$(PYTHON) setup.py install -f $(PYTHON_INSTALL) --install-scripts=$(BINDIR) > pysetup.log;\
	fi;\
	)

clean:
	/bin/rm -f *.o 90 *.mod  $(PROGS) molsurf.c
	/bin/rm -rf build pysetup.log

uninstall:
	-(for file in $(PROGS) $(PERLPROGS) $(PYTHONPROGS) $(SHELLPROGS) ; do \
		/bin/rm -fr $(BINDIR)/$${file} ; \
	done)

#  default rules for Fortran and C compilation:

.c.o:
	@echo "[ETC]  CC $<"
	$(VB)$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

include depend
