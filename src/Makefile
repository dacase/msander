#   top-level Makefile for msander

include ../config.h

install: $(INSTALLTYPE)

serial:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd msander && $(MAKE) $(LIBSANDER)
	@cd msander && $(MAKE) $(LIBSANDERLES)
	@cd pysander && $(MAKE) $(PYSANDER)
	@cd parmed && $(MAKE) $(PARMED)
	@cd xray && $(MAKE) xray
	@cd xtalutil && $(MAKE) xtalutil
	@cd addles && $(MAKE) addles
	@cd leap && $(MAKE) tleap
	@cd antechamber && $(MAKE) ac
	@cd ambpdb && $(MAKE) ambpdb
	@cd sqm && $(MAKE) sqm
	@cd moft && $(MAKE) $(MOFT)
#	@cd sff && $(MAKE) libsff
#	@cd nabc && $(MAKE) libnabc
#	@cp pdb4amber/pdb4amber $(BINDIR)
#	@cd gmml && $(MAKE) gmml

parallel:
	@cd msander && $(MAKE) msander
#	@cd msander && $(MAKE) msander.LES

cuda:
	@cd msander && $(MAKE) msander.cuda

mpicuda:
	@cd msander && $(MAKE) msander.cuda$(SFX)

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )
	-(cd xray/cuda && $(MAKE) clean )
	-(cd addles && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd ambpdb && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd moft && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nabc && $(MAKE) clean )
#	-(cd gmml && $(MAKE) clean )

uninstall: clean
	-(cd $(BINDIR) && /bin/rm -f msander msander.MPI \
        add_pdb add_xray teLeap tleap addles pdb4amber)
	-(cd $(LIBDIR) && /bin/rm -f liblapack.a libblas.a \
        libarpack.a libsqm.a librism.a /libsys.a \
        libsander.so libsanderles.so)
	/bin/rm -f $(INCDIR)/* ../config.h
	/bin/rm -rf $(BINDIR)/x $(LIBDIR)/python2.7

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)
	-(/bin/rm -f netcdf*.log )
	-(cd netcdf-c-4.7.2 && make uninstall)
	-(cd netcdf-f-4.5.2 && make uninstall)
	-(cd fftw-3.3 && make uninstall)

