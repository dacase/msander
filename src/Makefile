#   top-level Makefile for msander

include ../config.h

install: $(INSTALLTYPE)

full:
	@cd msander && $(MAKE) msander
	@cd msander && $(MAKE) msander.LES
	@cd xray && $(MAKE) xray
	@cd xtalutil && $(MAKE) xtalutil
	@cd addles && $(MAKE) addles
	@cd parmed && $(MAKE) $(PARMED)

serial:
	@cd msander && $(MAKE) msander
	@cd xray && $(MAKE) xray
	@cd xtalutil && $(MAKE) xtalutil

parallel:
	@cd msander && $(MAKE) msander.MPI
#	@cd msander && $(MAKE) msander.LES.MPI

cuda:
	@cd msander && $(MAKE) msander.cuda

mpicuda:
	@cd msander && $(MAKE) msander.cuda.MPI

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )
	-(cd xray/cuda && $(MAKE) clean )
	-(cd addles && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )

uninstall: clean
	-(cd $(BINDIR) && /bin/rm -f msander msander.MPI \
        add_pdb add_xray addles)
	-(cd $(LIBDIR) && /bin/rm -f liblapack.a libblas.a \
        libarpack.a libsqm.a librism.a /libsys.a )
	/bin/rm -f $(INCDIR)/* ../config.h
	/bin/rm -rf $(BINDIR)/x $(LIBDIR)/python3.*

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) logs config.h)
	-(/bin/rm -f netcdf*.log )
	-(cd netcdf-c-4.7.2 && make uninstall)
	-(cd netcdf-fortran-4.5.2 && make uninstall)
	-(cd fftw-3.3 && make uninstall)

