#   top-level Makefile for msander

include ../config.h

install:
	cd msander && $(MAKE) msander
	cd xray && $(MAKE) xray

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd xray && $(MAKE) clean )

uninstall: clean
	-(cd .. && /bin/rm -f $(BINDIR)/msander $(BINDIR)/msander.MPI \
	             $(BINDIR)/add_pdb $(BINDIR)/add_xray \
                 $(LIBDIR)/liblapack \
                 $(LIBDIR)/libblas \
                 $(LIBDIR)/libarpack \
                 $(LIBDIR)/libsqm \
                 $(LIBDIR)/librism \
                 $(LIBDIR)/libsys \
                 $(INCDIR)/* ../config.h)

distclean: clean
	-(cd .. && /bin/rm -rf $(BINDIR) $(LIBDIR) $(INCDIR) )
	-(/bin/rm -f netcdf*.log )
	-(cd netcdf-c-4.7.2 && make uninstall)
	-(cd netcdf-f-4.5.2 && make uninstall)
	-(cd fftw-3.3 && make uninstall)

