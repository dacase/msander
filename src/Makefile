#   top-level Makefile for msander

include config.h

install:
	cd msander && $(MAKE) msander

clean:
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd msander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )

uninstall: clean
	-(/bin/rm -f $(BINDIR)/msander $(BINDIR)/msander.MPI \
                 $(LIBDIR)/liblapack \
                 $(LIBDIR)/libblas \
                 $(LIBDIR)/libarpack \
                 $(LIBDIR)/libsqm \
                 $(LIBDIR)/librism \
                 $(LIBDIR)/libsys \
                 $(INCDIR)/* config.h)

