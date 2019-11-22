#   top-level Makefile for msander

include config.h

install:
	cd msander && $(MAKE) msander
	@(if [ -z "$(SFX)" ] ; then cd xray && $(MAKE) xray; fi );

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
	-(/bin/rm -f $(BINDIR)/msander $(BINDIR)/msander.MPI \
	             $(BINDIR)/add_pdb $(BINDIR)/add_xray \
                 $(LIBDIR)/liblapack \
                 $(LIBDIR)/libblas \
                 $(LIBDIR)/libarpack \
                 $(LIBDIR)/libsqm \
                 $(LIBDIR)/librism \
                 $(LIBDIR)/libsys \
                 $(INCDIR)/* config.h)

