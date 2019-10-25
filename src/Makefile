#   top-level Makefile for sandermini

include config.h

install:
	cd sander && $(MAKE) install

clean:
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd sander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )

uninstall: clean
	-(/bin/rm -rf $(BINDIR)/* $(LIBDIR)/* $(INCDIR)/* config.h)

