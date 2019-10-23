#   top-level Makefile for sandermini

include config.h

install:
	@echo "Starting installation of ${AMBERTOOLS} serial at `date`".
	@(if [ -d arpack ] ; then cd arpack && $(MAKE) libarpack; fi );
	@(if [ -s rism/Makefile ] ; then cd rism && $(MAKE) librism; fi )
	@(if [ -d sqm ] ; then cd sqm && $(MAKE) libsqm; fi )
	@(if [ -s sander/Makefile ] ; then cd sander && $(MAKE) install; fi )

	@echo ""
	@echo "Installation of ${AMBERTOOLS} is complete at `date`."
	@echo ""

clean:
	-(cd arpack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd sander && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )

uninstall: clean
	-(cd arpack && $(MAKE) uninstall )
	-(cd sqm && $(MAKE) uninstall )
	-(cd rism && $(MAKE) uninstall )
	-(cd sander && $(MAKE) uninstall )
	-(cd lib && $(MAKE) uninstall )
	-(/bin/rm -rf $(BINDIR)/* $(LIBDIR) $(INCDIR) config.h)

