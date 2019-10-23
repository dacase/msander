#   top-level Makefile for ambermini

include config.h

install:
	@echo "Starting installation of ${AMBERTOOLS} serial at `date`".
	@(if [ -d arpack ] ; then cd arpack && $(MAKE) install; fi );
	@(if [ -d lapack ] ; then cd lapack && $(MAKE) $(LB); fi )
	@(if [ -d blas ] ; then cd blas && $(MAKE) $(LB); fi )
	@(if [ -s rism/Makefile ] ; then cd rism && $(MAKE) install; fi )
	@(if [ -d sqm ] ; then cd sqm && $(MAKE) install; fi )
	@(if [ -s sander/Makefile ] ; then cd sander && $(MAKE) serial_notest; fi )
	@(if [ -d addles ] ; then cd addles && $(MAKE) install; fi )
#	@(if [ -d etc ] ; then cd etc && $(MAKE) ambpdb$(SFX); fi )

	@echo ""
	@echo "Installation of ${AMBERTOOLS} serial is complete at `date`."
	@echo ""

parallel:
	@echo "Starting installation of ${AMBERTOOLS} parallel at `date`".
	@(if [ -s rism/Makefile ] ; then cd rism && $(MAKE) parallel; fi )
	@(if [ -d sqm ] ; then cd sqm && $(MAKE) parallel; fi )
	@(if [ -d sander ] ; then cd sander && $(MAKE) parallel; fi )
	@echo ""
	@echo "Installation of ${AMBERTOOLS} parallel is complete at `date`."
	@echo ""

clean:
	-(cd arpack && $(MAKE) clean )
	-(cd blas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd sander && $(MAKE) clean )
	-(cd addles && $(MAKE) clean )
#	-(cd etc && $(MAKE) clean )
#
	@echo "skipping cleaning for netcdf, fftw3 and xblas"
#	-(cd fftw-3.3.8 && $(MAKE) clean )
#	-(cd netcdf-4.6.1 && $(MAKE) clean )
#	-(cd netcdf-fortran-4.4.4 && $(MAKE) clean )
#	-(cd xblas && $(MAKE) clean )

uninstall: clean
	-(cd arpack && $(MAKE) uninstall )
	-(cd blas && $(MAKE) uninstall )
	-(cd lapack && $(MAKE) uninstall )
	-(cd sqm && $(MAKE) uninstall )
	-(cd sander && $(MAKE) uninstall )
	-(cd addles && $(MAKE) uninstall )
	-(cd rism && $(MAKE) uninstall )
	-(cd etc && $(MAKE) uninstall )
#	-(cd fftw-3.3.8 && $(MAKE) distclean )
#	-(cd netcdf-4.6.1 && $(MAKE) distclean )
#	-(cd netcdf-fortran-4.4.4 && $(MAKE) distclean )
#	-(cd xblas && $(MAKE) distclean )
	-(/bin/rm -rf $(BINDIR)/* $(LIBDIR) $(INCDIR))
#	-(/bin/rm -f fftw3_config.log netcdf_config.log netcdf-fortran_config.log)

