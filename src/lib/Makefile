SHELL=/bin/sh
VPATH=
PREC=DPREC
include ../config.h
.SUFFIXES: .F .F90 .c .o

OBJ=  \
sys.o wallclock.o veclib.o amopen.o mexit.o nxtsec.o rgroup.o \
rfree.o matinv.o

#  default rules for Fortran and C compilation:

.F.o:   $<
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

.F90.o:   $<
	@echo "[LIB]  FC $<"
	$(VB)$(FC) $(FPPFLAGS) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

.c.o:
	@echo "[LIB]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

libsys:  $(LIBDIR)/libsys.a

$(LIBDIR)/libsys.a:  $(OBJ)
	@echo "[LIB]  AR $@"
	$(VB)$(AR) $@ $(OBJ)  > /dev/null
	$(VB)$(RANLIB) $@

clean:
	/bin/rm -f *.o  *.a  *.mod

uninstall:
	/bin/rm -f $(LIBDIR)/libsys.a

