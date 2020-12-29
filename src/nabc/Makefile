include ../../config.h

CUSTOMBUILDFLAGS = -I../sff

.SUFFIXES:
.SUFFIXES:  .c .o 

.c.o:
	@echo "[NABC]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(NETCDFINC) -o $@ $<

LIBNABOBJS=	\
	axis2frame.o	\
	bonded_atoms.o \
	database.o  \
	deform_ener.o \
    fd_helix.o \
	frameutil.o	\
	getdist.o	\
	jacobi.o	\
	match.o		\
	matop.o		\
	molio.o		\
	molsurf.o   \
	molutil.o	\
	pointutil.o	\
	prm.o \
	regexp.o \
	reslib.o	\
	select_atoms.o	\
	traceback.o \
    wc_basepair.o \
    wc_complement.o 

#===========================================================================

libnabc: $(LIBDIR)/libnabc.a

$(LIBDIR)/libnabc.a: $(LIBNABOBJS)
	@echo "[NABC]  AR $@"
	$(VB)$(AR) $@ $(LIBNABOBJS) > /dev/null
	$(VB)$(RANLIB) $@
	$(VB)cp nabc.h $(INCDIR)

lex.dg_options.c: dg_options.l
	@echo "[NABC]  $(LEX) $<"
	$(VB)$(LEX) -t dg_options.l | sed -e 's/{stdout}/{DONTchangeSTDOUThere}/' -e 's/yy/dgo/g' -e 's/stdout/nabout/' -e 's/DONTchangeSTDOUThere/stdout/' > lex.dg_options.c

#===========================================================================

clean:
	-rm -f *.o

uninstall: clean
	-rm -f  $(LIBDIR)/libnabc.a  $(INCDIR)/nabc.h
	-rm -f *.o

#===========================================================================

axis2frame.o:	\
	nabc.h
binpos.o:  \
	memutil.h
bonded_atoms.o: \
	nabc.h
boundsutil.o:	\
	nabc.h  \
	memutil.h \
    traceback.h \
	errormsg.h
database.o:  \
	database.h
dump.o:		\
	nabc.h
errormsg.o:	\
	errormsg.h
frameutil.o:	\
	nabc.h	\
	memutil.h
jacobi.o:	\
	memutil.h
matop.o: \
	matop.h
match.o:	\
	nabc.h \
	select_atoms.h 
molio.o:	\
	nabc.h	\
	errormsg.h
molsurf.o:  \
	nabc.h   \
	molsurf.h
molutil.o:	\
	nabc.h	\
	errormsg.h	\
	memutil.h
pointutil.o: \
	nabc.h
prm.o: \
	nabc.h
reslib.o: \
	nabc.h \
	errormsg.h \
	memutil.h  \
	database.h
select_atoms.o:	\
	nabc.h
mme_init.o: \
	nabc.h 
traceback.o: \
	errormsg.h

