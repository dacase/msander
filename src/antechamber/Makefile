# Makefile for antechamber

include ../config.h

ACDOCTOROBJS = acdoctor.o eprintf.o
AM1BCCOBJS = am1bcc.o eprintf.o
ANTECHAMBEROBJS = antechamber.o eprintf.o
ATOMTYPEOBJS = atomtype.o eprintf.o
BONDTYPEOBJS = bondtype.o eprintf.o ncsu-penalties.o
CHARMMGENOBJS = charmmgen.o eprintf.o
CRDGROWOBJS = crdgrow.o eprintf.o
DATABASEOBJS = database.o eprintf.o
DELPHIGENOBJS = delphigen.o
ESPGENOBJS = espgen.o eprintf.o
MATCHOBJS = eprintf.o match.o
MATCH_ATOMNAMEOBJS = match_atomname.o eprintf.o
PARMCALOBJS = eprintf.o parmcal.o
PARMCHK2OBJS = parmchk2.o eprintf.o
PREPGENOBJS = eprintf.o prepgen.o
RESIDUEGENOBJS = residuegen.o eprintf.o
RESPGENOBJS = respgen.o eprintf.o
TRANSLATEOBJS = translate.o eprintf.o
TOP2FFOBJS = top2ff.o eprintf.o
TOP2MOL2OBJS = top2mol2.o eprintf.o

PROGS = \
	$(BINDIR)/x/am1bcc$(SFX) \
	$(BINDIR)/x/antechamber$(SFX) \
	$(BINDIR)/x/atomtype$(SFX) \
	$(BINDIR)/x/bondtype$(SFX) \
	$(BINDIR)/x/espgen$(SFX) \
	$(BINDIR)/x/parmchk2$(SFX) \
	$(BINDIR)/x/prepgen$(SFX) \
	$(BINDIR)/x/residuegen$(SFX) \
	$(BINDIR)/x/respgen$(SFX) \

TESTS = test-penalties test-int-part test-combination

#############################################################################

antechamber: $(PROGS) dispatch

dispatch:
	SFX=$(SFX) BINDIR=$(BINDIR) /bin/sh dispatch.sh

clean:
	/bin/rm -f *.o $(TESTS)

uninstall:
	/bin/rm -f $(PROGS)

#############################################################################

$(LIBDIR)/libcifparse.a:
	cd ../cifparse && $(MAKE) install

.c.o:
	@echo "[ANTECHAMBER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) -o $@ $<

#_____________________________________________________________________________

$(BINDIR)/x/acdoctor$(SFX): $(ACDOCTOROBJS) $(LIBDIR)/libcifparse.a
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(ACDOCTOROBJS) \
	    -L$(LIBDIR) -lcifparse $(LM) -o $@

$(BINDIR)/x/am1bcc$(SFX): $(AM1BCCOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(AM1BCCOBJS)  $(LM) -o $@

$(BINDIR)/x/antechamber$(SFX): $(ANTECHAMBEROBJS) $(LIBDIR)/libcifparse.a
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(ANTECHAMBEROBJS) \
	    -L$(LIBDIR) -lcifparse $(LM) -o $@

$(BINDIR)/x/atomtype$(SFX): $(ATOMTYPEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(ATOMTYPEOBJS) $(LM) -o $@

$(BINDIR)/x/bondtype$(SFX): $(BONDTYPEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(BONDTYPEOBJS) $(LM) -o $@

$(BINDIR)/x/charmmgen$(SFX): $(CHARMMGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(CHARMMGENOBJS) $(LM) -o $@

$(BINDIR)/x/crdgrow$(SFX): $(CRDGROWOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(CRDGROWOBJS) $(LM) -o $@

$(BINDIR)/x/database$(SFX): $(DATABASEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(DATABASEOBJS) $(LM) -o $@

$(BINDIR)/x/delphigen$(SFX): $(DELPHIGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(DELPHIGENOBJS) -o $@

$(BINDIR)/x/espgen$(SFX): $(ESPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(ESPGENOBJS) $(LM) -o $@

$(BINDIR)/x/match$(SFX): $(MATCHOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(MATCHOBJS) $(LM) -o $@

$(BINDIR)/x/match_atomname$(SFX): $(MATCH_ATOMNAMEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(MATCH_ATOMNAMEOBJS) $(LM) -o $@

$(BINDIR)/x/parmcal$(SFX): $(PARMCALOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(PARMCALOBJS) $(LM) -o $@

$(BINDIR)/x/parmchk2$(SFX): $(PARMCHK2OBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(PARMCHK2OBJS) $(LM) -o $@

$(BINDIR)/x/parmjoin$(SFX): parmjoin.o
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) parmjoin.o -o $@

$(BINDIR)/x/prepgen$(SFX): $(PREPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(PREPGENOBJS) $(LM) -o $@

$(BINDIR)/x/residuegen$(SFX): $(RESIDUEGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(RESIDUEGENOBJS) $(LM) -o $@

$(BINDIR)/x/respgen$(SFX): $(RESPGENOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(RESPGENOBJS) $(LM) -o $@

$(BINDIR)/x/translate$(SFX): $(TRANSLATEOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(TRANSLATEOBJS) $(LM) -o $@

$(BINDIR)/x/top2ff$(SFX): $(TOP2FFOBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(TOP2FFOBJS) $(LM) -o $@

$(BINDIR)/x/top2mol2$(SFX): $(TOP2MOL2OBJS)
	@echo "[ANTECHAMBER]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(LDFLAGS) $(TOP2MOL2OBJS) $(LM) -o $@

#_____________________________________________________________________________
# Tests

test-penalties: ncsu-penalties.o test-penalties.o
	$(CC) $(CFLAGS) -o test-penalties $(LDFLAGS) ncsu-penalties.o test-penalties.o

test-int-part: ncsu-penalties.h ncsu-penalties.c
	$(CC)  $(CNOOPTFLAGS) $(CFLAGS) -o test-int-part ncsu-penalties.c $(LDFLAGS) \
   -DENABLE_INT_PART_MAIN=1

test-combination: ncsu-penalties.h ncsu-penalties.c
	$(CC) $(CNOOPTFLAGS) $(CFLAGS) -o test-combination ncsu-penalties.c $(LDFLAGS) \
   -DENABLE_COMBINATION_MAIN=1

#_____________________________________________________________________________
# Dependencies

CFILES= ac.c acdoctor.c alc.c am1bcc.c antechamber.c aromatic.c \
atomtype.c bondtype.c charge.c charmm.c charmmgen.c checkmolecule.c \
common.c crdgrow.c csd.c database.c delphigen.c divcrt.c divout.c \
eprintf.c equatom.c equatom2.c espgen.c fileformat.c gcrt.c gesp.c \
gout.c gzmat.c hin.c jcrt.c jout.c jzmat.c lsfit.c match.c match_atomname.c \
mdl.c mem.c mmcif.c mol2.c mopcrt.c mopint.c mopout.c ncsu-penalties.c \
parmcal.c parmchk2.c parmjoin.c pdb.c poltop.c prep.c prepgen.c \
residuegen.c respgen.c ring.c rotate.c rst.c smiles.c sqmcrt.c \
sqmout.c test-penalties.c top2ff.c top2mol2.c translate.c utility.c \

depend: $(CFILES)
	@echo "[ANTECHAMBER]  make depend"
	$(VB)$(CC) -MM $^ > $@

include depend
