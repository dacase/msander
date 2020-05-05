include ../config.h

LOCALFLAGS=-Wno-unused-result

.c.o:
	@echo "[XTALUTIL]  CC $<"
	@$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $< 

OBJS = crdmanip.o grid.o listRead.o matrix.o pdbRead.o QuatAlign.o ran2.o \
       topRead.o vector.o

install: $(BINDIR)/ChBox$(SFX) $(BINDIR)/UnitCell$(SFX) \
     $(BINDIR)/AddToBox$(SFX) $(BINDIR)/PropPDB$(SFX)

$(BINDIR)/AddToBox$(SFX): AddToBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/ChBox$(SFX): ChBox.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/CheckMD$(SFX): CheckMD.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/PropPDB$(SFX): PropPDB.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/UnitCell$(SFX): UnitCell.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

$(BINDIR)/pdb2rdf$(SFX): pdb2rdf.o $(OBJS)
	@echo "[XTALUTIL]  CC $@"
	@$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ -lm

uninstall:
	/bin/rm -f $(BINDIR)/AddToBox$(SFX) $(BINDIR)/ChBox$(SFX) \
       $(BINDIR)/PropPDB$(SFX) $(BINDIR)/UnitCell$(SFX)

clean:
	/bin/rm -f *.o

