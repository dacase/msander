include ../../../config.h

makeDIST_RST$(SFX):	makeDIST_RST.o
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o makeDIST_RST$(SFX) makeDIST_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

nef_to_RST$(SFX):	nef_to_RST.o $(LIBDIR)/libcifparse.a
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o nef_to_RST$(SFX) nef_to_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM) -L$(LIBDIR) -lcifparse

makeANG_RST$(SFX):	makeANG_RST.o
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o makeANG_RST$(SFX) makeANG_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM)

install: 	makeDIST_RST$(SFX) makeANG_RST$(SFX) nef_to_RST$(SFX)
	$(VB)mv makeDIST_RST$(SFX) makeANG_RST$(SFX) nef_to_RST$(SFX) $(BINDIR)
	$(VB)cp makeCHIR_RST $(BINDIR)
	$(VB)cp map.NEF-AMBER map.DG-AMBER $(DATDIR)

uninstall:
	/bin/rm -f $(BINDIR)/makeDIST_RST$(SFX) \
	           $(BINDIR)/makeANG_RST$(SFX) $(BINDIR)/nef_to_RST$(SFX)
	/bin/rm -f $(BINDIR)/makeCHIR_RST $(BINDIR)/makeDIP_RST.cyana
	/bin/rm -f $(DATDIR)/map.NEF-AMBER $(DATDIR)/map.DG-AMBER

clean:
	/bin/rm -f *.o 

makeDIST_RST.o:  global.h
nef_to_RST.o:  global.h

#  default rule for C compilation:

.c.o:
	@echo "[NMR_AUX]  CC $<"
	$(VB)$(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) -o $@ $<

$(LIBDIR)/libcifparse.a:
	cd ../../cifparse && $(MAKE) install
