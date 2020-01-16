include ../config.h

SRC = add_xray.F90 add_pdb.F90
OBJ = add_xray.o add_pdb.o

xray: $(BINDIR)/add_xray$(SFX) $(BINDIR)/add_pdb$(SFX) $(BINDIR)/XrayPrep

$(BINDIR)/add_xray$(SFX): add_xray.o 
	@echo "[XRAY]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LDFLAGS) $< $(LIBDIR)/libsys.a

$(BINDIR)/add_pdb$(SFX): add_pdb.o $(LIBS)
	@echo "[XRAY]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LDFLAGS) $< $(LIBDIR)/libsys.a

$(BINDIR)/XrayPrep: XrayPrep
	@echo "[XRAY]  install $@"
	$(VB)cp XrayPrep $(BINDIR)
	
clean:
	/bin/rm -f $(OBJ)

uninstall: clean
	-rm -f $(BINDIR)/add_xray$(SFX)
	-rm -f $(BINDIR)/add_pdb$(SFX)
	-rm -f $(BINDIR)/XrayPrep

add_xray.o:   add_xray.F90
	@echo "[XRAY]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

add_pdb.o:   add_pdb.F90
	@echo "[XRAY]  FC $<"
	$(VB)$(FC) -c $(FNOOPTFLAGS) $(FFLAGS) -o $@ $<

