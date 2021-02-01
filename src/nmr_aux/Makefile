include ../../config.h

install:  $(BINDIR)/nef_to_RST$(SFX)

$(BINDIR)/nef_to_RST$(SFX):	nef_to_RST.o $(LIBDIR)/libcifparse.a
	@echo "[NMR_AUX]  CC $@"
	$(VB)$(CC) $(CFLAGS) $(AMBERCFLAGS) -o $@ nef_to_RST.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LM) -L$(LIBDIR) -lcifparse

clean:
	/bin/rm -f *.o 

nef_to_RST.o:  global.h

#  default rule for C compilation:

.c.o:
	@echo "[NMR_AUX]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) -o $@ $<

$(LIBDIR)/libcifparse.a:
	cd ../cifparse && $(MAKE) $@
