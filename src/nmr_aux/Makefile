include ../../config.h

nmr_aux::	
	cd prepare_input && $(MAKE) clean && $(MAKE) install
	cp senergy sviol sviol2 $(BINDIR)

clean::	
	cd prepare_input && $(MAKE) clean 

uninstall:
	-(cd prepare_input && $(MAKE) uninstall)
	-rm -f $(BINDIR)/senergy
	-rm -f $(BINDIR)/sviol
	-rm -f $(BINDIR)/sviol2
