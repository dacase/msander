#  just redirect things to lower-level Makefiles

install::
	cd src && $(MAKE) install

full::
	cd src && $(MAKE) full

test::
	cd test && $(MAKE) test

clean::
	cd src && $(MAKE) clean

uninstall::
	cd src && $(MAKE) uninstall

distclean::
	cd src && $(MAKE) distclean
