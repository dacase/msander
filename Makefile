#  just redirect things to lower-level Makefiles

install::
	cd src && make install

install.ap::
	cd src && make -f Makefile.ap install

test::
	cd test && make test

clean::
	cd src && make clean

uninstall::
	cd src && make uninstall

distclean::
	cd src && make distclean
