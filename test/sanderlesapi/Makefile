include $(AMBERHOME)/config.h

test:
	(if [ -f "$(LIBDIR)/libsanderles$(SHARED_SUFFIX)" ]; then \
	   $(MAKE) test_all; \
	 else \
	   $(MAKE) test_none; \
	 fi;)

test_all: Fortran C CPP Python

# This rule should only be followed if the sander API was not installed.
test_none:
	@echo ""
	@echo "The SANDER LES API was not installed. Skipping SANDER LES API tests."
	@echo ""

Fortran: testf
	@echo "Running testf"
	@(./testf && /bin/rm -f testf) || echo "Program error"

testf::
	@echo "Testing the Fortran LES API"
	$(FC) $(FFLAGS) -I$(AMBERHOME)/include -L$(AMBERHOME)/lib $(LDFLAGS) test.F90 -o testf -lsanderles

C: testc
	@echo "Running testc"
	@(./testc && /bin/rm -f testc) || echo "Program error"

testc::
	@echo "Testing the C LES API"
	$(CC) $(CFLAGS) -I$(AMBERHOME)/include -L$(AMBERHOME)/lib $(LDFLAGS) test.c -o testc -lsanderles

CPP: testcpp
	@echo "Running testcpp"
	@(./testcpp && /bin/rm -f testcpp) || echo "Program error"

testcpp::
	@echo "Testing the C++ LES API"
	$(CXX) $(CFLAGS) -I$(AMBERHOME)/include -L$(AMBERHOME)/lib $(LDFLAGS) test.cpp -o testcpp -lsanderles

Python:
	@echo "Testing the Python LES API"
	@echo "$(PYTHON) test.py"
	@$(PYTHON) test.py || echo "Program error"

clean:
	/bin/rm -f testf testc testcpp

.SUFFIXES:
