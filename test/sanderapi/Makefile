include $(AMBERHOME)/config.h

test:
	(if [ -f "$(LIBDIR)/libsander$(SHARED_SUFFIX)" ]; then \
	   $(MAKE) test_all; \
	 else \
	   $(MAKE) test_none; \
	 fi;)

test_all: test.parm7 Fortran Fortran2 C CPP Python clean

test.parm7:
	@echo "Creating the bad prmtop test.parm7"
	@# Delete one of the critical flags in the prmtop with sed as if by:
	@#$(BINDIR)/parmed -Oen -i prepare.pmd > /dev/null
	@cat ../gb7_trx/prmtop_an                                 \
	    | sed '/%FLAG ATOM_NAME/,/%FLAG/p'                    \
            | sed '/%FLAG ATOM_NAME/N ; /^\(.*\)\n\1$//!P ; D'    \
            | sed '/%FLAG ATOM_NAME/,/%FLAG/d'                    \
            > ./test.parm7
	@#1. duplicate the %FLAG after %FLAG ATOM_NAME
	@#2. fix side effect from step 1:  unduplicate %FLAG ATOM_NAME
	@#3. delete between %FLAG ATOM_NAME and the %FLAG after %FLAG ATOM_NAME

Fortran: testf test.parm7
	@echo "Running testf"
	@(./testf && /bin/rm -f testf) || echo "Program error"

Fortran2: testf2 test.parm7
	@echo "Running testf2"
	@(./testf2 && /bin/rm -f testf2) || echo "Program error"

# This rule should only be followed if the sander API was not installed.
test_none:
	@echo ""
	@echo "The SANDER API was not installed. Skipping SANDER API tests."
	@echo ""

testf:: $(LIBDIR)/libsander$(SHARED_SUFFIX)
	@echo "Testing the Fortran API"
	$(FC) $(FFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LDFLAGS) test.F90 -o testf -lsander

testf2::
	@echo "Testing the Fortran file-less API"
	$(FC) $(FFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LDFLAGS) test2.F90 -o testf2 -lsander

C: testc test.parm7
	@echo "Running testc"
	@(./testc && /bin/rm -f testc) || echo "Program error"

testc::
	@echo "Testing the C API"
	$(CC) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LDFLAGS) test.c -o testc -lsander

CPP: testcpp test.parm7
	@echo "Running testcpp"
	@(./testcpp && /bin/rm -f testcpp) || echo "Program error"

testcpp::
	@echo "Testing the C++ API"
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LDFLAGS) test.cpp -o testcpp -lsander

Python: test.parm7
	@echo "Testing the Python API"
	@echo "$(PYTHON) test.py"
	@$(PYTHON) test.py || echo "Program error"

clean:
	/bin/rm -f test.parm7 testf testf2 testc testcpp

.SUFFIXES:

.NOTPARALLEL: Fortran Fortran2 C CPP Python clean
