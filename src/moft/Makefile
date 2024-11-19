include ../../config.h

HEADERS= \
   metaFFT/ConvCorr.hpp          \
   metaFFT/fftwmain.hpp          \
   metaFFT/fftwwrap.hpp          \
   metatwist/Density.hpp         \
   metatwist/boosthead.hpp       \
   metatwist/lutilities.hpp      \
   metatwist/readccp4.hpp        \
   metatwist/readdx.hpp          \
   metatwist/scattering.hpp      \
   utilities/boosthead.hpp       \
   utilities/interpolation.hpp   \
   utilities/lapack.hpp          \
   utilities/rotations.hpp       \
   utilities/statistics.hpp      \
   utilities/transrot.hpp        \
   utilities/utilities.hpp      

BOOSTLIB=-lboost_program_options -lboost_iostreams -lboost_regex -lboost_system -lboost_timer -lboost_chrono -lboost_filesystem -lboost_graph

install: $(BINDIR)/metatwist

metatwist.o:  metatwist.cpp $(HEADERS)
	@echo "[MoFT]  CXX $<"
	$(VB)$(CXX) $(CXXFLAGS) -c -std=c++0x -I. -I$(INCDIR) -o $@ $<

$(BINDIR)/metatwist:  metatwist.o
	@echo "[MoFT]  CXX $@"
	$(VB)$(CXX) -std=c++0x -o $@ $< -L$(LIBDIR) $(BOOSTLIB) -lfftw3

clean:
	/bin/rm -f *.o
