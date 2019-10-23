SHELL=/bin/sh
include ../config.h


# ------- rism3d information: ----------------------------------------------

LIBOBJ = quaternion.o array_util.o rism_util.o rism_report_c.o \
	rism3d_grid_c.o rism3d_closure_c.o \
	rism3d_kh_c.o rism3d_hnc_c.o rism3d_psen_c.o \
	rism3d_c.o rism3d_potential_c.o rism3d_restart.o rism3d_tree_cp.o \
	rism3d_xyzv.o rism3d_opendx.o rism3d_ccp4.o \
	rism3d_solvent_c.o rism3d_solute_c.o \
	pubfft.o rism3d_fft.o rism_parm.o rism_io.o \
	mdiis_orig_c.o mdiis_blas_c.o mdiis_blas2_c.o mdiis_c.o \
	fce_c.o erfcfun.o safemem.o blend.o rism_timer_c.o constants.o \
	getopts_c.o bspline.o  fftw3.o dsi.o dcpval.o \
	dqagie.o dqelg.o dqk15i.o dqpsrt.o d1mach.o i1mach.o \
	xermsg.o xerprn.o xersve.o xerhlt.o xercnt.o xgetua.o fdump.o j4save.o

librism: $(LIBDIR)/librism.a

$(LIBDIR)/librism.a: $(LIBOBJ)
	@echo "[RISM]  AR $@"
	$(VB)$(AR) $(LIBDIR)/librism.a $(LIBOBJ) > /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/librism.a

install: librism 

uninstall:
	/bin/rm -f $(LIBDIR)/librism.a

clean:
	-/bin/rm -f *.o  pmepot/*.o *.mod

depend::
	./makedepend > depend

include depend

# DO NOT DELETE
