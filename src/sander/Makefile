
include ../config.h

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI

SHELL=/bin/sh
VPATH=

# Sander Fortran source files are free format
# Add also search paths for .mod files
LOCALFLAGS = -I../sqm -I../rism 

NFE_OBJECTS = \
    nfe-abmd-hooks.o nfe-bbmd-ctxt.o nfe-bbmd-hooks.o  \
    nfe-colvar.o nfe-colvar-math.o nfe-colvar-type.o nfe-colvar-utils.o \
    nfe-constants.o nfe-cv-ANGLE.o nfe-cv-COM_ANGLE.o \
    nfe-cv-COM_DISTANCE.o nfe-cv-COM_TORSION.o nfe-cv-COS_OF_DIHEDRAL.o \
    nfe-cv-DISTANCE.o nfe-cv-HANDEDNESS.o nfe-cv-LCOD.o \
    nfe-cv-MULTI_RMSD.o nfe-cv-N_OF_BONDS.o nfe-cv-N_OF_STRUCTURES.o \
    nfe-cv-R_OF_GYRATION.o nfe-cv-TORSION.o  \
    nfe-pmd-hooks.o nfe-rmsd.o nfe-sander-hooks.o nfe-sander-proxy.o \
    nfe-stsm-hooks.o nfe-cv-PCA.o nfe-read-pca.o nfe-cv-SIN_OF_DIHEDRAL.o \
    nfe-cv-PAIR_DIHEDRAL.o nfe-cv-PATTERN_DIHEDRAL.o nfe-cv-DF_COM_DISTANCE.o \
    nfe-smd-hooks.o nfe-umbrella.o nfe-utils.o nfe-value.o mt19937.o

AMOEBAOBJ= \
    amoeba_valence.o amoeba_multipoles.o amoeba_recip.o amoeba_interface.o\
    amoeba_direct.o amoeba_mdin.o amoeba_adjust.o amoeba_self.o\
    amoeba_vdw.o amoeba_induced.o amoeba_runmd.o

QMOBJ= \
    abfqmmm_module.o qm_mm.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o \
    qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o

QMAPIOBJ = \
    abfqmmm_module.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o qm2_extern_nw_module.o \
    qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o

MMOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    multisander.o sander.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o rdparm.o sander_lib.o \
    mdread.o locmem.o runmd.o relaxmd.o quench.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
    set.o dynlib.o mdfil.o nmlsrc.o \
    force.o ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o pol_gauss.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o erfcfun.o sinr_t.o md_scheme.o \
    constantph.o constante.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o\
    $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o amd.o scaledMD.o barostats.o random.o les.o music_module.o

FULLPIMDOBJ = pimd_vars.o pimd_force.o pimd_init.o cmd_vars.o cmd_matrix.o

PARTPIMDOBJ = pimd_vars.o pimd_force.LES.o pimd_init.LES.o cmd_vars.o cmd_matrix.LES.o

LSCIVROBJ = lscivr_vars.o lsc_init.o lsc_xp.o

LESOBJ= memory_module.o stack.o file_io_dat.o \
    constants.o state.o nose_hoover.o nose_hoover_init_LES.o \
    trace.o rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o \
    relax_mat.o nmr.LES.o multitmd.o \
    multisander.LES.o sander.LES.o  cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o nonbond_list.LES.o force.LES.o rdparm.LES.o \
    mdread.LES.o locmem.LES.o runmd.LES.o relaxmd.o quench.o degcnt.LES.o decnvh.o \
    fastwt.o parallel.LES.o shake.o rattle.o ene.LES.o \
    dynlib.LES.o nmlsrc.o \
    ew_force.LES.o ew_setup.LES.o ew_bspline.o\
    ew_fft.o short_ene.LES.o pol_gauss.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    sander_lib.o random.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    new_time.o extra_pts.LES.o thermo_int.o matinv.o assert.o \
    erfcfun.o sinr_t.o md_scheme.LES.o \
    constantph.o constante.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(SANDER_RISM_INTERFACE) \
    les.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o barostats.o \
    qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o qm2_extern_nw_module.o \
    qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_genmpi_module.o \
    qmmm_adaptive_module.o abfqmmm_module.o prn_qmmm_dipole.o \
    AmberNetcdf.o bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o mdwrit.LES.o \
    lmod.LES.o trajene.o ew_box.LES.o minrit.LES.o set.LES.o qm_mm.LES.o amd.o scaledMD.o \
    music_module.o

XRAY_OBJS = \
   xray_fftpack.o xray_FFT.o xray_fourier.o xray_globals.o xray_interface.o \
   xray_real_space.o xray_reciprocal_space.o xray_utils.o

.PHONY: dummy

msander: $(BINDIR)/msander$(SFX)

msander.les: $(BINDIR)/msander.LES$(SFX)

dummy:

#---------------------------------------------------------------------------
$(BINDIR)/msander$(SFX): libs $(MMOBJ) $(QMOBJ) $(FULLPIMDOBJ) \
         $(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) $(PLUMED_DEPENDENCIES) \
         dummy
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@  $(MMOBJ) $(QMOBJ) \
	     $(FULLPIMDOBJ) $(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) \
	     -L$(LIBDIR) $(NETCDFLIBF) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/msander.LES$(SFX): libs $(LESOBJ) $(PARTPIMDOBJ) \
         $(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) \
         $(PLUMED_DEPENDENCIES) dummy
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(FFLAGS) -o $@ $(LESOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
	     $(XRAY_OBJS) $(NFE_OBJECTS) -L$(LIBDIR) $(NETCDFLIBF) \
	     $(FLIBS_RISMSANDER) $(FLIBSF) $(LDFLAGS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
libs:
	cd ../arpack && $(MAKE) libarpack
	cd ../rism && $(MAKE) librism
	cd ../sqm && $(MAKE) libsqm
	cd ../lib && $(MAKE) libsys
    
Plumed.o: Plumed.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
	        -o Plumed.o Plumed.c

rand2.o:  ../sff/rand2.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o rand2.o ../sff/rand2.c

lmodC.o:  ../sff/lmodC.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o lmodC.o ../sff/lmodC.c

xminC.o:  ../sff/xminC.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o xminC.o ../sff/xminC.c

clean:
	/bin/rm -f sander.PUPIL$(SFX)
	/bin/rm -f sander.APBS$(SFX)
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *.d
	/bin/rm -f *nbflag
	/bin/rm -f *__genmod.f90  # intel compiler-generated interface modules

depend::
	./makedepend > depend

#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
