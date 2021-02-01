test:
	@./run_tests.sh

test.all: test.general test.rism test.xray test.tleap test.antechamber

test.general:
	-(cd dhfr && ./Run.dhfr)
	-(cd dhfr && ./Run.dhfr.xmin)
	-(cd gb8_trx && ./Run.trxox)
	-(cd gb8_trx && ./Run.trxox_md)
	-(cd tip4p && ./Run.tip4p)
	-(cd tip4p && ./Run.tip4p_mcbar)
	-(cd middle-scheme/MD_Constr_Langevin_water && ./Run.LGV_LF)
	-(cd middle-scheme/MD_Unconstr_Langevin_water && ./Run.LGV_LF)
	-(cd middle-scheme/QMMM_Constr_ALA_Methol && Run.qmmm)
	-(cd softcore && ./Run_sc)
	-(cd qmmm_DFTB/lysine_DFTB && ./Run.lysine_md)
	-(cd qmmm2/crambin_PM6 && ./Run.crambin_md)

test.rism:
	-(cd rism3d/4lzt && ./Run.4lzt.pme)
	-(cd rism3d/2igd && ./Run.2igd.pme)
	-(cd xray/480d && ./Run.480d_pse2)

test.xray:
	-(cd xray/2wcu && ./Run.2wcu)
	-(cd xray/3onk && ./Run.3onk_ml)
	-(cd xray/3onk && ./Run.3onk_mlopt)
	-(cd xray/480d && ./Run.480d_xray)
#	-(cd xray/6o2h && ./Run.6o2h_min)
	-(cd cryoem/3j9d && ./Run.3j9d)
	-(cd cryoem/emap && ./Run.emap)

test.tleap:
	-(cd tleap/OL3 && ./Run.OL3.sh)
	-(cd tleap/BondAugmentations && ./Run.1p7e)

test.antechamber:
	-(cd antechamber/sustiva && Run.sustiva)
	-(cd antechamber/sustiva && Run.sustiva_g2)

test.cif:
	-(cd nmr && ./Run.nef)

clean::
	/bin/rm -f TEST_FAILURES.diff
