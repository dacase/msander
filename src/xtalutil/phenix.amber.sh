#!/bin/sh

#   standard amber refinement

if [ "$#" -ne 3 ]; then
   echo "Usage:  phenix.cdl.sh <pdbfile> <mtzfile> <serial-no>"
   exit 1
fi

cat <<EOF > amber.eff
refinement {
  input {
    xray_data {
      outliers_rejection = True
      r_free_flags {
        generate = True
      }
    }
  }
  output {
    prefix = "amber"
    serial = $3
    write_eff_file = False
    write_geo_file = False
    write_def_file = False
    write_model_cif_file = False
    write_map_coefficients = False
    export_final_f_model = False
  }
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls *occupancies group_anomalous
  }
  target_weights {
    optimize_xyz_weight = False
    fix_wxc = 1
    wc = 1.6667
  } 
  main {
    nqh_flips = True
    number_of_macro_cycles = 10
    target = auto *ml mlhl ml_sad ls mli
    use_experimental_phases = False
    scattering_table = wk1995 *it1992 n_gaussian electron neutron
  }
  hydrogens {
    refine = individual *riding Auto
  }
  pdb_interpretation {
    c_beta_restraints = False
  }
  mask {
    ignore_hydrogens = False
  }
  structure_factors_and_gradients_accuracy {
    algorithm = *fft direct
  }
  gui {
    skip_rsr = True
    skip_kinemage = True
  }
  amber {
    use_amber = True
    topology_file_name = "4amber_$1.prmtop"
    coordinate_file_name = "4amber_$1.rst7"
    order_file_name = "4amber_$1.order"
    wxc_factor = 0.2
    restraint_wt = 0
    restraintmask = ""
    reference_file_name = ""
    bellymask = ""
    netcdf_trajectory_file_name = ""
    print_amber_energies = True
  }
}
EOF

phenix.refine  4phenix_$1.pdb  $2  amber.eff --overwrite | tee amber$3.log
/bin/mv amber$3.log amber_00$3.log

/bin/rm -f amber.eff
