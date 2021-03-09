#!/bin/sh

#   standard amber refinement

if [ "$#" -lt 3 ]; then
   echo "Usage:  phenix.amber.sh <pdbfile> <mtzfile> <id> <cif-files>"
   exit 1
fi

cat <<EOF > ${3}_002.eff
refinement {
  input {
    xray_data {
      outliers_rejection = True
      r_free_flags {
        generate = False
      }
    }
  }
  output {
    prefix = "$3"
    serial = 2
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
    sites {
    }
    adp {
      individual {
         anisotropic = not element H
      }
    }
  }
  target_weights {
    optimize_xyz_weight = False
    fix_wxc = 1
    wc = 3.333
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
    ignore_hydrogens = True
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
    topology_file_name = "../XrayPrep_27uc/alt5.5_uc.parm7"
    coordinate_file_name = "../XrayPrep_27uc/alt5.5_uc.rst7"
    order_file_name = "../XrayPrep_27uc/alt5.5_uc.order"
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

phenix.refine  $1  $2  $4  ${3}_002.eff --overwrite > $3.amber.log

diff ${3}_002.log $3.amber.log | awk 'NF==8 && $2!="Amber" {print $2}' \
     > ${3}_002.dat

/bin/mv $3.amber.log ${3}_002.log


