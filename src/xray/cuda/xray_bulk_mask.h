#ifndef AMBER_XRAY_BULK_MASK_H
#define AMBER_XRAY_BULK_MASK_H

#include "xray_common.h"

extern "C" {
void pmemd_xray_bulk_mask_init_gpu(
    int *grid_dim,
    double *reciprocal_norms,
    int n_atom, double *mask_cutoffs,
    int n_grid_points, int *mask_bs_grid,
    double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
    double shrink_r, int n_hkl, complex_double* f_mask, int* hkl);
void pmemd_xray_bulk_mask_update_f_bulk(int n_atom, const double *frac);
void pmemd_xray_bulk_mask_finalize_gpu();
}

#endif // AMBER_XRAY_BULK_MASK_H
