#include "BulkMaskGPU.h"
#include "xray_non_bulk.h"
#include <cassert>
#include <memory>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace {

std::unique_ptr<xray::BulkMask>& bulk_mask_instance(){
  static std::unique_ptr<xray::BulkMask> mask;
  return mask;
}

} // namespace

extern "C" void
pmemd_xray_bulk_mask_init_gpu(
  int* grid_dim,
  double* reciprocal_norms,
  int n_atom, double* mask_cutoffs,
  int n_grid_points, int* mask_bs_grid,
  double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
  double shrink_r, int n_hkl, complex_double* f_mask, int* hkl
) {
  assert(n_grid_points == grid_dim[0] * grid_dim[1] * grid_dim[2]);
  auto cell = xray::UnitCell(a, b, c, alpha_deg, beta_deg, gamma_deg);
  bulk_mask_instance().reset(new xray::BulkMaskGPU{});
  bulk_mask_instance()->init(cell, n_atom, grid_dim, reciprocal_norms, mask_cutoffs,
             mask_bs_grid, shrink_r, n_hkl, f_mask, hkl);
}

extern "C" void pmemd_xray_bulk_mask_update_f_bulk(int n_atom, const double* frac) {
  assert(bulk_mask_instance());

  auto t1 = std::chrono::high_resolution_clock::now();
  bulk_mask_instance()->update_grid(n_atom, frac);
  auto t2 = std::chrono::high_resolution_clock::now();
  bulk_mask_instance()->calc_f_bulk();
  auto t3 = std::chrono::high_resolution_clock::now();
  #ifndef NDEBUG
  {
    long dt_us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cerr << std::setw(32) << "update_grid: " << std::setw(5) << dt_us << " us" << std::endl;
  }
  {
    long dt_us = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
    std::cerr << std::setw(32) << "calc_f_bulk: " << std::setw(5) << dt_us << " us" << std::endl;
  }
  #endif
}

extern "C" void pmemd_xray_bulk_mask_finalize_gpu() {
  bulk_mask_instance().reset();
}

