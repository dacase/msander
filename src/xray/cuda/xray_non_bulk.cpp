#include "xray_non_bulk.h"
#include "NonBulkCPU.h"
#include "NonBulkGPU.h"
#include <cassert>
#include <memory>
#include <chrono>
#include <iostream>
#include <iomanip>

namespace {
  std::unique_ptr<xray::NonBulk> non_bulk;
}

extern "C" void pmemd_xray_non_bulk_init_gpu(
  int n_hkl,
  const int* hkl,
  complex_double* f_non_bulk,
  const double* mSS4,
  int n_atom,
  const double* b_factor,
  int n_scatter_types,
  const int* scatter_type_index,
  const double* atomic_scatter_factor) {
  non_bulk = std::unique_ptr<xray::NonBulk>(
    new xray::NonBulkGPU<xray::NonBulkKernelVersion::StraightForward, xray::KernelPrecision::Single>(
      n_hkl, hkl, f_non_bulk, mSS4, n_atom, b_factor, n_scatter_types, scatter_type_index, atomic_scatter_factor
   ));
}

extern "C" void pmemd_xray_non_bulk_calc_f_non_bulk_gpu(
  int n_atoms,
  const double* frac_xyz
) {
  assert(non_bulk);

  auto t1 = std::chrono::high_resolution_clock::now();
  non_bulk->calc_f_non_bulk(n_atoms, frac_xyz);

  auto t2 = std::chrono::high_resolution_clock::now();
  #ifndef NDEBUG
  {
    long dt_us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cerr << std::setw(32) << "calc_f_non_bulk: " << std::setw(5) << dt_us << " us" << std::endl;
  }
  #endif
}


extern "C" void pmemd_xray_non_bulk_finalize_gpu() {
  assert(non_bulk);
  non_bulk = {};
}
