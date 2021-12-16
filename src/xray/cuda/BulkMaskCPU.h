#ifndef AMBER_XRAY_BULK_MASK_CPU_H
#define AMBER_XRAY_BULK_MASK_CPU_H

#include <array>
#include "BulkMask.h"
#include <vector>

namespace xray {

  class BulkMaskCPU : public BulkMask {
  public:
    BulkMaskCPU() noexcept = default;

    void init(xray::UnitCell unit_cell, int n_atom, int* mask_grid_size,
              double* reciprocal_norms, double* mask_cutoffs, int* mask_bs_grid, double shrink_r,
              int n_hkl, complex_double* f_mask, int* hkl) override;

    void update_grid(int n_atom, const double* frac) override;
    void calc_f_bulk() override;

  private:
    void shrink();

    [[nodiscard]] int to_index(int i, int j, int k) const;
  };
}

#endif // AMBER_XRAY_BULK_MASK_CPU_H
