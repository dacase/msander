#ifndef AMBER_XRAY_BULK_MASK_H
#define AMBER_XRAY_BULK_MASK_H

#include "xray_unit_cell.h"
#include "xray_common.h"
#include <vector>
#include <utility>

namespace xray {

  struct GridPoint {
    int x, y, z;
  };

  class BulkMask {
  public:
    BulkMask() noexcept = default;

    virtual ~BulkMask() = default;

    virtual void init(xray::UnitCell unit_cell, int n_atom, int* mask_grid_size,
                      double* reciprocal_norms, double* mask_cutoffs, int* mask_bs_grid, double shrink_r,
                      int n_hkl, complex_double* f_mask, int* hkl
    );

    virtual void update_grid(int n_atom, const double* frac) = 0;
    virtual void calc_f_bulk() = 0;

    [[nodiscard]] static std::pair<int, int> calc_grid_bounds(double x, double dx, int grid_size) noexcept;

  protected:
    xray::UnitCell m_unit_cell;
    int m_n_atom = 0;
    std::array<int, 3> m_grid_dim = {0, 0, 0};
    std::array<double, 3> reciprocal_norms = {0.0, 0.0, 0.0};
    double* mask_cutoffs = nullptr;
    int* mask_grid = nullptr;
    int* m_hkl = nullptr;
    complex_double* m_f_mask = nullptr;
    double shrink_r = 0;
    std::vector<GridPoint> m_shrink_neighbours;
    std::vector<int> m_hkl_grid_index;

  private:
    bool initialized = false;

    void init_shrink_neighbours(const double* reciprocal_norms, double shrink_r);
  };

}

#endif // AMBER_XRAY_BULK_MASK_H
