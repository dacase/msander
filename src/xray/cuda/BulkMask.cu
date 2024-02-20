#include "BulkMask.h"
#include <cassert>
#include <tuple>


namespace {
  int hkl_to_grid_index(int h, int k, int l, int nx, int ny, int nz) {
    // assert(-(nx - 1) / 2 <= h && h <= (nx - 1) / 2);
    // assert(-(ny - 1) / 2 <= k && k <= (ny - 1) / 2);
    // assert(0 <= l && l <= nz / 2 + 1);

    if (h < 0) {
      h += nx;
    }

    if (k < 0) {
      k += ny;
    }

    nz = nz / 2 + 1;
    return h * ny * nz + k * nz + l;
  }
}

void xray::BulkMask::init(xray::UnitCell unit_cell, int n_atom,
                          int* grid_dim, double* reciprocal_norms,
                          double* mask_cutoffs, int* mask_bs_grid, double shrink_r,
                          int n_hkl, complex_double* f_mask, int* hkl) {
  assert(!initialized);
  this->m_unit_cell = unit_cell;
  this->m_n_atom = n_atom;
  for (int i = 0; i < 3; i++) {
    this->m_grid_dim[i] = grid_dim[i];
    this->reciprocal_norms[i] = reciprocal_norms[i];
  }
  this->mask_cutoffs = mask_cutoffs;
  this->mask_grid = mask_bs_grid;
  this->m_f_mask = f_mask;
  this->m_hkl = hkl;
  this->shrink_r = shrink_r;

  init_shrink_neighbours(reciprocal_norms, shrink_r);

  m_hkl_grid_index.resize(n_hkl);
  for (int i = 0; i < n_hkl; ++i) {
    int h = hkl[i * 3];
    int k = hkl[i * 3 + 1];
    int l = hkl[i * 3 + 2];
    m_hkl_grid_index[i] = hkl_to_grid_index(h, k, l, grid_dim[0], grid_dim[1], grid_dim[2]);
  }

  this->initialized = true;
}

void xray::BulkMask::init_shrink_neighbours(const double* reciprocal_norms, double shrink_r) {
  int x_low, x_high;
  int y_low, y_high;
  int z_low, z_high;

  std::tie(x_low, x_high) = calc_grid_bounds(0, shrink_r * reciprocal_norms[0], m_grid_dim[0]);
  std::tie(y_low, y_high) = calc_grid_bounds(0, shrink_r * reciprocal_norms[1], m_grid_dim[1]);
  std::tie(z_low, z_high) = calc_grid_bounds(0, shrink_r * reciprocal_norms[2], m_grid_dim[2]);
  double cutoff_sq = shrink_r * shrink_r;
  for (int di = x_low; di <= x_high; ++di) {
    double dx = static_cast<double>(di) / m_grid_dim[0];
    for (int dj = y_low; dj <= y_high; ++dj) {
      double dy = static_cast<double>(dj) / m_grid_dim[1];
      for (int dk = z_low; dk <= z_high; ++dk) {
        double dz = static_cast<double>(dk) / m_grid_dim[2];
        double dist_sq = m_unit_cell.to_squared_orth_norm(dx, dy, dz);
        if (dist_sq < cutoff_sq) {
          m_shrink_neighbours.push_back(GridPoint{di, dj, dk});
        }
      }
    }
  }
}

std::pair<int, int> xray::BulkMask::calc_grid_bounds(double x, double dx, int grid_size) noexcept {
  return std::make_pair(static_cast<int>(std::floor((x - dx) * grid_size)),
                        static_cast<int>(std::ceil((x + dx) * grid_size)));
}
