#include "BulkMaskCPU.h"
#include <array>
#include <cassert>
#include <cmath>
#include <tuple>
#include <vector>
#include <fftw3.h>

namespace {

  [[nodiscard]] inline double wrap_frac(double x) noexcept {
    if (x < -0.5) {
      return x - std::floor(x);
    } else if (x >= 0.5) {
      return x - std::ceil(x);
    } else {
      return x;
    }
  }

  [[nodiscard]] inline int wrap_index(int i, int n) noexcept {
    if (i <= 0) {
      return i - ((i - n) / n) * n;
    } else if (i > n) {
      return i - ((i - 1) / n) * n;
    } else {
      return i;
    }
  }

} // namespace

void xray::BulkMaskCPU::init(xray::UnitCell unit_cell, int n_atom, int* mask_grid_size, double* reciprocal_norms,
                             double* mask_cutoffs, int* mask_bs_grid, double shrink_r,
                             int n_hkl, complex_double* f_mask, int* hkl) {
  BulkMask::init(unit_cell, n_atom, mask_grid_size, reciprocal_norms, mask_cutoffs, mask_bs_grid, shrink_r,
                 n_hkl, f_mask, hkl);
}

void xray::BulkMaskCPU::update_grid(int n_atom, const double* frac) {
  assert(n_atom == this->m_n_atom);

  const int grid_size = m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2];
  std::fill(mask_grid, mask_grid + grid_size, 1);

  for (int tid = 0; tid <= n_atom; ++tid) {
    const double cutoff = mask_cutoffs[tid];
    const double cutoff_sq = cutoff * cutoff;

    const double atomX = frac[tid * 3];
    const double atomY = frac[tid * 3 + 1];
    const double atomZ = frac[tid * 3 + 2];

    int x_low, x_high;
    int y_low, y_high;
    int z_low, z_high;
    std::tie(x_low, x_high) = calc_grid_bounds(
      atomX, cutoff * reciprocal_norms[0], m_grid_dim[0]);
    std::tie(y_low, y_high) = calc_grid_bounds(
      atomY, cutoff * reciprocal_norms[1], m_grid_dim[1]);
    std::tie(z_low, z_high) = calc_grid_bounds(
      atomZ, cutoff * reciprocal_norms[2], m_grid_dim[2]);

    for (int i = x_low; i <= x_high; ++i) {
      double dx = wrap_frac(atomX - static_cast<double>(i) / m_grid_dim[0]);
      int mdi = wrap_index(i, m_grid_dim[0]) - 1;
      for (int j = y_low; j <= y_high; ++j) {
        double dy =
          wrap_frac(atomY - static_cast<double>(j) / m_grid_dim[1]);
        int mdj = wrap_index(j, m_grid_dim[1]) - 1;
        for (int k = z_low; k <= z_high; ++k) {
          double dz =
            wrap_frac(atomZ - static_cast<double>(k) / m_grid_dim[2]);
          int mdk = wrap_index(k, m_grid_dim[2]) - 1;
          int index = to_index(mdi, mdj, mdk);
          if (mask_grid[index] != 0) {
            double dist_sq = m_unit_cell.to_squared_orth_norm(dx, dy, dz);

            if (dist_sq < cutoff_sq) {
              mask_grid[index] = 0;
            }
          }
        }
      }
    }
  }
  shrink();
}

void xray::BulkMaskCPU::calc_f_bulk() {
  const int grid_size = m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2];

  std::vector<double> in;
  std::vector<std::complex<double>> out;

  // https://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html
  out.resize((m_grid_dim[2] / 2 + 1) * m_grid_dim[1] * m_grid_dim[0]);
  in.resize(grid_size);
  std::copy(mask_grid, mask_grid + grid_size, in.begin());

  // https://fftw.org/doc/Reversing-array-dimensions.html
  auto plan = fftw_plan_dft_r2c_3d(
    m_grid_dim[0], m_grid_dim[1], m_grid_dim[2],
    in.data(), reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE
  );

  fftw_execute_dft_r2c(plan, in.data(), reinterpret_cast<fftw_complex*>(out.data()));

  double factor = m_unit_cell.get_volume() / grid_size;
  for (int i = 0; i < m_hkl_grid_index.size(); ++i) {
    reinterpret_cast<std::complex<double> &>(m_f_mask[i]) = std::conj(out[m_hkl_grid_index[i]]) * factor;
  }

  fftw_destroy_plan(plan);
}


void xray::BulkMaskCPU::shrink() {

  const int grid_size = m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2];
  std::vector<int> not_shrank(mask_grid, mask_grid + grid_size);

  for (int i0 = 0; i0 < m_grid_dim[0]; ++i0) {
    for (int j0 = 0; j0 < m_grid_dim[1]; ++j0) {
      for (int k0 = 0; k0 < m_grid_dim[2]; ++k0) {
        if (not_shrank[to_index(i0, j0, k0)] == 1) {
          for (auto &d: m_shrink_neighbours) {
            int i = (i0 + d.x + m_grid_dim[0]) % m_grid_dim[0];
            int j = (j0 + d.y + m_grid_dim[1]) % m_grid_dim[1];
            int k = (k0 + d.z + m_grid_dim[2]) % m_grid_dim[2];
            mask_grid[to_index(i, j, k)] = 1;
          }
        }
      }
    }
  }
}

int xray::BulkMaskCPU::to_index(int i, int j, int k) const {
  assert(0 <= i && i < m_grid_dim[0]);
  assert(0 <= j && j < m_grid_dim[1]);
  assert(0 <= k && k < m_grid_dim[2]);
  return k + m_grid_dim[2] * (j + i * m_grid_dim[1]);
}
