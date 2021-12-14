#include "BulkMaskGPU.h"
#include <thrust/copy.h>
#include <thrust/device_delete.h>
#include <thrust/device_new.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include <cufft.h>

namespace {

  __device__ float wrap_frac(float x) noexcept {
    if (x < -0.5) {
      return x - floor(x);
    } else if (x >= 0.5) {
      return x - ceil(x);
    } else {
      return x;
    }
  }

  __device__ inline int wrap_index(int i, int n) noexcept {
    if (i <= 0) {
      return i - ((i - n) / n) * n;
    } else if (i > n) {
      return i - ((i - 1) / n) * n;
    } else {
      return i;
    }
  }

  __device__ void calc_grid_bounds(double x, double dx, int grid_size, int &low, int &high) noexcept {
    low = static_cast<int>(floor((x - dx) * grid_size));
    high = static_cast<int>(ceil((x + dx) * grid_size));
  }

  __device__ int to_index(int i, int j, int k, int /* nx */, int ny, int nz) {
    return k + nz * (j + i * ny);
  }

} // namespace

__global__ void cut_atom_vicinity(int n_atoms, int* grid, double norm_x,
                                  double norm_y, double norm_z, int grid_x,
                                  int grid_y, int grid_z, const double* frac,
                                  const double* mask_cutoffs, xray::Sym33* m) {

  const int atom_index = blockIdx.x;

  if (atom_index < n_atoms) {

    const double cutoff = mask_cutoffs[atom_index];
    const double cutoff_sq = cutoff * cutoff;

    const float atomX = frac[atom_index * 3];
    const float atomY = frac[atom_index * 3 + 1];
    const float atomZ = frac[atom_index * 3 + 2];
    const float
        xx = m->xx,
        yy = m->yy,
        zz = m->zz,
        xy = m->xy,
        xz = m->xz,
        yz = m->yz;

    int x_low, x_high;
    int y_low, y_high;
    int z_low, z_high;

    calc_grid_bounds(atomX, cutoff * norm_x, grid_x, x_low, x_high);
    calc_grid_bounds(atomY, cutoff * norm_y, grid_y, y_low, y_high);
    calc_grid_bounds(atomZ, cutoff * norm_z, grid_z, z_low, z_high);

    for (int i = x_low + static_cast<int>(threadIdx.x); i <= x_high; i += static_cast<int>(blockDim.x)) {
      const int mdi = wrap_index(i, grid_x) - 1;
      float dx = wrap_frac(atomX - static_cast<float>(i) / grid_x);
      for (int j = y_low + static_cast<int>(threadIdx.y); j <= y_high; j += static_cast<int>(blockDim.y)) {
        const int mdj = wrap_index(j, grid_y) - 1;
        float dy = wrap_frac(atomY - static_cast<float>(j) / grid_y);
        for (int k = z_low + static_cast<int>(threadIdx.z); k <= z_high; k += static_cast<int>(blockDim.z)) {
          const int mdk = wrap_index(k, grid_z) - 1;
          const int index = mdk + grid_z * (mdj + mdi * grid_y);
          if (grid[index] == 1) {
            float dz = wrap_frac(atomZ - static_cast<float>(k) / grid_z);

            float dist_sq = xx * dx * dx + yy * dy * dy + zz * dz * dz +
                            2 * (xy * dx * dy + xz * dx * dz + yz * dy * dz);

            if (dist_sq < cutoff_sq) {
              grid[index] = 0;
            }
          }
        }
      }
    }
  }
}

__global__ void
shrink_kernel(int* not_shrank, int* mask_grid, int nx, int ny, int nz, xray::GridPoint* d_neigh, int n_neighbours) {
  const int i0 = blockIdx.x;
  const int j0 = blockIdx.y;
  const int k0 = blockIdx.z;
  if (not_shrank[to_index(i0, j0, k0, nx, ny, nz)] == 1) {
    for (int t = threadIdx.x; t < n_neighbours; t += blockDim.x) {
      xray::GridPoint &d = d_neigh[t];
      int i = (i0 + d.x + nx) % nx;
      int j = (j0 + d.y + ny) % ny;
      int k = (k0 + d.z + nz) % nz;
      mask_grid[to_index(i, j, k, nx, ny, nz)] = 1;
    }
  }
}

xray::BulkMaskGPU::~BulkMaskGPU() {
  if (m_dev_metric_tensor) {
    thrust::device_delete(m_dev_metric_tensor);
  }
  cufftDestroy(m_plan);
}

void xray::BulkMaskGPU::init(xray::UnitCell unit_cell, int n_atom,
                             int* mask_grid_size, double* reciprocal_norms,
                             double* mask_cutoffs, int* mask_bs_grid, double shrink_r,
                             int n_hkl, complex_double* f_mask, int* hkl) {
  BulkMask::init(unit_cell, n_atom, mask_grid_size, reciprocal_norms,
                 mask_cutoffs, mask_bs_grid, shrink_r,
                 n_hkl, f_mask, hkl);
  int n_grid_points = mask_grid_size[0] * mask_grid_size[1] * mask_grid_size[2];
  m_dev_mask_cutoffs.resize(n_atom);
  m_dev_frac.resize(n_atom * 3);
  m_dev_metric_tensor = thrust::device_new<Sym33>();
  *m_dev_metric_tensor = unit_cell.metric_tensor();
  m_dev_mask_grid.resize(n_grid_points);
  thrust::copy(mask_cutoffs, mask_cutoffs + n_atom, m_dev_mask_cutoffs.begin());

  m_dev_mask_grid_not_shrank.resize(n_grid_points);
  m_dev_shrink_neighbours.resize(m_shrink_neighbours.size());
  thrust::copy(m_shrink_neighbours.begin(), m_shrink_neighbours.end(), m_dev_shrink_neighbours.begin());

  m_dev_fft_in = thrust::device_vector<double>(m_dev_mask_grid.size());
  m_dev_fft_out = thrust::device_vector<cufftDoubleComplex>((m_grid_dim[2] / 2 + 1) * m_grid_dim[1] * m_grid_dim[0]);
  m_fft_out = std::vector<cufftDoubleComplex>(m_dev_fft_out.size());

  {
      auto result = cufftPlan3d(&m_plan, m_grid_dim[0], m_grid_dim[1], m_grid_dim[2], CUFFT_D2Z);
      assert(result == CUFFT_SUCCESS);
  }
}

void xray::BulkMaskGPU::update_grid(int n_atom, const double* frac) {
  assert(m_dev_mask_cutoffs.size() == n_atom);
  assert(m_dev_frac.size() == n_atom * 3);
  thrust::fill(m_dev_mask_grid.begin(), m_dev_mask_grid.end(), 1);
  thrust::copy(frac, frac + n_atom * 3, m_dev_frac.begin());

  dim3 threadsPerBlock(4, 4, 4);
  dim3 numBlocks(n_atom);

  const int grid_x = m_grid_dim[0];
  const int grid_y = m_grid_dim[1];
  const int grid_z = m_grid_dim[2];

  cut_atom_vicinity<<<numBlocks, threadsPerBlock>>>(
    n_atom, m_dev_mask_grid.data().get(), reciprocal_norms[0],
    reciprocal_norms[1], reciprocal_norms[2], grid_x, grid_y, grid_z,
    m_dev_frac.data().get(), m_dev_mask_cutoffs.data().get(),
    m_dev_metric_tensor.get());

  shrink();

  thrust::copy(m_dev_mask_grid.begin(), m_dev_mask_grid.end(), this->mask_grid);
}

void xray::BulkMaskGPU::shrink() {
  thrust::copy(m_dev_mask_grid.begin(), m_dev_mask_grid.end(), this->m_dev_mask_grid_not_shrank.begin());
  dim3 numBlocks(m_grid_dim[0], m_grid_dim[1], m_grid_dim[2]);
  dim3 threadsPerBlock(64);
  //int* not_shrank, int* mask_grid, int nx, int ny, int nz, xray::GridPoint* d_neigh, int n_neighbours
  shrink_kernel<<<numBlocks, threadsPerBlock>>>(
    m_dev_mask_grid_not_shrank.data().get(),
    m_dev_mask_grid.data().get(),
    m_grid_dim[0], m_grid_dim[1], m_grid_dim[2],
    m_dev_shrink_neighbours.data().get(),
    static_cast<int>(m_dev_shrink_neighbours.size())
  );
}

void xray::BulkMaskGPU::calc_f_bulk() {
  const int grid_size = m_grid_dim[0] * m_grid_dim[1] * m_grid_dim[2];

  thrust::copy(m_dev_mask_grid.begin(), m_dev_mask_grid.end(), m_dev_fft_in.begin());

  {
    auto result = cufftExecD2Z(m_plan, m_dev_fft_in.data().get(), m_dev_fft_out.data().get());
    assert(result == CUFFT_SUCCESS);
  }

  thrust::copy(m_dev_fft_out.begin(), m_dev_fft_out.end(), m_fft_out.begin());

  double factor = m_unit_cell.get_volume() / grid_size;
  for (int i = 0; i < m_hkl_grid_index.size(); ++i) {
    reinterpret_cast<std::complex<double> &>(m_f_mask[i]) =
      std::conj(reinterpret_cast<std::complex<double> &>(m_fft_out[m_hkl_grid_index[i]])) * factor;
  }
}
