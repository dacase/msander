#ifndef AMBER_XRAY_BULK_MASK_GPU_H
#define AMBER_XRAY_BULK_MASK_GPU_H

#include "BulkMask.h"
#include <thrust/device_vector.h>
#include <cufft.h>

namespace xray {

class BulkMaskGPU : public BulkMask {
public:
  BulkMaskGPU() = default;
  ~BulkMaskGPU() override;
  void init(xray::UnitCell unit_cell, int n_atom, int *mask_grid_size,
            double *reciprocal_norms, double *mask_cutoffs,
            int *mask_bs_grid, double shrink_r,
            int n_hkl, complex_double* f_mask, int* hkl) override;
  void update_grid(int n_atom, const double *frac) override;

  void calc_f_bulk() override;

private:
  void shrink();
  thrust::device_vector<int> m_dev_mask_grid;
  thrust::device_vector<int> m_dev_mask_grid_not_shrank;
  thrust::device_vector<GridPoint> m_dev_shrink_neighbours;
  thrust::device_vector<double> m_dev_mask_cutoffs;
  thrust::device_vector<double> m_dev_frac;
  thrust::device_ptr<Sym33> m_dev_metric_tensor;

  thrust::device_vector<double> m_dev_fft_in;
  thrust::device_vector<cufftDoubleComplex> m_dev_fft_out;

  std::vector<cufftDoubleComplex> m_fft_out;

  cufftHandle m_plan;
};
} // namespace xray

#endif // AMBER_XRAY_BULK_MASK_GPU_H
