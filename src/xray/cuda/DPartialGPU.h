#ifndef AMBER_DPARTIAL_GPU_H
#define AMBER_DPARTIAL_GPU_H

#include "DPartial.h"
#include <thrust/device_vector.h>
#include <thrust/complex.h>

namespace xray {
  class DPartialGPU : public DPartial {
  public:
    DPartialGPU(
      int n_hkl,
      const int* hkl,
      const double* m_mss4,
      std::complex<double>* f_calc,
      const double* abs_f_calc,
      int n_atom,
      const double* atom_b_factor,
      const int* atom_scatter_type,
      int n_scatter_types,
      const double* atomic_scatter_factor
    );

    ~DPartialGPU() override = default;

    void calc_d_target_d_frac(
      int n_atom,
      const double* frac,
      int n_hkl,
      const double* d_target_d_abs_f_calc,
      double* d_target_d_frac /* result variable */
    ) override;

  private:
    thrust::device_vector<double> m_dev_frac_by_2_pi;
    thrust::device_vector<double> m_dev_d_target_d_abs_f_calc_by_2_pi;
    thrust::device_vector<double> m_dev_atomic_scatter_factor;

    thrust::device_vector<double> m_dev_abs_f_calc;
    thrust::device_vector<double> m_dev_f_calc_phase;
    thrust::device_vector<int> m_dev_hkl;
    thrust::device_vector<double> m_dev_mss4;
    thrust::device_vector<double> m_dev_atom_b_factor;
    thrust::device_vector<int> m_dev_atom_scatter_type;
    thrust::device_vector<double> m_dev_d_target_d_frac;
    std::vector<double> m_f_calc_phase;

  };
}

#endif //AMBER_DPARTIAL_GPU_H
