#include <cassert>

#include "xray/DPartialGPU.h"
#include <thrust/device_vector.h>
#include <thrust/complex.h>
#include <thrust/functional.h>
#include <thrust/transform.h>

namespace {
  template<int blockDimX>
  __global__
  void calc_d_target_d_frac_kernel(
    int n_atom,
    const double* frac_by_2_pi,
    const double* b_factor,
    int n_hkl,
    const int* hkl,
    const double* mSS4,
    const double* atomic_scatter_factor,
    const int* scatter_type_index,
    const double* f_calc_phase,
    const double* abs_f_calc,
    const double* d_target_d_abs_f_calc_by_2_pi,
    double* d_target_d_frac) {

    const int tid = threadIdx.x;
    const int i = blockIdx.x;
    __shared__ double term_x[blockDimX];
    __shared__ double term_y[blockDimX];
    __shared__ double term_z[blockDimX];
    term_x[tid] = 0;
    term_y[tid] = 0;
    term_z[tid] = 0;

    if (i < n_atom) {

      const int offset = (scatter_type_index[i] - 1) * n_hkl;
      const double atom_b_factor = b_factor[i];

      for (int i_hkl = tid; i_hkl < n_hkl; i_hkl += blockDimX) {
        if (abs_f_calc[i_hkl] < 1e-3) {
          continue;
        }

        double phase = -(
          hkl[i_hkl * 3 + 0] * frac_by_2_pi[i * 3 + 0] +
          hkl[i_hkl * 3 + 1] * frac_by_2_pi[i * 3 + 1] +
          hkl[i_hkl * 3 + 2] * frac_by_2_pi[i * 3 + 2]
        );

        double scatter_factor = atomic_scatter_factor[offset + i_hkl];

        double f = scatter_factor * exp(mSS4[i_hkl] * atom_b_factor);
        double tmp = f * d_target_d_abs_f_calc_by_2_pi[i_hkl] * sin(phase + f_calc_phase[i_hkl]);

        term_x[tid] += hkl[i_hkl * 3 + 0] * tmp;
        term_y[tid] += hkl[i_hkl * 3 + 1] * tmp;
        term_z[tid] += hkl[i_hkl * 3 + 2] * tmp;
      }

      __syncthreads();

      for (int s = blockDimX / 2; s > 0; s >>= 1) {
        if (tid < s) {
          term_x[tid] += term_x[tid + s];
          term_y[tid] += term_y[tid + s];
          term_z[tid] += term_z[tid + s];
        }
        __syncthreads();
      }

      if (tid == 0) {
        d_target_d_frac[i * 3 + 0] = term_x[0];
        d_target_d_frac[i * 3 + 1] = term_y[0];
        d_target_d_frac[i * 3 + 2] = term_z[0];
      }
    }
  }

  template<typename T>
  void multiply_by(T factor, thrust::device_vector<T> &dev_vec) {
    using namespace thrust::placeholders;
    thrust::transform(dev_vec.begin(), dev_vec.end(), dev_vec.begin(), factor * thrust::placeholders::_1);
  }
}

xray::DPartialGPU::DPartialGPU(int n_hkl, const int* hkl, const double* mss4, std::complex<double>* f_calc,
                               const double* abs_f_calc, int n_atom, const double* atom_b_factor,
                               const int* atom_scatter_type, int n_scatter_types, const double* atomic_scatter_factor)
  : xray::DPartial(n_hkl, hkl, mss4, f_calc, abs_f_calc, n_atom, atom_b_factor, atom_scatter_type,
                   n_scatter_types, atomic_scatter_factor) {

  m_dev_frac_by_2_pi = thrust::device_vector<double>(n_atom * 3);
  m_dev_d_target_d_abs_f_calc_by_2_pi = thrust::device_vector<double>(n_hkl);
  m_dev_atomic_scatter_factor = thrust::device_vector<double>(m_atomic_scatter_factor,
                                                              m_atomic_scatter_factor + m_n_scatter_types * n_hkl);
  m_dev_abs_f_calc = thrust::device_vector<double>(n_hkl);
  m_dev_f_calc_phase = thrust::device_vector<double>(n_hkl);
  m_dev_hkl = thrust::device_vector<int>(m_hkl, m_hkl + n_hkl * 3);
  m_dev_mss4 = thrust::device_vector<double>(m_mss4, m_mss4 + n_hkl);
  m_dev_atom_b_factor = thrust::device_vector<double>(m_atom_b_factor, m_atom_b_factor + n_atom);
  m_dev_atom_scatter_type = thrust::device_vector<int>(m_atom_scatter_type, m_atom_scatter_type + n_atom);
  m_dev_d_target_d_frac = thrust::device_vector<double>(n_atom * 3);
  m_f_calc_phase = std::vector<double>(n_hkl);
}


void xray::DPartialGPU::calc_d_target_d_frac(
  int n_atom,
  const double* frac,
  int n_hkl,
  const double* d_target_d_abs_f_calc,
  double* d_target_d_frac) {

  assert(n_atom == m_n_atom);
  assert(n_hkl == m_n_hkl);

  thrust::copy(frac, frac + n_atom * 3, m_dev_frac_by_2_pi.begin());
  thrust::copy(m_abs_f_calc, m_abs_f_calc + n_hkl, m_dev_abs_f_calc.begin());
  multiply_by(2 * M_PI, m_dev_frac_by_2_pi);
  thrust::copy(d_target_d_abs_f_calc, d_target_d_abs_f_calc + n_hkl, m_dev_d_target_d_abs_f_calc_by_2_pi.begin());
  multiply_by(2 * M_PI, m_dev_d_target_d_abs_f_calc_by_2_pi);

  for (int i = 0; i < n_hkl; i++) {
    m_f_calc_phase[i] = arg(m_f_calc[i]);
  }

  thrust::copy(m_f_calc_phase.begin(), m_f_calc_phase.end(), m_dev_f_calc_phase.begin());

  dim3 numBlocks(n_atom);
  const int block_size = 64;
  dim3 threadsPerBlock(block_size);

  calc_d_target_d_frac_kernel<block_size>
  <<<numBlocks, threadsPerBlock>>>(
    n_atom,
    m_dev_frac_by_2_pi.data().get(),
    m_dev_atom_b_factor.data().get(),
    m_n_hkl,
    m_dev_hkl.data().get(),
    m_dev_mss4.data().get(),
    m_dev_atomic_scatter_factor.data().get(),
    m_dev_atom_scatter_type.data().get(),
    m_dev_f_calc_phase.data().get(),
    m_dev_abs_f_calc.data().get(),
    m_dev_d_target_d_abs_f_calc_by_2_pi.data().get(),
    m_dev_d_target_d_frac.data().get()
  );

  thrust::copy(m_dev_d_target_d_frac.begin(), m_dev_d_target_d_frac.end(), d_target_d_frac);
}