#ifndef AMBER_XRAY_COMMON_H
#define AMBER_XRAY_COMMON_H
#include <thrust/complex.h>

using complex_double = thrust::complex<double>;

namespace xray {
  enum class KernelPrecision {
    Single,
    Double,
  };

  template<KernelPrecision>
  struct KernelConfig;

  template<>
  struct KernelConfig<KernelPrecision::Single>{
    using FloatType = float;
  };

  template<>
  struct KernelConfig<KernelPrecision::Double>{
    using FloatType = double;
  };
}

#define CUDA_PRECISION Double

#endif //AMBER_XRAY_COMMON_H
