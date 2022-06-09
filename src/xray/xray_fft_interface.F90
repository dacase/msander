module xray_fft_interface_module

#define PMEMD_XRAY_FFT_BACKEND_IS_MKL 1

#ifdef PMEMD_XRAY_FFT_BACKEND_IS_FFTW
  use xray_fft_interface_impl_fftw_module
#else
#ifdef PMEMD_XRAY_FFT_BACKEND_IS_MKL
  use xray_fft_interface_impl_mkl_module
#else
#ifdef PMEMD_XRAY_FFT_BACKEND_IS_NONE
  use xray_fft_interface_impl_none_module
#else
#ifdef PMEMD_XRAY_FFT_BACKEND_IS_FFTPACK
!  use xray_fft_interface_impl_fftpack_module
#else
#error "One of PMEMD_XRAY_FFT_BACKEND_IS_* must be set"
#endif
#endif
#endif
#endif
  
  implicit none
  private
  
  public :: dft_i2c_3d

  
end module xray_fft_interface_module
