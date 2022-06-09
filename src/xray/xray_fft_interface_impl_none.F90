module xray_fft_interface_impl_none_module
  
  implicit none
  private
  
  public :: dft_i2c_3d

contains

subroutine dft_i2c_3d(dim, data_in, data_out)
  use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
  implicit none
  integer, intent(in) :: dim(3)
  integer, allocatable, intent(in) :: data_in(:)
  complex(8), allocatable, intent(inout) :: data_out(:)
  
  write(stderr, "(A)") "ERROR: xray module requires FFT library."
  write(stderr, "(A)") "This version of pmemd compiled without CPU FFT support."
  write(stderr, "(A)") "Reconfigure with PMEMD_XRAY_FFT_CPU_BACKEND set to MKL or FFTW, e.g.:"
  write(stderr, "(A)") "    cmake ... -DPMEMD_XRAY_FFT_CPU_BACKEND=MKL ... "
  write(stderr, "(A)") "Alternatively, you can use `pmemd.cuda`"
  call exit(1)
end subroutine dft_i2c_3d

  
end module xray_fft_interface_impl_none_module