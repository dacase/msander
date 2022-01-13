! This module relies on FFTW3 Interface to Intel® oneAPI Math Kernel Library
! therefore implementation is 1-to-1 copy of xray_fft_interface_impl_fftw_module
module xray_fft_interface_impl_mkl_module
  
  implicit none
  private
  
  public :: dft_i2c_3d

contains
  
  subroutine dft_i2c_3d(dim, data_in, data_out)
  
    implicit none
    include "fftw3.f"
    
    integer, intent(in) :: dim(3)
    integer, allocatable, intent(in) :: data_in(:)
    complex(8), allocatable, intent(inout) :: data_out(:)
    
    double precision :: mask_bs_grid_3d(dim(3), dim(2), dim(1))
    double complex ::   mask_bs_grid_3d_fft(dim(3)/2 + 1, dim(2), dim(1))
    integer(8) :: plan_forward(1)
    
    mask_bs_grid_3d = reshape(data_in, [dim(3) , dim(2), dim(1)])
    call dfftw_plan_dft_r2c_3d(plan_forward, dim(3), &
        dim(2), dim(1), mask_bs_grid_3d, &
        mask_bs_grid_3d_fft, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan_forward, mask_bs_grid_3d, &
        mask_bs_grid_3d_fft)
    data_out = reshape(mask_bs_grid_3d_fft, &
        [(dim(3) / 2 + 1) * dim(2) * dim(1)])
    call dfftw_destroy_plan(plan_forward)

  end subroutine dft_i2c_3d
  
end module xray_fft_interface_impl_mkl_module