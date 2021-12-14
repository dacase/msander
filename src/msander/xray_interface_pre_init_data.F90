! Contains intermediate data that exist between master_setup() and xray_init() calls
!
module xray_interface_pre_init_data

   use xray_pure_utils, only: real_kind

   implicit none

   public

   integer, parameter :: n_scatter_coeffs = 5
   integer, save :: n_scatter_types
   ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
   real(real_kind), allocatable, save :: scatter_coefficients(:,:,:)

contains

   subroutine clean_up()
      if (allocated(scatter_coefficients)) deallocate(scatter_coefficients)
   end subroutine clean_up

end module xray_interface_pre_init_data
