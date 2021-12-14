module xray_scaling_impl_gpu_module
  
  use xray_pure_utils, only : real_kind
  use xray_scaling_impl_cpu_module, only: optimize_scale_factors, combine, rescale

  implicit none
  
  ! Public module interface
  public :: combine
  public :: finalize
  public :: init
  public :: optimize_scale_factors
  public :: rescale

contains

  subroutine init(resolution, num_work_flags, hkl, max_resolution_bins, n_reflections_in_worst_resolution_bin, min_bin_size)
    use xray_scaling_impl_cpu_module, only : cpu_init => init
    implicit none
    real(real_kind), intent(in) :: resolution(:)
    integer, intent(in) :: num_work_flags
    integer, intent(in) :: hkl(3, size(resolution))

    ! Default values set to match defaults of set.log_binning in cctbx
    integer, intent(in), optional :: max_resolution_bins
    integer, intent(in), optional :: n_reflections_in_worst_resolution_bin
    integer, intent(in), optional :: min_bin_size

    call cpu_init(resolution, num_work_flags, hkl, max_resolution_bins, n_reflections_in_worst_resolution_bin, min_bin_size)
    call gpu_init()
  end subroutine init

  subroutine finalize()
    use xray_scaling_impl_cpu_module, only : cpu_finalize => finalize
    implicit none
    call gpu_finalize()
    call cpu_finalize()
  end subroutine finalize

  subroutine gpu_init()
    implicit none
  end subroutine gpu_init

  subroutine gpu_finalize()
    implicit none
    ! TODO: deallocate smth
  end subroutine gpu_finalize

end module xray_scaling_impl_gpu_module