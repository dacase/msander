! This module MUST be directly used ONLY by xray_scaling_impl_* modules
module xray_scaling_data_module

    use xray_pure_utils, only : real_kind

    implicit none

    public

    integer, save :: n_resolution_bins  ! Number of scaling resolution bins
    integer, save, allocatable :: hkl_scale_resolution_bin(:) ! Index of scale resolution bin per hkl
    integer, save, allocatable :: work_scale_resolution_bin_start(:) ! First inidex of scale resolution bin in WORK set
    integer, save, allocatable :: work_scale_resolution_bin_size(:)  ! Size of scale resolution bin in WORK set
    integer, save, allocatable :: free_scale_resolution_bin_start(:) ! First inidex of scale resolution bin in FREE set
    integer, save, allocatable :: free_scale_resolution_bin_size(:)  ! Size of scale resolution bin in FREE set

    integer, save :: n_work_high_resolution ! Number of WORK reflexes in high (=good) resolution bin

    integer, save :: n_work ! Number of work flags, partition point of mask arrays into [work|free] set

    real(real_kind), save, allocatable :: k_bulk(:)
    real(real_kind), save, allocatable :: k_iso(:)
    real(real_kind), save, allocatable :: k_iso_exp(:)
    real(real_kind), save, allocatable :: k_aniso(:)

    real(real_kind), save :: MUcryst_inv(6, 6)

end module xray_scaling_data_module