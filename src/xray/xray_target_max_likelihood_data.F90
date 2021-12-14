! This module MUST be directly used ONLY by xray_target_max_likelihood_impl_* modules
module xray_target_max_likelihood_data_module

    use xray_pure_utils, only : real_kind

    implicit none

    public

    integer, save :: n_work

    real(real_kind), save, allocatable :: ml_alpha(:)
    real(real_kind), save, allocatable :: ml_beta(:)

    real(real_kind), save, allocatable :: resolution_bin_average_abs_Fobs_pow_2(:)
    real(real_kind), save, allocatable :: resolution_bin_average_abs_Fobs_pow_4(:)

    integer, save :: num_ml_resolution_bins
    integer, allocatable, save :: hkl_index_to_ml_resolution_bin(:)

    integer, allocatable, save :: ml_bin_free_flag_start_index(:) ! First reflex index in resolution bin.
    integer, allocatable, save :: ml_bin_free_flag_count(:)       ! Number of reflexes in ml resolution bin.
    
    integer, save :: ml_update_period

end module xray_target_max_likelihood_data_module