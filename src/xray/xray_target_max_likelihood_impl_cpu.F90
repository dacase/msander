#include "../include/assert.fh"

module xray_target_max_likelihood_impl_cpu_module
  
  use xray_pure_utils, only : real_kind
  use xray_contracts_module
  use xray_target_max_likelihood_data_module

  implicit none

  private

  public :: init
  public :: calc_partial_d_target_d_absFcalc
  public :: finalize

contains

  subroutine init(resolution, abs_Fobs_free, ml_update_period_)
    use xray_pure_utils, only: create_equisized_bins, assign_resolution_bin_indices, set_start_size_from_bin_index
    implicit none
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind), intent(in) :: abs_Fobs_free(:)
    integer, intent(in) :: ml_update_period_
    real(real_kind), parameter :: resolution_tolerance = 1e-5
    integer :: i, bin_start, bin_size
    integer :: num_free_flags
    integer :: num_resolution_bins
    integer, allocatable :: free_hkl_scale_resolution_bin(:)
    
    ml_update_period = ml_update_period_
    num_free_flags = size(abs_Fobs_free)
    n_work = size(resolution) - num_free_flags
    num_resolution_bins = estimate_num_ml_resolution_bins(num_free_flags)
    num_ml_resolution_bins = num_resolution_bins
    
    ! assert num_free_flags > 0
    ! assert num_work_flags > 0
    ! assert is_sorted(resolution(1:num_work_flags))
    ! assert is_sorted(resolution(num_work_flags+1:))

    allocate(ml_alpha(n_work))
    allocate(ml_beta(n_work))
    allocate(hkl_index_to_ml_resolution_bin(n_work))
    
    allocate(resolution_bin_average_abs_Fobs_pow_2(num_resolution_bins))
    allocate(resolution_bin_average_abs_Fobs_pow_4(num_resolution_bins))

    allocate(ml_bin_free_flag_start_index(num_resolution_bins))
    allocate(ml_bin_free_flag_count(num_resolution_bins))
    
    call create_equisized_bins( &
      num_free_flags, num_resolution_bins, &
      ml_bin_free_flag_start_index, ml_bin_free_flag_count &
    )
    
    call assign_resolution_bin_indices( &
        resolution(:n_work), &
        resolution(n_work + ml_bin_free_flag_start_index + ml_bin_free_flag_count - 1) + resolution_tolerance, &
        hkl_index_to_ml_resolution_bin &
    )
    
    allocate(free_hkl_scale_resolution_bin(num_free_flags))
    
    call assign_resolution_bin_indices(&
        resolution(n_work + 1:), &
        resolution(n_work + ml_bin_free_flag_start_index + ml_bin_free_flag_count - 1) + resolution_tolerance, &
        free_hkl_scale_resolution_bin &
    )

    call set_start_size_from_bin_index( &
        free_hkl_scale_resolution_bin, &
        num_resolution_bins, &
        ml_bin_free_flag_start_index, &
        ml_bin_free_flag_count &
    )
    
    deallocate(free_hkl_scale_resolution_bin)
    
    do i = 1, num_resolution_bins
      bin_start = ml_bin_free_flag_start_index(i)
      bin_size = ml_bin_free_flag_count(i)

      resolution_bin_average_abs_Fobs_pow_2(i) = sum(abs_Fobs_free(bin_start:bin_start + bin_size - 1) ** 2) / bin_size
      resolution_bin_average_abs_Fobs_pow_4(i) = sum(abs_Fobs_free(bin_start:bin_start + bin_size - 1) ** 4) / bin_size
    end do
    
  end subroutine init
  
  subroutine finalize()
    
    if(allocated(ml_alpha)) deallocate(ml_alpha)
    if(allocated(ml_beta)) deallocate(ml_beta)
    
    if(allocated(resolution_bin_average_abs_Fobs_pow_2)) &
       deallocate(resolution_bin_average_abs_Fobs_pow_2)
    if(allocated(resolution_bin_average_abs_Fobs_pow_4)) &
       deallocate(resolution_bin_average_abs_Fobs_pow_4)
    
    if(allocated(hkl_index_to_ml_resolution_bin)) &
       deallocate(hkl_index_to_ml_resolution_bin)
    if(allocated(ml_bin_free_flag_start_index)) &
       deallocate(ml_bin_free_flag_start_index)
    if(allocated(ml_bin_free_flag_count)) &
       deallocate(ml_bin_free_flag_count)
    
  end subroutine finalize
  
  
  ! -------------------------------------------------------------------------
  ! Compute derivative of |Fobs|:|Fcalc| Maximum Likelihood (ML) target
  ! with respect to |Fcalc|.
  ! -------------------------------------------------------------------------
  subroutine calc_partial_d_target_d_absFcalc(current_step, absFobs, absFcalc, deriv, xray_energy)
    use xray_pure_utils, only: i1_over_i0, ln_of_i0
    implicit none
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    real(real_kind), intent(in) :: absFcalc(size(absFobs))
    real(real_kind), intent(out), optional :: deriv(size(absFobs))
    real(real_kind), intent(out), optional :: xray_energy
    
    real(real_kind), parameter :: epsilon = 1.0
    
    call estimate_alpha_beta(   &
        current_step, &
        absFobs(n_work + 1:),  &
        absFcalc(n_work + 1:)  &
        )
    
    if (present(deriv)) then
      deriv(n_work + 1:) = 0._real_kind
      deriv(:n_work) = &
          2 * ml_alpha / ml_beta / epsilon * ( &
              ml_alpha * absFcalc(:n_work) &
                  - i1_over_i0( &
                  2 * ml_alpha / ml_beta / epsilon * absFobs(:n_work) * absFcalc(:n_work) &
                  ) *  absFobs(:n_work)  &
              )
    end if
    
    if (present(xray_energy)) then
      xray_energy = &
        sum( &
            - log(2 * absFobs(:n_work) / ml_beta / epsilon) &
                + absFobs(:n_work) ** 2 / ml_beta / epsilon &
                + ml_alpha ** 2 * absFcalc(:n_work) ** 2 / ml_beta / epsilon &
                - ln_of_i0(2 * ml_alpha / ml_beta / epsilon * absFobs(:n_work) * absFcalc(:n_work)) &
            )
    end if
  
  end subroutine calc_partial_d_target_d_absFcalc
  
  pure function estimate_num_ml_resolution_bins(num_free_flags) result(result)
    integer, intent(in) :: num_free_flags
    integer :: result
    real(real_kind), parameter :: max_reflections_per_bin = 140
    
    result = max(1, int(num_free_flags / max_reflections_per_bin + 0.5 ))
  end function estimate_num_ml_resolution_bins
  
  !--------------------------------------------------------------------------------------------
  ! estimate_alpha_beta:  estimate alpha and beta parameters in resolution bins,
  !                       as described in https://doi.org/10.1107/S010876739500688X,
  !                       Appendix A:
  !                       eq. 29     - estimate_t_optimal() estimates root of function G(t)
  !                       eq. 30, 31 - calc_alpha_beta_in_bins() calculates alpha and beta from t
  !--------------------------------------------------------------------------------------------
  subroutine estimate_alpha_beta(current_step, abs_Fobs, abs_Fcalc)
    use xray_pure_utils, only: estimate_t_optimal, calc_bin_alpha_beta, smooth_resolution_bins
    implicit none
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: abs_Fobs(:)
    real(real_kind), intent(in) :: abs_Fcalc(size(abs_Fobs))
    real(real_kind) :: t_optimal(num_ml_resolution_bins)
    real(real_kind) :: alpha_in_bins(num_ml_resolution_bins)
    real(real_kind) :: beta_in_bins(num_ml_resolution_bins)
    real(real_kind) :: A_in_bins(num_ml_resolution_bins)
    real(real_kind) :: p_in_bins_i
    integer :: i, start, count
    
    ! Precondition
    ASSERT(size(abs_Fobs) == size(abs_Fcalc))
    ASSERT(size(abs_Fobs) == sum(ml_bin_free_flag_count))
    
    if (mod(current_step, ml_update_period) /= 0) then
      return
    end if
    
    do i = 1, num_ml_resolution_bins
      start = ml_bin_free_flag_start_index(i) ! index in FREE arrays
      count = ml_bin_free_flag_count(i)

      ! assert resolution_bin_average_abs_Fobs_pow_2(i) == sum(abs_Fobs(start:start + count - 1) ** 2) / count
      ! assert resolution_bin_average_abs_Fobs_pow_4(i) == sum(abs_Fobs(start:start + count - 1) ** 4) / count
      
      A_in_bins(i) = sum(abs_Fcalc(start:start + count - 1) ** 2) / count
      p_in_bins_i = sum(abs_Fcalc(start:start + count - 1) ** 4) / count
      
      t_optimal(i) = estimate_t_optimal( &
        abs_Fobs(start:start + count - 1), &
        abs_Fcalc(start:start + count - 1), &
        resolution_bin_average_abs_Fobs_pow_2(i), &
        resolution_bin_average_abs_Fobs_pow_4(i), &
        A_in_bins(i), &
        p_in_bins_i &
      )
      
    end do
    
    call smooth_resolution_bins(t_optimal)
    
    do i = 1, num_ml_resolution_bins
      call calc_bin_alpha_beta(                   &
        t_optimal(i),                             &
        A_in_bins(i),                             &
        resolution_bin_average_abs_Fobs_pow_2(i), &
        alpha_in_bins(i),                         &
        beta_in_bins(i)                           &
      )
    end do
    
    call smooth_resolution_bins(alpha_in_bins)
    call smooth_resolution_bins(beta_in_bins)
    
    ml_alpha = alpha_in_bins(hkl_index_to_ml_resolution_bin)
    ml_beta = beta_in_bins(hkl_index_to_ml_resolution_bin)
  
  end subroutine estimate_alpha_beta

end module xray_target_max_likelihood_impl_cpu_module
