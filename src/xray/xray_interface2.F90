module xray_interface2_module
  
  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind
  
  implicit none
  
  private
  
  public :: finalize
  public :: get_r_factors
  public :: get_f_calc
  public :: init
  public :: calc_force
  
  ! Globals
  
  real(real_kind), allocatable, save:: abs_Fobs(:)
  real(real_kind), allocatable, save:: abs_Fcalc(:)
  real(real_kind), allocatable, save:: resolution(:)
  real(real_kind), allocatable, save:: mSS4(:) ! minus s^2 divided by 4

contains
  
  subroutine init(target, bulk_model, hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
      &   atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk, &
      &   atom_atomic_number, mask_update_period, scale_update_period, target_meta_update_period,&
      &   k_sol, b_sol)
    use xray_interface2_data_module, only : init_data => init
    use xray_pure_utils, only : index_partition, index_sort, calc_resolution
    use constants_xray, only : set_xray_num_threads
    
    implicit none

    character(len=*), intent(in) :: target
    character(len=*), intent(in) :: bulk_model
    
    integer, intent(in) :: hkl(:, :)
    complex(real_kind), intent(in) :: Fobs(:)
    real(real_kind), intent(in) :: sigma_Fobs(:)
    logical, intent(in) :: work_flag(:)
    class(unit_cell_t), intent(in) :: unit_cell
    real(real_kind), intent(in) :: scatter_coefficients(:, :, :) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
    real(real_kind), intent(in) :: atom_b_factor(:)
    real(real_kind), intent(in) :: atom_occupancy(:)
    integer, intent(in) :: atom_scatter_type(:)
    logical, intent(in) :: atom_is_not_bulk(:)
    integer, intent(in) :: atom_atomic_number(:)
    integer, intent(in) :: mask_update_period
    integer, intent(in) :: scale_update_period
    integer, intent(in) :: target_meta_update_period
    real(real_kind), intent(in) :: k_sol
    real(real_kind), intent(in) :: b_sol

    call check_precondition(size(hkl, 1) == 3)
    call check_precondition(size(hkl, 2) == size(Fobs))
    call check_precondition(size(hkl, 2) == size(sigma_Fobs))
    call check_precondition(size(hkl, 2) == size(work_flag))
    call check_precondition(size(atom_b_factor) == size(atom_occupancy))
    call check_precondition(size(atom_b_factor) == size(atom_scatter_type))
    call check_precondition(size(atom_b_factor) == size(atom_atomic_number))
    call check_precondition(size(atom_b_factor) == size(atom_is_not_bulk))
    call check_precondition(minval(atom_b_factor, atom_is_not_bulk) >= 0)
    call check_precondition(all(atom_occupancy <= 1.0))
    call check_precondition(all(atom_occupancy >= 0.0))
    call check_precondition(minval(atom_scatter_type) >= 1)
    call check_precondition(maxval(atom_scatter_type) <= size(scatter_coefficients, 3))
    
    call set_xray_num_threads()

    call init_data(hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
        &   atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk)
  
    call init_submodules(target, bulk_model, atom_atomic_number, mask_update_period, scale_update_period, target_meta_update_period, k_sol, b_sol)

  end subroutine init
  
  subroutine calc_force(xyz, current_step, xray_weight, force, energy)
    use xray_interface2_data_module
    use xray_target_module, only: calc_partial_d_target_d_absFcalc
    use xray_non_bulk_module, only: calc_f_non_bulk, get_f_non_bulk
    use xray_bulk_model_module, only: add_bulk_contribution_and_rescale, get_f_scale
    use xray_dpartial_module, only: calc_partial_d_target_d_frac
    implicit none
    real(real_kind), intent(in) :: xyz(:, :)
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: xray_weight
    real(real_kind), intent(out) :: force(:, :)
    real(real_kind), intent(out) :: energy
    real(real_kind), allocatable :: d_target_d_absFcalc(:)
    real(real_kind), allocatable :: frac(:, :)
    real(real_kind), allocatable :: grad_xyz(:, :)
#include "../msander/def_time.h"

    call check_precondition(size(xyz, 1) == 3)
    call check_precondition(size(xyz, 2) == n_atom)
    call check_precondition(size(force, 1) == 3)
    call check_precondition(size(force, 2) == n_atom)
    
    allocate(grad_xyz(3, size(non_bulk_atom_indices)))
    
    frac = modulo(unit_cell%to_frac(xyz(:, non_bulk_atom_indices)), 1.0_real_kind)
    
    call check_assertion(all(frac <= 1))
    call check_assertion(all(frac >= 0))

    call timer_start(TIME_IHKL)
    call calc_f_non_bulk(frac)
    Fcalc = get_f_non_bulk()
    call timer_stop(TIME_IHKL)

    call add_bulk_contribution_and_rescale(&
        frac, &
        current_step, &
        abs_Fobs, Fcalc, &
        mSS4, hkl &
    )

    abs_Fcalc(:) = abs(Fcalc(:))
    
    allocate(d_target_d_absFcalc(size(Fcalc)))
    call calc_partial_d_target_d_absFcalc(current_step, &
            abs_Fobs, abs_Fcalc, &
            deriv=d_target_d_absFcalc, xray_energy=energy &
    )

    energy = xray_weight * energy
    call timer_start(TIME_DHKL)
    grad_xyz = xray_weight * unit_cell%to_orth_derivative( &
        & calc_partial_d_target_d_frac(frac, get_f_scale(size(abs_Fobs)), d_target_d_absFcalc) &
        &)
    call check_assertion(size(grad_xyz, 2) == size(non_bulk_atom_indices))

    force(:,non_bulk_atom_indices) = force(:,non_bulk_atom_indices) - grad_xyz
    call timer_stop(TIME_DHKL)

  end subroutine calc_force
  
  subroutine get_r_factors(r_work, r_free)
    use xray_interface2_data_module, only: n_work
    use xray_pure_utils, only: calc_r_factor
    implicit none
    real(real_kind), intent(out) :: r_work
    real(real_kind), intent(out) :: r_free
    r_work = calc_r_factor(abs_Fobs(:n_work), abs_Fcalc(:n_work))
    r_free = calc_r_factor(abs_Fobs(n_work+1:), abs_Fcalc(n_work+1:))
  end subroutine get_r_factors

  function get_f_calc() result(result)
    use xray_interface2_data_module, only: Fcalc, hkl_io_order
    implicit none
    complex(real_kind), allocatable :: result(:)
    result = Fcalc(hkl_io_order)
  end function get_f_calc
  
  subroutine finalize()
    use xray_interface2_data_module, only : finalize_data => finalize
    implicit none
    call finalize_submodules()
    call finalize_data()
  end subroutine finalize
  
  ! Private procedures
  
  subroutine init_submodules(target, bulk_model, atom_atomic_number, mask_update_period, scale_update_period, target_meta_update_period, k_sol, b_sol)
    use xray_interface2_data_module
    
    use xray_atomic_scatter_factor_module, only: init_atomic_scatter_factor => init
    use xray_bulk_model_module, only: init_bulk => init
    use xray_non_bulk_module, only: init_non_bulk => init
    use xray_scaling_module, only: init_scaling => init
    use xray_target_module, only: init_target => init
    use xray_dpartial_module, only: init_dpartial => init
    
    implicit none
    
    character(len=*), intent(in) :: target
    character(len=*), intent(in) :: bulk_model
    integer, intent(in) :: atom_atomic_number(:)
    integer, intent(in) :: mask_update_period
    integer, intent(in) :: scale_update_period
    integer, intent(in) :: target_meta_update_period
    real(real_kind), intent(in) :: k_sol
    real(real_kind), intent(in) :: b_sol
    
    real(real_kind), allocatable :: s2(:)
    
    s2 = unit_cell%get_s2(hkl)
    resolution = 1 / sqrt(s2)
    mSS4 = - s2 / 4
    
    abs_Fobs = abs(Fobs)
    allocate(abs_Fcalc(n_hkl))
    
    call init_target(target, resolution, n_work, abs_Fobs, sigma_Fobs, target_meta_update_period)
    call init_bulk(bulk_model, mask_update_period, scale_update_period, minval(resolution), hkl, unit_cell, atom_atomic_number(non_bulk_atom_indices), k_sol, b_sol)
    call init_scaling(resolution, n_work, hkl)
    call init_atomic_scatter_factor(mSS4, scatter_coefficients)
    call init_non_bulk(hkl, mSS4, atom_b_factor(non_bulk_atom_indices), atom_scatter_type(non_bulk_atom_indices), atom_occupancy(non_bulk_atom_indices))
    call init_dpartial(hkl, mss4, Fcalc, abs_Fcalc, atom_b_factor(non_bulk_atom_indices), atom_scatter_type(non_bulk_atom_indices))

  end subroutine init_submodules
  
  subroutine finalize_submodules()
    use xray_atomic_scatter_factor_module, only: finalize_atomic_scatter_factor => finalize
    use xray_bulk_model_module, only: finalize_bulk => finalize
    use xray_non_bulk_module, only: finalize_non_bulk => finalize
    use xray_scaling_module, only: finalize_scaling => finalize
    use xray_target_module, only: finalize_target => finalize
    use xray_dpartial_module, only: finalize_dpartial => finalize
  
    implicit none
    
    call finalize_dpartial()
    call finalize_non_bulk()
    call finalize_atomic_scatter_factor()
    call finalize_scaling()
    call finalize_bulk()
    call finalize_target()
    
    if(allocated(abs_Fobs)) deallocate(abs_Fobs)
    if(allocated(abs_Fcalc)) deallocate(abs_Fcalc)
    if(allocated(resolution)) deallocate(resolution)
    if(allocated(mSS4)) deallocate(mSS4)
  end subroutine finalize_submodules
  
  
end module xray_interface2_module
