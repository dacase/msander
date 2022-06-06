#include "../include/assert.fh"

module xray_bulk_model_simple_module
  
  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind
  
  implicit none
  private
  
  public :: add_bulk_contribution_and_rescale
  public :: finalize
  public :: get_f_scale
  public :: init
  
  integer, save :: mask_update_period = 50 ! in steps
  integer, save :: scale_update_period = 50 ! in steps

  real(real_kind) :: k_overall
  real(real_kind) :: k_sol
  real(real_kind) :: b_sol

contains
  
    subroutine init(k_sol_, b_sol_, mask_update_period_, scale_update_period_, &
        resolution_high, hkl, unit_cell, atm_atomicnumber, &
        solvent_mask_adjustment, solvent_mask_probe_radius)
    use xray_bulk_mask_module, only : init_mask => init
    implicit none
    
    real(real_kind), intent(in) :: k_sol_
    real(real_kind), intent(in) :: b_sol_
    integer, intent(in) :: mask_update_period_
    integer, intent(in) :: scale_update_period_
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    real(real_kind), intent(in) :: solvent_mask_adjustment
    real(real_kind), intent(in) :: solvent_mask_probe_radius
    
    ASSERT(mask_update_period_ > 0)
    ASSERT(scale_update_period_ > 0)
    ASSERT(mod(scale_update_period_, mask_update_period_) == 0)
    
    mask_update_period = mask_update_period_
    scale_update_period = scale_update_period_
    k_sol = k_sol_
    b_sol = b_sol_
    
    call init_mask(resolution_high, hkl, unit_cell, atm_atomicnumber, &
        & solvent_mask_adjustment, solvent_mask_probe_radius)
  
  end subroutine init
  
  subroutine finalize()
    use xray_bulk_mask_module, only : finalize_mask => finalize
    implicit none
    call finalize_mask()
  end subroutine finalize
  
  subroutine add_bulk_contribution_and_rescale(frac, current_step, absFobs, &
        Fcalc, mSS4, Fuser)
    use xray_pure_utils, only : calc_k_overall, calc_k_overallc
    use xray_bulk_mask_module, only : update_f_bulk
    use xray_bulk_mask_data_module, only : f_mask
    use xray_target_module, only : target_function_id
    use xray_interface2_data_module, only : Fobs, n_work
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc
    real(real_kind), intent(in) :: mSS4(:)
    complex(real_kind), allocatable, intent(in) :: Fuser(:)
    
    ASSERT(size(frac, 1) == 3)
    
    if (mod(current_step, mask_update_period) == 0) then
      call update_f_bulk(frac, Fuser)
    end if
    
    Fcalc = Fcalc + f_mask * k_sol * exp(b_sol * mSS4)

    if(mod(current_step, scale_update_period) == 0) then
      if( target_function_id == 1 ) then
         k_overall = calc_k_overallc(Fobs, Fcalc, n_work)
      else
         k_overall = calc_k_overall(absFobs, abs(Fcalc))
      endif
    end if

    Fcalc = k_overall * Fcalc

  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl) result(result)
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    result = k_overall
  end function get_f_scale
  
end module xray_bulk_model_simple_module
