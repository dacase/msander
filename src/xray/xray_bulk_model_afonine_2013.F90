module xray_bulk_model_afonine_2013_module
  
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

contains
  
  subroutine init(mask_update_period_, scale_update_period_, &
        resolution_high, hkl, unit_cell, atm_atomicnumber, &
        solvent_mask_adjustment, solvent_mask_probe_radius)
    use xray_bulk_mask_module, only : init_mask => init
    implicit none
    
    integer, intent(in) :: mask_update_period_
    integer, intent(in) :: scale_update_period_
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    real(real_kind), intent(in) :: solvent_mask_adjustment
    real(real_kind), intent(in) :: solvent_mask_probe_radius
    
    call check_precondition(mask_update_period_ > 0)
    call check_precondition(scale_update_period_ > 0)
    call check_precondition(mod(scale_update_period_, mask_update_period_) == 0)
    
    mask_update_period = mask_update_period_
    scale_update_period = scale_update_period_
    
    call init_mask(resolution_high, hkl, unit_cell, atm_atomicnumber, &
          solvent_mask_adjustment, solvent_mask_probe_radius)
  
  end subroutine init
  
  subroutine finalize()
    use xray_bulk_mask_module, only : finalize_mask => finalize
    implicit none
    call finalize_mask()
  end subroutine finalize
  
  subroutine add_bulk_contribution_and_rescale(frac, current_step, absFobs, Fcalc, mSS4, hkl)
    use xray_scaling_module, only : rescale, combine, optimize_scale_factors
    use xray_bulk_mask_module, only : update_f_bulk
    use xray_bulk_mask_data_module, only : f_mask
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc
    real(real_kind), intent(in) :: mSS4(:)
    integer, intent(in) :: hkl(:, :)

    call check_precondition(size(frac, 1) == 3)

    if (mod(current_step, mask_update_period) == 0) then
      call update_f_bulk(frac)
    end if
    
    if(mod(current_step, scale_update_period) == 0) then
      call optimize_scale_factors(absFobs, Fcalc, f_mask, mSS4, hkl)
    end if
    
    Fcalc = combine(Fcalc, f_mask)
    Fcalc = rescale(Fcalc)

  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl) result(result)
    use xray_scaling_module, only: scaling_get_f_scale => get_f_scale
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    result = scaling_get_f_scale(n_hkl)
  end function get_f_scale
  
end module xray_bulk_model_afonine_2013_module
