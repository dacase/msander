module xray_bulk_model_user_module
  
  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind
  
  implicit none
  private
  
  public :: add_bulk_contribution_and_rescale
  public :: finalize
  public :: get_f_scale
  public :: init
  
  integer, save :: scale_update_period = 50 ! in steps
  real(real_kind) :: k_overall

contains
  
  subroutine init(k_sol, b_sol, scale_update_period_, resolution_high, hkl, unit_cell, atm_atomicnumber)
    use xray_bulk_mask_module, only : init_mask => init
    implicit none
    
    real(real_kind), intent(in) :: k_sol
    real(real_kind), intent(in) :: b_sol
    integer, intent(in) :: scale_update_period_
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    
    call check_precondition(scale_update_period_ > 0)
    scale_update_period = scale_update_period_
  
  end subroutine init
  
  subroutine finalize()
    implicit none
  end subroutine finalize
  
  subroutine add_bulk_contribution_and_rescale(current_step, absFobs, Fcalc, Fuser)
    use xray_pure_utils, only : calc_k_overall
    use xray_bulk_mask_module, only : update_f_bulk
    use xray_bulk_mask_data_module, only : f_mask
    implicit none
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc
    complex(real_kind), intent(in) :: Fuser(size(absFobs))
    
    Fcalc = Fcalc + Fuser

    if(mod(current_step, scale_update_period) == 0) then
      k_overall = calc_k_overall(absFobs, abs(Fcalc))
    end if

    Fcalc = k_overall * Fcalc

  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl) result(result)
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    result = k_overall
  end function get_f_scale
  
end module xray_bulk_model_user_module
