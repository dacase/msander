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
  real(real_kind) :: k_overall, k_sol, b_sol

contains
  
  subroutine init(scale_update_period_, k_sol_, b_sol_)
    use xray_bulk_mask_module, only : init_mask => init
    implicit none
    integer, intent(in) :: scale_update_period_
    real(real_kind), intent(in) :: k_sol_, b_sol_
    
    call check_precondition(scale_update_period_ > 0)
    scale_update_period = scale_update_period_
    k_sol = k_sol_
    b_sol = b_sol_
  end subroutine init
  
  subroutine finalize()
    implicit none
  end subroutine finalize
  
  subroutine add_bulk_contribution_and_rescale(current_step, absFobs, Fcalc, &
         Fuser, mSS4, hkl)
    use xray_pure_utils, only : calc_k_overall
    use xray_interface2_data_module, only : new_order
    use xray_scaling_module, only : rescale, combine, optimize_scale_factors
    implicit none
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc
    complex(real_kind), intent(in) :: Fuser(size(absFobs))
    real(real_kind), intent(in) :: mSS4(:)
    integer, intent(in) :: hkl(:,:)

    complex(real_kind) :: Fuser_new(size(absFobs))

#if 0
    Fcalc = Fcalc + Fuser(new_order) * k_sol * exp(b_sol * mSS4)

    if(mod(current_step, scale_update_period) == 0) then
      k_overall = calc_k_overall(absFobs, abs(Fcalc))
    end if

    Fcalc = k_overall * Fcalc
#else
    Fuser_new = Fuser(new_order)
    if(mod(current_step, scale_update_period) == 0) then
      call optimize_scale_factors(absFobs, Fcalc, Fuser_new, mSS4, hkl)
    end if

    Fcalc = combine(Fcalc, Fuser_new)
    Fcalc = rescale(Fcalc)
#endif

  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl) result(result)
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    result = k_overall
  end function get_f_scale
  
end module xray_bulk_model_user_module
