#include "../include/assert.fh"

module xray_bulk_model_none_module

  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind

  implicit none
  private

  public :: add_bulk_contribution_and_rescale
  public :: finalize
  public :: get_f_scale
  public :: init

  integer, save :: scale_update_period ! in steps
  real(real_kind) :: k_overall

contains

  subroutine init(scale_update_period_)
    implicit none
    integer, intent(in) :: scale_update_period_
    ASSERT(scale_update_period_ > 0)

    scale_update_period = scale_update_period_
  end subroutine init

  subroutine finalize()
  end subroutine finalize

  subroutine add_bulk_contribution_and_rescale(current_step, absFobs, Fcalc)
    use xray_pure_utils, only : calc_k_overall, calc_k_overallc
    use xray_target_module, only : target_function_id
    use xray_interface2_data_module, only : n_work, Fobs
    implicit none
    integer, intent(in) :: current_step
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc

    if(mod(current_step, scale_update_period) == 0) then
      if( target_function_id == 1 ) then
         k_overall = calc_k_overallc(Fobs, Fcalc, n_work)
      else
         k_overall = calc_k_overall(absFobs, abs(Fcalc))
      endif
      write(6,'(a, e14.7)') '| setting isotropic scaling to ', k_overall
    end if

    Fcalc = k_overall * Fcalc

  end subroutine add_bulk_contribution_and_rescale
  
  function get_f_scale(n_hkl)  result(result)
    implicit none
    integer, intent(in) :: n_hkl
    real(real_kind) :: result(n_hkl)
    result = k_overall
  end function get_f_scale

end module xray_bulk_model_none_module
