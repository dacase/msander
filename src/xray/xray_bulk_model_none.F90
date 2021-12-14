module xray_bulk_model_none_module

  use xray_contracts_module
  use xray_unit_cell_module
  use xray_pure_utils, only : real_kind

  implicit none
  private

  public :: add_bulk_contribution_and_rescale
  public :: finalize
  public :: init

  integer, save :: scale_update_period ! in steps
  real(real_kind) :: k_overall

contains

  subroutine init(scale_update_period_)
    implicit none
    integer, intent(in) :: scale_update_period_
    call check_precondition(scale_update_period_ > 0)

    scale_update_period = scale_update_period_
  end subroutine init

  subroutine finalize()
  end subroutine finalize

  subroutine add_bulk_contribution_and_rescale(absFobs, Fcalc)
    use xray_pure_utils, only : calc_k_overall
    implicit none
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(inout) :: Fcalc(size(absFobs)) !< input: Fcalc=Fprot, output Fcalc=Fcalc

    integer, save :: nstep = 0

    if(mod(nstep, scale_update_period) == 0) then
      k_overall = calc_k_overall(absFobs, abs(Fcalc))
    end if

    Fcalc = k_overall * Fcalc

    nstep = nstep + 1
  end subroutine add_bulk_contribution_and_rescale

end module xray_bulk_model_none_module