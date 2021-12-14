module xray_target_least_squares_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind

    implicit none

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_absFcalc
    public :: finalize

    real(real_kind), save :: norm_scale
    integer, save :: n_work

contains

    subroutine init(abs_Fobs_work)
        implicit none
        real(real_kind), intent(in) :: abs_Fobs_work(:)

        n_work = size(abs_Fobs_work)
        norm_scale = 1.0 / sum(abs_Fobs_work ** 2)

    end subroutine init

    subroutine finalize()
    end subroutine finalize

    ! -------------------------------------------------------------------------
    ! This routine computes the force gradient on Fcalc as a harmonic
    ! restraint on the magnitudes of Fobs and Fcalc
    ! -------------------------------------------------------------------------
    subroutine calc_partial_d_target_d_absFcalc(absFobs, absFcalc, weight, deriv, xray_energy)
        implicit none
        real(real_kind), intent(in) :: absFobs(:)
        real(real_kind), intent(in) :: absFcalc(size(absFobs))
        real(real_kind), intent(in), optional :: weight(size(absFobs))
        real(real_kind), intent(out), optional :: deriv(size(absFobs))
        real(real_kind), intent(out), optional :: xray_energy

        if (present(deriv)) then
            deriv(n_work + 1:) = 0 ! no force for things unselected here
            deriv(:n_work) = 2 * (absFcalc(:n_work) - absFobs(:n_work)) * norm_scale
        endif
        if (present(xray_energy))then
            xray_energy = sum((absFobs(:n_work) - absFcalc(:n_work))**2) * norm_scale
        end if
    end subroutine calc_partial_d_target_d_absFcalc

end module xray_target_least_squares_module