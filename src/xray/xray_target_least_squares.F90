module xray_target_least_squares_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_interface2_data_module, only: n_work, sigma_Fobs

    implicit none

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_absFcalc
    public :: finalize

contains

    subroutine init(abs_Fobs_work)
        implicit none
        real(real_kind), intent(in) :: abs_Fobs_work(:)

        ! convert sigma_Fobs to weights:
        sigma_Fobs(:) = 1.d0 / ( 2.0 * sigma_Fobs(:)**2 )

    end subroutine init

    subroutine finalize()
    end subroutine finalize

    ! -------------------------------------------------------------------------
    ! This routine computes the force gradient on Fcalc as a harmonic
    ! restraint on the magnitudes of Fobs and Fcalc
    ! -------------------------------------------------------------------------
    subroutine calc_partial_d_target_d_absFcalc(absFobs, absFcalc, deriv, xray_energy)
        implicit none
        real(real_kind), intent(in) :: absFobs(:)
        real(real_kind), intent(in) :: absFcalc(size(absFobs))
        real(real_kind), intent(out), optional :: deriv(size(absFobs))
        real(real_kind), intent(out), optional :: xray_energy

        if (present(deriv)) then
            deriv(n_work + 1:) = 0 ! no force for things unselected here
            deriv(:n_work) = 2.d0 * sigma_Fobs(:n_work) * &
               (absFcalc(:n_work) - absFobs(:n_work))
        endif
        if (present(xray_energy))then
            xray_energy = sum( sigma_Fobs(:n_work) * &
               (absFobs(:n_work) - absFcalc(:n_work))**2)
        end if
    end subroutine calc_partial_d_target_d_absFcalc

end module xray_target_least_squares_module
