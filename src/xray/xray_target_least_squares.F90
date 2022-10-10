module xray_target_least_squares_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_interface2_data_module, only: n_work, sigma_Fobs, penalty

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
        use xray_globals_module, only : ls_r3, ls_r4
        implicit none
        real(real_kind), intent(in) :: absFobs(:)
        real(real_kind), intent(in) :: absFcalc(size(absFobs))
        real(real_kind), intent(out), optional :: deriv(size(absFobs))
        real(real_kind), intent(out), optional :: xray_energy

        real(real_kind) :: e,rij,r1,dif1,r2,df,r3,r4,dif
        integer :: i

#if 1   /* optimize R_work, based on disnrg()  */
        deriv(n_work + 1:) = 0 ! no force for things unselected here
        penalty(n_work + 1:) = 0 ! no penalty for things unselected here
        xray_energy = 0

        ! disnrg() -like parameters
        r1 = -ls_r4
        r2 = -ls_r3
        r3 =  ls_r3
        r4 =  ls_r4

        do i = 1,n_work
           rij = absFcalc(i) - absFobs(i)
           if (rij < r1) then
               dif1 = r1-r2
               df = 2.0d0 * dif1
               e = df * (rij-r1) + dif1*dif1
            else if (rij < r2) then
               dif = rij - r2
               df = 2.0d0 * dif
               e = dif*dif
            else if (rij <= r3) then
               df = 0.d0
               e = 0.0d0
            else if (rij < r4) then
               dif = rij - r3
               df = 2.0d0 * dif
               e = dif*dif
            else
               dif1 = r4-r3
               df = 2.0d0 * dif1
               e = df * (rij-r4) + dif1*dif1
            end if
            deriv(i) = sigma_Fobs(i)*df
            penalty(i) = sigma_Fobs(i)*e
            xray_energy = xray_energy + e
        end do

#else  /* weighted least-squares  */
        if (present(deriv)) then
            deriv(n_work + 1:) = 0 ! no force for things unselected here
            deriv(:n_work) = 2.d0 * sigma_Fobs(:n_work) * &
               (absFcalc(:n_work) - absFobs(:n_work))
        endif
        if (present(xray_energy))then
            xray_energy = sum( sigma_Fobs(:n_work) * &
               (absFobs(:n_work) - absFcalc(:n_work))**2)
        end if
#endif

    end subroutine calc_partial_d_target_d_absFcalc

end module xray_target_least_squares_module
