module xray_target_vector_least_squares_impl_cpu_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_target_vector_least_squares_data_module

    implicit none

    private

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_absFcalc
    public :: finalize

contains

    subroutine init(abs_Fobs, sig_Fobs)
        use xray_interface2_data_module, only:  Fobs
        implicit none
        real(real_kind), intent(in) :: abs_Fobs(:)
        real(real_kind), intent(in) :: sig_Fobs(:)

        integer :: i, num_hkl
        real(real_kind) :: phi
        num_hkl = size(abs_Fobs)
        allocate(Fobs(num_hkl))
        do i = 1, num_hkl
            !  sigFobs() here is assumed to be really phi()
            phi = sig_Fobs(i) * 0.0174532925d0
            Fobs(i) = cmplx(abs_Fobs(i) * cos(phi), abs_Fobs(i) * sin(phi), real_kind)
        end do

        norm_scale = 1 / sum(abs_Fobs ** 2)

    end subroutine init

    subroutine finalize()
        deallocate(Fobs)
    end subroutine finalize

    ! This routine computes the force gradient on Fcalc as a harmonic
    ! restraint on the vector (complex) difference between Fcalc and
    ! Fobs
    subroutine calc_partial_d_target_d_Fcalc(Fcalc, deriv, xray_energy)
        implicit none
        complex(real_kind), intent(in) :: Fcalc(:)
        complex(real_kind), intent(out) :: deriv(:)
        real(real_kind), intent(out) :: xray_energy
        complex(real_kind) :: vecdif(size(Fcalc))

        real(real_kind) :: Fcalc_scale

        Fcalc_scale = sum(real(Fobs * conjg(Fcalc))) / sum(abs(Fcalc)**2)

        vecdif(:) = Fobs(:) - Fcalc(:)
        xray_energy = norm_scale * sum(vecdif(:) * conjg(vecdif(:)))
        deriv(:) = - norm_scale * 2 * Fcalc_scale * vecdif(:)

    end subroutine calc_partial_d_target_d_Fcalc

end module xray_target_vector_least_squares_impl_cpu_module
