module xray_target_vector_least_squares_impl_cpu_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_interface2_data_module, only : n_hkl, Fobs
    use xray_target_vector_least_squares_data_module

    implicit none

    private

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_Fcalc
    public :: finalize

contains

    subroutine init(abs_Fobs)
        implicit none
        real(real_kind), intent(in) :: abs_Fobs(:)

        allocate(derivc(n_hkl))
        norm_scale = 1 / sum(abs_Fobs ** 2)

    end subroutine init

    subroutine finalize()
    end subroutine finalize

    ! This routine computes the force gradient on Fcalc as a harmonic
    ! restraint on the vector (complex) difference between Fcalc and
    ! Fobs

    ! deriv is treated as a complex variable, but this is not really
    ! complex analysis!  Rather, deriv%re is the derivative of the target
    ! with respect to Fcalc%re, and deriv%im is the derivative of the 
    ! target with respect to Fcalc%im.

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
