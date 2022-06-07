module xray_target_vector_least_squares_impl_gpu_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind
    use xray_target_vector_least_squares_data_module

    implicit none

    private

    ! Public module interface
    public :: init
    public :: calc_partial_d_target_d_Fcalc
    public :: finalize

contains

    subroutine init(abs_Fobs, sig_Fobs)
        use xray_target_vector_least_squares_impl_cpu_module, only : cpu_init => init
        implicit none
        real(real_kind), intent(in) :: abs_Fobs(:)
        real(real_kind), intent(in) :: sig_Fobs(:)
        call cpu_init(abs_Fobs)
        call gpu_init()
    end subroutine init

    subroutine finalize()
        use xray_target_vector_least_squares_impl_cpu_module, only : cpu_finalize => finalize
        call gpu_finalize()
        call cpu_finalize()
    end subroutine finalize

    subroutine calc_partial_d_target_d_Fcalc(Fcalc, deriv, xray_energy)
        implicit none
        complex(real_kind), intent(in) :: Fcalc(:)
        complex(real_kind), intent(out) :: deriv(:)
        real(real_kind), intent(out) :: xray_energy
    end subroutine calc_partial_d_target_d_Fcalc

    subroutine gpu_init()
        implicit none

        ! real(real_kind), allocatable :: r_Fobs(:), i_Fobs(:)
        ! r_Fobs(:) = real(Fobs(:))
        ! i_Fobs(:) = aimag(Fobs(:))
        ! TODO: copy r_Fobs(:), i_Fobs(:) to gpu
    end subroutine gpu_init

    subroutine gpu_finalize()
        implicit none
        ! TODO: dealloc gpu
    end subroutine gpu_finalize

end module xray_target_vector_least_squares_impl_gpu_module
