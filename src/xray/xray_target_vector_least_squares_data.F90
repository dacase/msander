module xray_target_vector_least_squares_data_module

    use xray_pure_utils, only : real_kind

    implicit none

    public

    complex(real_kind), allocatable, save :: Fobs(:)
    real(real_kind), save :: norm_scale

end module xray_target_vector_least_squares_data_module