module xray_target_module

    use xray_contracts_module
    use xray_pure_utils, only: real_kind
    implicit none

    private

    public :: init
    public :: finalize
    public :: calc_partial_d_target_d_absFcalc

    ! Enumeration
    !   0 -- Least Squares
    !   1 -- Vector Least Squares
    !   2 -- Max Likelihood
    integer :: target_function_id
    integer, parameter :: least_squares_id = 0
    integer, parameter :: vector_least_squares_id = 1
    integer, parameter :: max_likehood_id = 2

    character(len = *), parameter :: least_squares_name = "ls"
    character(len = *), parameter :: vector_least_squares_name = "vls"
    character(len = *), parameter :: max_likehood_name = "ml"

contains

    subroutine init(target, resolution, num_work_flags, absFobs, sigFobs, meta_update_period)
        use xray_target_least_squares_module, only : ls_init => init
        use xray_target_vector_least_squares_module, only : vls_init => init
        use xray_target_max_likelihood_module, only : ml_init => init
        implicit none
        character(len = *), intent(in) :: target
        real(real_kind), intent(in) :: resolution(:)
        integer, intent(in) :: num_work_flags
        real(real_kind), intent(in) :: absFobs(:)
        real(real_kind), intent(in) :: sigFobs(:)
        integer, intent(in) :: meta_update_period
        target_function_id = target_function_name_to_id(target)

        select case(target_function_id)
        case (least_squares_id)
            call ls_init(absFobs(:num_work_flags))
        case(vector_least_squares_id)
            call vls_init(absFobs, sigFobs)
        case (max_likehood_id)
            call ml_init(resolution, absFobs(num_work_flags + 1:), meta_update_period)
        end select
    end subroutine init

    subroutine finalize()
        use xray_target_least_squares_module, only : ls_finalize => finalize
        use xray_target_vector_least_squares_module, only : vls_finalize => finalize
        use xray_target_max_likelihood_module, only : ml_finalize => finalize
        implicit none
        select case(target_function_id)
        case (least_squares_id)
            call ls_finalize()
        case (vector_least_squares_id)
            call vls_finalize()
        case (max_likehood_id)
            call ml_finalize()
        end select
    end subroutine finalize

    subroutine calc_partial_d_target_d_absFcalc(current_step, absFobs, absFcalc, deriv, xray_energy)
        use xray_target_least_squares_module, only : ls_partial => calc_partial_d_target_d_absFcalc
        use xray_target_max_likelihood_module, only : ml_partial => calc_partial_d_target_d_absFcalc
        implicit none
        integer, intent(in) :: current_step
        real(real_kind), intent(in) :: absFobs(:)
        real(real_kind), intent(in) :: absFcalc(size(absFobs))
        real(real_kind), intent(out) :: deriv(size(absFobs))
        real(real_kind), intent(out), optional :: xray_energy

        logical, save :: first_call = .TRUE.

        call check_requirement(.not. first_call .or. current_step == 0, &
            & "First call of `xray_target_module::calc_partial_d_target_d_absFcalc(...)` &
            & must be made with current_step=0")
        first_call = .FALSE.

        select case(target_function_id)
        case (least_squares_id)
            call ls_partial(absFobs, absFcalc, deriv=deriv, xray_energy=xray_energy)
        case (vector_least_squares_id)
            ! FIXME
        case (max_likehood_id)
            call ml_partial(current_step, absFobs, absFcalc, deriv, xray_energy)
        end select
    end subroutine calc_partial_d_target_d_absFcalc

    function target_function_name_to_id(target) result(f_id)
        implicit none
        character(len = *), intent(in) :: target
        integer :: f_id
        select case(trim(target))
        case (least_squares_name)
            f_id = least_squares_id
        case (vector_least_squares_name)
            f_id = vector_least_squares_id
        case (max_likehood_name)
            f_id = max_likehood_id
        case default
            call check_requirement(.FALSE., "Unknown target function name: '" // trim(target) // "'")
        end select
    end function target_function_name_to_id

end module xray_target_module