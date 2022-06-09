module xray_interface_impl_gpu_module

    use xray_contracts_module
    use xray_pure_utils, only: real_kind

    ! Reuse functions from CPU module
    use xray_interface_impl_cpu_module, only: &
            & xray_read_mdin, &
            & xray_read_parm, &
            & xray_write_md_state, &
            & xray_write_min_state, &
            & xray_write_options

    implicit none

    private

    public :: finalize
    public :: init
    public :: xray_get_derivative
    public :: xray_read_mdin
    public :: xray_read_parm
    public :: xray_write_md_state
    public :: xray_write_min_state
    public :: xray_write_options

contains

    subroutine init()
        use xray_interface_impl_cpu_module, only : cpu_init => init
        implicit none
        call cpu_init()
        call gpu_init()
    end subroutine init


    subroutine finalize()
        use xray_interface_impl_cpu_module, only : cpu_finalize => finalize
        call gpu_finalize()
        call cpu_finalize()
    end subroutine finalize


    subroutine xray_get_derivative(xyz, force, current_step, xray_e)
    
        use xray_interface_impl_cpu_module, only:  cpu_xray_get_derivative => xray_get_derivative
        use xray_globals_module, only: xray_active, xray_energy

        implicit none
        real(real_kind), intent(inout) :: xyz(:, :)
        real(real_kind), intent(inout) :: force(:,:)
        integer, intent(in) :: current_step
        real(real_kind), intent(out) :: xray_e
        
        call check_precondition(size(xyz, 1) == 3)
        call check_precondition(size(xyz, 2) == size(force, 2))
        call check_precondition(size(force, 1) == 3)
        !      call check_precondition(size(force, 2) == n_atom)
        
        if (.not. xray_active) then
            xray_e = 0
            xray_energy = 0
            return
        end if
 
        ! TODO: sync mss4, unit cell
        ! call gpu_download_crd(xyz)
        ! call gpu_download_frc(force)
        call cpu_xray_get_derivative(xyz, force, current_step, xray_e)
        ! call gpu_upload_frc(force)
        
    end subroutine xray_get_derivative


    subroutine gpu_init()
    end subroutine gpu_init


    subroutine gpu_finalize()
    end subroutine gpu_finalize

end module xray_interface_impl_gpu_module
