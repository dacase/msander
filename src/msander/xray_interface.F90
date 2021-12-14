module xray_interface_module

#ifdef CUDA
    use xray_interface_impl_gpu_module
#else
    use xray_interface_impl_cpu_module
#endif

    implicit none

    private

    public :: finalize
    public :: init
    public :: init_clean_up
    public :: xray_get_derivative
    public :: xray_read_mdin
    public :: xray_read_parm
    public :: xray_write_md_state
    public :: xray_write_min_state
    public :: xray_write_options

contains

    subroutine init_clean_up()
        use xray_interface_pre_init_data, only: pre_init_clean_up => clean_up
        call pre_init_clean_up()
    end subroutine init_clean_up

end module xray_interface_module