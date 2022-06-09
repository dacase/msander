module xray_non_bulk_module

#ifdef CUDA
    use xray_non_bulk_impl_gpu_module
#else
    use xray_non_bulk_impl_cpu_module
#endif

    implicit none

    private

    public :: init
    public :: finalize
    public :: calc_f_non_bulk
    public :: get_f_non_bulk

end module xray_non_bulk_module