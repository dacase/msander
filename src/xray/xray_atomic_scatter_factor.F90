module xray_atomic_scatter_factor_module

#ifdef CUDA
    use xray_atomic_scatter_factor_impl_gpu_module
#else
    use xray_atomic_scatter_factor_impl_cpu_module
#endif
    implicit none

    private

    public :: init
    public :: finalize
    public :: atomic_scatter_factor

end module xray_atomic_scatter_factor_module