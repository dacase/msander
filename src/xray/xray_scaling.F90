module xray_scaling_module

#ifdef CUDA
    use xray_scaling_impl_gpu_module
#else
    use xray_scaling_impl_cpu_module
#endif
    implicit none

    ! Public module interface
    public :: init
    public :: optimize_scale_factors
    public :: combine
    public :: get_f_scale
    public :: rescale
    public :: finalize

end module xray_scaling_module