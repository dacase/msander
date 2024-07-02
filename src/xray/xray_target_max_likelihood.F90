module xray_target_max_likelihood_module

    use xray_target_max_likelihood_impl_cpu_module

    implicit none

    private

    public :: init
    public :: calc_partial_d_target_d_absFcalc
    public :: calc_xray_energy
    public :: finalize

end module xray_target_max_likelihood_module
