#include "../include/assert.fh"

module xray_atomic_scatter_factor_impl_gpu_module

    use xray_contracts_module
    use xray_pure_utils, only : real_kind

    use xray_atomic_scatter_factor_impl_cpu_module, only: atomic_scatter_factor

    implicit none

    private

    public :: init
    public :: finalize
    public :: atomic_scatter_factor

contains
    
    subroutine init(mSS4, scatter_coefficients)
        use xray_atomic_scatter_factor_impl_cpu_module, only : cpu_init => init
        implicit none
        real(real_kind), intent(in) :: mSS4(:)
        real(real_kind), intent(in) :: scatter_coefficients(:, :, :)
        
        integer :: ihkl, i
        integer :: n_hkl, n_scatter_types
        
        ASSERT(size(scatter_coefficients, 1) == 2)

        call cpu_init(mSS4, scatter_coefficients)
        call gpu_init()
    end subroutine init

    subroutine finalize()
        use xray_atomic_scatter_factor_impl_cpu_module, only : cpu_finalize => finalize
        call gpu_finalize()
        call cpu_finalize()
    end subroutine finalize

    subroutine gpu_init()
    end subroutine gpu_init

    subroutine gpu_finalize()
    end subroutine gpu_finalize

end module xray_atomic_scatter_factor_impl_gpu_module
