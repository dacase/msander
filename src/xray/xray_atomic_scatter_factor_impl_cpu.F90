#include "../include/assert.fh"

module xray_atomic_scatter_factor_impl_cpu_module

    use xray_contracts_module
    use xray_pure_utils, only: real_kind
    
    implicit none

    private

    public :: init
    public :: finalize
    public :: atomic_scatter_factor

    real(real_kind), allocatable, save :: atomic_scatter_factor(:, :) ! (reflex number, scatter type)

    real(real_kind), parameter :: m_twopi = 2 * 3.1415926535897932384626433832795d0 ! FIXME: use constants module

contains

    subroutine init(mSS4, scatter_coefficients)
        implicit none
        real(real_kind), intent(in) :: mSS4(:)
        real(real_kind), intent(in) :: scatter_coefficients(:, :, :)
        
        integer :: ihkl, i
        integer :: n_hkl, n_scatter_types
        
        ASSERT(size(scatter_coefficients, 1) == 2)
        
        n_hkl = size(mSS4)
        n_scatter_types = size(scatter_coefficients, 3)

        allocate(atomic_scatter_factor(n_hkl, n_scatter_types))
        
        do ihkl = 1, n_hkl
            do i = 1, n_scatter_types
                atomic_scatter_factor(ihkl, i) = &
                        atom_scatter_factor_mss4(scatter_coefficients(:, :, i), mSS4(ihkl))
            end do
        enddo
    end subroutine init

    subroutine finalize()
        if (allocated(atomic_scatter_factor)) deallocate(atomic_scatter_factor)
    end subroutine finalize

    ! Private

    function atom_scatter_factor_mss4(coeffs, mss4) result(result)
        implicit none
        real(real_kind), intent(in) :: coeffs(:, :)
        real(real_kind), intent(in) :: mss4
        real(real_kind) :: result
        
        integer :: n_scatter_coeffs
        
        ASSERT(size(coeffs, 1) == 2)
        
        n_scatter_coeffs = size(coeffs, 2)

        result = coeffs(1, n_scatter_coeffs) + &
            sum(coeffs(1, 1:n_scatter_coeffs - 1) &
            * exp(mss4 * coeffs(2, 1:n_scatter_coeffs - 1)))
    end function atom_scatter_factor_mss4

end module xray_atomic_scatter_factor_impl_cpu_module
