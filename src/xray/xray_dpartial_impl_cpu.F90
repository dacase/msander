#include "../include/assert.fh"

module xray_dpartial_impl_cpu_module
  
  use xray_contracts_module
  use xray_dpartial_data_module
  use xray_pure_utils, only : real_kind
  use constants_xray, only : xray_num_threads
  
  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: calc_partial_d_vls_d_frac
  public :: finalize
  public :: init

contains
  
  function calc_partial_d_target_d_frac(frac, f_scale, &
        d_target_d_abs_Fcalc) result(d_target_d_frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    use xray_pure_utils, only : PI
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: f_scale(:)
    real(real_kind), intent(in) :: d_target_d_abs_Fcalc(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))
    real(real_kind) :: hkl_v(3)
    real(real_kind) :: phase
    complex(real_kind) :: f
    integer :: i
    integer :: ihkl
    
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_b_factor))
    ASSERT(size(frac, 2) == size(atom_scatter_type))
    ASSERT(size(f_scale) == size(hkl, 2))
    ASSERT(size(d_target_d_abs_Fcalc) == size(hkl, 2))
    
    ASSERT(all(abs_Fcalc >= 0))
    ASSERT(all(mSS4 <= 0))
    
    d_target_d_frac = 0

!$omp parallel do private(i,ihkl,hkl_v,phase,f) num_threads(xray_num_threads)
    do i = 1, size(frac, 2)
      do ihkl = 1, size(hkl, 2)
        
        if (abs_Fcalc(ihkl) < 1e-3) then
          ! Note: when Fcalc is approximately zero the phase is undefined,
          ! so no force can be determined even if the energy is high. (Similar
          ! to a linear bond angle.)
          cycle
        end if
        
        ! hkl-vector by 2pi
        hkl_v = hkl(:, ihkl) * 2 * PI
        
        phase = -sum(hkl_v * frac(:, i))
        ! f_n(s)          = atomic_scatter_factor(ihkl, atom_scatter_type(iatom))
        ! exp(-B_n*s^2/4) = exp(mSS4(ihkl) * atom_b_factor(iatom))
        f = atomic_scatter_factor(ihkl, atom_scatter_type(i)) &
            * exp(mSS4(ihkl) * atom_b_factor(i))
        
        f = f * cmplx(cos(phase), sin(phase), real_kind)
        ! iatom's term of F^protein_calc (S1)
        
        d_target_d_frac(:, i) = d_target_d_frac(:, i) &
              + atom_occupancy(i) * f_scale(ihkl) * hkl_v(:) * &
                aimag(f * Fcalc(ihkl)) * &
                d_target_d_abs_Fcalc(ihkl) / abs_Fcalc(ihkl)
      end do
    end do
!$omp end parallel do
  
  end function calc_partial_d_target_d_frac
  
  function calc_partial_d_vls_d_frac(frac, f_scale) &
         result(d_target_d_frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    use xray_target_vector_least_squares_data_module, only: derivc
    use xray_pure_utils, only : PI
    use xray_interface2_data_module, only : n_work
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: f_scale(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))
    real(real_kind) :: hkl_v(3)
    real(real_kind) :: f, phase
    integer :: ihkl, i
    
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_b_factor))
    ASSERT(size(frac, 2) == size(atom_scatter_type))
    ASSERT(size(f_scale) == size(hkl, 2))
    
    ASSERT(all(abs_Fcalc >= 0))
    ASSERT(all(mSS4 <= 0))
    
    d_target_d_frac = 0

!$omp parallel do private(i,ihkl,hkl_v,phase,f) num_threads(xray_num_threads)
    do i = 1, size(frac, 2)
      do ihkl = 1, n_work

        if (abs_Fcalc(ihkl) < 1e-3) then
          ! Note: when Fcalc is approximately zero the phase is undefined,
          ! so no force can be determined even if the energy is high. (Similar
          ! to a linear bond angle.)
          cycle
        end if
        
        ! hkl-vector by 2pi
        hkl_v = hkl(:, ihkl) * 2 * PI
        
        phase = -sum(hkl_v * frac(:, i))
        f = atomic_scatter_factor(ihkl, atom_scatter_type(i)) &
            * exp(mSS4(ihkl) * atom_b_factor(i)) &
            * ( sin(phase) * real(derivc(ihkl)) + cos(phase) * aimag(derivc(ihkl)) )
        
        d_target_d_frac(:, i) = d_target_d_frac(:, i) &
            + atom_occupancy(i) * f_scale(ihkl) * hkl_v(:) * f

      end do
    end do
!$omp end parallel do
  
  end function calc_partial_d_vls_d_frac
  
  
  subroutine init(hkl_, mss4_, Fcalc_, abs_Fcalc_, atom_b_factor_,  &
        atom_occupancy_, atom_scatter_type_)
    implicit none
    integer, target, intent(in) :: hkl_(:, :)
    real(real_kind), target, intent(in) :: mSS4_(:)
    complex(real_kind), target, intent(in) :: Fcalc_(:)
    real(real_kind), target, intent(in) :: abs_Fcalc_(:)
    real(real_kind), intent(in) :: atom_b_factor_(:)
    real(real_kind), intent(in) :: atom_occupancy_(:)
    integer, intent(in) :: atom_scatter_type_(:)
    
    ASSERT(size(hkl_, 1) == 3)
    ASSERT(size(mSS4_) == size(hkl_, 2))
    ASSERT(size(abs_Fcalc_) == size(hkl_, 2))
    ASSERT(size(Fcalc_) == size(hkl_, 2))
    
    ASSERT(size(atom_scatter_type_) == size(atom_b_factor_))
    
    hkl => hkl_
    mSS4 => mss4_
    Fcalc => Fcalc_
    abs_Fcalc => abs_Fcalc_
    atom_b_factor = atom_b_factor_
    atom_occupancy = atom_occupancy_
    atom_scatter_type = atom_scatter_type_
  
  end subroutine init
  
  subroutine finalize()
    hkl => null()
    mSS4 => null()
    Fcalc => null()
    abs_Fcalc => null()
    !  dac: These are also deallocated in xray_interface2_data.F90
    !  deallocate(atom_b_factor)
    ! deallocate(atom_scatter_type)
  end subroutine finalize

end module xray_dpartial_impl_cpu_module
