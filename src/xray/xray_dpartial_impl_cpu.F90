module xray_dpartial_impl_cpu_module
  
  use xray_contracts_module
  use xray_dpartial_data_module
  use xray_pure_utils, only : real_kind
  
  implicit none
  private
  
  public :: calc_partial_d_target_d_frac
  public :: finalize
  public :: init

contains
  
  function calc_partial_d_target_d_frac(frac, d_target_d_abs_Fcalc) result(d_target_d_frac)
    use xray_atomic_scatter_factor_module, only : atomic_scatter_factor
    use xray_pure_utils, only : PI
    use constants_xray, only : omp_num_threads
    use xray_scaling_data_module, only : k_scale
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind), intent(in) :: d_target_d_abs_Fcalc(:)
    real(real_kind) :: d_target_d_frac(3, size(frac, 2))
    real(real_kind) :: hkl_v(3)
    real(real_kind) :: f, phase
    real(real_kind) :: prefac(size(hkl,2))
    integer :: i
    integer :: ihkl
    
    call check_precondition(size(frac, 1) == 3)
    call check_precondition(size(frac, 2) == size(atom_b_factor))
    call check_precondition(size(frac, 2) == size(atom_scatter_type))
    call check_precondition(size(d_target_d_abs_Fcalc) == size(hkl, 2))
    
    call check_precondition(all(abs_Fcalc >= 0))
    call check_precondition(all(mSS4 <= 0))
    
    d_target_d_frac = 0
    prefac(:) = k_scale(:) * d_target_d_abs_Fcalc(:) / abs(Fcalc(:))
    
!$omp parallel do private(i,ihkl,hkl_v,phase,f) num_threads(omp_num_threads)
    do i = 1, size(frac, 2)
      do ihkl = 1, size(hkl, 2)
        
#if 0
        if (abs_Fcalc_ihkl < 1e-3) then
          ! Note: when Fcalc is approximately zero the phase is undefined,
          ! so no force can be determined even if the energy is high. (Similar
          ! to a linear bond angle.)
          cycle
        end if
#endif
        
        ! hkl-vector by 2pi
        hkl_v = hkl(:, ihkl) * 2 * PI
        
        phase = -sum(hkl_v * frac(:, i))
        ! f_n(s)          = atomic_scatter_factor(ihkl, atom_scatter_type(iatom))
        ! exp(-B_n*s^2/4) = exp(mSS4(ihkl) * atom_b_factor(iatom))
        f = prefac(ihkl) * atomic_scatter_factor(ihkl, atom_scatter_type(i)) &
            * exp(mSS4(ihkl) * atom_b_factor(i)) &
            * ( sin(phase)*Fcalc(ihkl)%re + cos(phase)*Fcalc(ihkl)%im )
            
        ! iatom's term of F^protein_calc (S1)
        
        d_target_d_frac(:, i) = d_target_d_frac(:, i) + hkl_v(:) * f
                
      end do
    end do
!$omp end parallel do
  
  end function calc_partial_d_target_d_frac
  
  
  subroutine init(hkl_, mss4_, Fcalc_, abs_Fcalc_, atom_b_factor_, atom_scatter_type_)
    implicit none
    integer, target, intent(in) :: hkl_(:, :)
    real(real_kind), target, intent(in) :: mSS4_(:)
    complex(real_kind), target, intent(in) :: Fcalc_(:)
    real(real_kind), target, intent(in) :: abs_Fcalc_(:)
    real(real_kind), intent(in) :: atom_b_factor_(:)
    integer, intent(in) :: atom_scatter_type_(:)
    
    call check_precondition(size(hkl_, 1) == 3)
    call check_precondition(size(mSS4_) == size(hkl_, 2))
    call check_precondition(size(abs_Fcalc_) == size(hkl_, 2))
    call check_precondition(size(Fcalc_) == size(hkl_, 2))
    
    call check_precondition(size(atom_scatter_type_) == size(atom_b_factor_))
    
    hkl => hkl_
    mSS4 => mss4_
    Fcalc => Fcalc_
    abs_Fcalc => abs_Fcalc_
    atom_b_factor = atom_b_factor_
    atom_scatter_type = atom_scatter_type_
  
  end subroutine init
  
  subroutine finalize()
    hkl => null()
    mSS4 => null()
    Fcalc => null()
    abs_Fcalc => null()
    deallocate(atom_b_factor)
    deallocate(atom_scatter_type)
  end subroutine finalize

end module xray_dpartial_impl_cpu_module
