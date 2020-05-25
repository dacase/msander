! <compile=optimized>
#include "../include/assert.fh"
module xray_fourier_module
! This module has the fourier (non-FFT) routines for calculating X-ray forces.
! Fourier (reciprocal) space values are stored as list of H,K,L values rather
! than a 3D array. Normally, only H,K,L indices with an observed value are
! saved.
!
! These routines all pass F90-style dimensioned arrays.
!
!  SUBROUTINES:
!
!  fourier_Fcalc     --   Caclulate structure factors using a direct
!                         method.
!
!  get_residual      --   Compute standard r_work/r_free statistics
!
!  scale_Fcalc       --   Carry out scaling (just anisotropic for now)
!
!  get_solvent_contribution -- bulk_solvent mask or other models
!
!  fourier_dTarget_dXYZBQ  --  Calculate the derivative of the xray restraint
!                         energy with respect to coordinate, B, or occupancy.
!                         Uses chain rule to combine dTarget/dF (passed
!                         from that routine in array dF) and dF/dXYZ or
!                         dF/dB or dF/dQ (computed in this routine).
!
!  dTargetLS_dF       --  Calculate a structure-factor restraint force
!                         from the scalar difference of |Fobs| and |Fcalc|.
!                         This uses a simple least-squares target function.
!
!  dTargetV_dF        --  Calculate a structure-factor restraint force
!                         from the vector (complex) difference of Fobs 
!                         and Fcalc.  (For use when phases are available,
!                         as in cryoEM.
!
!  dTargetML_dF       --  Use the phenix maximum-likelihood target function.
!
! FUNCTIONS:
!
! atom_scatter_factor_mss4 --  Calculate the atomic scatter (f) at a specific 
!                              resolution, given the Gaussian coefficients, 
!                              and a modified resolution, defined as -S*S/4.0, 
!                              where S is the reciprocal resolution. (mss4 
!                              stands for "Minus S Squared over 4") See comments
!                              at the beginning of fourier_Fcalc()

   use xray_globals_module
   use constants, only: M_TWOPI => TWOPI
   implicit none
#ifdef MPI
#include "parallel.h"
#else
      integer :: mytaskid = 0
#endif
      real(real_kind), allocatable, save  :: atomic_scatter_factor(:,:)

   !-------------------------------------------------------------------
contains

   !-------------------------------------------------------------------
   ! Caution: Some literature uses S to represent S^2
   !
   ! (for d*, see p. 93 of Glusker, Lewis, Rossi)
   ! S == d* = sqrt(sum(HKL * orth_to_frac)^2) = sqrt(-4*mSS4)
   ! mSS4 = -S*S/4
   ! mSS4 is more relevant to the formulas used than S.

   subroutine fourier_Fcalc(num_hkl,hkl,Fcalc,mSS4, &
      num_atoms,xyz,tempFactor,scatter_type_index, occupancy)

      implicit none
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl(3,num_hkl)
      complex(real_kind), intent(out) :: Fcalc(num_hkl)
      real(real_kind), intent(in) :: mSS4(num_hkl)
      integer, intent(in) :: num_atoms
      real(real_kind), intent(in), target :: xyz(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      ! index into the scatter - type coeffs table for each atom:
      integer, intent(in), target :: scatter_type_index(num_atoms) 
      real(real_kind), intent(in), target, optional :: occupancy(num_atoms)

      ! locals
      integer :: ihkl, i, ier, ierr
      ! automatic
      real(real_kind) :: f(num_atoms), angle(num_atoms)
      double precision :: time0, time1
      logical, save :: first=.true.

      call wallclock( time0 )

      if( first ) then
         ! set up reflection partitioning for MPI
#ifdef MPI
         ihkl1 = mytaskid*num_hkl/numtasks + 1
         ihkl2 = (mytaskid + 1) * num_hkl/numtasks
         if(mytaskid == numtasks - 1) ihkl2 = num_hkl
#else
         ihkl1 = 1
         ihkl2 = num_hkl
#endif

         allocate(atomic_scatter_factor(num_hkl,num_scatter_types), stat=ierr)
         REQUIRE(ierr==0)
         do ihkl=ihkl1,ihkl2
            do i = 1,num_scatter_types
               atomic_scatter_factor(ihkl,i) = &
               atom_scatter_factor_mss4(scatter_coefficients(:,:,i),mSS4(ihkl))
            end do
         enddo
         first = .false. 
      endif

#ifdef MPI
      Fcalc(:) = 0._rk_   ! needed since we will do an allreduce later
#endif

      ! special kludge to just get bulk_solvent factors:
      ! Fcalc(:) = 0._rk_
      ! return

!$omp parallel do private(ihkl,i,f,angle)  
      do ihkl = ihkl1, ihkl2

         ! Fhkl = SUM( fj * exp(2 * M_PI * i * (h * xj + k * yj + l * zj)) ),
         !      j = 1,num_selected_atoms
         ! where:
         !    The sum is versus j, over all selected atoms
         !    fj is the atomic scatter for atom j:   atomic_scatter_factor(j)
         !    h,k,l are from the hkl list for index ihkl:   hkl(1:3,ihkl)
         !    x,y,z are coordinates for atom j:   xyz(1:3,j)
         !        xyz(:) may be a reduced list.
         !
         ! Rather than using a complex exponential where the real part is 
         ! always zero, this is optimized to calculate sin and cosine parts,
         ! then convert to a complex number
         ! after the A and B components are summed over all selected atoms.
         ! This can be written as:
         !
         ! Ahkl = SUM( fj * cos(2 * M_PI * (h * xj + k * yj + l * zj)) ), 
         ! Bhkl = SUM( fj * sin(2 * M_PI * (h * xj + k * yj + l * zj)) ),
         !    j = 1,num_selected_atoms

         f(:) = exp( mSS4(ihkl) * tempFactor(:) ) &
              * atomic_scatter_factor(ihkl,scatter_type_index(:))
         if (present(occupancy)) then
            f(:) = f(:)*occupancy(:)
         endif
         angle(:) = matmul(M_TWOPI * hkl(1:3,ihkl),xyz(1:3,:))

         Fcalc(ihkl) = cmplx( sum(f(:) * cos(angle(:))), &
              sum(f(:) * sin(angle(:))), rk_ )

      end do
!$omp end parallel do
      call wallclock( time1 )
      ! write(6,'(a,f8.3)') '| ihkl loop time: ', time1 - time0

   end subroutine fourier_Fcalc

   ! -------------------------------------------------------------------------
   ! Combine dFcalc with respect to dXYZ, dB and/or dQ, 
   !    with dTarget/dFcalc (passed in as dF) to get dTarget/DXYZ,
   !    dTarget/dB or dTarget/dQ.

   ! Note the here XYZ is in the fractional coordinate system; conversion
   ! to Cartesian coordinates is done by the calling program in 
   ! xray_get_derivative() in xray_interface.F90

   ! Note Q is short for occupancy; B is short for tempFactor (aka B - factor).
   ! Small molecule software uses a U parameter in place of B - factor.
   ! U has a physical meaning, but B - factor is a more natural fit to Fouriers.

   !     isotropic B - factor = 8 * pi ** 2 * isotropic - U
   !     isotropic U = [U(1,1) + U(2,2) + U(3,3)]/3.0

   subroutine fourier_dTarget_dXYZBQ( num_hkl,hkl,dF,mSS4, & 
            num_atoms,xyz,tempFactor,scatter_type_index, &
            occupancy, dxyz,d_occupancy,d_tempFactor )
      implicit none

      ! reciprocal space arrays:
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl(3,num_hkl)
      complex(real_kind), intent(in) :: dF(num_hkl)
      real(real_kind), intent(in) :: mSS4(num_hkl)

      ! coordinate arrays:
      integer, intent(in) :: num_atoms
      real(real_kind), intent(in) :: xyz(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      integer, intent(in) :: scatter_type_index(num_atoms)
      real(real_kind), intent(in), optional :: occupancy(num_atoms)

      ! output derivatives:
      real(real_kind), intent(out), optional :: dxyz(3,num_atoms)
      real(real_kind), intent(out), optional :: d_occupancy(num_atoms)
      real(real_kind), intent(out), optional :: d_tempFactor(num_atoms)
      !real(real_kind), intent(out), optional :: d_aniso_Bij(6,num_atoms)

      ! locals
      integer :: ihkl, iatom, i
      real(real_kind) :: dhkl(3)
      complex(real_kind) :: f
      real(real_kind) :: phase
      double precision time0, time1
      logical, save :: first=.true.

      if (present(dxyz)) dxyz(:,:) = 0._rk_
      if (present(d_tempFactor)) d_tempFactor(:) = 0._rk_

      call wallclock( time0 )

      ! TODO: does it hurt to have if statements inside the double loop?
#ifdef MPI
      REFLECTION: do ihkl = ihkl1,ihkl2
         ATOM: do iatom = 1,num_atoms
#else
!$omp parallel do private(ihkl,dhkl,iatom,phase,f) 
      ATOM: do iatom = 1,num_atoms
         REFLECTION: do ihkl = ihkl1,ihkl2
#endif

            dhkl = hkl(:,ihkl) * M_TWOPI ! * symmop...

            phase = sum( dhkl * xyz(:,iatom) )
            f = atomic_scatter_factor(ihkl, scatter_type_index(iatom)) &
                  * exp(mSS4(ihkl) * tempFactor(iatom)) 
   
            f = f * cmplx(cos(phase),sin(phase), rk_)

#if 0
            if (present(d_occupancy)) then
               d_occupancy(iatom) = d_occupancy(iatom) + &
                 real(f) * real(dF(ihkl)) + aimag(f) * aimag(dF(ihkl))
            end if

            if (present(occupancy)) then
               f = f * occupancy(iatom)
            end if

            if (present(d_tempFactor)) then
               d_tempFactor(iatom) = d_tempFactor(iatom) &
                  + ( real(f) * real(dF(ihkl)) + aimag(f) * aimag(dF(ihkl)) ) &
                  * mSS4(ihkl)
            end if
#endif

            ! if (present(dxyz)) then
               dxyz(:,iatom) = dxyz(:,iatom) - dhkl(:) * &
                   ( aimag(f) * real(dF(ihkl)) - real(f) * aimag(dF(ihkl)) )
            ! end if

#ifdef MPI
         end do ATOM
      end do REFLECTION
#else
         end do REFLECTION
      end do ATOM
!$omp end parallel do
#endif
      call wallclock( time1 )
      ! write(6,'(a,f8.3)') '| dhkl loop time: ', time1 - time0
      return

   end subroutine fourier_dTarget_dXYZBQ

   ! -------------------------------------------------------------------------
   ! R - factor = sum(abs(Fobs - Fcalc)) / sum(Fobs)
   ! -------------------------------------------------------------------------
   subroutine get_residual (num_hkl,abs_Fobs,abs_Fcalc,residual,selected)
      implicit none
      integer, intent(in) :: num_hkl
      real(real_kind), intent(in) :: abs_Fobs(num_hkl), abs_Fcalc(num_hkl)
      real(real_kind), intent(out) :: residual
      integer, intent(in), optional :: selected(num_hkl)
      real(real_kind) :: denom
   
      if (present(selected)) then
         denom = sum(abs_Fobs, selected/=0)
         if( denom > 0._rk_ ) then
            residual = sum( abs(abs_Fobs - abs_Fcalc), selected/=0) / denom
         else
            residual = 0._rk_
         endif
      else
         residual = sum (abs( abs_Fobs - abs_Fcalc) ) / sum(abs_Fobs)
      end if
      return
   end subroutine get_residual

   subroutine get_solvent_contribution(nstep,crd)
      use bulk_solvent_mod, only : f_mask, k_mask, grid_bulk_solvent, &
          shrink_bulk_solvent, fft_bs_mask, mask_bs_grid_t_c, &
          hkl_indexing_bs_mask, mask_cell_params, mask_grid_size
      implicit none
      integer, intent(in) :: nstep
      real(real_kind), intent(in) :: crd(3*num_atoms)

      integer :: i

      if (has_f_solvent == 0 .and. mod(nstep, mask_update_frequency) == 0) then
         call grid_bulk_solvent(num_atoms, crd)
         call shrink_bulk_solvent()
         call fft_bs_mask()

         do i=1,num_hkl
            f_mask(i) = conjg(mask_bs_grid_t_c(hkl_indexing_bs_mask(i) + 1)) * &
                        mask_cell_params(16) / mask_grid_size(4)
         end do
         if (mytaskid == 0 ) &
           write(6,'(a,5e14.6)') '| updating bulk solvent: ', &
               Fcalc(1),k_mask(1),f_mask(1)
      endif
      if( bulk_solvent_model == 'simple') &
        Fcalc(:) = Fcalc(:) + k_mask(:)*f_mask(:)
      return

   end subroutine get_solvent_contribution

   subroutine scale_Fcalc(nstep)
      use ml_mod, only : b_vector_base, NRF_work, NRF_work_sq, &
           h_sq, k_sq, l_sq, hk, kl, hl, MUcryst_inv
      implicit none
      integer, intent(in) :: nstep

      real(real_kind) :: sum_fo_fc, sum_fc_fc
      real(real_kind) :: b(7), Uaniso(7)

      if (mod(nstep, scale_update_frequency) == 0) then
#if 1
         ! isotropic scaling for Fcalc, with same general notation:
         NRF_work_sq = NRF_work * NRF_work
         b_vector_base = log(abs_Fobs(1:NRF_work) / abs(Fcalc(1:NRF_work))) &
                      / NRF_work_sq
         b(1) = sum(b_vector_base)
         k_scale(:) = exp(b(1))
         if(mytaskid==0) &
           write(6,'(a,f10.5)') '| updating   isotropic scaling: ', k_scale(1)
#else
         ! isotropic scaling from dtargetLS:

         sum_fo_fc = sum(abs_Fobs(1:NRF_work) * abs(Fcalc(1:NRF_work)))
         sum_fc_fc = sum(abs(Fcalc(1:NRF_work)) ** 2)
         Fcalc_scale = sum_fo_fc / sum_fc_fc
         if (mytaskid == 0 ) &
            write(6,'(a,f12.5)') '| simple   isotropic scaling: ',Fcalc_scale
         k_scale(:) = Fcalc_scale

         ! anisotropic scaling for Fcalc:
         NRF_work_sq = NRF_work * NRF_work
         b_vector_base = log(abs_Fobs(1:NRF_work) / abs(Fcalc(1:NRF_work))) &
                      / NRF_work_sq
         b(1) = sum(b_vector_base)
         b(2) = sum(b_vector_base * h_sq(1:NRF_work))
         b(3) = sum(b_vector_base * k_sq(1:NRF_work))
         b(4) = sum(b_vector_base * l_sq(1:NRF_work))
         b(5) = sum(b_vector_base * hk(1:NRF_work))
         b(6) = sum(b_vector_base * hl(1:NRF_work))
         b(7) = sum(b_vector_base * kl(1:NRF_work))

         Uaniso = matmul(MUcryst_inv, b)
         k_scale = exp(Uaniso(1) + Uaniso(2)*h_sq + Uaniso(3)*k_sq + &
             Uaniso(4)*l_sq + Uaniso(5)*hk + Uaniso(6)*hl + Uaniso(7)*kl)
         if (mytaskid == 0 ) then
           write(6,'(a)') '| updating anisotropic scaling: '
           write(6,'(a,7f10.5)') '|     ', Uaniso
           write(6,'(a,2f10.5)')  '|     ', exp(b(1))
           write(6,'(a,7f10.5)')  '|     ', k_scale(1:7)
         endif
#endif
      endif

      Fcalc = Fcalc * k_scale  !either using newly computed or stored k_scale
      return
   end subroutine scale_Fcalc

   ! -------------------------------------------------------------------------
   ! This routine computes the force gradient on Fcalc as a harmonic
   ! restraint on the magnitudes of Fobs and Fcalc
   ! -------------------------------------------------------------------------
   subroutine dTargetLS_dF(crd,weight,selected,deriv,xray_energy)
      implicit none
      real(real_kind), intent(in) :: crd(3*num_atoms)
      real(real_kind), intent(in), optional :: weight (num_hkl)
      integer, intent(in), optional :: selected(num_hkl)
      complex(real_kind), intent(out), optional :: deriv(num_hkl)
      real(real_kind), intent(out), optional :: xray_energy

      real(real_kind) :: sum_fo_fc, sum_fo_fo, sum_fc_fc
      real(real_kind) :: abs_Fcalc(num_hkl)
      real(real_kind), parameter :: F_EPSILON = 1.0e-20_rk_
      integer, save :: nstep=0
      integer :: i
      logical :: iso_scale=.true. ! no input varible yet, just change the code

      if( bulk_solvent_model .eq. 'simple' ) &
         call get_solvent_contribution(nstep, crd)

      if( iso_scale ) then
         abs_Fcalc(:) = abs(Fcalc(:))
         if( mod(nstep,scale_update_frequency) == 0 ) then
            if (present(selected)) then
               sum_fo_fc = sum(abs_Fobs * abs_Fcalc,selected/=0)
               sum_fo_fo = sum(abs_Fobs ** 2,selected/=0)
               sum_fc_fc = sum(abs_Fcalc ** 2,selected/=0)
            else
               sum_fo_fc = sum(abs_Fobs * abs_Fcalc)
               sum_fo_fo = sum(abs_Fobs ** 2)
               sum_fc_fc = sum(abs_Fcalc ** 2)
            end if
            Fcalc_scale = sum_fo_fc / sum_fc_fc
            norm_scale = 1.0_rk_  / sum_fo_fo
            if (mytaskid == 0 ) &
               write(6,'(a,f12.5,e12.5)') '| updating isotropic scaling: ', &
                   Fcalc_scale,norm_scale
         endif
         Fcalc(:) = Fcalc_scale * Fcalc(:)
      else
         Fcalc_scale = 1._rk_
      endif

      nstep = nstep + 1
      abs_Fcalc(:) = abs(Fcalc(:))

      ! Note: when Fcalc is approximately zero the phase is undefined, 
      ! so no force can be determined even if the energy is high. (Similar 
      ! to a linear bond angle.)

      ! We get the phase from Fcalc, so it needs to be divided by abs(Fcalc)
      ! to become a unit vector.
      !
      ! deriv = Fcalc/abs(Fcalc) * K * ( Fobs - abs(Fcalc) )
      !      ....plus the terms arising from differentiating Fcalc_scale with
      !      respect to Fcalc, which are ignored here, since Fcalc_scale is
      !      generally not updated on every step.

      if (present(deriv)) then
         deriv(:) = 0._rk_  ! so no force for things unselected here
         if (present(selected)) then
            where (selected/=0 .and. abs_Fcalc > 1.d-3)
               deriv(:) = - 2.0_rk_ * f_weight(:) * Fcalc(:) * norm_scale * &
                  ( abs_Fobs(:) - abs_Fcalc(:) ) *  &
                  ( Fcalc_scale/abs_Fcalc(:) )
            end where
         else ! no selected
            ! where( abs_Fcalc > 1.d-3 )
               deriv(:) = - 2.0_rk_ * f_weight(:) * Fcalc(:) * norm_scale * &
                  ( abs_Fobs(:) - abs_Fcalc(:) ) *  &
                  ( Fcalc_scale/abs_Fcalc(:) ) 
            ! end where
         end if
      end if

      if (present(selected)) then
         xray_energy = norm_scale * &
              sum(f_weight*(abs_Fobs - abs_Fcalc)**2, selected/=0)
      else
         xray_energy = norm_scale * sum(f_weight*(abs_Fobs - abs_Fcalc)**2)
      end if

   end subroutine dTargetLS_dF

   ! This routine computes the force gradient on Fcalc as a harmonic
   ! restraint on the vector (complex) difference between Fcalc and
   ! Fobs

   subroutine dTargetV_dF(crd,deriv,residual,xray_energy)
      implicit none
      real(real_kind), intent(in) :: crd(3*num_atoms)
      complex(real_kind), intent(out) :: deriv(:)
      real(real_kind), intent(out) :: residual
      real(real_kind), intent(out) :: xray_energy

      complex(real_kind) :: vecdif(num_hkl)
      integer, save :: nstep=0

      ! step 1: get fcalc, including solvent mask contribution:
      !  (atomic part already done in fourier_Fcalc)
      ! assume for now that cryoEM won't need a solvent contribution
      ! call get_solvent_contribution(nstep, crd)

      ! step 2: isotropic scaling:
      if (scale_update_frequency > 0 ) then
         if (mod(nstep,scale_update_frequency) == 0) then
            Fcalc_scale = sum( real(Fobs*conjg(Fcalc)) ) / sum(abs(Fcalc)**2)
            norm_scale = 1.0_rk_ / sum(abs_Fobs ** 2)
            if (mytaskid == 0 ) &
              write(6,'(a,f12.5,e12.5)') '| updating isotropic scaling: ',&
                   Fcalc_scale, norm_scale
         endif
      else
         if( nstep==0 ) norm_scale = 1.0_rk_ / sum(abs_Fobs ** 2)
         Fcalc_scale = 1.0_rk_
      endif

      Fcalc(:) = Fcalc_scale * Fcalc(:)
      nstep = nstep + 1

      vecdif(:) = Fobs(:) - Fcalc(:)
      xray_energy = norm_scale * sum( vecdif(:)*conjg(vecdif(:)) )
      deriv(:) = - norm_scale * 2._rk_ * Fcalc_scale * vecdif(:)
      residual = sum (abs(vecdif)) / sum(abs_Fobs)

   end subroutine dTargetV_dF

   ! This routine computes the force gradient on Fcalc, using the
   ! phenix maximum likelihood function

   subroutine dTargetML_dF(crd,deriv,xray_energy)
      use ml_mod, only : b_vector_base, &
           alpha_array, beta_array, delta_array, NRF_work, &
           i1_over_i0, ln_of_i0, estimate_alpha_beta, NRF, &
           init_scales, optimize_k_scale_k_mask, k_iso, k_iso_exp, &
           k_aniso
      use bulk_solvent_mod, only: f_mask, k_mask
      implicit none
      real(real_kind), intent(in) :: crd(3*num_atoms)
      complex(real_kind), intent(out) :: deriv(num_hkl)
      real(real_kind), intent(out) :: xray_energy

      real(real_kind) :: abs_Fcalc(num_hkl)
      integer, save :: nstep = 0
      integer :: i
      double precision :: eterm1, eterm2, x

      if( bulk_solvent_model .eq. 'simple' ) then
         ! step 1: get fcalc, including solvent mask contribution:
         !  (atomic part already done in fourier_Fcalc, and passed in here.)
         call get_solvent_contribution(nstep, crd)
         call scale_Fcalc( nstep )

      else if( bulk_solvent_model .eq. 'opt' ) then
         if (mod(nstep, mask_update_frequency) == 0) then
           call get_solvent_contribution(nstep, crd)
           call init_scales()
           call optimize_k_scale_k_mask()
           k_scale = k_iso * k_iso_exp * k_aniso
         endif
         Fcalc = k_scale * (Fcalc + k_mask * f_mask)

      else if( bulk_solvent_model .eq. 'none' ) then
         call scale_Fcalc( nstep )
         
      else
         write(6,*) 'Bad value for bulk_solvent_model: ', trim(bulk_solvent_model)
      endif
      abs_Fcalc(:) = abs(Fcalc(:))

      ! step 3: get ml parameters

      if (mod(nstep, ml_update_frequency) == 0) then
        if (mytaskid == 0 ) &
           write(6,'(a)') '| updating alpha and beta'
        call estimate_alpha_beta()
        delta_array = alpha_array / beta_array
      endif

      ! ML target function
#if 0
      xray_energy = 0._rk_
      do i=1,NRF_work
         x = 2.d0*delta_array(i)*abs_Fcalc(i)*abs_Fobs(i)
         eterm1 = alpha_array(i)*delta_array(i)*abs_Fcalc(i)*abs_Fcalc(i)
         eterm2 = ln_of_i0(x)
         xray_energy = xray_energy + eterm1 - eterm2
      enddo
#else
      ! Oleg new version, Dec. 2019:
      xray_energy = sum(abs_Fobs(1:NRF_work)**2 / beta_array(1:NRF_work) - &
            log(2*abs_Fobs(1:NRF_work)/beta_array(1:NRF_work)) +  &
            alpha_array(1:NRF_work) * delta_array(1:NRF_work) * &
            abs_Fcalc(1:NRF_work) * abs_Fcalc(1:NRF_work) - &
            ln_of_i0(2 * delta_array(1:NRF_work) * abs_Fcalc(1:NRF_work) * &
            abs_Fobs(1:NRF_work)))
#endif

      nstep = nstep + 1

      ! step 4: put dTargetML/dF into deriv(:)
      do i=1,NRF_work
         x = 2.d0*delta_array(i)*abs_Fcalc(i)*abs_Fobs(i)
         deriv(i) = k_scale(i)*Fcalc(i)*2.d0*delta_array(i) &
            * ( alpha_array(i)*abs_Fcalc(i) - i1_over_i0(x)*abs_Fobs(i) ) &
            / abs_Fcalc(i)
      end do
      do i=NRF_work+1, NRF
         deriv(i) = 0.d0
      end do
      return
   end subroutine dTargetML_dF

   function atom_scatter_factor_mss4(coeffs,mss4) result(sfac)
      real(real_kind) :: sfac
      real(real_kind), intent(in) :: coeffs(2,scatter_ncoeffs), mss4
      sfac = coeffs(1,scatter_ncoeffs) + &
             sum( coeffs(1,1:scatter_ncoeffs-1) &
             * exp(mss4*coeffs(2,1:scatter_ncoeffs-1)))
   end function atom_scatter_factor_mss4

   !dac addition:
   subroutine get_mss4(num_hkl,hkl_index,mSS4)
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl_index(3,num_hkl)
      real(real_kind), intent(inout) :: mss4(:)
      integer :: ihkl, h,k,l
      real(real_kind) :: a,b,c,alpha,beta,gamma,V,S2
      real(real_kind) :: sina,cosa,sinb,cosb,sing,cosg
      real(real_kind) :: astar,bstar,cstar,cosas,cosbs,cosgs

      a = unit_cell(1)
      b = unit_cell(2)
      c = unit_cell(3)
      alpha = 3.1415926536d0*unit_cell(4)/180.
      beta = 3.1415926536d0*unit_cell(5)/180.
      gamma = 3.1415926536d0*unit_cell(6)/180.
      sina = sin(alpha)
      cosa = cos(alpha)
      sinb = sin(beta)
      cosb = cos(beta)
      sing = sin(gamma)
      cosg = cos(gamma)

      V = a*b*c*sqrt( 1.d0 - cosa**2 - cosb**2 -cosg**2 &
          + 2.d0*cosa*cosb*cosg )
      astar = b*c*sina/V
      bstar = a*c*sinb/V
      cstar = b*a*sing/V
      cosas = (cosb*cosg - cosa)/(sinb*sing)
      cosbs = (cosa*cosg - cosb)/(sina*sing)
      cosgs = (cosb*cosa - cosg)/(sinb*sina)

      ! work from p. 93 of Glusker, Lewis, Rossi:

      do ihkl=1,num_hkl
        h = hkl_index(1,ihkl)
        k = hkl_index(2,ihkl)
        l = hkl_index(3,ihkl)
        S2 = (h*astar)**2 + (k*bstar)**2 + (l*cstar)**2  &
           + 2*k*l*bstar*cstar*cosas + 2*l*h*cstar*astar*cosbs  &
           + 2*h*k*astar*bstar*cosgs
        mSS4(ihkl) = -S2/4.
      end do
   end subroutine get_mss4

end module xray_fourier_module
