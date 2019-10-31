!<compile=optimized>

#include "../include/dprec.fh"

!> Closure super class for 3D-RISM.  Closure sub-classes (i.e.,
!! actual closure implementations) are registered here.  Subroutine
!! calls then call the appropriate subroutine of the subclass
!! interface.
!!
!! This is an explicit implementation of class inheritance. See
!! V. K. Decyk, C. D. Norton, B. K. Szymanski.  How to express C++
!! concepts in Fortran 90. Scientific Programming. 6, 363-390 (1997).
!!
!! Some closure independent properties are calculated within this
!! class. Uvv is the site-site potential.
!!
!! In gerneral, this class is MPI aware only through the rism3d_grid
!! class (it knows about the total system size and its own piece of
!! it).  It does not know about processes and does not perform
!! reductions (all thermodynamic quantities are distributed and each
!! process has only the contribution of the local slab).
module rism3d_closure_c
  use rism3d_potential_c
  use rism3d_grid_c
  use rism3d_kh_c
  use rism3d_hnc_c
  use rism3d_psen_c
  use rism_report_c
  use safemem
  implicit none

  type rism3d_closure
     !> Currenly selected closure.
     character(len=4) :: type
     !> Kovalenko-Hirata closure.
     type(rism3d_kh), pointer :: kh => NULL()
     !> Hypernetted chain equation closure.
     type(rism3d_hnc), pointer :: hnc => NULL()
     !> Partial series expansion of order n (PSE-n) closure.
     type(rism3d_psen), pointer :: psen => NULL()
     !> Electric potential object.
     type(rism3d_potential), pointer :: potential => NULL()
     !> Box grid, stored in potential object.
     type(rism3d_grid), pointer :: grid => NULL()
     !> Solvent, stored in potential object.
     type(rism3d_solvent), pointer :: solvent => NULL()
     !> Solute, stored in potential object.
     type(rism3d_solute), pointer :: solute => NULL()
  end type rism3d_closure

  private lennardJonesForcePeriodic, particleMeshEwaldForce
  
contains

  !> Creates a new closure object of the requested type.
  !! @param[in,out] this the closure object
  !! @param[in] pot rism3d_potential object.  Must be initialized.
  !! @param[in] type one of 'KH', 'HNC', 'PSEn', where 'n' is the
  !!          order of the PSE-n closure
  subroutine rism3d_closure_new(this,type,pot)
    use rism_util, only : caseup
    implicit none
    type(rism3d_closure), intent(inout) :: this
    type(rism3d_potential), target, intent(in)  :: pot
    character(len=*), intent(in) :: type
    integer :: order, iostat
    this%potential => pot
    this%grid => this%potential%grid
    this%solvent => this%potential%solvent
    this%solute => this%potential%solute
    this%type = trim(type)
    call caseup(this%type)
    if (this%type .eq. "KH") then
       allocate(this%kh)
       call rism3d_kh_new(this%kh,this%potential)
    else if (index(this%type,"PSE") == 1) then
       read(this%type(4:),*, iostat=iostat) order
       if (iostat/=0)&
          call rism_report_error("'"//trim(this%type)//"' not a valid closure")
       allocate(this%psen)
       call rism3d_psen_new(this%psen,this%potential,order)
    else if (trim(this%type) .eq. "HNC") then
       allocate(this%hnc)
       call rism3d_hnc_new(this%hnc,this%potential)
    else
       call rism_report_error("'"//trim(this%type)//"' not a valid closure")
    end if
  end subroutine rism3d_closure_new

  !> Returns a identifier string for the closure type.
  !! @param[in] this The closure object.
  function rism3d_closure_type(this) result(type)
    implicit none
    type(rism3d_closure), intent(in) :: this
    character(len=4) :: type
    type=this%type
  end function rism3d_closure_type


  !> Calculates Guv from Uuv, Huv, and Cuv using the associated
  !! closure.
  !! @param[in] this The closure object.
  !! @param[in] guv Site-site pair correlation function.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  subroutine rism3d_closure_guv(this,guv, huv, cuv)
    implicit none
    type(rism3d_closure), intent(inout) :: this
    _REAL_, intent(out) :: guv(:,:)
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    if (associated(this%kh)) then
       call rism3d_kh_guv(this%kh,guv,huv,cuv)
    else if (associated(this%psen)) then
       call rism3d_psen_guv(this%psen,guv,huv,cuv)
    else if (associated(this%hnc)) then
       call rism3d_hnc_guv(this%hnc,guv,huv,cuv)
    end if
  end subroutine rism3d_closure_guv

  !> Calculates the excess chemical potential in kT for each site.
  !! @param[in,out] this The closure object.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  function rism3d_closure_excessChemicalPotential(this, huv, cuv) &
           result(excessChemicalPotential)
    implicit none
    type(rism3d_closure), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:), cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%solvent%numAtomTypes)
    if (associated(this%kh)) then
       excessChemicalPotential = rism3d_kh_excessChemicalPotential(this%kh,huv,cuv)
    else if (associated(this%psen)) then
       excessChemicalPotential = rism3d_psen_excessChemicalPotential(this%psen,huv,cuv)
    else if (associated(this%hnc)) then
       excessChemicalPotential = rism3d_hnc_excessChemicalPotential(this%hnc,huv,cuv)
    end if
  end function rism3d_closure_excessChemicalPotential

  !!Calculate the total solvation interaction energy: de = density sum g*u for
  !!each solvent site.  I.e., the direct intection potential energy of
  !!solute and solvent and not the total solvation energy (see solvationEnergy).
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!OUT:
  !!    the contribution of each solvent site to the total
  !!    solute-solvent potential energy [kT]
  function rism3d_closure_solvPotEne (this,guv) result(ene)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer ::  igk, iv, ix,iy,iz
    _REAL_ ::  ene(this%solvent%numAtomTypes)
    ene = 0.d0
#ifdef RISM_DEBUG
    write(6,*)"EXENER"
#endif /*RISM_DEBUG*/

    !!FIX - can we avoid the loops and use BLAS array operations?
    do iv=1,this%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
                ene(iv) = ene(iv) &
                   + rism3d_closure_solvPotEne_ijk (this,guv,(/ix,iy,iz/),iv)
             end do
          end do
       end do
    end do
    !!endfix
    ene = ene*this%grid%voxelVolume
  end function rism3d_closure_solvPotEne


  !!Calculate the total solvation interaction energy: de = density sum g*u
  !!for the request grid point and solvent site.  I.e., the direct
  !!intection potential energy of solute and solvent and not the total
  !!solvation energy (see solvationEnergy).
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    guv  :: site-site pair distribution function
  !!    ijk  : 3d-grid index
  !!OUT:
  !!    the contribution of each solvent site to the total
  !!    solute-solvent potential energy [kT]
  function rism3d_closure_solvPotEne_ijk (this,guv,ijk,iv) result(ene)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    integer, intent(in) :: ijk(3), iv
    integer ::  igk, ix,iy,iz
    _REAL_ ::  ene
    ix=ijk(1)
    iy=ijk(2)
    iz=ijk(3)
#ifdef MPI
    igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) &
             + (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
    igk = ix + (iy-1)*this%grid%localDimsR(1) &
             + (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif
    ene = guv(igk,iv) * this%potential%uuv(ix,iy,iz,iv) &
                      * this%potential%solvent%density(iv)
  end function rism3d_closure_solvPotEne_ijk

  !> Calculates the partial molar volume.
  !! @param[in] this rism3d object with calculated solution.
  !! @param[in] cuv Site-site direct correlation function.
  !! @return The calculated partial molar volume. [A^3]
  function rism3d_closure_partialMolarVolume(this, cuv) result(partialMolarVolume)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_ :: partialMolarVolume
    _REAL_ :: rCuv(this%solvent%numAtomTypes),rCuv_dT(this%solvent%numAtomTypes)
    integer ::  ix, iy, iz, iv, ierr
    rcuv = rism3d_closure_DCFintegral(this,cuv)
    partialMolarVolume = -this%solvent%xikt &
         * sum(rcuv* this%solvent%density)
    if (this%grid%mpirank == 0) then
       partialMolarVolume = partialMolarVolume +this%solvent%xikt
    end if
  end function rism3d_closure_partialMolarVolume

  !> Calculate the excess number of particles about the solute
  !! compared to the bulk solvation.  No attempt is made to account
  !! for excluded volume.  This is also-known-as a molar preferential
  !! interaction parameter.
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] guv site-site pair distribution function.
  !! @return The number of excess particles for each solvent site.
  function rism3d_closure_excessParticles (this, guv) result(num)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: num(this%solvent%numAtomTypes)
    num = this%solvent%density * rism3d_closure_kirkwoodBuff(this, guv)
  end function rism3d_closure_excessParticles

  !> Calculate the Kirkwood-Buff integral for the solute. This is the
  !! all space integral of h_{uv}.
  !!
  !! J. G. Kirkwood; F. P. Buff. J. Chem. Phys. 1951, 19, 774-777
  !! @param[in] this rism3d object with computed solution.
  !! @param[in] guv Site-site pair distribution function.
  !! @return Kirkwood-Buff integral for each solvent site.
  function rism3d_closure_kirkwoodBuff (this, guv) result(H)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: guv(:,:)
    _REAL_ :: H(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    
    H = 0
    do iv = 1, this%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
#ifdef MPI
                igk = ix + (iy - 1) * (this%grid%localDimsR(1) + 2) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                                    * (this%grid%localDimsR(1) + 2)
#else
                igk = ix + (iy - 1) * this%grid%localDimsR(1) &
                         + (iz - 1) * this%grid%localDimsR(2) &
                         * this%grid%localDimsR(1)
#endif
                H(iv) = H(iv) + (guv(igk, iv) - 1d0)
             end do
          end do
       end do
    end do
    H = H * this%grid%voxelVolume
  end function rism3d_closure_kirkwoodBuff

  !!Calculates the direct correlation function integral for the solute. This is the
  !!all space integral of cuv.
  !!
  !! IN:
  !!    this :: rism3d object with computed solution
  !!    cuv  :: site-site pair direct correlation function
  !!OUT:
  !!    DCF integeral for each solvent site
  function rism3d_closure_DCFintegral (this,cuv) result(C)
    implicit none
    type(rism3d_closure) :: this
    _REAL_, intent(in) :: cuv(:,:,:,:)
    _REAL_ ::  C(this%solvent%numAtomTypes)
    integer :: iv, ix, iy, iz, igk
    do iv=1,this%solvent%numAtomTypes
       C(iv) = sum(cuv(:,:,:,iv))
    end do
    C = C*this%grid%voxelVolume
  end function rism3d_closure_DCFintegral

  !> Calculates the forces on the solute contributed by the solvent
  !! according to 3D-RISM. In fact, this subroutine calls the
  !! appropriate subroutines to calculate this.
  subroutine rism3d_closure_force(this, ff, guv, periodicPotential)
    implicit none
    type(rism3d_closure):: this !< Closure object with computed solution.
    _REAL_, intent(out) :: ff(3,this%solute%numAtoms) !< 3D-RISM forces [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair distribution function.
    character(len=*), intent(in) :: periodicPotential !< Label for
    ! periodic potential, else empty string.
    integer :: atom
    integer :: atomRange
    integer :: k, iat

    ff = 0
    call lennardJonesForcePeriodic(this%potential, ff, guv)
    call particleMeshEwaldForce(this%potential, ff, guv)

  end subroutine rism3d_closure_force

  !>Frees memory and resets object state
  !!IN:
  !!   this : the closure object
  subroutine rism3d_closure_destroy(this)
    use safemem
    implicit none
    type(rism3d_closure), intent(inout) :: this
    if (associated(this%kh)) then
       call rism3d_kh_destroy(this%kh)
       deallocate(this%kh)
    end if
    if (associated(this%psen)) then
       call rism3d_psen_destroy(this%psen)
       deallocate(this%psen)
    end if
    if (associated(this%hnc)) then
       call rism3d_hnc_destroy(this%hnc)
       deallocate(this%hnc)
    end if
    nullify(this%potential)
    nullify(this%grid)
    nullify(this%solute)
    nullify(this%solvent)
  end subroutine rism3d_closure_destroy


  !!                         PRIVATE

! Numerical, finite difference way to check force computation.
  subroutine rism3d_checkForceNumDeriv(this, ff, guv)
    use constants, only : PI, KB
    use rism_util, only : checksum
    implicit none
    type(rism3d_closure), intent(in) :: this !> potential object.
    _REAL_, intent(in) :: ff(:, :) !> Force on each atom.
    _REAL_, intent(in) :: guv(:, :) !< Site-site pair correlation function.
    integer :: iu,id
    _REAL_  :: delta 
    _REAL_  :: ene(this%solvent%numAtomTypes)
    _REAL_  :: ffn (3,this%potential%solute%numAtoms)
    delta = 1e-4
    ene = 0.d0
    ffn = 0.d0
    write (6,*) "Entering force debug procedure."
    write (6,*) "-------------------------------"
    call flush(6)
    do iu = 1, this%potential%solute%numAtoms
       do id=1,3
          this%potential%solute%position(id, iu) = this%potential%solute%position(id, iu) - delta
          call rism3d_potential_calc(this%potential)
          ene = rism3d_closure_solvPotEne(this,guv)
          
          this%potential%solute%position(id, iu) = this%potential%solute%position(id, iu) + 2.0*delta
          call rism3d_potential_calc(this%potential)
          ene = ene - rism3d_closure_solvPotEne(this,guv)
          ffn(id,iu) = sum(ene)/(2d0*delta)

          this%potential%solute%position(id, iu) = this%potential%solute%position(id, iu) - delta
      end do
      write (6,*) "Force - Numerical:  ", iu, ffn(:,iu) * KB * this%solvent%temperature
      write (6,*) "Force - Analytical: ", iu, ff(:,iu) * KB * this%solvent%temperature
      write (6,*) "Force diff (A-F):   ", iu, ( ff(:,iu) - ffn(:,iu) ) * KB * this%solvent%temperature

      write (6,*) "Force dev angle:    ", iu, (180d0/PI) * acos(dot_product(ffn(:,iu), ff(:,iu)) &
           / sqrt( dot_product(ff(:,iu),ff(:,iu) ) * dot_product(ffn(:,iu),ffn(:,iu) )))
      write (6,*) "Relative error:     ", sqrt(dot_product(ffn(:,iu)-ff(:,iu),ffn(:,iu)-ff(:,iu))&
          /dot_product(ffn(:,iu),ffn(:,iu)))

      call flush(6)  
     
    end do
    write (6,*) "-------------------------------"
    write (6,*) "Exiting force debug procedure."
  end subroutine rism3d_checkForceNumDeriv


  !> Periodic solute-solvent 12-6 Lennard-Jones force 
  !! subject to the minimum image convention.
  subroutine lennardJonesForcePeriodic(this, ff, guv)
    use constants, only : PI, KB
    use rism_util, only : checksum
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_potential), intent(in) :: this !> potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !> Force on each atom.
    _REAL_, intent(in) :: guv(:, :) !< Site-site pair correlation function.
    
    ! Amount to offset the grid in the z-axis (when using MPI).
    _REAL_ :: offset
    ! Grid, solute, slovent, and dimension indices.
    integer :: igx, igy, igz, ig, iu, iv, id
    ! Grid point position.
    _REAL_ :: rx(3), ry(3), rz(3)
    ! Distance of grid point from solute.
    _REAL_ :: sd2
    ! Base term in LJ equation (ratio of sigma and distance).
    _REAL_ :: ljBaseTerm
    ! Minimum distance to prevent divide by zero.
    _REAL_, parameter :: minDistance = 0.002d0
    _REAL_, parameter :: minDistanceSquared = minDistance**2
    _REAL_ :: solutePosition(3)
    ! Lennard-Jones force between two particles.
    _REAL_ :: dUlj_dr
    _REAL_ :: debugenergy, time0, time1
#ifdef OPENMP
    integer :: numtasks
    character(len=5) omp_num_threads

    !  Following should not be necessary, but ifort doesn't seem to
    !  do the right thing with the OMP_NUM_THREADS environment variable here
    call get_environment_variable('OMP_NUM_THREADS', omp_num_threads)
    read( omp_num_threads, * ) numtasks
#endif
    call wallclock(time0)

!$omp parallel do private (rx,ry,rz,solutePosition,sd2,dUlj_dr,ljBaseTerm, &
!$omp&   igx,igy,igz,iu,iv,ig) num_threads(numtasks)

    do iu = 1, this%solute%numAtoms
       do igz = 1, this%grid%localDimsR(3)
          rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
          do igy = 1, this%grid%localDimsR(2)
             ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
             do igx = 1, this%grid%localDimsR(1)
                rx = (igx - 1) * this%grid%voxelVectorsR(1, :)

#ifdef MPI
                ig = 1 + (igx - 1) + (igy -1) * (this%grid%globalDimsR(1) + 2) &
                       + (igz - 1) * this%grid%globalDimsR(2) &
                             * (this%grid%globalDimsR(1) + 2)
#else
                ig = 1 + (igx - 1) + (igy - 1) * this%grid%globalDimsR(1) &
                       + (igz - 1) * this%grid%globalDimsR(2) &
                                   * this%grid%globalDimsR(1)
#endif
                
                ! Vector from solute atom to grid point, in real space.
                solutePosition = rx + ry + rz - this%solute%position(:, iu)
                
                ! Adjust distance to be to the closest solute atom
                ! image, which may be in an adjacent cell.
                solutePosition = minimumImage(this, solutePosition)
   
                sd2 = dot_product(solutePosition, solutePosition)
                sd2 = max (minDistanceSquared, sd2)

                if (sd2 < this%cutoff2) then
                   dUlj_dr = 0
                   do iv = 1, this%solvent%numAtomTypes
                      ljBaseTerm = sd2 / this%ljSigmaUV(iu, iv)**2
                      ljBaseTerm = 1d0 / ljBaseTerm**3
                      dUlj_dr = dUlj_dr + this%ljEpsilonUV(iu, iv) &
                           * ljBaseTerm * (ljBaseTerm - 1.d0) &
                            * this%solvent%density(iv) * guv(ig, iv)
                   end do
                   ff(:, iu) = ff(:, iu) &
                        - dUlj_dr * solutePosition(:) / sd2
                end if
             end do
          end do
       end do
    end do
!$omp end parallel do
    ff = ff * this%grid%voxelVolume *12d0
    call wallclock(time1)
  end subroutine lennardJonesForcePeriodic

! Computes the PME reciprocal, long range part of electrostatic forces exerted 
! by the solvent charge density onto the solute.
  subroutine particleMeshEwaldForce (this,ff,guv)
    use, intrinsic :: iso_c_binding
    use bspline
    use constants, only : pi, KB
    use FFTW3
    use rism_util, only: r2c_pointer
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif

    type(rism3d_potential), intent(inout) :: this !< potential object.
    _REAL_, intent(inout) :: ff(3, this%solute%numAtoms) !< Force array [kT/A].
    _REAL_, intent(in) :: guv(:,:) !< Site-site pair correlation function.

    _REAL_ :: smear, zeta, t
    ! Order of b-spline interpolation.
    integer, parameter :: splineOrder = 6
    !logical, parameter :: DEBUG = .false.

    integer :: i, j, k
    integer :: ix, iy, iz, ixyz
    integer :: igx, igy, igz
    integer :: igk, igk_fort, igk_cpp, igk_fort_global, igk_cpp_global
    integer :: id
    integer :: iu, iv, ig1


    _REAL_, pointer :: kxi(:,:) => NULL()
    _REAL_, pointer :: kyi(:,:) => NULL()
    _REAL_, pointer :: kzi(:,:) => NULL()
    _REAL_ :: k2
    _REAL_ :: waveVector(3)

    integer :: gridPoints(splineOrder, 3)
    integer :: lgx, lgy, lgz

    _REAL_ :: byz, bxyz

    _REAL_ :: weights(splineOrder, 3), potGrad(3)
    _REAL_ :: weightDerivs(splineOrder, 3), weightYZ(3)
    _REAL_, pointer :: bsplineFourierCoeffX(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffY(:) => NULL()
    _REAL_, pointer :: bsplineFourierCoeffZ(:) => NULL()
    _REAL_, pointer :: gaussianFourierCoeff(:) => NULL()
    _REAL_, pointer :: kernel(:) => NULL()

    type(C_PTR) :: rho_r_cptr, rho_c_cptr, rho_final_cptr
    real(C_DOUBLE), pointer :: rho_r(:,:,:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: rho_c(:,:,:) => NULL()
    real(C_DOUBLE), pointer :: rho_final(:,:,:) => NULL()
    integer(C_INTPTR_T) :: localPtsK
    real(C_DOUBLE), pointer :: outr_1d(:) => NULL()
    integer(C_INTPTR_T) :: L, M, N
    integer(C_INTPTR_T) :: local_N, local_k_offset

    type(C_PTR) :: planfwd = C_NULL_PTR, planbwd = C_NULL_PTR

    _REAL_, pointer :: rhoR(:) => NULL()
    _REAL_, pointer :: rhoK(:) => NULL()

    _REAL_ :: reciprocalPos(3)
    _REAL_ :: chargeCorrection

    integer :: gridDimX_k

    character(len=120) :: filename, suffix
    integer :: ierr
    _REAL_  :: debug_energy
    _REAL_  :: rhoijk

    _REAL_ :: factor, smear2
    _REAL_ :: rx(3), ry(3), rz(3)
    _REAL_, parameter :: minDistance = 0.002
    _REAL_, parameter :: minDistance2 = minDistance**2
    _REAL_ :: gridPoint(3), solutePosition(3)
    _REAL_ :: sd, sd2

    _REAL_ :: qall

    _REAL_ :: time0, time1
#ifdef OPENMP
    integer :: numtasks
    character(len=50) omp_threads

    !  Following should not be necessary, but ifort doesn't seem to
    !  do the right thing with the OMP_NUM_THREADS environment variable here
    call get_environment_variable('OMP_NUM_THREADS', omp_threads)
    read( omp_threads, * ) numtasks
#endif
    call wallclock(time0)


    smear = this%chargeSmear
    zeta = smear * smear
    t = -0.25 * zeta

    factor = 2d0 / (this%chargeSmear * sqrt(pi))
    smear2 = this%chargeSmear**2

    L = this%grid%localDimsR(1)
    M = this%grid%localDimsR(2)
    N = this%grid%globalDimsR(3)


    gridDimX_k = this%grid%localDimsR(1) / 2 + 1
    kxi => safemem_realloc(kxi,  gridDimX_k, 3, .false.)
    kyi => safemem_realloc(kyi,  this%grid%localDimsR(2), 3, .false.)
    kzi => safemem_realloc(kzi,  this%grid%localDimsR(3), 3, .false.)


    bsplineFourierCoeffX => safemem_realloc(bsplineFourierCoeffX, this%grid%globalDimsR(1), .false.)
    bsplineFourierCoeffY => safemem_realloc(bsplineFourierCoeffY, this%grid%globalDimsR(2), .false.)
    bsplineFourierCoeffZ => safemem_realloc(bsplineFourierCoeffZ, this%grid%globalDimsR(3), .false.)

    gaussianFourierCoeff => safemem_realloc(gaussianFourierCoeff, &
        gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

#if defined(MPI)
    localPtsK = fftw_mpi_local_size_3d(N, M, L / 2 + 1, &
         this%grid%mpicomm, local_N, local_k_offset)
#else
    localPtsK = gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3)
    local_N = this%grid%localDimsR(3)
    local_k_offset = 0
#endif

    rho_r_cptr = fftw_alloc_real(2 * localPtsK)
    rho_c_cptr = fftw_alloc_complex(localPtsK)
#if defined(MPI)
    call c_f_pointer(rho_r_cptr, rho_r, [2*(L/2+1), M, local_N])
    call c_f_pointer(rho_r_cptr, outr_1d, [2*(L/2+1) * M * local_N])
    call c_f_pointer(rho_c_cptr, rho_c, [L/2+1, M, local_N])
    if (this%grid%mpirank == 0) then
       rho_final_cptr = fftw_alloc_real(L * M * N)
       call c_f_pointer(rho_final_cptr, rho_final, [L, M, N])
    end if
#else
    call c_f_pointer(rho_r_cptr, rho_r, [L, M, local_N])
    call c_f_pointer(rho_r_cptr, outr_1d, [L * M * local_N])
    call c_f_pointer(rho_c_cptr, rho_c, [L/2+1, M, local_N])
#endif

   ! Angular wave numbers.
   do igk = 0, gridDimX_k - 1
      kxi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(1,:) * igk
   end do
   do igk = 0, this%grid%localDimsR(2) - 1
      kyi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(2,:) &
           * merge(igk, -(this%grid%localDimsR(2) - igk), &
           igk < this%grid%localDimsR(2) / 2 + 1)
   end do
   do igk = 0, this%grid%localDimsR(3) - 1
      kzi(igk + 1,:) = 2 * pi * this%grid%unitCellVectorsK(3,:) &
           * merge(igk + this%grid%offsetR(3), &
           -(this%grid%globalDimsR(3) - (igk + this%grid%offsetR(3))), &
           (igk + this%grid%offsetR(3)) < this%grid%globalDimsR(3) / 2 + 1)
   end do

    kernel => safemem_realloc(kernel, &
         gridDimX_k * this%grid%localDimsR(2) * this%grid%localDimsR(3), .false.)

    ! Compute the discrete Fourier transform coefficients of the b-spline.
    call cardinal_bspline(merge(0d0, 0.5d0, mod(splineOrder, 2) == 0), &
         splineOrder, weights(:,1))

    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(1), splineOrder, &
         weights(:,1), bsplineFourierCoeffX, .true.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(2), splineOrder, &
         weights(:,1), bsplineFourierCoeffY, .false.)
    call cardinal_bspline_Fourier_coefficients( &
         this%grid%globalDimsR(3), splineOrder, &
         weights(:,1), bsplineFourierCoeffZ, .false.)

    ixyz = 1
    do iz = 0, this%grid%localDimsR(3) - 1
       do iy = 0, this%grid%localDimsR(2) - 1
          byz = bsplineFourierCoeffY(iy + 1) * bsplineFourierCoeffZ(iz + 1 + this%grid%offsetR(3))
          do ix = 0, gridDimX_k - 1

             bxyz = bsplineFourierCoeffX(ix + 1) * byz
             waveVector = kxi(ix+1,:) + kyi(iy+1,:) + kzi(iz+1,:)

             k2 = dot_product(waveVector, waveVector)
             !k2 = this%grid%waveVectors2(igk)
             kernel(ixyz) = (4 * pi / k2) * exp(k2 * t) / bxyz
             ixyz = ixyz + 1
          end do
       end do
    end do

    ! Remove the k = 0 term (tinfoil boundary conditions).
    if (this%grid%offsetR(3) == 0) then
       kernel(1) = 0d0
    end if

    ! Allocate the FFT.
#if defined(MPI)
    planfwd = fftw_mpi_plan_dft_r2c_3d(N, M, L, &
         rho_r, rho_c, this%grid%mpicomm, &
         FFTW_ESTIMATE)

    planbwd = fftw_mpi_plan_dft_c2r_3d(N, M, L, &
         rho_c, rho_r, this%grid%mpicomm, &
         FFTW_ESTIMATE)
#else
    planfwd = fftw_plan_dft_r2c_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         rho_r, rho_c, FFTW_ESTIMATE)
    planbwd = fftw_plan_dft_c2r_3d( &
         this%grid%globalDimsR(3), &
         this%grid%globalDimsR(2), &
         this%grid%globalDimsR(1), &
         rho_c, rho_r, FFTW_ESTIMATE)
#endif /*defined(MPI)*/

    ! (1) get solvent charge density, rho(r) on the grid
    rho_r = 0d0
    rho_c = 0d0

    do iv = 1, this%solvent%numAtomTypes
       do igz = 1, this%grid%localDimsR(3)
          do igy = 1, this%grid%localDimsR(2)
             do igx = 1, this%grid%localDimsR(1)
                ig1 = igx + (igy-1) * this%grid%localDimsR(1) + &
                     (igz - 1) * this%grid%localDimsR(2) * this%grid%localDimsR(1)
#if defined(MPI)
                   igk = igx + (igy - 1) * (this%grid%localDimsR(1) + 2) &
                        + (igz - 1) * this%grid%localDimsR(2) * (this%grid%localDimsR(1) + 2)
#else
                   igk = ig1
#endif /*defined(MPI)*/
                   rho_r(igx,igy,igz) = rho_r(igx,igy,igz) &
                     + this%solvent%charge(iv) &
                     * guv(igk,iv) * this%solvent%density(iv)
             end do
          end do
       end do
    end do

    rho_r = rho_r * this%grid%voxelVolume

    !
    ! short range part of PME.
    !

!$omp parallel do private (rx,ry,rz,solutePosition,sd2,sd, &
!$omp&   igx,igy,igz,iu) num_threads(numtasks)

    do iu =1, this%solute%numAtoms
        do igz = 1, this%grid%localDimsR(3)
           rz = (igz - 1 + this%grid%offsetR(3)) * this%grid%voxelVectorsR(3, :)
           do igy = 1, this%grid%localDimsR(2)
              ry = (igy - 1) * this%grid%voxelVectorsR(2, :)
              do igx = 1, this%grid%localDimsR(1)
                 rx = (igx - 1) * this%grid%voxelVectorsR(1, :)
                 solutePosition = rx + ry + rz - this%solute%position(:, iu)
                 solutePosition = minimumImage(this, solutePosition)
                 sd2 = dot_product(solutePosition, solutePosition)
                 sd2 = max(minDistance2,sd2)
                 if (sd2 < this%cutoff2) then
                    sd = sqrt(sd2)
                    ff(:,iu) = ff(:,iu) - rho_r(igx,igy,igz) * this%solute%charge(iu) * &
                          (( factor * exp(-sd2 / smear2) &
                          + erfc(sd / this%chargeSmear) &
                          / sd )/ sd2) * solutePosition
                 end if 
              end do
           end do
        end do
    end do
!$omp end parallel do

    ! (2) FT [ rho(r) ]
    ! Convert the charge density into reciprocal space.

#if defined(MPI)
    call fftw_mpi_execute_dft_r2c(planfwd, rho_r, rho_c)
#else
    call fftw_execute_dft_r2c(planfwd, rho_r, rho_c)
#endif

    ! Convolution.
    do igz = 0, local_N - 1
       do igy = 0, M - 1
          do igx = 0, (L / 2 + 1) - 1
             igk_fort = 1 + igx + (igy + igz * M) * (L / 2 + 1)
             rho_c(igx + 1, igy + 1, igz + 1) &
                  = rho_c(igx + 1, igy + 1, igz + 1) * kernel(igk_fort)
          end do
       end do
    end do

    ! Evaluate the recip-space potential at the grid points.
    rho_r = 0
#if defined(MPI)
    call fftw_mpi_execute_dft_c2r(planbwd, rho_c, rho_r)
#else
    call fftw_execute_dft_c2r(planbwd, rho_c, rho_r)
#endif
    rho_r = rho_r / this%grid%boxVolume

! Here is where interpolation/differentiation starts.

    do iu = 1, this%solute%numAtoms
       ! Convert Cartesian position to reciprocal space by projecting
       ! to reciprocal unit cell vectors.
       do id = 1, 3
          reciprocalPos(id) = dot_product(this%solute%position(:,iu), &
               this%grid%unitCellVectorsK(id, :))
       end do

       ! Set the box length to unity since reciprocal space already
       ! divides positions by the box length.
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(1), 1d0, &
            this%grid%globalDimsR(1), &
            splineOrder, gridPoints(:,1), weights(:,1), weightDerivs(:,1))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(2), 1d0, &
            this%grid%globalDimsR(2), &
            splineOrder, gridPoints(:,2), weights(:,2), weightDerivs(:,2))
       call cardinal_bspline_periodic_grid( &
            reciprocalPos(3), 1d0, &
            this%grid%globalDimsR(3), &
            splineOrder, gridPoints(:,3), weights(:,3), weightDerivs(:,3))

       potGrad = 0

       do k = 1, splineOrder
#if defined(MPI)
          if (    (gridPoints(k,3) < this%grid%offsetR(3)) &
             .or. (gridPoints(k,3) + 1 > this%grid%offsetR(3) &
                                       + this%grid%localDimsR(3))) cycle
#endif
          do i = 1, splineOrder
             do j = 1, splineOrder

                rhoijk = rho_R(1 + gridPoints(i,1), &
                               1 + gridPoints(j,2), &
                               1 + gridPoints(k,3) - this%grid%offsetR(3))
               
                potGrad(1) = potGrad(1) + weights(k,3) &
                    * weightDerivs(i,1) * weights(j,2) * rhoijk 
                potGrad(2) = potGrad(2) + weights(k,3) &
                    * weights(i,1)      * weightDerivs(j,2) * rhoijk 
                potGrad(3) = potGrad(3) + weightDerivs(k,3) &
                    * weights(i,1)      * weights(j,2) * rhoijk

            end do
          end do
       end do
       
       ! Apply chain rule from converting b-spline derivatives from
       ! reciprocal space to real space positions.
       potGrad = this%grid%globalDimsR * this%solute%charge(iu) * potGrad

       do id = 1, 3
          ff(id, iu) = ff(id, iu) - &
               dot_product(this%grid%unitCellVectorsK(:,id), potGrad)
       end do
    end do

    ! Deallocate the FFT plans.
    call fftw_destroy_plan(planfwd)
    call fftw_destroy_plan(planbwd)

    if (safemem_dealloc(kxi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kyi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(kzi) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if

    if (safemem_dealloc(bsplineFourierCoeffX) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffY) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(bsplineFourierCoeffZ) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    if (safemem_dealloc(gaussianFourierCoeff) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    call fftw_free(rho_r_cptr)
    call fftw_free(rho_c_cptr)
#if defined(MPI)
    if (this%grid%mpirank == 0) then
       call fftw_free(rho_final_cptr)
    end if
#endif
    if (safemem_dealloc(kernel) /= 0) then
       call rism_report_error("uvParticleMeshEwaldPotential: Failed to deallocate arrays.")
    end if
    call wallclock(time1)
    ! write(6,*) 'short-range PME force: ', time1 - time0

  end subroutine particleMeshEwaldForce 

end module rism3d_closure_c
