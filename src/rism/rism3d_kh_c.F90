!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999);
!ibid. 112:10391-10417 (2000).
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010).

#include "../include/dprec.fh"

!> Kovalenko-Hirata closure class for 1D-RISM.
!! Kovalenko, A.; Hirata, F. J. Chem. Phys. 1999, 110, 10095–10112.
module rism3d_kh_c
  use rism3d_potential_c
  use rism3d_grid_c
  !the KH type
  type rism3d_kh
     type(rism3d_potential),pointer :: pot => NULL()
     !grid : points to grid in potential object
     type(rism3d_grid),pointer :: grid => NULL()
  end type rism3d_kh

  public rism3d_kh_new, rism3d_kh_destroy!, rism3d_kh_guv
  
contains
  

  !> Initializes the KH closure.
  !! @param[in,out] this KH object.
  subroutine rism3d_kh_new(this, pot)
    implicit none
    type(rism3d_kh), intent(inout) :: this
    type(rism3d_potential), target, intent(in) :: pot
    this%pot => pot
    this%grid => this%pot%grid
  end subroutine rism3d_kh_new


  !> Calculates Guv from Uuv, Huv, and Cuv using the KH closure.
  !! @param[in] this KH closure object.
  !! @param[out] guv Site-site pair correlation function.
  !! @param[in] huv Site-site total correlation function.
  !! @param[in] cuv Site-site direct correlation function.
  subroutine rism3d_kh_guv(this, guv, huv, cuv)
    implicit none
    type(rism3d_kh), intent(in) :: this
    _REAL_, intent(out) :: guv(:,:)
    _REAL_, intent(in) :: huv(:,:), cuv(:,:,:,:)
    integer :: iv, ir, ix, iy, iz, ig 
    _REAL_ :: exponent

    do iv = 1,this%pot%solvent%numAtomTypes
       do iz = 1, this%grid%localDimsR(3)
          do iy = 1, this%grid%localDimsR(2)
             do ix = 1, this%grid%localDimsR(1)
#if defined(MPI)
                ig = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                    (iz-1)*(this%grid%localDimsR(1)+2)*this%grid%localDimsR(2)
#else
                ig = ix + (iy-1)*this%grid%localDimsR(1) + &
                    (iz-1)*this%grid%localDimsR(1)*this%grid%localDimsR(2)
#endif /*defined(MPI)*/
                ! exponent = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv) &
                !      - this%pot%tcfBackgroundChargeCorrection(iv)
                ! if (exponent >= 0d0) then
                !    guv(ig,iv) = 1d0 + exponent + this%pot%tcfBackgroundChargeCorrection(iv)
                ! else
                !    guv(ig,iv) = exp(exponent) + this%pot%tcfBackgroundChargeCorrection(iv)
                ! end if
                !FIXME: Below is the non-periodic case.
                exponent = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                if (exponent >= 0d0) then
                   guv(ig,iv) = 1d0 + exponent
                else
                   guv(ig,iv) = exp(exponent)
                end if
             end do
          end do
       end do
    end do
  end subroutine rism3d_kh_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_kh_excessChemicalPotential(this, huv, cuv) result(excessChemicalPotential)
    implicit none
    type(rism3d_kh), intent(in) :: this
    _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
    _REAL_ :: excessChemicalPotential(this%pot%solvent%numAtomTypes)
    _REAL_ :: tuv
    integer :: ix, iy, iz, iv, igk
    excessChemicalPotential = 0.d0
    do iv=1,this%pot%solvent%numAtomTypes
       do iz=1,this%grid%localDimsR(3)
          do iy=1,this%grid%localDimsR(2)
             do ix=1,this%grid%localDimsR(1)
#if defined(MPI)
                igk = ix + (iy-1)*(this%grid%localDimsR(1)+2) + &
                   (iz-1)*this%grid%localDimsR(2)*(this%grid%localDimsR(1)+2)
#else
                igk = ix + (iy-1)*this%grid%localDimsR(1) + &
                    (iz-1)*this%grid%localDimsR(2)*this%grid%localDimsR(1)
#endif /*defined(MPI)*/
                tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                if (huv(igk,iv) > 0d0) then
                   excessChemicalPotential(iv) = excessChemicalPotential(iv) &
                        - (1.d0 + 0.5d0 * huv(igk,iv)) * cuv(ix,iy,iz,iv)
                else
                   excessChemicalPotential(iv) = excessChemicalPotential(iv) &
                        + 0.5d0 * huv(igk,iv) * tuv - cuv(ix,iy,iz,iv)
                end if

             end do
          end do
       end do
       excessChemicalPotential(iv) =  this%pot%solvent%density(iv) * &
            excessChemicalPotential(iv) * this%grid%voxelVolume
    enddo
  end function rism3d_kh_excessChemicalPotential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the KH closure
!!!IN:
!!!   this : KH object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_kh_destroy(this)
    implicit none
    type(rism3d_kh), intent(inout) :: this
    nullify(this%pot)
    nullify(this%grid)
  end subroutine rism3d_kh_destroy
end module rism3d_kh_c
