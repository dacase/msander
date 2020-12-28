! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"
#define VACUUM3 603

! simplified egb() routines: no LES, no MPI.  Hope this may help
!   with other optimizations, like openmp

module genborn

_REAL_, private, dimension(:), allocatable :: r2x,rjx,vectmp1,vectmp2, &
                             vectmp3,vectmp4,vectmp5,sumdeijda,psi
integer, private, dimension(:), allocatable :: temp_jj,k_vals,j_vals, neckidx
logical, private, dimension(:), allocatable :: skipv


public allocate_gb, deallocate_gb, egb, igb7_init

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ allocates scratch space for egb()
subroutine allocate_gb( natom,ncopy )

   implicit none
   integer, intent(in) :: natom, ncopy
   integer ier

   allocate( r2x(natom), rjx(natom), vectmp1(natom), vectmp2(natom),   &
             vectmp3(natom), vectmp4(natom), vectmp5(natom), &
             sumdeijda(natom), psi(natom), temp_jj(natom),  &
             skipv(0:natom),k_vals(natom),j_vals(natom), neckidx(natom), &
             stat = ier )

   REQUIRE( ier == 0 )
   return

end subroutine allocate_gb

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ deallocates scratch space for egb()
subroutine deallocate_gb( )

   implicit none
   integer ier

   ! assume that if r2x is allocated then all are allocated
   if ( allocated( r2x ) ) then
      deallocate( skipv, neckidx, j_vals, k_vals, temp_jj, psi, sumdeijda, &
            vectmp5, vectmp4, vectmp3, vectmp2, vectmp1, rjx, r2x, stat = ier )
      REQUIRE( ier == 0 )
   else
      REQUIRE( .false. )  ! cannot deallocate un-allocated array
   end if
   return

end subroutine deallocate_gb

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize table of indexes for GB neck lookup table.
subroutine igb7_init(natom, rborn)

  implicit none
  integer, intent(in) :: natom

  integer i

  _REAL_, intent(in) :: rborn(*)

  do i=1,natom
     neckidx(i) = nint((rborn(i)-1.0d0)*20d0)
     if (neckidx(i) < 0 .or. neckidx(i) > 20) then
        write (6,*) "Atom ",i," has radius ",rborn(i), &
                    " outside of allowed range"
        write (6,*) "of 1.0 to 2.0 angstroms for igb=7 or 8. Regenerate &
                    &prmtop file with bondi radii."
        call mexit(6,1)
     end if
  end do

end subroutine igb7_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ handles generalized Born functionality, plus reg. nonbon, plus surface area
subroutine egb(x,f,rborn,fs,reff,onereff,charge,iac,ico,numex, &
      natex,dcharge,cut,ntypes,natom,natbel, &
      epol,eelt,evdw,esurf,dvdl,vdwrad,ineighbor,p1,p2,p3,p4, &
      ncopy &
#ifndef LES
      ,gbvalpha,gbvbeta,gbvgamma &
#endif
      )
      ! gbvalpha,gbvbeta,gbvgamma arrays for GB
      ! put all gbalpha, gbbeta, gbgamma for each atom in these 3 arrays
      ! (look at mdread.f for details)

   !--------------------------------------------------------------------------

   !     Compute nonbonded interactions with a generalized Born model,
   !     getting the "effective" Born radii via the approximate pairwise method
   !     Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
   !     (1996).  Aside from the scaling of the radii, this is the same
   !     approach developed in Schaefer and Froemmel, JMB 216:1045 (1990).

   !     The input coordinates are in the "x" array, and the forces in "f"
   !     get updated; energy components are returned in "epol", "eelt" and
   !     "evdw".

   !     Input parameters for the generalized Born model are "rborn(i)", the
   !     intrinsic dielectric radius of atom "i", and "fs(i)", which is
   !     set (in routine mdread) to (rborn(i) - offset)*si.

   !     Input parameters for the "gas-phase" electrostatic energies are
   !     the charges, in the "charge()" array.

   !     Input parameters for the van der Waals terms are "cn1()" and "cn2()",
   !     containing LJ 12-6 parameters, and "asol" and "bsol" containing
   !     LJ 12-10 parameters.  (The latter are not used in 1994 are more
   !     forcefields.)  The "iac" and "ico" arrays are used to point into
   !     these matrices of coefficients.

   !     The "numex" and "natex" arrays are used to find "excluded" pairs of
   !     atoms, for which gas-phase electrostatics and LJ terms are skipped;
   !     note that GB terms are computed for all pairs of atoms.

   !     If gbsa=1, then an approximate surface-area dependent term is
   !     computed, with the resulting energy placed into "esurf".  The
   !     algorithm is from J. Weiser, P.S. Shenkin, and W.C. Still,
   !     "Approximate atomic sufraces from linear combinations of pariwise
   !     overlaps (LCPO)", J. Computat. Chem. 20:217 (1999).

   !     The code also supports a multiple-time-step facility:

   !       pairs closer than sqrt(cut_inner) are evaluated every nrespai steps;
   !         "   between sqrt(cut_inner) and sqrt(cut) are evaluated
   !                        every nrespa steps
   !         "   beyond sqrt(cut) are ignored

   !       the forces arising from the derivatives of the GB terms with respect
   !          to the effective Born radii are evaulated every nrespa steps

   !       the surface-area dependent term is evaluated every nrespa steps

   !       the effective radii are only updated every nrespai steps

   !     (Be careful with the above: what seems to work is dt=0.001,
   !     nrespai=2, nrespa=4; anything beyond this seems dangerous.)

   !     Written 1999-2000, primarily by D.A. Case, with help from C. Brooks,
   !       T. Simonson, R. Sinkovits  and V. Tsui.  The LCPO implementation
   !       was written by V. Tsui.

   !     Vectorization and optimization 1999-2000, primarily by C. P. Sosa,
   !       T. Hewitt, and D. A. Case.  Work presented at CUG Fall of 2000.
   !--------------------------------------------------------------------------

   use icosasurf, only : icosa_init, icosa_sphere_approx
   use qmmm_module, only : qmmm_nml,qmmm_struct,qm2_struct
   use parms, only: cn1,cn2
   use constants, only: zero, one, two, three, four, five, six, seven, &
                        eight, nine, ten, eleven, twelve, half, third, &
                        fourth, eighth, pi, fourpi, alpb_alpha, &
                        AMBER_ELECTROSTATIC
   use crg_reloc, only: ifcr, cr_add_dcdr_factor
   use commandline_module, only: cpein_specified

   implicit none

#  include "../include/md.h"
#  include "def_time.h"
   _REAL_ temp7

   _REAL_ deel
   logical onstep,onstepi,oncpstep,oncestep,doeel, dovdw
   _REAL_ x,f,rborn,fs,reff,onereff,charge,cut, &
         epol,eelt,evdw,esurf,vdwrad,p1,p2,p3,p4, &
         dcharge
   _REAL_ totsasa,extdieli,intdieli, &
         lastxj,lastyj,lastzj,xi,yi,zi,ri,four_ri,ri1i,xij,yij,zij, &
         dij1i,r2,dij,sj,sj2,frespa,si,sumaij,sumajk, &
         sumaijajk,sumdaijddijdxi,sumdaijddijdyi,sumdaijddijdzi, &
         sumdaijddijdxiajk,sumdaijddijdyiajk,sumdaijddijdziajk, &
         xj,yj,zj,rij,tmpaij,aij,daijddij,daijddijdxj, daijddijdyj, &
         daijddijdzj,sumajk2,sumdajkddjkdxj,sumdajkddjkdyj, &
         sumdajkddjkdzj,p3p4aij,xk,yk,zk,rjk2,djk1i,rjk,vdw2dif, &
         tmpajk,ajk,dajkddjk,dajkddjkdxj,dajkddjkdyj,dajkddjkdzj, &
         daidxj,daidyj,daidzj,ai,daidxi,daidyi,daidzi,qi,dumx, &
         dumy,dumz,de,rj,temp1,fgbi,rinv,r2inv,qiqj,fgbk,expmkf, &
         dl,e,temp4,temp5,temp6,eel,r6inv,f6,f12, &
         dedx,dedy,dedz,qi2h,dij2i,datmp,dij3i, &
         qid2h,dvdl,thi,thi2, self_e,reff_j
!  _REAL_ intdiele_inv !variable for gas phase calculations (future fix by Dan Parkin)
    _REAL_ :: alpb_beta, one_Arad_beta
    _REAL_ :: gbvalpha(*),gbvbeta(*),gbvgamma(*) !add gbvalpha,gbvbeta,gbvgamma

   integer count,count2,icount,ineighbor(*),max_count, iminus
   integer iac,ico,numex,natex,ntypes,natom,natbel,ncopy
   integer i,j,k,kk1,maxi,num_j_vals,jjj,count2_fin,num_k_vals, &
           iexcl,iaci,jexcl,jexcl_last,jjv,ic,kk
   integer j3
   _REAL_ f_x,f_y,f_z,f_xi,f_yi,f_zi
   _REAL_ dumbo, tmpsd, rborn_i, psi_i
   _REAL_ gba_i, gbb_i, gbg_i ! temporary variables for GB (calculating GB force)

   ! Variables for QMMM specific loops
   integer qm_temp_count
!Locals for link atoms
   _REAL_ :: forcemod(3)
   integer :: lnk_no, mm_no, qm_no


   ! variables needed for icosa surface area calculation
   integer ineighborpt

   ! variables needed for smooth integration cutoff in Reff:
   _REAL_ rgbmax2, rgbmax1i, rgbmax2i,rgbmaxpsmax2
   !     _REAL_  datmp2

   _REAL_ onekappa !1/kappa

   ! Scratch variables used for calculating neck correction
   _REAL_ mdist,mdist2,mdist3,mdist5,mdist6

#include "gbneck.h"

   dimension x(*),f(*),rborn(*),charge(*),iac(*), &
         ico(*),numex(*),natex(*),fs(*),reff(natom), onereff(natom), &
         vdwrad(*),p1(*),p2(*),p3(*),p4(*), &
         dcharge(*)

   integer ncopy2

   !   FGB taylor coefficients follow
   !   from A to H :
   !   1/3 , 2/5 , 3/7 , 4/9 , 5/11
   !   4/3 , 12/5 , 24/7 , 40/9 , 60/11

   _REAL_  ta
   _REAL_  tb
   _REAL_  tc
   _REAL_  td
   _REAL_  tdd
   _REAL_  te
   _REAL_  tf
   _REAL_  tg
   _REAL_  th
   _REAL_  thh
   parameter ( ta   = third )
   parameter ( tb   = two / five )
   parameter ( tc   = three / seven )
   parameter ( td   = four / nine )
   parameter ( tdd  = five / eleven )
   parameter ( te   = four / three )
   parameter ( tf   = twelve / five )
   parameter ( tg   = three * eight / seven )
   parameter ( th   = five * eight / nine )
   parameter ( thh  = three * four * five / eleven )

   qid2h = 0.d0
   dl = 0.d0
   f6 = 0.d0
   f12 = 0.d0

!============ QMMM ==========================================================
!  If ifqnt == True and igb /=0 and igb /= 6 then we will do EGB with QMMM.

!  In this case the charge array should currently be zero for the QM atoms.
!  This corresponds to qmgb = 0 where we do only EGB(MM-MM) -- skipping
!  QM-MM and QM-QM interactions.

!  If qmgb==1 then we need to copy back the qm resp charges from
!  qm_resp_charges, so we can do EGB on everything. When done we need to
!  ensure we zero the QM charge array again.

!  Before doing the non-bonded electrostatics:
   if (qmmm_nml%ifqnt) then
      if (qmmm_nml%qmgb == 1) then

         call qmmm_restore_mm_charges(qmmm_struct%nquant,qmmm_struct%qm_resp_charges,charge, &
                                      qmmm_struct%scaled_mm_charges, qmmm_struct%iqmatoms, &
                                      qmmm_nml%chg_lambda,qmmm_struct%nlink,qmmm_struct%link_pairs, &
                                      qmmm_struct%mm_link_pair_resp_charges, &
                                      .false.)

         ! Unlike ewald where we don't explicitly skip QM-MM interactions,
         ! in GB we do explicitly skip them so it doesn't matter that the
         ! QM charges are not set to zero in the charge array.

      else if (qmmm_nml%qmgb > 1 ) then

        ! In this case we need to fill the charge array with the Mulliken
        ! charges from the SCF.  Note if igb==2 we skip the radii calculation
        ! below as this has already been done before the call to qm_mm.
        do qm_temp_count = 1, qmmm_struct%nquant_nlink
           charge(qmmm_struct%iqmatoms(qm_temp_count)) = &
                   qm2_struct%scf_mchg(qm_temp_count)*AMBER_ELECTROSTATIC
        end do
      end if

      !We need to replace the MM link pair coordinates with
      !the MM link atom coordinates.
      call adj_mm_link_pair_crd(x)

      !We also need to zero the nlink part of dxyzqm so we can accumulate link atom forces in this routine.
      qmmm_struct%dxyzqm(1:3,qmmm_struct%nquant+1:qmmm_struct%nquant_nlink) = zero

   end if
!============ END QMMM =======================================================

   epol = zero
   eelt = zero
   evdw = zero
   esurf = zero
   totsasa = zero
   thi2 = zero
   onekappa = zero
   count2_fin = 0
   onstep = mod(irespa,nrespa) == 0
   onstepi = mod(irespa,nrespai) == 0
   if( .not.onstepi ) return

   oncpstep = ((icnstph == 1 .or. (icnste == 1 .and. cpein_specified)) .and. &
              mod(irespa, ntcnstph) == 0) .or. (icnstph == 2 .or. (icnste == 2 &
              .and. cpein_specified))

   oncestep = ((icnste == 1 .and. mod(irespa, ntcnste) == 0) .or. icnste == 2) &
              .and. .not. cpein_specified

!     variable for gas phase calculations (future fix by Dan Parkin)
!  intdiele_inv = one/intdiel

   if(alpb == 1) then
!     Sigalov Onufriev ALPB (epsilon-dependent GB):
      alpb_beta = alpb_alpha*(intdiel/extdiel)
      extdieli = one/(extdiel*(one + alpb_beta))
      intdieli = one/(intdiel*(one + alpb_beta))
      one_Arad_beta = alpb_beta/Arad
      if (kappa/=zero) onekappa = one/kappa
   else
   !  Standard Still's GB - alpb=0
      extdieli = one/extdiel
      intdieli = one/intdiel
      one_Arad_beta = zero
   end if

   ! Smooth "cut-off" in calculating GB effective radii.
   ! Implementd by Andreas Svrcek-Seiler and Alexey Onufriev.
   ! The integration over solute is performed up to rgbmax and includes
   ! parts of spheres; that is an atom is not just "in" or "out", as
   ! with standard non-bonded cut.  As a result, calclated effective
   ! radii are less than rgbmax. This saves time, and there is no
   ! discontinuity in dReff/drij.

   ! Only the case rgbmax > 5*max(sij) = 5*fsmax ~ 9A is handled; this is
   ! enforced in mdread().  Smaller values would not make much physical
   ! sense anyway.

   rgbmax2 = rgbmax*rgbmax
   rgbmax1i = one/rgbmax
   rgbmax2i = rgbmax1i*rgbmax1i
   rgbmaxpsmax2 = (rgbmax+fsmax)**2

   ncopy2 = ncopy

   !---------------------------------------------------------------------------
   !      Step 1: loop over pairs of atoms to compute the effective Born radii.
   !---------------------------------------------------------------------------
   !
   ! The effective Born radii are now calculated via a call at the
   ! beginning of force.
   !
   !---------------------------------------------------------------------------

   iexcl = 1 !moved to outside the index loop from original location in "step 2"

   maxi = natom
   if(natbel > 0) maxi = natbel

   !--------------------------------------------------------------------------
   !
   !     Step 2: loop over all pairs of atoms, computing the gas-phase
   !             electrostatic energies, the LJ terms, and the off-diagonal
   !             GB terms.  Also accumulate the derivatives of these off-
   !             diagonal terms with respect to the inverse effective radii,
   !             sumdeijda(k) will hold  sum over i,j>i ( deij/dak ),  where
   !             "ak" is the inverse of the effective radius for atom "k".
   !
   !             Update the forces with the negative derivatives of the
   !             gas-phase terms, plus the derivatives of the explicit
   !             distance dependence in Fgb, i.e. the derivatives of the
   !             GB energy terms assuming that the effective radii are
   !             constant.
   !
   !--------------------------------------------------------------------------

   sumdeijda(1:natom) = zero

   call timer_start(TIME_GBFRC)

   !       Note: this code assumes that the belly atoms are the first natbel
   !             atoms...this is checked in mdread.

   do i=1,maxi

      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      qi = charge(i)
      ri = reff(i)
      four_ri = four*reff(i)
      iaci = ntypes * (iac(i) - 1)
      jexcl = iexcl
      jexcl_last = iexcl + numex(i) -1
      dumx = zero
      dumy = zero
      dumz = zero

#if defined(LES)
      icnum = cnum(i)
      nrg_vdw_tmp = zero
      nrg_ele_tmp = zero
      nrg_egb_tmp = zero
#endif

      !         -- check the exclusion list for eel and vdw:

      do k=i+1,natom
         skipv(k) = .false.
      end do
      do jjv=jexcl,jexcl_last
         skipv(natex(jjv))=.true.
      end do

      ! QMMM:  We have 2 lists here:
      !        atom_mask which is natom long and .true. for each qm atom
      !        skipv which is natom long and is .true. if atom v should
      !             be skipped for this atom

      !Step 1 - is current atom i a QM atom - NOT link atoms?

      if (qmmm_nml%ifqnt) then
         if ( qmmm_struct%atom_mask(i) ) then  !yes it is
           !step 2 - exclude all other QM atoms
           do qm_temp_count = 1, qmmm_struct%nquant
              skipv(qmmm_struct%iqmatoms(qm_temp_count)) = .true.
              !The VDW between the link atom and the MM atoms is not
              !calculated. Instead MM link pair atom with Real QM atoms is calculated.
           end do
         !Is atom i an MM link pair atom?
         else if ( qmmm_struct%mm_link_mask(i) ) then !yes it is
           !Exclude all interactions (eel and VDW) with this
           !atom. We will do the MM link pair VDW terms manually
           !at the end of this routine.
           skipv(i+1:natom) = .true. !We can't just cycle here because we want the GB interactions.
         end if
      end if

      ! END QMMM

      icount = 0
      do j=i+1,natom

         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         r2 = xij*xij + yij*yij + zij*zij

         if( r2 <= cut .and. (onstep .or. r2 <= cut_inner)) then
! set for non-LES case so we can refer to reff_j not reff(j)
! this is needed with LES since each atom j has multiple reff

             reff_j = reff(j)

             icount = icount + 1
             temp_jj(icount) = j
             r2x(icount) = r2

! carlos changed this, we are STILL INSIDE LOOP OVER K for multiple LES reff_j
! so we need to use reff_j not reff(j) (since LES may have multiple reff_j for atom j)

             rjx(icount) = reff_j

        end if !r2 <= cut
      end do

      if( igb/=6 ) then

         vectmp1(1:icount) = four_ri*rjx(1:icount)

         call vdinv( icount, vectmp1, vectmp1 ) !Invert things
         vectmp1(1:icount) = -r2x(1:icount)*vectmp1(1:icount)
         call vdexp( icount, vectmp1, vectmp1 )
              ! [ends up with Exp(-rij^2/[4*ai*aj])]
         vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount)*ri*vectmp1(1:icount)
              ! [ends up with fij]

         ! CARLOS: LES
         ! we now have two sets of vectors: one is 1 to last atom, the other
         ! is 1 to icount
         ! the atom one is for vdw/ele and for excl, the other is for gb offdiag
         ! which do we loop over? unless we have a pointer we need to loop
         ! over the bigger one, have a pointer j in it, and pull excl out
         ! of that. but how to know to skip the nonbonds for all after the
         ! first in a loop over a pair? maybe we should have the excl loop
         ! done for icount, not atoms.
         ! SKIPV PART OK, DO THE DISTANCE PART SO WE DON'T HAVE TO INV ALL RIJ
         ! FOR NOW IT'S OK, JUST SLOW

         call vdinvsqrt( icount, vectmp3, vectmp2 ) !vectmp2 = 1/fij

         if( kappa /= zero ) then
            call vdinv( icount, vectmp2, vectmp3 )
            vectmp3(1:icount) = -kappa*vectmp3(1:icount)
            call vdexp( icount, vectmp3, vectmp4 ) !exp(-kappa*fij)
         end if

      end if
      call vdinvsqrt( icount, r2x, vectmp5 ) !1/rij

      ! vectmp1 = Exp(-rij^2/[4*ai*aj])
      ! vectmp2 = 1/fij
      ! vectmp3 = -kappa*fij - if kappa/=zero, otherwise = fij
      ! vectmp4 = exp(-kappa*fij)
      ! vectmp5 = 1/rij - note with qmmm this contains the
                          !distance to mm_link pairs not to QM link atoms.

  !---- Start first outer loop ----
!!! !dir$ ivdep

      do k=1,icount
         j = temp_jj(k)

         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         r2 = r2x(k)
         qiqj = qi * charge(j)
         if( igb/=6 ) then

           if( kappa == zero ) then
              fgbk = zero
              expmkf = extdieli
           else
              expmkf = vectmp4(k)*extdieli


              fgbk = vectmp3(k)*expmkf !-kappa*fij*exp(-kappa*fij)/Eout
              if(alpb == 1) then ! Sigalov Onufriev ALPB:
                 fgbk = fgbk+(fgbk*one_Arad_beta*(-vectmp3(k)*onekappa))

                 ! (-kappa*fij*exp(-kappa*fij)(1 + fij*ab/A)/Eout)*(1/fij+ab/A)
                 ! Note: -vectmp3(k)*onekappa = fij
              end if
           end if
           dl = intdieli - expmkf

           fgbi = vectmp2(k)  !1/fij

           temp1 = -dl*(fgbi + one_Arad_beta)
           e = qiqj*temp1
           if ( ifcr /= 0 ) then
              call cr_add_dcdr_factor( i, temp1*charge(j) )
              call cr_add_dcdr_factor( j, temp1*qi )
           end if
           epol = epol + e

           if (oncpstep .or. oncestep) then
              dvdl = dvdl - (dl*fgbi*dcharge(i)*dcharge(j)+e)
           end if
           temp4 = fgbi*fgbi*fgbi !1/fij^3

           !   [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]

           temp6 = -qiqj*temp4*(dl + fgbk)

           ! -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
           ! exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]

           temp1 = vectmp1(k) !Exp(-rij^2/[4*ai*aj])
           de = temp6*(one - fourth*temp1)

           rj = rjx(k)
           temp5 = half*temp1*temp6*(ri*rj + fourth*r2)

           sumdeijda(i) = sumdeijda(i) + ri*temp5
           sumdeijda(j) = sumdeijda(j) + rj*temp5

         !    -- skip exclusions for remaining terms:
         else
           de = zero
         end if !if( igb/=6 )

         if( .not. skipv(j) ) then

            !   -- gas-phase Coulomb energy:

            ! Note: With QMMM we don't need to explicitly exclude QM-MM
            !       interactions here since the QM charges should
            !       be zero.  Unless qmgb/=0 in which case they are not
            !       zero so we need to exclude them.
            ! We also need to exclude the QM link atom interactions.

            doeel = .true.
            dovdw = .true.
            !  check if i or j is a quantum atom:
            if ( qmmm_nml%ifqnt ) then
               if ( qmmm_struct%atom_mask(i) .OR. qmmm_struct%atom_mask(j) ) then
                  doeel = .false.
               end if

               ! if i or j is an MM link pair atom we also need to exclude it.
               ! We also need to skip the VDW term if either is a mm link
               ! pair atom
               if (qmmm_struct%mm_link_mask(i) .or. qmmm_struct%mm_link_mask(j)) then
                 doeel = .false.
                 dovdw = .false.
               end if
            end if

            !we can use the cached values.
            rinv = vectmp5(k) !1/rij
            r2inv = rinv*rinv

            if( doeel ) then

               ! A future fix implemented by Dan Parkin from Waseda Uni.
               ! intdieli is not 1/Ein when alpb=1
               ! intdiel               : Ein
               ! intdieli (alpb==1)    : one/(intdiel*(one + alpb_beta))
               ! intdieli (alpb/=1)    : one/intdiel
               ! intdiele_inv (always) : one/intdiel [newly added variable]
               temp1 = intdieli*rinv !comment out to implement Dan Parkin's fix
!              temp1 = intdiele_inv*rinv !add to implement Dan Parkin's fix
               eel = qiqj*temp1
               if ( ifcr /= 0 ) then
                  call cr_add_dcdr_factor( i, temp1*charge(j) )
                  call cr_add_dcdr_factor( j, temp1*qi )
               end if
               deel=eel*r2inv
            else
               eel = zero
               deel = zero
            end if

            eelt = eelt + eel
            de = de + deel
            if (oncpstep .or. oncestep) then
               dvdl = dvdl + (intdieli*rinv*dcharge(i) &
                     *dcharge(j) - eel)
            end if

            !    -- van der Waals energy:

            ic = ico( iaci + iac(j) )
            if( ic > 0 .and. dovdw ) then
               !                                    6-12 potential:
               r6inv = r2inv*r2inv*r2inv
               f6 = cn2(ic)*r6inv
               f12 = cn1(ic)*(r6inv*r6inv)
               evdw = evdw + (f12 - f6)
               de = de + (twelve*f12 - six*f6)*r2inv

            end if  ! ( ic > 0 )
         end if  ! ( .not. skipv(j) )

         !    -- derivatives:
         if( onstep .and. r2 > cut_inner ) then
            de = de*nrespa
         else
            de = de*nrespai
         end if
         dedx = de * xij
         dedy = de * yij
         dedz = de * zij
         dumx = dumx + dedx
         dumy = dumy + dedy
         dumz = dumz + dedz

         if (qmmm_nml%ifqnt) then
           !Currently if j is a link atom then
           !dedxy,y,z are force on the link atom, not on the MM link pair.
           !we have to use the chain rule to put the forces back onto the MM link
           !pair atom.
           !For the moment we just accumulate the forces on the end of dxyzqm
           !this should have been zeroed at the beginning of this routine.
           if (qmmm_struct%mm_link_mask(j)) then
             !j is a link atom, instead of adding it to the main force
             !array, accumulate it in dxyzqm
             !Find the link id for this j
             do qm_temp_count = 1,qmmm_struct%nlink
               if (qmmm_struct%link_pairs(1,qm_temp_count) == j) exit
             end do
             qm_temp_count = qm_temp_count + qmmm_struct%nquant
             qmmm_struct%dxyzqm(1,qm_temp_count) = qmmm_struct%dxyzqm(1,qm_temp_count) + dedx
             qmmm_struct%dxyzqm(2,qm_temp_count) = qmmm_struct%dxyzqm(2,qm_temp_count) + dedy
             qmmm_struct%dxyzqm(3,qm_temp_count) = qmmm_struct%dxyzqm(3,qm_temp_count) + dedz
           else
             !Not a link atom, can just add to main force array.
             f(3*j-2) = f(3*j-2) - dedx
             f(3*j-1) = f(3*j-1) - dedy
             f(3*j  ) = f(3*j  ) - dedz
           end if
         else
           f(3*j-2) = f(3*j-2) - dedx
           f(3*j-1) = f(3*j-1) - dedy
           f(3*j  ) = f(3*j  ) - dedz
         end if
      end do !k=1,icount

  !---- End first outer loop ----

      if (qmmm_nml%ifqnt) then
        if (qmmm_struct%mm_link_mask(i)) then
          do qm_temp_count = 1,qmmm_struct%nlink
            if (qmmm_struct%link_pairs(1,qm_temp_count) == i) exit
          end do
          qm_temp_count = qm_temp_count + qmmm_struct%nquant
          qmmm_struct%dxyzqm(1,qm_temp_count) = qmmm_struct%dxyzqm(1,qm_temp_count) - dumx
          qmmm_struct%dxyzqm(2,qm_temp_count) = qmmm_struct%dxyzqm(2,qm_temp_count) - dumy
          qmmm_struct%dxyzqm(3,qm_temp_count) = qmmm_struct%dxyzqm(3,qm_temp_count) - dumz
        else
          !Not a link atom, can just add to main force array.
          f(3*i-2) = f(3*i-2) + dumx
          f(3*i-1) = f(3*i-1) + dumy
          f(3*i  ) = f(3*i  ) + dumz
        end if
      else
        f(3*i-2) = f(3*i-2) + dumx
        f(3*i-1) = f(3*i-1) + dumy
        f(3*i  ) = f(3*i  ) + dumz
      end if
      iexcl = iexcl + numex(i)
   end do  !  i=1,maxi
   call timer_stop(TIME_GBFRC)

   if( igb==6 ) goto VACUUM3

   !--------------------------------------------------------------------------
   !
   !    Step 3:  Finally, do the reduction over the sumdeijda terms:, adding
   !             into the forces those terms that involve derivatives of
   !             the GB terms (including the diagonal or "self" terms) with
   !             respect to the effective radii.  This is done by computing
   !             the vector dai/dxj, and using the chain rule with the
   !             previously-computed sumdeijda vector.
   !
   !             Also, compute a surface-area dependent term if igbsa=1
   !
   !             Do these terms only at "nrespa" multiple-time step intervals;
   !             (when igb=2 or 5, one may need to do this at every step)
   !
   !--------------------------------------------------------------------------

   if( onstep ) then
      count=0
      frespa = nrespa
      call timer_start(TIME_GBRAD2)

      !  -- diagonal epol term, plus off-diag derivs wrt alpha == reff^-1:

      do i=1,maxi

         f_xi = zero
         f_yi = zero
         f_zi = zero
         qi = charge(i)

         expmkf = exp( -kappa * reff(i) )*extdieli
         dl = intdieli - expmkf
         qi2h = half*qi*qi
         qid2h = qi2h * dl

         temp1 = (onereff(i) + one_Arad_beta)
         self_e = qid2h*temp1
         if ( ifcr /= 0 ) &
            call cr_add_dcdr_factor( i, -qi*temp1*dl )

         temp7 = -sumdeijda(i) + qid2h - &
               kappa*qi2h*expmkf*reff(i)*(one+one_Arad_beta*reff(i))

         epol = epol - self_e

         if(oncpstep .or. oncestep) then
            dvdl = dvdl - (half*dl*dcharge(i)*dcharge(i)-qid2h)*onereff(i)
         end if

         ! qi2h = qiqj/2
         ! qid2h = qiqj/2*[1/Ein - exp[-kappa*effbornrad]/Eout]
         ! temp7 = ... + qiqj/2*[1/Ein - exp[-kappa*effbornrad]/Eout]
         !         -kappa*qiqj/2*exp[-kappa*effbornrad]/Eout*effbornrad
         ! temp7 without the -sumdeijda part is the diagonal gradient.

         ! carlos: moved temp7 calculation up into the LES region above

         xi = x(3*i-2)
         yi = x(3*i-1)
         zi = x(3*i)
         ri = rborn(i)-offset
         ri1i = one/ri
         iaci = ntypes * (iac(i) - 1)

         if( igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) then

            !  --- new onufriev: we have to later scale values by a
            !      alpha,beta,gamma -dependent factor:

            rborn_i = rborn(i)
            ri = rborn_i - offset

            psi_i = psi(i)
            gba_i = gbvalpha(i) ! take gbalpha_i,gbbeta_i,gbgamma_i
            gbb_i = gbvbeta(i)
            gbg_i = gbvgamma(i)
            thi   = tanh( (gba_i + gbg_i*psi_i*psi_i - gbb_i*psi_i )*psi_i )
            thi2  = (gba_i + three*gbg_i*psi_i*psi_i - &
                    two*gbb_i*psi_i )*(one - thi*thi)*ri/rborn_i
         end if

         icount = 0
         do j=1,natom
            if( i /= j ) then
               xij = xi - x(3*j-2)
               yij = yi - x(3*j-1)
               zij = zi - x(3*j  )
               r2 = xij*xij + yij*yij + zij*zij
               if( r2 <= rgbmaxpsmax2 ) then
                  ! pairlist contains only atoms within rgbmax + safety margin
                    icount = icount + 1
                     temp_jj(icount) = j
                     r2x(icount) = r2
              end if ! r2 <= rgbmaxpsmax2
            end if !i/=j
         end do
         call vdinvsqrt( icount, r2x, vectmp1 )

         kk1 = 0
         do k=1,icount
            j = temp_jj(k)
            r2 = r2x(k)
            sj =  fs(j)

            dij1i = vectmp1(k)
            dij = r2*dij1i
            sj2 = sj * sj

            if( dij <= four*sj ) then
              kk1 = kk1 + 1
              vectmp3(kk1) = dij + sj
              if( dij > ri+sj ) then
                 vectmp2(kk1) = r2 - sj2
                 vectmp4(kk1) = dij - sj
              else if ( dij > abs(ri-sj) ) then
                 vectmp2(kk1) = dij + sj
                 vectmp4(kk1) = ri
              else if ( ri < sj ) then
                 vectmp2(kk1) = r2 - sj2
                 vectmp4(kk1) = sj - dij
              else
                 vectmp2(kk1) = one
                 vectmp4(kk1) = one
              end if
            end if !dij <= four*sj
         end do

         call vdinv( kk1, vectmp2, vectmp2 )
         call vdinv( kk1, vectmp3, vectmp3 )
         vectmp4(1:kk1) = vectmp4(1:kk1)*vectmp3(1:kk1)
         call vdln( kk1, vectmp4, vectmp4 )

         kk1 = 0
         do k=1,icount
            j = temp_jj(k)
            j3 = 3*j
            r2 = r2x(k)
            xij = xi - x(j3-2)
            yij = yi - x(j3-1)
            zij = zi - x(j3  )

            dij1i = vectmp1(k)
            dij = r2*dij1i
            sj =  fs(j)
            if (dij <= rgbmax +sj) then
              sj2 = sj * sj

              !           datmp will hold (1/r)(dai/dr):

              dij2i = dij1i*dij1i
              dij3i = dij2i*dij1i

              if (dij > rgbmax - sj ) then

                 temp1 = 1.0d0/(dij-sj)
                 datmp = eighth * dij3i * ((r2 + sj2) * &
                         (temp1*temp1 - rgbmax2i) - two * log(rgbmax*temp1))

              else if( dij > four*sj ) then

                 tmpsd = sj2*dij2i
                 dumbo = te+tmpsd* (tf+tmpsd* (tg+tmpsd* (th+tmpsd* thh)))
                 datmp = tmpsd*sj*dij2i*dij2i*dumbo

                 !     ---check accuracy of above Taylor series:
                 !      kk1 = kk1 + 1
                 !      datmp2 = vectmp2(kk1)*sj*(-Half*dij2i + vectmp2(kk1)) +
                 !              Fourth*dij3i*vectmp4(kk1)
                 !      if( abs( datmp/datmp2 - 1.d0 ) .gt. 0.00001) &
                 !             write(6,*) i,j, datmp, datmp2

              else if( dij > ri+sj ) then

                 kk1 = kk1 + 1
                 datmp = vectmp2(kk1)*sj*(-half*dij2i + vectmp2(kk1)) + &
                       fourth*dij3i*vectmp4(kk1)

              else if ( dij > abs(ri-sj) ) then
                 kk1 = kk1 + 1
                 datmp = -fourth*(-half*(r2 - ri*ri + sj2)*dij3i*ri1i*ri1i &
                       + dij1i*vectmp2(kk1)*(vectmp2(kk1) - dij1i) &
                       - dij3i*vectmp4(kk1) )

              else if ( ri < sj ) then
                 kk1 = kk1 + 1
                 datmp = -half*(sj*dij2i*vectmp2(kk1) &
                       - two*sj*vectmp2(kk1)*vectmp2(kk1) &
                       - half*dij3i*vectmp4(kk1) )

              else
                 kk1 = kk1 + 1
                 datmp = zero
              end if  ! ( dij > 4.d0*sj )

              if( (igb == 7 .or. igb ==8) .and. dij < rborn(i) +rborn(j) + GBNECKCUT) then

                 ! (no changes needed for LES)

                 ! Derivative of neck with respect to dij is:
                 !                     5
                 !              9 mdist
                 !   (2 mdist + --------) neckMaxVal gbneckscale
                 !                 5
                 ! -(------------------------)
                 !                        6
                 !             2   3 mdist  2
                 !   (1 + mdist  + --------)
                 !                    10

                 mdist = dij - neckMaxPos(neckidx(i),neckidx(j))
                 mdist2 = mdist * mdist
                 mdist3 = mdist2 * mdist
                 mdist5 = mdist2 * mdist3
                 mdist6 = mdist3 * mdist3

                 ! temp1 will be divisor of above fraction * dij
                 ! (datmp is deriv * 1/r)

                 temp1 = 1 + mdist2 + (three/ten)*mdist6
                 temp1 = temp1 * temp1 * dij

                 ! (Note "+" means subtracting derivative, since above
                 !     expression has leading "-")

                 datmp = datmp + ((2 * mdist + (nine/five) * mdist5) &
                       *neckMaxVal(neckidx(i),neckidx(j))*gbneckscale)/temp1
              end if !f( (igb == 7 .or. igb == 8) .and. dij < rborn(i) +rborn(j) + GBNECKCUT)

              datmp = -datmp*frespa*temp7

              if( igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) datmp = datmp*thi2
              f_x = xij*datmp
              f_y = yij*datmp
              f_z = zij*datmp
              if (qmmm_nml%ifqnt) then
                if (qmmm_struct%mm_link_mask(j)) then
                  do qm_temp_count = 1,qmmm_struct%nlink
                    if (qmmm_struct%link_pairs(1,qm_temp_count) == j) exit
                  end do
                  qm_temp_count = qm_temp_count + qmmm_struct%nquant
                  qmmm_struct%dxyzqm(1,qm_temp_count) = qmmm_struct%dxyzqm(1,qm_temp_count) - f_x
                  qmmm_struct%dxyzqm(2,qm_temp_count) = qmmm_struct%dxyzqm(2,qm_temp_count) - f_y
                  qmmm_struct%dxyzqm(3,qm_temp_count) = qmmm_struct%dxyzqm(3,qm_temp_count) - f_z
                else
                  !Not a link atom, can just add to main force array.
                  f(j3-2) = f(j3-2) + f_x
                  f(j3-1) = f(j3-1) + f_y
                  f(j3  ) = f(j3  ) + f_z
                end if
              else
                f(j3-2) = f(j3-2) + f_x
                f(j3-1) = f(j3-1) + f_y
                f(j3  ) = f(j3  ) + f_z
              end if

              f_xi = f_xi - f_x
              f_yi = f_yi - f_y
              f_zi = f_zi - f_z

           end if ! (dij <= rgbmax +sj)
         end do  !  k=1,icount

         if (qmmm_nml%ifqnt) then
           if (qmmm_struct%mm_link_mask(i)) then
             do qm_temp_count = 1,qmmm_struct%nlink
               if (qmmm_struct%link_pairs(1,qm_temp_count) == i) exit
             end do
             qm_temp_count = qm_temp_count + qmmm_struct%nquant
             qmmm_struct%dxyzqm(1,qm_temp_count) = qmmm_struct%dxyzqm(1,qm_temp_count) - f_xi
             qmmm_struct%dxyzqm(2,qm_temp_count) = qmmm_struct%dxyzqm(2,qm_temp_count) - f_yi
             qmmm_struct%dxyzqm(3,qm_temp_count) = qmmm_struct%dxyzqm(3,qm_temp_count) - f_zi
           else
             !Not a link atom, can just add to main force array.
             f(3*i-2) = f(3*i-2) + f_xi
             f(3*i-1) = f(3*i-1) + f_yi
             f(3*i  ) = f(3*i  ) + f_zi
           end if
         else
           f(3*i-2) = f(3*i-2) + f_xi
           f(3*i-1) = f(3*i-1) + f_yi
           f(3*i  ) = f(3*i  ) + f_zi
         end if

         !  --- Define neighbor list ineighbor for calc of LCPO areas ---

         if ( gbsa > 0 ) then
            do k=1,icount
               j = temp_jj(k)
               dij = sqrt(r2x(k))
               if ((vdwrad(i) + vdwrad(j)) > dij ) then

                  ! Consider all atoms for icosa, only non-H's for LCPO:

                  if ( (gbsa == 2) .or. &
                       (vdwrad(i) > 2.5).and.(vdwrad(j) > 2.5) ) then
                     count=count+1
                     ineighbor(count) = j
                  end if
               end if
            end do
            count=count+1
            ineighbor(count)=0
         end if

         if ( gbsa == 2 ) then

            !  --- Calc surface area with icosasurf algo;
            !  does not provide forces ==> only for single point calculations

            if(i == 1) then
               ineighborpt = 1
               call icosa_init(2, 3, zero)
            end if

          !  if(i == 1) then
          !     ineighborpt = 1
          !     call icosa_init(2, 3, zero)
          !  end if
            totsasa = totsasa + icosa_sphere_approx(i,x, &
                   vdwrad,ineighborpt,ineighbor)

         end if  !  ( gbsa == 2 )

      end do   ! end loop over atom i

      call timer_stop(TIME_GBRAD2)

      if ( gbsa == 1 ) then

         !       --- calculate surface area by computing LCPO (Still) over
         !           all atoms ---

#        include "gbsa.h"

      end if  !  ( gbsa == 1 )
      esurf = surften*totsasa
   end if  !  onstep


VACUUM3 &

!======== QMMM ==========
   if (qmmm_nml%ifqnt) then
     call timer_start(TIME_QMMM)
     if (qmmm_nml%qmgb /= 0) then

        ! We filled the main charge array with charges for QM atoms; make sure
        ! we zero it again so the 1-4's which are done outside of EGB will be
        ! correctly skipped for QM-MM on the next step. Note here we don't have the routine
        ! save the charges again since depending on the GB option they may have been filled
        ! with either RESP charges or mulliken charges.
        call qm_zero_charges(charge,qmmm_struct%scaled_mm_charges,.false.)
        if (qmmm_struct%nlink > 0 ) then
           !Don't save the charges here since they could be the resp charges.
           call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink,qmmm_struct%link_pairs,charge, &
                                              qmmm_struct%scaled_mm_charges,.false.)
        end if
     end if

     call timer_start(TIME_QMMMCOLLATEF)
     !We need to restore the MM link pair coordinates and then
     !Use the chain rule to put the link pair forces back onto the
     !QM and MM link pairs.
     call rst_mm_link_pair_crd(x)
     do i=1,qmmm_struct%nlink
       mm_no = 3*qmmm_struct%link_pairs(1,i)-2  !location of atom in x array
       lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom bound to link atom
       qm_no = 3*qmmm_struct%iqmatoms(lnk_no)-2
       !Note this routine uses the flink in the form -flink.
       call distribute_lnk_f(forcemod,qmmm_struct%dxyzqm(1:3,qmmm_struct%nquant+i),x(mm_no), &
                             x(qm_no),qmmm_nml%lnk_dis)

       !NOTE: forces are reversed in QM calc with respect to amber force array
       !so we subtract forcemod from MM atom and add it to QM atom.
       j = (qmmm_struct%link_pairs(1,i)-1)*3 !Natom number of MM link pair.
       !MM atom's new force = FMM(x,y,z) - FORCEMOD(x,y,z)
       f(j+1) = f(j+1) - forcemod(1)
       f(j+2) = f(j+2) - forcemod(2)
       f(j+3) = f(j+3) - forcemod(3)

       j = (qmmm_struct%iqmatoms(lnk_no)-1)*3
       !QM atom's new force = FQM(x,y,z) - Flink(x,y,z) + FORCEMOD(x,y,z)
       !Note QM forces should be subtracted from sander F array to leave total force.
       f(j+1) = f(j+1) - qmmm_struct%dxyzqm(1,qmmm_struct%nquant+i) + forcemod(1)
       f(j+2) = f(j+2) - qmmm_struct%dxyzqm(2,qmmm_struct%nquant+i) + forcemod(2)
       f(j+3) = f(j+3) - qmmm_struct%dxyzqm(3,qmmm_struct%nquant+i) + forcemod(3)
     end do
     call timer_stop_start(TIME_QMMMCOLLATEF,TIME_QMMMENERGY)
     call timer_start(TIME_QMMMGBENERGY)
     !Finally we need to calculate the VDW terms
     !for the MM link pair atoms as they were skipped above.

!------------ START MM LINK PAIR VDW -------------------
     do i = 1, qmmm_struct%nlink
       skipv(1:natom) = .false.
       jexcl = 1
       mm_no = qmmm_struct%link_pairs(1,i)
       iminus = mm_no-1
       do j = 1, iminus
         !As we loop over each atom check if this atom is supposed to exclude this MM
         !link pair atom.
         jexcl_last = jexcl + numex(j) -1
         do jjv = jexcl, jexcl_last
           iexcl = natex(jjv)
           !Should this atom exclude the MM atom?
           if (iexcl == mm_no) skipv(j) = .true.
         end do
         jexcl = jexcl + numex(j)
       end do
       jexcl_last = jexcl + numex(mm_no) -1
       do jjv = jexcl, jexcl_last
         skipv(natex(jjv)) = .true.
       end do
       !We need to avoid double counting so we should exclude ourself and
       !all link atoms of lower number
       !Note this assumes the link atom list is numerically sorted...
       do k = 1, i !To i ensures we skip the self interaction.
         skipv(qmmm_struct%link_pairs(1,k)) = .true.
       end do
       xi = x(3*mm_no-2)
       yi = x(3*mm_no-1)
       zi = x(3*mm_no  )

       iaci = ntypes * (iac(mm_no) - 1)
       dumx = zero
       dumy = zero
       dumz = zero
       de = zero
       do j=1,maxi !1 to natom

         !Do all atoms that are not excluded and
         !are within the cutoff for this i. Skip other
         !MM interactions that are less than i to avoid double
         !counting.
         if ( .not. skipv(j) ) then
           xij = xi - x(3*j-2)
           yij = yi - x(3*j-1)
           zij = zi - x(3*j  )
           r2 = xij*xij + yij*yij + zij*zij
           if ( r2 <= cut ) then
             r2inv = one/r2
             ic = ico( iaci + iac(j) )
             if (ic > 0) then
               !6-12 potential:
               r6inv = r2inv*r2inv*r2inv
               f6 = cn2(ic)*r6inv
               f12 = cn1(ic)*(r6inv*r6inv)
               evdw = evdw + (f12 - f6)
               de = (twelve*f12 - six*f6)*r2inv
             end if !if ic>0
             de = de*nrespai
             dedx = de * xij
             dedy = de * yij
             dedz = de * zij
             dumx = dumx + dedx
             dumy = dumy + dedy
             dumz = dumz + dedz
             f(3*j-2) = f(3*j-2) - dedx
             f(3*j-1) = f(3*j-1) - dedy
             f(3*j  ) = f(3*j  ) - dedz
           end if !if r2<=cut
         end if !if not skipv(i)
       end do

       f(3*mm_no-2) = f(3*mm_no-2) + dumx
       f(3*mm_no-1) = f(3*mm_no-1) + dumy
       f(3*mm_no  ) = f(3*mm_no  ) + dumz
     end do ! i = 1, qmmm_struct%nlink
     call timer_stop(TIME_QMMMGBENERGY)
     call timer_stop(TIME_QMMMENERGY)
!------------ END MM LINK PAIR VDW ---------------------
     call timer_stop(TIME_QMMM)
   end if !if ifqnt.
!======== END QMMM ==========

   return
end subroutine egb

subroutine egb_calc_radii(igb,natom,x,fs,reff,onereff,fsmax,rgbmax, &
                          rborn, offset,rbornstat, &
                          rbave,rbfluct,rbmax,rbmin, gbneckscale, ncopy, rdt, &
                          gbalpha,gbbeta,gbgamma )

! Calculates effective GB radii and puts result in reff;
! onereff contains 1.0d0/reff.

  use constants, only : zero, eighth, fourth, third, half, &
                        one, two, three, four, five, seven, eight, nine, &
                        eleven, twelve, thirtieth
  implicit none

#include "gbneck.h"

   ! Scratch variables used for calculating neck correction
   _REAL_ mdist,mdist2,mdist3,mdist6,neck

  integer, intent(in) :: igb, natom, rbornstat, ncopy
  _REAL_, intent(in) :: x(3*natom), fs(natom), rdt

  _REAL_, intent(out) :: reff(natom), onereff(natom)

  _REAL_, intent(in) :: fsmax, rgbmax, gbneckscale
  _REAL_, intent(in) :: rborn(natom),offset

  _REAL_, intent(in) :: gbalpha(natom), gbbeta(natom), gbgamma(natom)
  _REAL_ :: reff_i
  integer :: ncopy2
  _REAL_, intent(out) :: rbave(natom), rbfluct(natom), rbmax(natom),rbmin(natom)

! Local:
  integer :: i, icount, iplus, kk1, kk2, k, j, vecend
  _REAL_ :: xi,yi,zi,xij, yij, zij,r2,ri,ri1i,rj1i,si,si2,sj, sj2
  _REAL_ :: dij1i, dij2i, dij, rj, uij, tmpsd, dumbo, theta
  _REAL_ :: rgbmax1i, rgbmax2i, rgbmaxpsmax2
  _REAL_ :: temp3,temp4

   !   FGB taylor coefficients follow
   !   from A to D :
   !   1/3 , 2/5 , 3/7 , 4/9 , 5/11
   _REAL_  ta
   _REAL_  tb
   _REAL_  tc
   _REAL_  td
   _REAL_  tdd
   parameter ( ta   = third )
   parameter ( tb   = two / five )
   parameter ( tc   = three / seven )
   parameter ( td   = four / nine )
   parameter ( tdd  = five / eleven )

   rgbmax1i = one/rgbmax
   rgbmax2i = rgbmax1i*rgbmax1i
   rgbmaxpsmax2 = (rgbmax+fsmax)**2

   onereff(1:natom) = zero
   ncopy2 = ncopy ! BUGBUG -- only here to suppress warnings

   do i=1,natom
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i)
      reff_i = onereff(i)

      ri = rborn(i)-offset
      ri1i = one/ri
      si = fs(i)
      si2 = si*si

      !  Here, reff_i will sum the contributions to the inverse effective
      !  radius from all of the atoms surrounding atom "i"; later the
      !  inverse of its own intrinsic radius will be added in

      icount = 0
      iplus=i+1

      do j=iplus,natom

         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )
         r2 = xij*xij + yij*yij + zij*zij

         if( r2 <= rgbmaxpsmax2 ) then
           icount = icount + 1
           temp_jj(icount) = j
           r2x(icount) = r2
         end if

      end do

      call vdinvsqrt( icount, r2x, vectmp1 )

      kk1 = 0
      kk2 = 0
!!!      !dir$ ivdep
      do k = 1, icount

         j = temp_jj(k)
         r2 = r2x(k)
         sj =  fs(j)

!        don't fill the remaining vectmp arrays if atoms don't see each other:
         dij1i = vectmp1(k) !1/sqrt(r^2)
         dij = r2*dij1i != rij
         if (dij <= rgbmax+si .or. dij <= rgbmax+sj) then
           rj = rborn(j) - offset

           if( dij <= four*sj) then
              kk1 = kk1 + 1
              vectmp2(kk1) = dij + sj
              if( dij > ri+sj ) then
                 vectmp4(kk1) = dij - sj
              else if ( dij > abs(ri-sj) ) then
                 vectmp4(kk1) = ri
              else if ( ri < sj ) then
                 vectmp4(kk1) = sj - dij
              else
                 vectmp4(kk1) = one
              end if
           end if

           if( dij <= four*si) then
              kk2 = kk2 + 1
              vectmp3(kk2) = dij + si
              if( dij > rj+si) then
                 vectmp5(kk2) = dij - si
              else if ( dij > abs(rj-si) ) then
                 vectmp5(kk2) = rj
              else if ( rj < si ) then
                 vectmp5(kk2) = si - dij
              else
                 vectmp5(kk2) = one
              end if
           end if
        end if ! (dij <= rgbmax+si .or. dij <= rgbmax+sj)
      end do  !  k = 1, icount

      call vdinv( kk1, vectmp2, vectmp2 )
      call vdinv( kk2, vectmp3, vectmp3 )
      vectmp4(1:kk1) = vectmp2(1:kk1)*vectmp4(1:kk1)
      vectmp5(1:kk2) = vectmp3(1:kk2)*vectmp5(1:kk2)
      call vdln( kk1, vectmp4, vectmp4 )
      call vdln( kk2, vectmp5, vectmp5 )

      kk1 = 0
      kk2 = 0
      do k = 1, icount

         j = temp_jj(k)
         r2 = r2x(k)
         rj = rborn(j) - offset
         rj1i = one/rj
         sj =  fs(j)

         sj2 = sj * sj

         xij = xi - x(3*j-2)
         yij = yi - x(3*j-1)
         zij = zi - x(3*j  )

         dij1i = vectmp1(k)
         dij = r2*dij1i

         temp3 = zero
         temp4 = zero

         if (dij <= rgbmax + sj) then

            if ((dij > rgbmax - sj)) then
               uij = 1.0d0/(dij -sj)

! carlos: store descreening contrib in temp3, this makes it easier to
! apply to multiple atoms for LES

               temp3 = temp3  - eighth * dij1i * (one + two * dij *uij + &
                     rgbmax2i * (r2 - four * rgbmax * dij - sj2) + &
                     two * log((dij-sj)*rgbmax1i))

             else if( dij > four*sj ) then

               dij2i = dij1i*dij1i
               tmpsd = sj2*dij2i
               dumbo = ta+tmpsd* (tb+tmpsd* (tc+tmpsd* (td+tmpsd* tdd)))

               temp3 = temp3 - tmpsd*sj*dij2i*dumbo

            !     ---following are from the Appendix of Schaefer and Froemmel,
            !        J. Mol. Biol. 216:1045-1066, 1990, divided by (4*Pi):

            else if( dij > ri+sj ) then

               kk1 = kk1 + 1
               temp3 = temp3 - half*( sj/(r2-sj2) + half*dij1i*vectmp4(kk1) )

            !-----------------------------------------------------------------

            else if ( dij > abs(ri-sj) ) then

               kk1 = kk1 + 1
               theta = half*ri1i*dij1i*(r2 + ri*ri -sj2)
               temp3 = temp3- fourth*( ri1i*(two-theta) &
                  - vectmp2(kk1) + dij1i*vectmp4(kk1) )

            !-----------------------------------------------------------------

            else if ( ri < sj ) then
               kk1 = kk1 + 1
               temp3 = temp3- half*sj/(r2-sj2) + ri1i    &
                  + fourth*dij1i*vectmp4(kk1)

            !-----------------------------------------------------------------

            else
               kk1 = kk1 + 1
            end if  ! ( dij > 4.d0*sj )

            if ( (igb == 7 .or. igb ==8) .and. dij < rborn(i) +rborn(j) + GBNECKCUT) then
               mdist = dij - neckMaxPos(neckidx(i),neckidx(j))
               mdist2 = mdist *mdist
               mdist3 = mdist2 * mdist
               mdist6 = mdist3 *mdist3
               neck = neckMaxVal(neckidx(i),neckidx(j))/ &
                       (one + mdist2 +0.3d0*mdist6)
               temp3 = temp3 - gbneckscale * neck
            end if

         end if

         !   --- Now the same thing, but swap i and j and use temp4 for descreening contrib:

         if (dij <= rgbmax +si) then

           if (dij > rgbmax - si) then
              uij = 1.0d0/(dij -si)
              temp4 = temp4 - eighth * dij1i * (one + two * dij *uij + &
                       rgbmax2i * (r2 - four * rgbmax * dij - si2) + &
                       two * log((dij-si)*rgbmax1i))
           else if( dij > four*si ) then
              dij2i = dij1i*dij1i
              tmpsd = si2*dij2i
              dumbo = ta+tmpsd* (tb+tmpsd* (tc+tmpsd* (td+tmpsd* tdd)))
              temp4 = temp4 - tmpsd*si*dij2i*dumbo
           else if( dij > rj+si ) then
              kk2 = kk2 + 1
              temp4 = temp4 - half*( si/(r2-si2) + &
                    half*dij1i*vectmp5(kk2) )
            !-----------------------------------------------------------------
           else if ( dij > abs(rj-si) ) then
              kk2 = kk2 + 1
              theta = half*rj1i*dij1i*(r2 + rj*rj -si2)
              temp4 = temp4 - fourth*( rj1i*(two-theta) &
                    - vectmp3(kk2) + dij1i*vectmp5(kk2) )
            !-----------------------------------------------------------------
           else if ( rj < si ) then
              kk2 = kk2 + 1
              temp4 = temp4 - half*si/(r2-si2) + rj1i &
                    + fourth*dij1i*vectmp5(kk2)
            !-----------------------------------------------------------------
           else
              kk2 = kk2 + 1
           end if  ! ( dij > 4.d0*si )

           if ((igb == 7 .or. igb ==8) .and. dij < rborn(j) +rborn(i) +GBNECKCUT) then
              mdist = dij - neckMaxPos(neckidx(j),neckidx(i))
              mdist2 = mdist * mdist
              mdist3 = mdist2 * mdist
              mdist6 = mdist3 * mdist3
              neck = neckMaxVal(neckidx(j),neckidx(i))/ &
                      (one + mdist2 +0.3d0*mdist6)
              temp4 = temp4 - gbneckscale *neck
           end if

        end if !(dij <= rgbmax +si)
! now add the calculated descreening component to onereff

            reff_i = reff_i + temp3
            onereff(j) = onereff(j) + temp4

      end do                    !  k = 1, icount


      ! we are ending the do-i-loop, reassign the scalar to the original array:

      onereff(i) = reff_i

   end do  !  i = 1,natom

! set end variable since LES arrays are longer than normal
! this will make the changes for LES code simpler to understand and modify

   vecend = natom

   if( igb == 2 .or. igb == 5 .or. igb == 7 .or. igb == 8) then

      ! --- apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:


      vectmp1(1:vecend) = rborn(1:vecend)
      vectmp2(1:vecend) = vectmp1(1:vecend)-offset
      call vdinv(vecend, vectmp1, vectmp1) !1.0d0/rborn

! use of onereff here means that we need changes for LES to accomodate multiple onereff per atom

      psi(1:vecend) = -vectmp2(1:vecend)*onereff(1:vecend)
      call vdinv(vecend, vectmp2, vectmp2) !1.0d0/(rborn-offset)

  !if not LES: use gbalpha, gbgamma, gbbeta arrays instead
      vectmp3(1:vecend) = ((gbalpha(1:vecend)+gbgamma(1:vecend)*psi(1:vecend)*   &
                            psi(1:vecend)-gbbeta(1:vecend)*psi(1:vecend))*psi(1:vecend))
      call vdtanh(vecend, vectmp3, vectmp3)

      ! vectmp1 = 1.0d0/rborn
      ! vectmp2 = 1.0d0/(rborn-offset)
      ! vectmp3 = tanh( (gbalpha + gbgamma*psi_i*psi_i - gbbeta*psi_i )*psi_i )

      onereff(1:vecend) = vectmp2(1:vecend) - &
                            (vectmp3(1:vecend)*vectmp1(1:vecend))
      do j =1,vecend
         if ( onereff(j) < zero ) onereff(j) = thirtieth
      end do
      call vdinv(vecend, onereff, reff)
   else

      !       "standard" GB, including the "diagonal" term here:

      vectmp1(1:vecend) = rborn(1:vecend)-offset

      call vdinv(vecend, vectmp1, vectmp1)
      onereff(1:vecend) = onereff(1:vecend) + vectmp1(1:vecend)
      call vdinv(vecend, onereff, reff)
   end if

   REQUIRE( rdt == 0 )

   if ( rbornstat == 1 ) then

      do k=1,natom
         i=k
         rbave(i) = rbave(i) + reff(i)
         rbfluct(i) = rbfluct(i) + reff(i)*reff(i)
         if ( rbmax(i) <= reff(i) ) rbmax(i) = reff(i)
         if ( rbmin(i) >= reff(i) ) rbmin(i) = reff(i)
      end do
   end if

   return

end subroutine egb_calc_radii

! checking if an atom belongs to nucleic or not
! make sure to add new DNA/RNA residue name to this list
! anyway not to hardcoded nucnamenum? (=60)
! tag: gbneck2nu
! TODO: use DNA/RNA some where?
subroutine isnucat(nucat,atomindex,nres,nucnamenum,ipres,lbres)
    implicit none
    integer, intent(in) :: atomindex,nres,nucnamenum,ipres(nres)
    integer :: j,k,nucat
    character(4) :: nucname(nucnamenum),lbres(nres)
    !nucnamenum = 60
    nucat = 0
    nucname = (/ &
     "A   ", & !1
     "A3  ", &
     "A5  ", &
     "AN  ", &
     "C   ", & !5
     "C3  ", &
     "C5  ", &
     "CN  ", &
     "DA  ", &
     "DA3 ", & !10
     "DA5 ", &
     "DAN ", &
     "DC  ", &
     "DC3 ", &
     "DC5 ", & !15
     "DCN ", &
     "DG  ", &
     "DG3 ", &
     "DG5 ", &
     "DGN ", & !20
     "DT  ", &
     "DT3 ", &
     "DT5 ", &
     "DTN ", &
     "G   ", & !25
     "G3  ", &
     "G5  ", &
     "8OG ", &
     "GN  ", &
     "U   ", & !30
     "U3  ", &
     "U5  ", &
     "UN  ", &
     "AP  ", &
     "DAP ", & !35
     "CP  ", &
     "AF2 ", &
     "AF5 ", &
     "AF3 ", &
     "GF2 ", & !40
     "GF3 ", &
     "GF5 ", &
     "CF2 ", &
     "CFZ ", &
     "CF3 ", & !45
     "CF5 ", &
     "UF2 ", &
     "UF5 ", &
     "UF3 ", &
     "UFT ", & !50
     "OMC ", &
     "OMG ", &
     "MC3 ", &
     "MC5 ", &
     "MG3 ", & !55
     "MG5 ", &
     "SIC ", &
     "NAM ", &
     "NM5 ", &
     "CAP "  /)!60

    do j=1,nres-1
        !find residue to which atom i belongs
        if (ipres(j)  <= atomindex .and. atomindex < ipres(j+1)) then
            do k=1,nucnamenum
                !if atom i belongs to nucleic, return 1
                if (lbres(j)(1:4) == nucname(k)(1:4)) then
                    nucat = 1
                end if
            end do
        end if
    end do
    if (atomindex >= ipres(nres)) then
       do k=1,nucnamenum
           !if atom i belongs to nucleic, return 1
           if (lbres(j)(1:4) == nucname(k)(1:4)) then
               nucat = 1
           end if
       end do
    end if
end subroutine isnucat

end module genborn
