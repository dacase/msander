! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

module sgld

! MODULE: sgld
! ================== SELF-GUIDED MOLECULAR/LANGEVIN DYNAMICS ====================
! Xiongwu Wu, 2022

#ifdef MPI
   use mpi
#endif

implicit none

!     Head file for the self-guided Langevin Dynamics simulation  
!
!      variables for SGLD simulation
!


! ... integers:
!
!
!    SGMD/SGLD applying range
!      ISGSTA      Begining atom index applying SGLD
!      ISGEND      Ending atom index applying SGLD
!
!
      integer,save:: ISGSTA,ISGEND
!
! ... floats:
!
!
!    SGMD/SGLD VARIABLES
!     SGFT    !  Guiding factor 
!     TSGAVG  !  Local average time, ps
!     TSGAVP  !  Convergence time, ps
!
!
!     SGAVG0  !  Local average remains
!     SGAVG1  !  Local average factor, SGAVG1=1-SGAVG0
!     SGAVP0  !  Convergence average remains
!     SGAVP1  !  Convergency average factor, SGAVP1=1-SGAVP0
!     GAMMAS  !  friction coefficient
!     TEMPLF  !  Low frequency motion temperature 
!     TEMPHG  !  High frequency temperature
!

      _REAL_  SGFT,SGFF,SGFG,TSGAVG,TSGSET,TSGAVP,GAMMAS, &
          SGAVG0,SGAVG1,SGAVP0,SGAVP1,SGFTI,SGFFI,SGFGI,EPOTLF,EPOTHF,EPOTLLF, &
          TEMPLF,TEMPHF, &
          SGWT, sgrndf, &
          MYSCALSG,SGSCALE,FSGLDG,PSGLDG,tsgfac

!  common block for parallel broadcast
      common/sgldr/SGFT,SGFF,SGFG,TSGAVG,TSGSET,TSGAVP,GAMMAS, &
      SGAVG0,SGAVG1,SGAVP0,SGAVP1,SGFTI,SGFFI,SGFGI,EPOTLF,EPOTHF,EPOTLLF, &
      TEMPLF,TEMPHF, &
      SGWT, sgrndf, &
      MYSCALSG,SGSCALE,FSGLDG,PSGLDG,tsgfac
!  Number of broadcasting variables
integer, parameter :: nsgld_real=26
       
!
! ... flags:
!
!
!     isgld   ! input control; positive values activate SGLD
!     TSGBOND ! Perform guiding force average over bonded structure
!     TSGLDGLE ! Perform SGLD-GLE to mantain canonical ensemble distribution
!
!
      integer, save :: isgld,nsgsize
      LOGICAL, save, private :: TSGLD,TSGLDGLE,TSGBOND
      LOGICAL, save :: TRXSGLD
      character(8), save:: sglabel
      
!*******************************************************************************
!
! ...allocatable arrays:
!
!     avgx1     ! local averages of position
!     avgx2     ! local averages of local averages of position
!     avgp     ! local averages of momentum
!     avgr     ! local average of random forces
!     sgfps   !  averages of force-momentum product
!     sgpps   ! averages of momentum-momentum product
!*******************************************************************************

      double precision, dimension(:,:), allocatable, save :: avgx1,avgx2,avgp,avgr
      double precision, dimension(:), allocatable, private,save :: sgfps,sgpps,sgmass

      type listdata_rec
      integer             :: offset
      integer             :: cnt
      end type listdata_rec
  

      type(listdata_rec), allocatable, save :: atm_sg_maskdata(:)
      integer, allocatable, save            :: atm_sg_mask(:)
    
   
contains


    SUBROUTINE PSGLD(atm_cnt,numex,natex,AMASS,crd,vel, rem)
!-----------------------------------------------------------------------
!     This routine performs initiation for the Self-Guided        
!       Langevin Dynamcs (SGLD) simulaiton                 
!
      use md_scheme, only: gamma_ln
      implicit none
#include "../include/md.h"
#include "../include/memory.h"
#ifdef MPI
#  include "parallel.h"
#endif
      INTEGER atm_cnt,numex(*),natex(*)
      _REAL_ AMASS(*),crd(3,*),vel(3,*)
      integer rem
      INTEGER I,I3,M,ierror,j,nsgsubi,idx_nbex,jatm
      _REAL_ AMASSI,XI3,VI3,ekin,ekinsg,GAMM,FACT1,FACT2
      logical is_langevin  ! Is this a Langevin dynamics simulation
!
      is_langevin = gamma_ln > 0.0d0
      TSGLD = (isgld >0)
      tsgbond = (nsgsize > 1)
      tsgldgle = (abs(sgfg) > 1.0d-6).and.is_langevin
      trxsgld = rem > 0 .and. tsgld
!  Check for invalid sgld setting
      IF(ISGSTA < 1)ISGSTA=1
      IF(ISGEND > NATOM .OR. ISGEND < 1)ISGEND=NATOM
      IF(TSGAVG.LT.DT)TSGAVG=DT
      IF(TSGAVP.LT.DT)TSGAVP=10.0D0*TSGAVG
      SGAVG1=DT/TSGAVG
      SGAVG0=1.0D0-SGAVG1
      SGAVP1=DT/TSGAVP
      SGAVP0=1.0D0-SGAVP1
      gammas=gamma_ln/20.455d0
      tsgfac=1.0d0/20.455d0/tsgavg
      sglabel="SGMD: "
      if(is_langevin)sglabel="SGLD: "
      if(sgft<-1.0d0.or.sgft>1.0d0)then
        if(sgff<-1.0d0 .or. sgff>1.0d0)then
          sgfti=0.0d0
          sgffi=0.0d0
        else
          sgffi=sgff
          sgfti=(1.0d0+SGFFI)*(1.0d0+SGFFI)-1.0d0/(1.0d0+SGFFI)
        endif
      else
        sgfti=sgft
        if(sgff<-1.0d0 .or. sgff>1.0d0)then
          PSGLDG=0.0d0
          if(ABS(SGFTI)>1.0d-8)then
            FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
            FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
            PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
          endif
          sgffi=psgldg
        else
          sgffi=sgff
        endif
      endif
      PSGLDG=0.0d0
      if(ABS(SGFTI)>1.0d-8)then
        FACT1=9.0d0-SQRT(81.0d0-12.0d0*SGFTI*SGFTI*SGFTI)
        FACT2=(ABS(FACT1)*1.5d0)**(1.0/3.0d0)
        PSGLDG=SIGN(1.0d0,FACT1)*(FACT2/3.0d0+SGFTI/FACT2)-1.0d0
      endif
      if(tsgldgle)then
        sgfgi=sgfg
        FSGLDG=SQRT(1.0d0-SGFGI)-1.0d0
      else
        sgfgi=0.0d0
        FSGLDG=0.0d0
      endif
#ifdef MPI
      if(mytaskid.eq.0)THEN
#endif
      write(6,910)isgsta,isgend,tsgavg
      if(is_langevin)then
        if(tsgldgle)then
          write(6,928)
          write(6,927)sgfgi,fsgldg
        else
          write(6,940)
        endif
        write(6,930)gamma_ln
      else
          write(6,941)
      endif
      write(6,925)sgfti,psgldg
      write(6,926)sgffi
      if(tsgbond)then
        if(nsgsize==2)then
          write(6,942)
        else
          write(6,943)
        endif
      endif      
      write(6,935)
#ifdef MPI
    ENDIF
#endif
      tsgset=temp0
      if(tsgset<1.0d-6)tsgset=300.0d0
      !     allocate working arrays
      allocate( avgx1(3,natom),avgx2(3,natom),  &
      avgp(3,natom ),avgr(3,natom ),sgfps(natom ),sgpps(natom ),&
                stat=ierror)
      REQUIRE( ierror == 0 )
     ! build bidireectional  exclusion lists
      if(tsgbond)then 
        allocate(atm_sg_maskdata(natom), &
                 atm_sg_mask(nnb*2), &
                 sgmass(natom),stat = ierror)
        call make_sgavg_mask_list(natom,nnb,numex, natex)
        REQUIRE( ierror == 0 )
      endif
  !    Initialize arrays
      gamm=sqrt(dt/tsgavg)
      ekin=0.0d0
      ekinsg=0.0d0
       DO I=1,NATOM
        AMASSI = AMASS(I)
        !IF((I>=ISGSTA).AND.(I<=ISGEND))THEN
        !ENDIF
        DO M=1,3
          xi3=crd(m,i)
          vi3=vel(m,i)
          ekin=ekin+amassi*vi3*vi3
          if(tsgbond)then 
            nsgsubi=atm_sg_maskdata(i)%cnt
            xi3=amassi*xi3
            vi3=amassi*vi3
            sgmass(i)=amassi
            do j=1,nsgsubi
              idx_nbex=atm_sg_maskdata(i)%offset + j 
              jatm=atm_sg_mask(idx_nbex)
              xi3=xi3+amass(jatm)*crd(m,jatm)
              vi3=vi3+amass(jatm)*vel(m,jatm)
              sgmass(i)=sgmass(i)+amass(jatm)
            enddo
            xi3=xi3/sgmass(i)
            vi3=vi3/sgmass(i)
          endif
          avgx1(m,i)=xi3
          avgx2(m,i)=avgx1(m,i)
          avgp(m,i)=0.0d0
          avgr(m,i)=0.0d0
          ekinsg=ekinsg+amassi*vi3*vi3
        END DO
        sgfps(i)=0.0d0
        sgpps(i)=amassi*0.001987*tsgset*sgavg1  
      END DO
      epotlf=2.0d10
      epotllf=2.0d10
      templf=tsgset*gamm*ekinsg/ekin
      temphf=0.0d0
      sgwt=0.0d0
      910   format("  _________________ SGLD parameters _________________"/  &
      "  Parameters for self-guided Langevin dynamics (SGLD) simulation"//  &
          "  Guiding range from ",i5,"  to ",i5 /  &
          "  Local averaging time: ",f10.4," ps ")
925   format("  sgfti: ",f8.4," psgldg: ",f8.4)
926   format("  sgffi: ",f8.4)
927   format("  momentum factor sgfgi= ",f8.4," random force factor fsgldg=",f8.4)
928   format("  SGLD-GLE method is used to mantain a canonical distribution. ")
940   format("  SGLDg  method is used to enhance conformational search. ")
941   format("  SGMDg  method is used to enhance conformational search. ")
942   format("  NSGSIZE=2, Guiding forces are averaged over 1-2,1-3 bonded structures" )
943   format("  NSGSIZE>2, Guiding forces are averaged over 1-2,1-3,1-4 bonded structures" )
930   format("  Collision frequency:",f8.2," /ps" )
935   format("  SGMD/SGLD output properties:"    /  &
             "  SGLABEL:  SGGAMMA TEMPLF  TEMPHF  EPOTLF EPOTHF EPOTLLF SGWT" /  &
             "         SGMD/SGLD weighting factor=exp(SGWT)"/  &
              " _______________________________________________________"/)
      RETURN
      END SUBROUTINE PSGLD


      subroutine make_sgavg_mask_list(atm_cnt, nnb, numex, natex)

     
        implicit none
      
      ! Formal arguments:
      
        integer               :: atm_cnt,nnb
        ! Excluded atom count for each atom.
        integer               :: numex(atm_cnt)
        ! Excluded atom concatenated list:
        integer               :: natex(nnb)
      
      ! Local variables:
      
        integer               :: atm_i, atm_j
        integer               :: lst_idx, sublst_idx, num_sublst
        integer               :: mask_idx
        integer               :: offset
        integer               :: total_excl
      
      ! Double the mask to deal with our list generator
      
      ! Pass 1: get pointers, check size
      
        lst_idx = 0
      
        atm_sg_maskdata(:)%cnt = 0       ! array assignment
      
        do atm_i = 1, atm_cnt - 1         ! last atom never has any...
          num_sublst = numex(atm_i)
          do sublst_idx = 1, num_sublst
            atm_j = natex(lst_idx + sublst_idx)
            if (atm_j .gt. 0 .or. nsgsize>2) then
              atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
              atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
            end if
          end do
          lst_idx = lst_idx + num_sublst
        end do
      
        total_excl = 0
      
        do atm_i = 1, atm_cnt
          total_excl = total_excl + atm_sg_maskdata(atm_i)%cnt
        end do
        !write(6,*)"total_excl, nnb: ",total_excl, nnb
        if (total_excl .gt. nnb*2) then
          write(6, '(a,a)') "SGBOND: ", &
               'The total number of sg substructure exceeds that stipulated by the'
          write(6, '(a,a)') "SGBOND: ", &
               'prmtop.  This is likely due to a very high density of added extra points.'
          write(6, '(a,a)') "SGBOND: ", &
               'Scale back the model detail, or contact the developers for a workaround.'
          call mexit(6, 1)
        end if
      
        offset = 0
      
        do atm_i = 1, atm_cnt
          atm_sg_maskdata(atm_i)%offset = offset
          offset = offset + atm_sg_maskdata(atm_i)%cnt
        end do
      
      ! Pass 2: fill mask array
      
        lst_idx = 0
      
        atm_sg_maskdata(:)%cnt = 0       ! array assignment
        
          do atm_i = 1, atm_cnt - 1
            num_sublst = numex(atm_i)
            do sublst_idx = 1, num_sublst
              atm_j = natex(lst_idx + sublst_idx)
              if (atm_j .gt. 0 .or. nsgsize>2) then
                if(atm_j==0)cycle
                atm_j=abs(atm_j)
                atm_sg_maskdata(atm_j)%cnt = atm_sg_maskdata(atm_j)%cnt + 1
                mask_idx = atm_sg_maskdata(atm_j)%offset + &
                           atm_sg_maskdata(atm_j)%cnt
                atm_sg_mask(mask_idx) = atm_i
      
                atm_sg_maskdata(atm_i)%cnt = atm_sg_maskdata(atm_i)%cnt + 1
                mask_idx = atm_sg_maskdata(atm_i)%offset + &
                           atm_sg_maskdata(atm_i)%cnt
                atm_sg_mask(mask_idx) = atm_j
      
              end if
            end do
            lst_idx = lst_idx + num_sublst
          end do
        return
      
      end subroutine make_sgavg_mask_list
      
      
    subroutine sg_fix_degree_count(sgsta_rndfp, sgend_rndfp, ndfmin, rndf)
!-----------------------------------------------------------------------
!     Correct the total number of degrees of freedom for a translatable COM,
!       and compute the number of degrees of freedom in the SGLD part.
!       The latter is mostly done by the caller to avoid passing the long
!       argument list needed by routine degcnt which also differs between
!       sander and pmemd.
!
      implicit none
      _REAL_, intent(in)    :: sgsta_rndfp, sgend_rndfp
      integer, intent(in)   :: ndfmin
      _REAL_, intent(inout) :: rndf

      sgrndf = sgend_rndfp - sgsta_rndfp
      return
      end subroutine sg_fix_degree_count



      SUBROUTINE SGLDW(NATOM,ISTART,IEND, &
             DTX,TEMP0,ENER,AMASS,WINV,crd,Frc,Vel)
!-----------------------------------------------------------------------
!     This routine perform SGLD integration        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(20)
# endif
#endif
      INTEGER NATOM,ISTART,IEND
      _REAL_ DTX,TEMP0
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),crd(3,*),Frc(3,*),Vel(3,*)
!
      INTEGER I,M,JSTA,JEND,j,nsgsubi,idx_nbex,jatm
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC,GAM,RSD,FLN
      _REAL_ EKIN,EKINSG,SGBETA
      double precision sumgam,sumfp,sumpp,sumgv,sumpv
      double precision sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      double precision xi3,x1i3,x2i3,vi3t,vi3,fi3
      PARAMETER (BOLTZ = 1.987192d-3)
!
      gam=gammas*dtx
      JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        sumgam=0.0d0
        EKIN=0.0D0
        EKINSG=0.0D0
        DO  I = 1,NATOM 
          AMASSI = AMASS(I)
          WFAC =  2.0D0*DTX*WINV(I)
          RSD = SQRT(2.D0*GAMMAS*BOLTZ*TEMP0*AMASSI/DTX)
          IF(I>=JSTA.AND.I<=JEND)THEN
            ! sggamma
          sggammai=-sgfps(i)/sgpps(i)
          sumgam=sumgam+sggammai
          sumfp=0.0d0
          sumpp=0.0d0
          sumgv=0.0d0
          sumpv=0.0d0
          DO  M = 1,3
!   Keep random number series the same as that in a single cpu simulation
            CALL GAUSS( 0.D0, RSD, FLN )
              ! avg(x)
              xi3=crd(m,i)
              if(tsgbond)then 
                nsgsubi=atm_sg_maskdata(i)%cnt
                xi3=amassi*xi3
                do j=1,nsgsubi
                  idx_nbex=atm_sg_maskdata(i)%offset + j 
                  jatm=atm_sg_mask(idx_nbex)
                  xi3=xi3+amass(jatm)*crd(m,jatm)
                enddo
                xi3=xi3/sgmass(i)
              endif
              x1i3=sgavg0*avgx1(m,i)+sgavg1*xi3
              avgx1(m,i)=x1i3
              ! avgavg(x)
              x2i3=sgavg0*avgx2(m,i)+sgavg1*x1i3
              avgx2(m,i)=x2i3
              ! avg(p)
              avgpi3=tsgfac*amassi*(xi3-x1i3)
              pi3t=(avgpi3-sgavg0*avgp(m,i))/sgavg1
              avgp(m,i)=avgpi3
              ! avg(f-avg(f))
              avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
              ! sum(avg(f-avg(f))avg(p))
              sumfp=sumfp+avgdfi3*avgpi3
              ! sum(avg(p)avg(p))
              sumpp=sumpp+avgpi3*avgpi3
              ! average random forces
              avgri3=sgavg0*avgr(m,i)+sgavg1*fln
              avgr(m,i)=avgri3
              ! guiding forces
              fsgpi=(sgfti*sggammai)*avgpi3
              fsgfi=sgffi*avgdfi3
              fsgi3=fsgpi+fsgfi
              fi3=frc(m,i)+fln+fsgi3+sgfgi*gammas*avgpi3+fsgldg*avgri3
              !fi3=frc(m,i)+fln
              frc(m,i)=fi3
             ! estimate velocity at t+dt/2
              vi3t=vel(m,i)+0.5d0*dtx/amassi*fi3
              ! sum(g*v)
              sumgv=sumgv+fsgi3*vi3t
              ! sum(p*v)
              sumpv=sumpv+amassi*vi3t*vi3t
              ekin=ekin+pi3t*pi3t/amassi
              ekinsg=ekinsg+avgpi3*avgpi3/amassi
          end do
            ! <(avg(f-avg(f))avg(v))>
            sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
            ! <(avg(p)avg(v))>
            sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
            ! energy conservation friction constant
            sgbeta=(2.0d0+gam)*sumgv/(2.0d0*sumpv-sumgv*dtx)
            !sgbeta=0.0d0
            fact=dtx*(gammas+sgbeta)
            do  m = 1,3
              fi3=frc(m,i)
              vi3t=((2.0d0-fact)*vel(m,i)+fi3*wfac)/(2.0d0+fact)
              vel(m,i)=vi3t
            end do
          ELSE 
               ! without guiding forces
            do  m = 1,3
              !   generate random number 
              call gauss( 0.d0, rsd, fln )
              FI3=FRC(m,I)+fln
              Frc(m,I)=FI3
              vi3t=((2.0d0-gam)*vel(m,i)+fi3*wfac)/(2.0d0+gam)
              vel(m,i)=vi3t
          END DO
        ENDIF
      END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=sumgam
          TEMP1(2)=EKIN
          TEMP1(3)=EKINSG
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,3,&
             MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          sumgam=TEMP1(1)
          EKIN=TEMP1(2)
          EKINSG=TEMP1(3)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,11, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          sumgam=TEMP2(1)
          EKIN=TEMP2(2)
          EKINSG=TEMP2(3)
# endif
        ENDIF
#endif
    ! Estimate low frequency temperatures
        TEMPI=EKIN/sgrndf/BOLTZ
        TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPHF=Tempi-TEMPLF
        sgscale=20.455d0*sumgam/(isgend-isgsta+1)
        ! update accumulators
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGLDW

        SUBROUTINE SGMDW(NATOM,ISTART,IEND, &
             DTX,ENER,AMASS,WINV,crd,Frc,Vel)
!-----------------------------------------------------------------------
!     This routine calculate guiding force using SGLD method 
!     for MD simulation        
!
      use state
      use random, only: GAUSS
      implicit none
#ifdef MPI
      integer ierr
# include "parallel.h"
      _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
      _REAL_ :: temp2(20)
# endif
#endif
      INTEGER NATOM,ISTART,IEND,NTP
      _REAL_ DTX
      type(state_rec) :: ener
      _REAL_ AMASS(*),WINV(*),crd(3,*),Frc(3,*),Vel(3,*)
!
      INTEGER JSTA,JEND,I,M,j,nsgsubi,idx_nbex,jatm
      _REAL_ BOLTZ,AMASSI,TEMPI
      _REAL_ FACT,WFAC
      _REAL_ EKIN,EKINSG,SGBETA
      double precision sumgam,sumfp,sumpp,sumgv,sumpv
      double precision sggammai,avgpi3,pi3t,avgdfi3,avgri3,fsgpi,fsgfi,fsgi3,frici
      double precision xi3,x1i3,x2i3,vi3t,vi3,fi3
      PARAMETER (BOLTZ = 1.987192d-3)
!
        JSTA=ISTART
        JEND=IEND
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        sumgam=0.0d0
        EKIN=0.0D0
        EKINSG=0.0D0
        DO  I = jsta,jend
          AMASSI = AMASS(I)
          WFAC =  DTX*0.5D0*WINV(I)
          sggammai=-sgfps(i)/sgpps(i)
          sumgam=sumgam+sggammai

          sumfp=0.0d0
          sumpp=0.0d0
          sumgv=0.0d0
          sumpv=0.0d0
          DO  M = 1,3
              ! avg(x)
            xi3=crd(m,i)
            if(tsgbond)then 
              nsgsubi=atm_sg_maskdata(i)%cnt
              xi3=amassi*xi3
              do j=1,nsgsubi
                idx_nbex=atm_sg_maskdata(i)%offset + j 
                jatm=atm_sg_mask(idx_nbex)
                xi3=xi3+amass(jatm)*crd(m,jatm)
              enddo
              xi3=xi3/sgmass(i)
            endif
            x1i3=sgavg0*avgx1(m,i)+sgavg1*xi3
            avgx1(m,i)=x1i3
            ! avgavg(x)
            x2i3=sgavg0*avgx2(m,i)+sgavg1*x1i3
            avgx2(m,i)=x2i3
            ! avg(p)
            avgpi3=tsgfac*amassi*(xi3-x1i3)
            pi3t=(avgpi3-sgavg0*avgp(m,i))/sgavg1
            avgp(m,i)=avgpi3
            ! avg(f-avg(f))
            avgdfi3=tsgfac*(pi3t-2.0d0*avgpi3+tsgfac*amassi*(x1i3-x2i3))
            ! sum(avg(f-avg(f))avg(p))
            sumfp=sumfp+avgdfi3*avgpi3
            ! sum(avg(p)avg(p))
            sumpp=sumpp+avgpi3*avgpi3
            ! guiding forces
            fsgpi=sgfti*sggammai*avgpi3
            fsgfi=sgffi*avgdfi3
            fsgi3=fsgpi+fsgfi
            fi3=frc(m,i)+fsgi3
            frc(m,i)=fi3
           ! estimate velocity at t+dt/2
            vi3t=vel(m,i)+fi3*wfac
            ! sum(g*v)
            sumgv=sumgv+fsgi3*vi3t
            ! sum(p*v)
            sumpv=sumpv+amassi*vi3t*vi3t
            ekin=ekin+pi3t*pi3t/amassi
            ekinsg=ekinsg+avgpi3*avgpi3/amassi
          enddo
            ! <(avg(f-avg(f))avg(v))>
          sgfps(i)=sgavp0*sgfps(i)+sgavp1*sumfp
          ! <(avg(p)avg(v))>
          sgpps(i)=sgavp0*sgpps(i)+sgavp1*sumpp
          ! energy conservation friction constant
          sgbeta=2.0d0*sumgv/(2.0d0*sumpv-sumgv*dtx)
          fact=sgbeta/(1.0d0+0.5d0*sgbeta*dtx)
          do  m = 1,3
            fi3=frc(m,i)
            vi3t = vel(m,i) + fi3*wfac
            frici=fact*amassi*vi3t
            frc(m,i)=fi3-frici
          end do
        END DO
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=sumgam
          TEMP1(2)=EKIN
          TEMP1(3)=EKINSG
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,3, &
            MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          sumgam=TEMP1(1)
          EKIN=TEMP1(2)
          EKINSG=TEMP1(3)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,3, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          sumgam=TEMP2(1)
          EKIN=TEMP2(2)
          EKINSG=TEMP2(3)
# endif
        ENDIF
#endif
    ! Estimate low frequency temperatures
        TEMPI=EKIN/sgrndf/BOLTZ
        TEMPLF=SGAVP0*TEMPLF+SGAVP1*EKINSG/sgrndf/BOLTZ
        TEMPHF=TEMPI-TEMPLF
        sgscale=20.455d0*sumgam/(isgend-isgsta+1)
        ! update accumulators
        CALL SGENERGY(ENER)
        RETURN
        END SUBROUTINE SGMDW
     
        SUBROUTINE SGENERGY(ENER)
!-----------------------------------------------------------------------
!     This routine set the ener fields of SGLD variables
!
      use state
      implicit none
      type(state_rec) :: ener
      _REAL_ BOLTZ
      PARAMETER (BOLTZ = 1.987192d-3)
      _REAL_ EPOTI
    ! Weighting accumulators
        EPOTI=ENER%POT%TOT
        IF(EPOTLF>1.0D10)THEN
          EPOTLF=EPOTI
          EPOTLLF=EPOTI
        ELSE
          EPOTLF=SGAVG0*EPOTLF+SGAVG1*EPOTI
          EPOTLLF=SGAVG0*EPOTLLF+SGAVG1*EPOTLF
        ENDIF
        EPOTHF=EPOTI-EPOTLF
        sgwt=(psgldg-sgffi)*(epotlf-epotllf)/(boltz*tsgset)
    ! Update ENER structure
        ENER%SGLD%SGSCALE=SGSCALE
        ENER%SGLD%TEMPLF=TEMPLF
        ENER%SGLD%TEMPHF=TEMPHF
        ENER%SGLD%EPOTLF=EPOTLF
        ENER%SGLD%EPOTHF=EPOTHF
        ENER%SGLD%EPOTLLF=EPOTLLF
        ENER%SGLD%SGWT=SGWT
       RETURN
       END SUBROUTINE SGENERGY

#ifdef MPI

!*********************************************************************
!               SUBROUTINE SGLD_EXCHG
!*********************************************************************
!  exchange all data in the SGLDR common block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sgld_exchg(irep)

   implicit none
#  include "parallel.h"

   integer, intent(in) :: irep
   integer  ierror, istatus(mpi_status_size)
   call mpi_sendrecv_replace(sgft,nsgld_real, mpi_double_precision, &
                   irep, 511, irep, 511, commmaster, istatus, ierror)
    !call mpi_barrier(commmaster, ierror)
   return
   end subroutine sgld_exchg

!*********************************************************************
!               SUBROUTINE REMD_SCALE_VELO
!*********************************************************************
! Scale velocities based on new temps after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine rxsgld_scale(stagid,nr,myscaling,amass,crd, vel)

   implicit none
#  include "parallel.h"

   integer, intent(in) :: stagid,nr
   _REAL_, intent(in) :: myscaling
   _REAL_,  intent(in) :: crd(3,nr)
   _REAL_,  intent(inout) :: vel(3,nr)
   _REAL_, dimension (*), intent(in)    :: amass
   integer ierror
   integer i,j,jsta,jend
   _REAL_ amassi,xi,x1i,x2i
   _REAL_ temp1(10),temp2(10)
!--------------------
         !if (sanderrank==0) then
         !   write (6,'(a,i4,2x,f8.3,a,2f8.3)') &
         !      "RXSGLD: stagid, scalsg  ",stagid,myscalsg,&
         !      " to match sgft,tempsg,: ",sgft,tempsg
         !endif
      call mpi_bcast(stagid,1,mpi_integer,0,commsander,ierror)

         ! All processes scale velocities.
         ! DAN ROE: This could potentially be divided up as in runmd
         !  since when there are mutiple threads per group each thread 
         !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
         if (sanderrank==0) then
            write (6,'(a,f8.3,a,f8.3)') &
               "| RXSGLD: scaling guiding properties by ",myscalsg,&
               " to match a new guiding factors sgfti, sgffi ",sgfti,sgffi
         endif
#endif
! ---=== Broadcast RXSGLD guiding effect ===---
      IF(SANDERSIZE > 1)call mpi_bcast(sgft,nsgld_real,mpi_double_precision,&
                                              0,commsander,ierror)
      if (myscaling > 0.0d0) then
        vel(:,:)=myscaling*vel(:,:)
      endif                         
      if (myscalsg > 0.0d0) then
        templf=myscalsg*myscalsg*templf
        avgp(:,:)=myscalsg*avgp(:,:)
        avgr(:,:)=myscalsg*avgr(:,:)
        sgfps(:)=myscalsg*sgfps(:)
        sgpps(:)=myscalsg*myscalsg*sgpps(:)
        jsta = iparpt(mytaskid) + 1
        jend = iparpt(mytaskid+1)
        IF(JSTA < ISGSTA)JSTA=ISGSTA
        IF(JEND > ISGEND)JEND=ISGEND
        do i = jsta,jend
            amassi=amass(i)
            do j=1,3
              xi=crd(j,i)
              x1i=avgx1(j,i)
              x2i=avgx2(j,i)
              avgx1(j,i)=xi+myscalsg*(x1i-xi)
              avgx2(j,i)=xi+myscalsg*(x2i-xi)
           enddo
         enddo
      endif
   return

end subroutine rxsgld_scale

!*********************************************************************
!               FUNCTION TEMPSGLOOKUP
!*********************************************************************
! lookup temp in templist and return its index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function tempsglookup(numreps,temp,sgft,sgff,temps,sgfts,sgffs)

   implicit none
#  include "parallel.h"

   integer numreps
   _REAL_, intent(in) :: temp,sgft,sgff
   _REAL_, dimension(numreps), intent(in) :: temps,sgfts,sgffs

   integer i
   
   tempsglookup=0
   do i = 1, numreps
      if(abs(temp-temps(i)) < 1.0d-6 &
      .and. abs(sgft-sgfts(i)) < 1.0d-6 &
      .and. abs(sgff-sgffs(i)) < 1.0d-6) then
         if(tempsglookup>0)then
            write (6,*) "================================"
            write (6,*) "Two replicas are the same: ",tempsglookup,i
            write (6,*) "================================"
            call mexit(6,1)
         endif
         tempsglookup = i
      end if
   end do
   return
end function tempsglookup

!*********************************************************************
!               FUNCTION STAGIDLOOKUP
!*********************************************************************
! lookup stagid in stagidlist and return its neighboring index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function stagidlookup(numreps,id,idtable)

   implicit none
#  include "parallel.h"

   integer numreps
   INTEGER, dimension(numreps), intent(in) :: idtable
   INTEGER, intent(in) :: id

   integer i
   
   stagidlookup=-1
   do i = 1, numreps
      if(id==idtable(i)) stagidlookup=i
   end do
   return
end function stagidlookup

!*********************************************************************
!               SUBROUTINE SORTTEMPSG
!*********************************************************************
! sort temp ascendingly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine sorttempsg(numreps,temps,psgldgs,sgffs)

   implicit none

#  include "parallel.h"

   integer numreps
   _REAL_, dimension(numreps), intent(inout) :: temps,psgldgs,sgffs
   _REAL_, dimension(numreps) :: tmp
   INTEGER, dimension(numreps) :: tmpid

   _REAL_ tempt
   integer i, j, ii

   do i = 1, numreps
     tmp(i)=1000000*temps(i)+100*(psgldgs(i)- sgffs(i))
     tmpid(i)=i
   enddo
   do i = 1, numreps
      do j = i + 1, numreps
         if(tmp(j) < tmp(i)) then
            tempt = tmp(i)
            tmp(i) = tmp(j)
            tmp(j) = tempt
            ii=tmpid(i)
            tmpid(i)=tmpid(j)
            tmpid(j)=ii
         end if
      end do
   end do
   do i = 1, numreps
     ii=tmpid(i)
     tmp(i)=temps(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     temps(i)=tmp(i)
     tmp(i)=psgldgs(ii)
   enddo
   do i = 1, numreps
     ii=tmpid(i)
     psgldgs(i)=tmp(i)
     tmp(i)=sgffs(ii)
   enddo
   do i = 1, numreps
     sgffs(i)=tmp(i)
   enddo
   return
end subroutine sorttempsg

#endif /* MPI */

end module sgld

