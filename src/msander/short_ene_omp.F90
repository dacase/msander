!<compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

!------------------------------------------------------------------------------
! get_nb_energy: the main routine for vdw, hbond, and direct space Ewald sum
!                computations. 
!
!------------------------------------------------------------------------------
subroutine get_nb_energy(iac, ico, ntypes, charge, cn1, cn2, cn6, force, &
                         numatoms, ipairs, ewaldcof, eedtbdns, eed_cub, &
                         eed_lin, maxnblst, eelt, evdw, ehb, dir_vir, eedvir, &
                         filter_cut, ee_type, eedmeth, dxdr, pol, pol2, cn3, &
                         cn4, cn5, epol, dipole, field, mpoltype)

  use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
                     myindexlo, myindexhi, numimg, tranvec, nucgrd
  use stack
#ifdef MPI /* SOFT CORE */
  use softcore, only : sc_ener
  use nblist, only : numsc
#endif
  use constants, only: zero, one, two, half, third, TWOPI, six, twelve
  use file_io_dat
#ifdef LES
  use les_data, only: cnum, lestyp, lestmp, lesfac, lfac, nlesty
#endif
  use decomp, only: decpr, decpair
  use nbips, only: teips, tvips, nnbips, rips2, ripsr, rips2r, rips6r, &
                   rips12r, aipse, aipsvc, aipsva, bipse, bipsvc, bipsva, &
                   pipsec, pipsvcc, pipsvac
#ifdef MPI /* SOFT CORE */
  use softcore, only: scalpha, scbeta, sigma6, foureps, sc_dvdl, &
                      isProcessV1, sc_ener, nsc, oneweight, weight0, &
                      weight1, sceeorder, sc_dvdl_ee
  use mbar, only: ifmbar, bar_i, bar_states, bar_lambda, bar_cont, do_mbar
#endif
  use crg_reloc, only: ifcr, cr_add_dcdr_factor
  use omp_lib

  implicit none
  character(kind=1, len=13) :: routine="get_nb_energy"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
#include "flocntrl.h"
#include "def_time.h"
   
  integer numatoms, maxnblst, mpoltype
  integer iac(*), ico(*), ntypes, ee_type, eedmeth
  _REAL_ charge(*), cn1(*), cn2(*), cn6(*)
  _REAL_ ewaldcof, eedtbdns, dxdr, eed_cub(*), eed_lin(2,*), dir_vir(3,3)
  integer ipairs(maxnblst)
  _REAL_ force(3,numatoms), eelt, epol, evdw, ehb
  _REAL_ eedvir, filter_cut, dipole(3,*), field(3,*), pol(*)
  _REAL_ pol2(*)
  _REAL_ cn3(*), cn4(*), cn5(*)
   
  integer index, numpack, i, k, ncell_lo, ncell_hi, ntot, nvdw, nhbnd
  _REAL_ xk, yk, zk
  integer ic,j,m,n,ind,iaci,inddel
  _REAL_ del, delrinv, delr12inv
  _REAL_ switch, d_switch_dx
  _REAL_ b0, b1, b2
  _REAL_ filter_cut2, xx
  _REAL_ comm1
  _REAL_ xktran(3,18)
  _REAL_ e3dx, e4dx
#ifdef TVDW
  _REAL_ r4, r6pinv
#endif
  _REAL_ ecur
  integer, parameter :: mask27 = 2**27 - 1
#ifdef LES
#  include "ew_cntrl.h"
#endif
  _REAL_ delx(3), delr, delr2, cgi, cgj, delr2inv, r6, f6, f12, df, &
         dfee, dx, x, dfx(3), dumx(3)
#ifdef MPI
  _REAL_ denom, denom_n, delr_n, switch_c, denom2, denom3, rfour
#endif

#include "../include/md.h"
  _REAL_ uips, uips2, uips4, uips2r, uips6r, uips12r
  _REAL_ pipse, dpipse, pvc, dvcu, pva, dvau
  integer itran

  _REAL_ time0, time1
   
  if (do_dir == 0) then
    return
  end if
  eelt = zero
  epol = zero
  evdw = zero
  ehb = zero
  eedvir = zero
#ifdef MPI /* SOFT CORE */
  sc_ener(7) = 0.0d0
  sc_ener(8) = 0.0d0  
  sc_ener(9) = 0.0d0
  sc_ener(10) = 0.d0
  sc_ener(11) = 0.d0
#endif
  dir_vir(1:3,1:3) = zero
  del = one / eedtbdns
  filter_cut2 = filter_cut * filter_cut
  numpack = 1

  call timer_start(TIME_SHORT_ENE)

!$omp parallel private( inddel, myindexlo, myindexhi, ncell_lo, ncell_hi, &
!$omp&  k,i,xk,yk,zk,ntot,nvdw,nhbnd,dumx,cgi,iaci,m,xktran, &
!$omp&  n,itran,j,delx,delr2,delrinv,x,  &
!$omp&  ind,dx,e3dx,e4dx,switch,d_switch_dx,b0,b1,cgj,comm1,ecur,  &
!$omp&  dfee,delr2inv,ic,r6,delr12inv,f6,f12,dfx)  &
!$omp&  firstprivate(numpack)

!$  inddel = (nucgrd-1) / OMP_GET_NUM_THREADS() + 1
!$  if (inddel == 0) inddel = 1
!$  myindexlo = 1 + OMP_GET_THREAD_NUM()*inddel
!$  myindexhi = myindexlo + inddel - 1
!$  if (OMP_GET_THREAD_NUM() == OMP_GET_NUM_THREADS()-1) myindexhi = nucgrd

  ! do index = myindexlo, myindexhi
  do index = 1,nucgrd
    if (numimg(index) <= 0) cycle
    ncell_lo = nlogrid(index)
    ncell_hi = nhigrid(index)
    do k = ncell_lo, ncell_hi
      i = bckptr(k)
#ifdef MPI /* SOFT CORE */
      ! SOFT CORE contribution in numsc
      ntot = numvdw(i) + numhbnd(i) + numsc(i)
#else
      ntot = numvdw(i) + numhbnd(i)
#endif
      nvdw = numvdw(i)
      nhbnd = numhbnd(i)
      if (ntot <= 0) cycle
!$    if( index < myindexlo .or. index > myindexhi ) then
!$       ! some other thread will do the work here:
!$       numpack = numpack + ntot
!$       cycle
!$    endif
      xk = imagcrds(1,k)
      yk = imagcrds(2,k)
      zk = imagcrds(3,k)
!==========  short ene routine in-lined below:===============================
  dumx(:) = zero
  cgi = charge(i)
  iaci = ntypes * (iac(i) - 1)
   
  do m = 1, 18
    xktran(1,m) = tranvec(1,m) - xk
    xktran(2,m) = tranvec(2,m) - yk
    xktran(3,m) = tranvec(3,m) - zk
  end do

#ifdef LES
  lestmp=nlesty*(lestyp(i)-1)
#endif
   
  !-----------------------------------------------
  ! Loop over the 12-6 LJ and electrostatic terms 
  !-----------------------------------------------
  do m = 1, nvdw+nhbnd

    n=ipairs(m+numpack-1)
    itran=ishft(n,-27)
    n = iand(n,mask27)
    j = bckptr(n)
    delx(:) = imagcrds(:,n) + xktran(:,itran)
    delr2 = delx(1)*delx(1) + delx(2)*delx(2) + delx(3)*delx(3)

    !  write(6,'(9i5)') index,k,i,j,OMP_GET_THREAD_NUM(),m,n,itran,numpack
 
    if ( delr2 >= filter_cut2 ) cycle

    delrinv = 1.d0/sqrt(delr2)
    x = dxdr*delrinv*delr2

    ! Cubic spline on switch:
    ind = int(eedtbdns*x)
    dx = x - ind*del
    ind = 4*ind
    e3dx = dx * eed_cub(3+ind)
    e4dx = dx * dx * eed_cub(4+ind)
    switch = eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
             (e3dx + e4dx*third)*half)
    d_switch_dx = eed_cub(2+ind) + e3dx + e4dx*half
       
    ! Tom Darden got the idea for B_l from Walter Smith's CCP5 article 1982
    ! Ewald for point multipoles.
    b0 = switch*delrinv
    b1 = b0 - d_switch_dx*dxdr
    cgj = charge(j)
#ifdef LES
    if (use_pme .ne. 0) then

      ! If we are using PME, then the correction for lfac will
      ! be done after the reciprocal space calculation is done,
      ! so no need for it here
      comm1 = cgi * cgj

      ! calculate contribution of dc/dr to force
      if (ifcr .ne. 0) then
        call cr_add_dcdr_factor(i, b0*cgj)
        call cr_add_dcdr_factor(j, b0*cgi)
      end if
    else
      lfac = lesfac(lestmp + lestyp(j))
      comm1 = cgi * cgj * lfac

      ! Calculate contribution of dc/dr to force
      if (ifcr .ne. 0) then
        b2 = b0*lfac
        call cr_add_dcdr_factor(i, b2*cgj)
        call cr_add_dcdr_factor(j, b2*cgi)
      end if
    end if
#else
    comm1 = cgi*cgj

    ! Calculate contribution of dc/dr to force
    if (ifcr .ne. 0) then
      call cr_add_dcdr_factor(i, b0*cgj)
      call cr_add_dcdr_factor(j, b0*cgi)
    end if
#endif
    ecur = comm1 * b0
!$omp critical
    eelt = eelt + ecur
!$omp end critical

#if 0
    ! Thermodynamic Integration decomposition
    if (decpr .and. idecomp > 0) then
      call decpair(2, i, j, ecur/(nstlim/ntpr))
    end if
#endif
    dfee = comm1*b1
#ifdef LES
#  include "ene_decomp.h"
#endif
    delr2inv = delrinv*delrinv
    dfee = dfee*delr2inv

    ! epilogue: 12-6 LF terms

    if (m<=nvdw) then
       ic = ico(iaci+iac(j))
       r6 = delr2inv*delr2inv*delr2inv
       delr12inv = r6 * r6
#ifdef LES 
       lfac=lesfac(lestmp+lestyp(j))
       f6 = cn2(ic)*r6*lfac
       f12 = cn1(ic)*delr12inv*lfac
#else
       f6 = cn2(ic)*r6
       f12 = cn1(ic)*delr12inv
#endif
!$omp critical
       evdw = evdw + f12 - f6
!$omp end critical
       dfee = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
    endif

    dfx(:) = delx(:)*dfee
    dumx(:) = dumx(:) + dfx(:)
!$omp critical
    force(:,j) = force(:,j) + dfx(:)
!$omp end critical
  end do 

  !  end of what is now a single loop over j that are in the list as
  !    being close to i

!   #include "soft_core.inc"   /* MPI code for soft-core electrostatics */
      
  ! Contribute forces on the ith particle, accumulated as forces
  ! from the ith particle were added onto all other particles.
!$omp critical
  force(:,i) = force(:,i) - dumx(:)
!$omp end critical

!==========  short ene routine in-lined above===============================
                        
      numpack = numpack + ntot
    end do  !  k = ncell_lo,ncell_hi
  end do
!$omp end parallel
  ! nd loop over this process's assigned atoms

  call timer_stop(TIME_SHORT_ENE)

  return

end subroutine get_nb_energy 

