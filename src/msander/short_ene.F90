!<compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

!------------------------------------------------------------------------------
! get_nb_energy: the main routine for vdw, hbond, and direct space Ewald sum
!                computations. 
!------------------------------------------------------------------------------
subroutine get_nb_energy(iac, ico, ntypes, charge, cn1, cn2, cn6, force, &
                         numatoms, ipairs, ewaldcof, eedtbdns, eed_cub, &
                         eed_lin, maxnblst, eelt, evdw, ehb, dir_vir, eedvir, &
                         filter_cut, ee_type, eedmeth, dxdr, pol, pol2, cn3, &
                         cn4, cn5, epol, dipole, field, mpoltype)

  use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
                     myindexlo, myindexhi, numimg, tranvec, nucgrd
  use stack
  use constants, only: zero, one, two, half, third, TWOPI, six, twelve
  use file_io_dat
#ifdef LES
  use les_data, only: cnum, lestyp, lestmp, lesfac, lfac, nlesty
#endif
  use nbips, only: teips, tvips, nnbips, rips2, ripsr, rips2r, rips6r, &
                   rips12r, aipse, aipsvc, aipsva, bipse, bipsvc, bipsva, &
                   pipsec, pipsvcc, pipsvac
  use omp_lib

  implicit none
  character(kind=1, len=13) :: routine="get_nb_energy"
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
  integer ic,j,m,n,ind,iaci,inddel
  _REAL_ del, delrinv, delr12inv
  _REAL_ b0
  _REAL_ filter_cut2
  _REAL_ comm1
  _REAL_ xktran(3,18)
  _REAL_ e3dx, e4dx, eeltl, evdwl
  _REAL_ forcel(3,numatoms)
#ifdef TVDW
  _REAL_ r4, r6pinv
#endif
  integer, parameter :: mask27 = 2**27 - 1
#ifdef LES
#  include "ew_cntrl.h"
#endif
  _REAL_ delx(3), delr, delr2, cgi, delr2inv, r6, f6, f12, df, &
         dfee, dx, x, dfx(3)

!$  character(len=30) omp_num_threads
!$  integer max_threads, ier

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
  dir_vir(1:3,1:3) = zero
  del = one / eedtbdns
  filter_cut2 = filter_cut * filter_cut
  numpack = 1

#ifdef OPENMP
  call get_environment_variable('OMP_NUM_THREADS', omp_num_threads, &
           status=ier)
  if( ier == 1 ) then
     max_threads = 1
  else
     read(omp_num_threads,*) max_threads
  endif
  max_threads = min( 8, max_threads )
#endif

  call timer_start(TIME_SHORT_ENE)

!$omp parallel private( inddel, myindexlo, myindexhi, index, &
!$omp&  k,i,ntot,nvdw,nhbnd,cgi,iaci,m,xktran, &
!$omp&  n,itran,j,delx,delr2,delrinv,x,  &
!$omp&  ind,dx,e3dx,e4dx,b0,comm1,  &
!$omp&  dfee,delr2inv,ic,r6,delr12inv,f6,f12,dfx,evdwl,eeltl,forcel)  &
!$omp&  firstprivate(numpack) num_threads(max_threads)

!$  inddel = (nucgrd-1) / OMP_GET_NUM_THREADS() + 1
!$  if (inddel == 0) inddel = 1
!$  myindexlo = 1 + OMP_GET_THREAD_NUM()*inddel
!$  myindexhi = myindexlo + inddel - 1
!$  if (OMP_GET_THREAD_NUM() == OMP_GET_NUM_THREADS()-1) myindexhi = nucgrd

  ! local versions of reducations: will be reduced in a single critical
  ! section after the end of the do loop
  evdwl = zero
  eeltl = zero
  forcel(:,:) = zero

  do index = 1,nucgrd
    if (numimg(index) <= 0) cycle
    do k = nlogrid(index), nhigrid(index)
      i = bckptr(k)
      ntot = numvdw(i) + numhbnd(i)
      if (ntot <= 0) cycle
!$    if( index < myindexlo .or. index > myindexhi ) then
!$       ! some other thread will do the work here; just need to update numpack
!$       numpack = numpack + ntot
!$       cycle
!$    endif
      nvdw = numvdw(i)
      nhbnd = numhbnd(i)

      !==========  short ene routine in-lined below:==========================
      cgi = charge(i)
      iaci = ntypes * (iac(i) - 1)
   
      xktran(1,:) = tranvec(1,:) - imagcrds(1,k)
      xktran(2,:) = tranvec(2,:) - imagcrds(2,k)
      xktran(3,:) = tranvec(3,:) - imagcrds(3,k)

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
        delr2 = dot_product( delx, delx )

        if ( delr2 >= filter_cut2 ) cycle

        delrinv = 1.d0/sqrt(delr2)
        x = dxdr*delrinv*delr2

        ! Cubic spline on switch:
        ind = int(eedtbdns*x)
        dx = x - ind*del
        ind = 4*ind
        e3dx = dx * eed_cub(3+ind)
        e4dx = dx * dx * eed_cub(4+ind)
        b0 = delrinv*(eed_cub(1+ind) + dx*(eed_cub(2+ind) + &
                 (e3dx + e4dx*third)*half))
#ifdef LES
        if (use_pme .ne. 0) then

          ! If we are using PME, then the correction for lfac will
          ! be done after the reciprocal space calculation is done,
          ! so no need for it here
          comm1 = cgi * charge(j)
        else
          lfac = lesfac(lestmp + lestyp(j))
          comm1 = cgi * charge(j) * lfac
        end if
#else
        comm1 = cgi*charge(j)
#endif
        dfee = comm1*(b0 - dxdr*(eed_cub(2+ind) + e3dx + e4dx*half))
        eeltl = eeltl + comm1*b0
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
           evdwl = evdwl + f12 - f6
           dfee = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
        endif

        dfx(:) = delx(:)*dfee
        forcel(:,j) = forcel(:,j) + dfx(:)
        forcel(:,i) = forcel(:,i) - dfx(:)
      end do 

      !  end of what is now a single loop over j that are in the list as
      !    being close to i

      ! Contribute forces on the ith particle, accumulated as forces
      ! from the ith particle were added onto all other particles.

      !==========  short ene routine in-lined above===========================
                        
      numpack = numpack + ntot
    end do  !  k = ncell_lo,ncell_hi
  end do
!$omp critical
  force(:,:) = force(:,:) + forcel(:,:)
  evdw = evdw + evdwl
  eelt = eelt + eeltl
!$omp end critical
!$omp end parallel

  call timer_stop(TIME_SHORT_ENE)

  return

end subroutine get_nb_energy 

