
  else if ( eedmeth == 2 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 2
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do

    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1,icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Linear lookup on switch:
      delrinv = one / sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      xx = eedtbdns*x + 1
      ind = int(xx)
      dx = xx - ind
      switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
      d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
      end if
#else
      comm1 = cgi*cgj
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End prologue loop

    ! regular epilogue:
#   include "ew_directe.h"

    ! Now loop over the 12-10 LJ terms for eedmeth = 2
    icount = 0
    do m = nvdw + 1, ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Linear lookup on switch:
      delrinv = one/sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      xx = eedtbdns*x + 1
      ind = int(xx)
      dx = xx - ind
      switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
      d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi*cgj
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
      end if
#else
      comm1 = cgi*cgj
#endif
      ecur = comm1*b0        
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End epilogue loop

#   include "ew_directe2.h"

  else if ( eedmeth == 3 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 3
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do
      
    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1,icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Explicit function call:
      delrinv = one/sqrt(delr2)
      delr = delr2*delrinv
      delr2inv = delrinv*delrinv
      x = dxdr * delr
      call get_ee_func(x, switch, d_switch_dx, ee_type)
      b0 = switch * delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi * cgj * lfac
      end if
#else
      comm1 = cgi * cgj
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new)=delr2inv
      cache_df(im_new)=dfee
    end do
    ! End electrostatic loop

    ! regular epilogue:
#   include "ew_directe.h"

    ! Now loop over the 12-10 LJ terms for eedmeth = 3
    icount = 0
    do m = nvdw + 1, ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Explicit function call:
      delrinv = one / sqrt(delr2)
      delr = delr2 * delrinv
      delr2inv = delrinv * delrinv
      x = dxdr * delr
      call get_ee_func(x, switch, d_switch_dx, ee_type)
      b0 = switch*delrinv
      b1 = (b0 - d_switch_dx*dxdr)*delr2inv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then

        ! If we are using PME, then the correction for lfac will
        ! be done after the reciprocal space calculation is done,
        ! so no need for it here.
        comm1 = cgi * cgj
      else
        lfac = lesfac(lestmp+lestyp(j))
        comm1 = cgi*cgj*lfac
      end if
#else
      comm1 = cgi * cgj
#endif
      ecur = comm1 * b0
      eelt = eelt + ecur
      dfee = comm1 * b1

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End epilogue loop

#   include "ew_directe2.h"
  else if (eedmeth == 4) then

    ! Loop over the 12-6 LJ terms for eedmeth = 4
    icount = 0
    do m = 1, nvdw
#     include "ew_directp.h"
    end do
      
    ! Calculation starts: loop over the data gathered in the temporary
    ! array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi*cgj*lfac*delrinv
      end if
#else
      b0 = cgi * cgj * delrinv
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = b0 * delr2inv
         
      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do
    ! End prologue loop

    ! regular epilogue:
#   include "ew_directe.h"

    ! Now loop over the 12-10 LJ terms for eedmeth = 4
    icount = 0
    do m = nvdw+1, nvdw+nhbnd
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount

      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delrinv
      end if
#else
      b0 = cgi * cgj * delrinv
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = b0 * delr2inv
      cache_r2(im_new)=delr2inv
      cache_df(im_new)=dfee
    end do

#   include "ew_directe2.h"

#ifdef MPI /* SOFT CORE */
    ! Now loop over the softcore (modified 12-6) LJ terms for eedmeth = 4
    icount = 0
    do m = nvdw+nhbnd+1,ntot
#     include "ew_directp.h"
    end do
      
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)

      ! Don't use a switch: straight Coulomb
      delrinv = one / sqrt(delr2)
      delr2inv = delrinv * delrinv
      cgj = charge(j)
#  if defined(LES) 
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delrinv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delrinv
      end if
#  else
      b0 = cgi*cgj*delrinv
#  endif /* LES */
      eelt = eelt + b0
      dfee = b0 * delr2inv
      cache_r2(im_new) = delr2

      ! Contrary to the 12-6 and 12-10 cases above, cache_r2 contains r^2 here
      cache_df(im_new) = dfee
    end do

    ! V1 uses ew_directe3.h, in which softcore atoms are treated as appearing,
    ! i.e. fully interacting at lambda=1 and 'soft' at small lambda
    ! V0 uses ew_directe4.h, in which softcore atoms are treated as vanishing, 
    ! i.e. fully interacting at lambda=0 and 'soft' at large lambda
    if (isProcessV1) then
#     include "ew_directe3.h"
    else
#     include "ew_directe4.h"
    end if
#endif /* MPI for SOFT CORE */

  else if ( eedmeth == 5 )then

    ! Loop over the 12-6 LJ terms for eedmeth = 5
    icount = 0
    do m = 1,nvdw
#     include "ew_directp.h"
    end do
      
    !  Calculation starts: loop over the data gathered in the temporary
    !  array caches.
    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Use a distance-dependent dielectric of 1/r:
      delr2inv = one / delr2
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delr2inv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delr2inv
      end if
#else
      b0 = cgi*cgj*delr2inv
#endif 
      ecur = b0
      eelt = eelt + ecur
      dfee = two * b0 * delr2inv

      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee
    end do

    ! regular epilogue:
#   include "ew_directe.h"

    !     Now loop over the 12-10 LJ terms for eedmeth = 5
    icount = 0
    do m = nvdw+1, ntot
#     include "ew_directp.h"
    end do

    ! SGI compiler directive to prevent compiler loop fusioning.
    !*$* NO FUSION
    do im_new = 1, icount
      j = cache_bckptr(im_new)
      delr2 = cache_r2(im_new)
         
      ! Use dielectric of 1/r:
      delr2inv = one / delr2
      cgj = charge(j)
#ifdef LES
      if (use_pme .ne. 0) then
        b0 = cgi * cgj * delr2inv
      else
        lfac = lesfac(lestmp+lestyp(j))
        b0 = cgi * cgj * lfac * delr2inv
      end if
#else
      b0 = cgi * cgj * delr2inv
#endif
      ecur = b0
      eelt = eelt + ecur
      dfee = two * b0 * delr2inv
      cache_r2(im_new) = delr2inv
      cache_df(im_new) = dfee

    end do
#   include "ew_directe2.h"
  else if (eedmeth == 6) then

    write(6,*) 'Error:  eedmeth=6 is not supported in msander'
    call mexit(6,1)

#   include "ew_directe2.h"
