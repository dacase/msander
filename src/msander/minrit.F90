! <compile=optimized>
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write restart data to file 'restrt' and possibly to a numbered restrt
subroutine minrit(n_force_calls,nrp,ntxo,x)
   
   use file_io_dat
   use binrestart, only: write_nc_restart

   implicit none
   
   ! Parameters
   integer n_force_calls,nrp
   _REAL_ :: x(*)
   
   ! Local variables
   character(len=89) restrt2
   character(len=12) num
   integer istart,iend,ntxo
   logical first
   save first
   data first/.true./
! DRR: This box.h include is only needed to get ntb for write_nc_restart
#  include "box.h"

   ! Netcdf restart
   if ( ntxo == 2) then
     call write_nc_restart(restrt,title,owrite,nrp,ntb,first,x,x,0.0d0,.false.&
#    ifdef MPI
                 , 0.0d0, 0, 0, (/ 0 /), (/ 0 /), (/ 0 /), 0, 0, 0, 0.0d0, 0.0d0 &
#    endif
                          )
     if (first) first=.false.
   else
   ! Standard formatted/unformatted restart   
      if (first) then
         call amopen(16,restrt,owrite,'F','W')
         first = .false.
      else
         call amopen(16,restrt,'O','F','W')
      end if
      call minri2(16,nrp,x)
      close(16)
   endif
   
   ! Ben Roberts: Added support for writing multiple
   ! restrt files, based on similar provisions in
   ! mdwrit.f.
      
   if (ntwr >= 0) return
   
   do iend=1,80
      if (restrt(iend:iend) <= ' ') goto 1
   end do
   1 continue
   iend = iend - 1
   
   ! Only do this if we have an actual step number
   ! to write.
   if (n_force_calls > 0) write(num,'(i12)') n_force_calls
   
   do istart=1,12
      if (num(istart:istart) /= ' ') goto 2
   end do
   2 continue
   write(restrt2, '(a,a,a)') restrt(1:iend), '_', num(istart:12)
   write(6,'(a,a)') ' writing ', restrt2
   
   ! Open a file on filehandle 17, with name restrt2
   if (ntxo == 2) then
      call write_nc_restart(restrt2,title,owrite,nrp,ntb,.true.,x,x,0.0d0,.false.&
#    ifdef MPI
                 , 0.0d0, 0, 0, (/ 0 /), (/ 0 /), (/ 0 /), 0, 0, 0, 0.0d0, 0.0d0 &
#    endif
                          )
      return
   else
      call amopen(17,restrt2,owrite,'F','W')
   end if
   
   ! Write coordinates to filehandle 17
   call minri2(17,nrp,x)
   close(17)
   
   return
   
   entry minwrit_reset
   first = .false.

end subroutine minrit

!------------------------------------------------------------------
! Ben Roberts: Incorporated coordinate-writing code into a new subroutine
! so it can be called for multiple files, not only "restrt". Enables
! writing of multiple restart files during a minimisation.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write restart data (coords only) to an arbitrary filehandle
subroutine minri2(nf,nrp,x)
   
   use nblist, only: a,b,c,alpha,beta,gamma
   use file_io_dat
   
   implicit none
   
   integer i,nf,nr,nrp,nr3
   _REAL_ x(*),tt
#  include "box.h"
   
   nr = nrp
   nr3 = 3*nr
   tt = 0.0d0
   
   write(nf,'(a80)') title
   if( nr < 100000 ) then
     write(nf,'(i5)') nr
   elseif ( nr < 1000000 ) then ! sander 7/8/9/10 large system format...
     write(nf,'(i6)') nr
   elseif ( nr < 10000000 ) then ! Sander 11 - 1 mil+ format
     write(nf,'(i7)') nr
   else
     write(nf,'(i8)') nr ! assume amber 11 VERY large system format. 10 mil+
   end if
   write(nf,'(6f12.7)') (x(i),i=1,nr3)
   if ( ntb /= 0 ) write(nf,'(6f12.7)') a,b,c,alpha,beta,gamma
   
   return
end subroutine minri2
