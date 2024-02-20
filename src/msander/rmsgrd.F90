! <compile=optimized>
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate the RMS gradient
!-----------------------------------------------------------------------
!     --- RMSGRD ---
!-----------------------------------------------------------------------

subroutine rmsgrd(forces, grms)
   
   use constants, only : ZERO
   
   implicit none
   
#include "../include/md.h"
#include "box.h"
#include "../include/memory.h"
#include "nmr.h"

   _REAL_, intent(in) :: forces(*)
   _REAL_, intent(out) :: grms
   
   ! the ddot function is external
   _REAL_ ddot
   
   _REAL_ :: dotprod
   
   integer :: numcomponents
   
   if (ibelly > 0) then
      numcomponents = natbel * 3 + iscale
   else
      numcomponents = nrp * 3 + iscale
   end if
   
   ! Initialise rmsgrad so it at least has a valid value
   grms = ZERO
   dotprod = ddot(3*nrp+iscale, forces, 1, forces, 1)
   if (numcomponents /= 0) grms = sqrt(dotprod / numcomponents)
   
end subroutine rmsgrd


