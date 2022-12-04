
!------------------------------------------------------------------------------
! mexit: Platform independent exit, designed to produce the same behavior and
!        return the same status to the OS no matter the architecture.  This
!        routine detects the specific platform and proceeds accordingly.
!
! Arguments:
!   output_unit:   this unit will be closed if its index is greater than zero
!                  and this is not an MPI slave process
!   status:        exit status (returned)
!------------------------------------------------------------------------------
subroutine mexit(output_unit, status)

#ifdef MPI
   use mpi
#endif
  use xray_cpu_module, only: xray_fini=>finalize
  implicit none
  integer output_unit
  integer status

#ifdef MPI
  integer ierr
#  include "parallel.h"
   
  ! We always need to call xray_fini() in order to properly shut down the GPU
  if( mytaskid .eq. 0 ) call xray_fini()

  ! Status .gt. 0 implies an error condition, therefore
  ! kill all the nodes.  mpi_abort on the world communicator
  ! should do this, but it does not on some implementations.
  ! Some MPI's have an environmental sledge hammer that kills
  ! every MPI process if one dies: mpiexec -kill   
  if (status /= 0) then
    call flush(output_unit)
    call mpi_abort(MPI_COMM_WORLD, status, ierr)
  else
    call mpi_finalize(ierr)
  endif
#else
  call xray_fini()
#endif

  if (output_unit > 0 .and. status/=0) then
    close(unit = output_unit)
  endif

  call exit(status)
end subroutine mexit 
