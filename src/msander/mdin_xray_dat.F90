!*******************************************************************************
!
! Module: mdin_xray_dat_mod
!
! Description: <TBS>
!
!*******************************************************************************

module mdin_xray_dat_mod

use file_io_dat_mod

  implicit none

! The mdin xray input options are grouped in the order of input, with the
! common blocks containing valid options and the private data containing
! any deprecated options.

  integer                  ixray, nsfc, itetherSF, itetherCRD
  double precision         wtSF, wtCRD

  common / mdin_xray_int /     ixray, nsfc, itetherSF, itetherCRD
  
  integer, parameter    :: mdin_xray_int_cnt = 4

  save  :: / mdin_xray_int /

  common / mdin_xray_dbl /     wtSF, wtCRD

  integer, parameter    :: mdin_xray_dbl_cnt = 2

  save  :: / mdin_xray_dbl /

  ! Amber masks
  character(256), public   :: sfrstmask, sfobsmask, sfInputFile, sfOutputFile, pdbRefFile, &
                              tethermask

  ! The active namelist:
  private       :: xray

  namelist /xray/    nsfc, itetherSF, itetherCRD, sfrstmask, sfobsmask, sfInputFile, &
                     sfOutputFile, pdbRefFile, tethermask
      
contains

!*******************************************************************************
!
! Subroutine:  init_mdin_xray_dat
!
! Description: <TBS>
!
!*******************************************************************************
subroutine init_mdin_xray_dat()

  use file_io_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

  integer :: ifind
  
  ixray        = 0
  nsfc         = 0
  itetherSF    = 0
  itetherCRD   = 0
  wtSF         = 0.0
  wtCRD        = 0.0
  sfrstmask    = ''
  sfobsmask    = ''
  tethermask   = ''
  sfInputFile  = ''
  sfOutputFile = ''
  pdbRefFile   = ''
  
  call nmlsrc('xray', mdin, ifind)

  if (ifind .ne. 0) then
    ixray = 1
    read(mdin, nml = xray)
  end if
  
  return

end subroutine init_mdin_xray_dat

!*******************************************************************************
!
! Subroutine:  validate_mdin_xray_dat
!
! Description: <TBS>
!
!*******************************************************************************
subroutine validate_mdin_xray_dat()

  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

  integer :: inerr

  inerr = 0
  if (nsfc .lt. 0) then
    write(mdout, '(a,i8)') 'Invalid frequency for structure factor &
                           &computation: ', nsfc
    inerr = 1
  end if

  ! Field any errors and bag out.
  if (inerr .eq. 1) then
   write(mdout, '(/,a)') ' Input errors occurred. Terminating execution.'
   call mexit(6, 1)
  end if
  
  return

end subroutine validate_mdin_xray_dat

!*******************************************************************************
!
! Subroutine:  print_mdin_xray_dat
!
! Description: <TBS>
!
!*******************************************************************************
subroutine print_mdin_xray_dat()

  implicit none

  if (ixray .eq. 1) then
    continue
  end if
 
  return

end subroutine print_mdin_xray_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_ctrl_dat
!
! Description: <TBS>
!
!*******************************************************************************
subroutine bcast_mdin_xray_dat()

  use parallel_dat_mod

  implicit none

  !call mpi_bcast(<variable names>, count, mpi_integer, 0, pmemd_comm, err_code_mpi)
  !call mpi_bcast(<variable names>, count, mpi_double_precision, 0, pmemd_comm, err_code_mpi)

  return

end subroutine bcast_mdin_xray_dat
#endif

end module mdin_xray_dat_mod
