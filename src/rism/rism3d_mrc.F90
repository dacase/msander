!<compile=optimized>

#include "../include/dprec.fh"

module rism3d_mrc
  use safemem
  use rism_report_c
  implicit none

contains

  !> Write volumetric data to a file in MRC 2014 format.  see
  !!
  !! https://www.ccpem.ac.uk/mrc_format/mrc2014.php
  !!
  !! When writing in parallel, each process must call this function
  !! with its local data. Data transfer is handled internally.  We
  !! assume decomposition in the z-axis.
  !!
  !! MRC format is also compatible with the older CCP4
  !! format. However, some readers may ignore MRC-specific information
  !! such as the origin offset.
  !!
  !! @param[in] file File name to write to.
  !! @param[in] data Data to write in a n(1)*n(2)*n(3) linear array.
  !! @param[in] grid Grid object.
  !! @param[in] solute Solute object.
  !! @param[in] o_rank (optional) MPI process rank.
  !! @param[in] o_nproc (optional) MPI number of processes.
  !! @param[in] o_comm (optional) MPI communicator.
  subroutine rism3d_mrc_map_write (file, data, grid, solute, o_rank, o_nproc, o_comm)
    use constants_rism, only: PI
    use rism_util, only : freeUnit, rmExPrec
    use rism3d_grid_c
    use rism3d_solute_c
    use iso_fortran_env
    implicit none
#if defined(MPI)
    include 'mpif.h'
    integer status(MPI_STATUS_SIZE)
#endif /*defined(MPI)*/
    character(len=*), intent(in) :: file
    type(rism3d_grid), intent(in) :: grid
    _REAL_, target, intent(in) :: data(grid%localDimsR(1), grid%localDimsR(2), grid%localDimsR(3)) !, centerOfMass(3)
    type(rism3d_solute), intent(in) :: solute
    integer, optional :: o_rank, o_nproc, o_comm
    
    integer :: rank = 0, nproc = 1, comm = 0
    integer :: i,j,k, irank, err, count
    integer, parameter :: dataperline = 3
    integer :: unit, iostat

#ifdef MPI
    _REAL_, pointer :: wrk_data(:, :, :) => NULL()
#endif /*MPI*/

    integer :: id
    _REAL_ :: minValue, maxValue, meanValue, rmsd
    logical, parameter :: bigEndian = ichar(transfer(1,'a')) == 0
    ! Up to 80-character long label describing file origin.
    character(len=*), parameter :: amberLabel = 'Amber 3D-RISM MRC map volumetric data.'
    
#ifdef RISM_DEBUG
    write(0, *) "writeMRC", rank
    call flush(6)
#endif /*RISM_DEBUG*/
    
    unit = freeUnit()
    if (present(o_rank)) rank = o_rank
    if (present(o_nproc)) nproc = o_nproc
    if (present(o_comm)) comm = o_comm
    if (rank == 0) then
       ! Unfortunately gfortran does not support form='BINARY' for
       ! compiler-independent binary output, but that probably won't
       ! be an issue for mrc readers.
       ! if (gfortran) then
       open(unit=unit, file=file, iostat=iostat, access='stream', status='replace', form='unformatted')
       ! else ! Intel Fortran
       !     open(unit=unit, file=file, iostat=iostat, access='sequential', status='replace', form='binary')
       ! end if
       if (iostat /= 0) then
          call rism_report_error("opening "//trim(file))
       end if
    end if
    ! Minimum, maximum, and mean density values.
#if defined(MPI)
    call MPI_REDUCE(minval(data), minValue, 1, MPI_DOUBLE, MPI_MIN, 0, comm, err)
    call MPI_REDUCE(maxval(data), maxValue, 1, MPI_DOUBLE, MPI_MAX, 0, comm, err)       
    call MPI_ALLREDUCE(sum(data), meanValue, 1, MPI_DOUBLE, MPI_SUM, comm, err)       
    call MPI_ALLREDUCE(size(data), count, 1, MPI_INTEGER, MPI_SUM, comm, err)
    meanValue = meanValue/count
    call MPI_REDUCE(sum((meanValue - data)**2), rmsd, 1, MPI_DOUBLE, MPI_SUM, 0, comm, err)
    rmsd = sqrt(rmsd/count)
#else
    minValue = minval(data)
    maxValue = maxval(data)
    meanValue = sum(data) / size(data)
    rmsd = sqrt(sum((meanValue - data)**2) / size(data))
#endif /*defined(MPI)*/
    if (rank == 0) then
       ! Write header.

       ! Number of columns, rows, and sections (fastest to slowest changing).
       ! NX, NY, NZ
       write(unit) int(grid%globalDimsR, int32)

       ! Since values are stored as reals, mode == 2.
       ! MODE
       write(unit) int(2, int32)

       ! There is no offset for column, row, or section.
       ! NXSTART, NYSTART, NZSTART
       write(unit) int((/ 0, 0, 0 /), int32)
       
       ! Number of intervals along X, Y, Z.
       ! MX, MY, MZ
       write(unit) int(grid%globalDimsR, int32)

       ! Cell dimensions (Angstroms).
       ! CELLA
       write(unit) real(grid%boxLength, real32)

       ! Cell angles (degrees).
       ! CELLB
       write(unit) real(grid%unitCellAngles * 180 / PI, real32)

       ! Map column, rows, sects to X, Y, Z (1, 2, 3).
       ! MAPC, MAPR, MAPs
       write(unit) int((/ 1, 2, 3 /), int32)

       ! rmsd = sqrt(rmsd / size(data))
       ! meanValue = totalValue / size(data)
       ! DMIN, DMAX, DMEAN
       write(unit) real(minValue, real32), real(maxValue, real32), &
                   real(meanValue, real32)

       ! Space group number.  We assume P 1.
       ! ISPG
       write(unit) int(1, int32)

       ! Number of bytes used for storing symmetry operators.
       ! In our case, none.
       ! NSYMBT
       write(unit) int(0, int32)

       ! extra space used for anything - 0 by default
       ! The format is a bit confusing here.  It indicates 25 words
       ! but it is the third and fourth are EXTTYP AND NVERSION
       ! EXTRA
       do i=1, 2
          write(unit) int(i, int32)
       end do

       ! code for the type of extended header. One of CCP4, MCRO, SERI, AGAR, FEI1, HDF5
       ! We don't use the extended header so it shouldn't matter
       ! EXTTYP
       write(unit) 'MRCO'

       ! version of the MRC format
       ! The version of the MRC format that the file adheres to, specified as a 32-bit integer and calculated as:
       ! Year * 10 + version within the year (base 0)
       ! For the current format change, the value would be 20140. 
       ! NVERSION
       write(unit) int(20140, int32)
       do i=1, 21
          write(unit) int(i, int32)
       end do

       ! The rest of EXTRA, words 29-49
       
#if 1
       ! phase origin (pixels) or origin of subvolume (A)
       ! ORIGIN
       ! dac: origin is always zero for periodic code:
       write(unit) real((/ 0., 0., 0. /), real32)
#else
       ! Here is Tyler's code from 2/24/23, which tries to find an optimum
       ! origin for visualization:

       ! For periodic solutes, there is no translation, but the
       ! best view should have the CoM of the solute somewhere in
       ! the grid. Translate the origin to choose the periodic
       ! image that contains the CoM.

       ! 1. Transform the CoM to fractional coordinates of the unit cell
       nCellVec = transpose( matmul( grid%unitCellVectorsK, &
                  (solute%centerOfMass, [3,1]) ) )
       ! 2. Round down to the nearest whole unit cell vector
       nCellVec = floor(nCellVec)
       ! 3. Set the origin to whatever multiple of each unit cell vectors
       write(unit) real( matmul(nCellVec, grid%unitCellVectorsR), real32)
#endif
       
       ! Character string 'MAP ' to identify file type.
       ! MAP
       write(unit) 'MAP '

       ! Machine stamp indicating endianness.
       ! MACHST
       if (bigEndian) then
          ! 0x11 0x11 0x00 0x00
          write(unit) int(z'11110000', int32)
       else
          ! 0x44 0x44 0x00 0x00 but 0x44 0x41 0x00 0x00 (z'00004144') is also ok
          write(unit) int(z'00004444', int32)
       end if
       
       ! RMS deviation of map from mean density.
       ! RMS
       write(unit) real(rmsd, real32)

       ! Number of labels being used.
       ! NLABL
       write(unit) int(1, int32)

       ! Ten 80-character labels.
       ! LABEL
       write(unit) amberLabel
       do id = 1, (9 * 80) + (80 - len(amberLabel))
          write(unit) int(0, int8)
       end do

       ! Symmetry records would go here, but we are not using any.
       
       ! Write volumetric data. This is in column-major format, so
       ! rank-0 always writes out first.
       write(unit) real(data, real32)
    end if
#if defined(MPI)
    ! only rank-0 needs temp data
    if(rank == 0) then
       wrk_data => safemem_realloc(wrk_data, ubound(data, 1), ubound(data, 2), ubound(data, 3), .false.)
    end if
    ! proceed through each process in order and transfer data to
    ! rank-0 for output.  This walks through the z-index.
    do i = 1, nproc-1
       if (rank == 0) then
          call mpi_recv(wrk_data, size(wrk_data), mpi_double, i, 0, comm, status, err)
          write(unit) real(wrk_data,real32)
       elseif(rank == i) then
          call mpi_send(data, size(data), mpi_double, 0, 0, comm, err)
       end if
    end do
    if (safemem_dealloc(wrk_data) /= 0) call rism_report_error("MRC_MAP_WRITE: wrk_data deallocation failed")
#endif /*defined(MPI)*/
  end subroutine rism3d_mrc_map_write
  
end module rism3d_mrc
