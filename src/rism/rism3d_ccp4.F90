#include "../include/dprec.fh"

!> Minimal support for binary format CCP4 volumetric data output.
!! The format is described here: http://www.ccp4.ac.uk/html/maplib.html#description
!! It involves a structured file header of 256 longwords, then
!! symmetry information, then the map stored a 3-dimensional array.
!! The header itself is broken into 56 longwords followed by ten
!! 80-character text labels.
!! Supports triclinic unit cells.
!! Does not work with MPI parallel.
module rism3d_ccp4
  use safemem
  use rism_report_c
  implicit none

contains

  !> Write volumetric data to a file in CCP4 format.  When writing in
  !! parallel, each process must call this function with its local
  !! data. Data transfer is handled internally.  We assume
  !! decomposition in the z-axis.
  !! @param[in] file File name to write to.
  !! @param[in] data Data to write in a n(1)*n(2)*n(3) linear array.
  !! @param[in] grid Grid object.
  !! @param[in] solute Solute object.
  !! @param[in] o_rank (optional) MPI process rank.
  !! @param[in] o_nproc (optional) MPI number of processes.
  !! @param[in] o_comm (optional) MPI communicator.
  subroutine rism3d_ccp4_map_write (file, data, grid, solute, o_rank, o_nproc, o_comm)
    use constants_rism, only: PI
    use rism_util, only : freeUnit, rmExPrec
    use rism3d_grid_c
    use rism3d_solute_c
    implicit none
    character(len=*), intent(in) :: file
    type(rism3d_grid), intent(in) :: grid
    _REAL_, target, intent(in) :: data(grid%localDimsR(1), &
        grid%localDimsR(2), grid%localDimsR(3))
    type(rism3d_solute), intent(in) :: solute
    integer, optional :: o_rank, o_nproc, o_comm
    
    integer :: rank = 0, nproc = 1, comm = 0
    integer :: i,j,k, irank, err, count, icount
    integer, parameter :: dataperline = 3
    integer :: unit, iostat

    integer :: id
    _REAL_ :: minValue, maxValue, meanValue, rmsd!, totalValue
    logical, parameter :: bigEndian = ichar(transfer(1,'a')) == 0
    ! Up to 80-character long label describing file origin.
    character(len=*), parameter :: amberLabel = &
        'Amber 3D-RISM CCP4 map volumetric data. Format revision A.'
    
    unit = freeUnit()
    ! for future MPI use:
    if (present(o_rank)) rank = o_rank
    if (present(o_nproc)) nproc = o_nproc
    if (present(o_comm)) comm = o_comm
    if (rank /= 0) return

    open(unit=unit, file=file, iostat=iostat, access='stream', &
         status='replace', form='unformatted')
    if (iostat /= 0) then
       call rism_report_error("opening "//trim(file))
    end if
    ! Write header.

    ! Number of columns, rows, and sections (fastest to slowest changing).
    write(unit) int(grid%localDimsR, 4)
    ! Since values are stored as reals, mode == 2.
    write(unit) int(2, 4)
    ! There is no offset for column, row, or section.
    write(unit) int((/ 0, 0, 0 /), 4)
    ! Number of intervals along X, Y, Z.
    write(unit) int(grid%localDimsR, 4)
    ! Cell dimensions (Angstroms).
    write(unit) real(grid%boxLength, 4)
    ! Cell angles (degrees).
    write(unit) real(grid%unitCellAngles * 180 / PI, 4)
    ! Map column, rows, sects to X, Y, Z (1, 2, 3).
    write(unit) int((/ 1, 2, 3 /), 4)
    ! Minimum, maximum, and mean density values.
    !TODO: This should be MPI-ified.
    minValue = minval(data)
    maxValue = maxval(data)
    meanValue = sum(data) / size(data)
    rmsd = sqrt(sum((meanValue - data)**2) / size(data))
    write(unit) real(minValue, 4), real(maxValue, 4), real(meanValue, 4)
    ! Space group number.  We assume P 1.
    write(unit) int(1, 4)
    ! Number of bytes used for storing symmetry operators.
    ! In our case, none.
    write(unit) int(0, 4)
    ! Flag for skew transformation, where 0 = none, 1 = foll.  Skew
    ! transformation is from standard orthogonal coordinate frame
    ! (as used for atoms) to orthogonal map frame, as:
    !          Xo(map) = S * (Xo(atoms) - t).
    write(unit) int(0, 4)
    ! Skew matrix S (in order S11, S12, S13, S21, etc.) if above
    ! flag is non-zero.
    write(unit) int((/ 0, 0, 0 /), 4)
    write(unit) int((/ 0, 0, 0 /), 4)
    write(unit) int((/ 0, 0, 0 /), 4)
    ! Skew translation t if above flag is non-zero.
    write(unit) int((/ 0, 0, 0 /), 4)
    ! Future use (zero default).
    do id = 1, 15
       write(unit) int(0, 4)
    end do
    ! Character string 'MAP ' to identify file type.
    write(unit) 'MAP '
    ! Machine stamp indicating endianness.
    if (bigEndian) then
       write(unit) real(z'11110000', 4)
    else
       write(unit) real(z'00004144', 4)
    end if
    ! RMS deviation of map from mean density.
    write(unit) real(rmsd, 4)
    ! Number of labels being used.
    write(unit) int(1, 4)
    ! Ten 80-character labels.
    write(unit) amberLabel
    do id = 1, (9 * 80) + (80 - len(amberLabel))
       write(unit) int(0, 1)
    end do

    ! Symmetry records would go here, but we are not using any.
    
    ! Write volumetric data.
    write(unit) real(data, 4)

  end subroutine rism3d_ccp4_map_write
  
end module rism3d_ccp4
