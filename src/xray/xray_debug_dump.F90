#include "../include/assert.fh"

module xray_debug_dump_module
  
  use xray_contracts_module
  use xray_unit_cell_module, only: unit_cell_t
  use xray_pure_utils, only : real_kind
  
  implicit none
  
  private
  
  public :: dump
  public :: load
  
contains
  
  subroutine dump(filename, &
      & hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
      & atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk, atom_atomic_number &
      &)
    
    implicit none
    
    character(len = *), intent(in) :: filename
    integer, intent(in) :: hkl(:, :)
    complex(real_kind), intent(in) :: Fobs(:)
    real(real_kind), intent(in) :: sigma_Fobs(:)
    logical, intent(in) :: work_flag(:)
    class(unit_cell_t), intent(in) :: unit_cell
    real(real_kind), intent(in) :: scatter_coefficients(:, :, :) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
    real(real_kind), intent(in) :: atom_b_factor(:)
    real(real_kind), intent(in) :: atom_occupancy(:)
    integer, intent(in) :: atom_scatter_type(:)
    logical, intent(in) :: atom_is_not_bulk(:)
    integer, intent(in) :: atom_atomic_number(:)
    
    integer :: unit_id = 19
    integer :: i, j

    ASSERT(size(hkl, 1) == 3)
    ASSERT(size(hkl, 2) == size(Fobs))
    ASSERT(size(hkl, 2) == size(sigma_Fobs))
    ASSERT(size(hkl, 2) == size(work_flag))
    ASSERT(size(atom_b_factor) == size(atom_occupancy))
    ASSERT(size(atom_b_factor) == size(atom_scatter_type))
    ASSERT(size(atom_b_factor) == size(atom_atomic_number))
    ASSERT(size(atom_b_factor) == size(atom_is_not_bulk))
    ASSERT(minval(atom_b_factor) > 0)
    ASSERT(all(atom_occupancy <= 1.0))
    ASSERT(all(atom_occupancy >= 0.0))
    ASSERT(minval(atom_scatter_type) >= 1)
    ASSERT(maxval(atom_scatter_type) <= size(scatter_coefficients, 3))
    
    
    open(unit = unit_id, FILE = filename, action = 'WRITE')
    write (unit_id, *) "! n_atom"
    write (unit_id, *) size(atom_b_factor)
    write (unit_id, *) "! i, atom_b_factor(i), atom_occupancy(i), atom_scatter_type(i), atom_is_not_bulk(i), atom_atomic_number(i)"
    do i = 1, size(atom_b_factor)
      write (unit_id, *) i, atom_b_factor(i), atom_occupancy(i), atom_scatter_type(i), atom_is_not_bulk(i), atom_atomic_number(i)
    end do
    write (unit_id, *) "! n_scatter_coefficients, n_scatter_types"
    write (unit_id, *) size(scatter_coefficients, 2), size(scatter_coefficients, 3)
    write (unit_id, *) "! i, j, scatter_coefficients(1:2, i, j)"
    do i = 1, size(scatter_coefficients, 2)
      do j = 1, size(scatter_coefficients, 3)
        write (unit_id, *) i, j, scatter_coefficients(1:2, i, j)
      end do
    end do
    write (unit_id, *) "! a, b, c, alpha, beta, gamma"
    write (unit_id, *) unit_cell%a(), unit_cell%b(), unit_cell%c(), unit_cell%alpha_deg(), unit_cell%beta_deg(), unit_cell%gamma_deg()
    write (unit_id, *) "! n_hkl"
    write (unit_id, *) size(hkl, 2)
    write (unit_id, *) "! i, h(i), k(i), l(i), Fobs(i), sigma_Fobs(i), is_work(i)"
    do i = 1, size(hkl, 2)
      write (unit_id, *) i, hkl(:, i), Fobs(i), sigma_Fobs(i), work_flag(i)
    end do
    close(unit_id)
  end subroutine dump
  
  subroutine load(filename, &
      & hkl, Fobs, sigma_Fobs, work_flag, unit_cell, scatter_coefficients, &
      & atom_b_factor, atom_occupancy, atom_scatter_type, atom_is_not_bulk, atom_atomic_number &
      &)
    character(len = *), intent(in) :: filename
    integer, intent(out), allocatable :: hkl(:, :)
    complex(real_kind), intent(out), allocatable :: Fobs(:)
    real(real_kind), intent(out), allocatable :: sigma_Fobs(:)
    logical, intent(out), allocatable :: work_flag(:)
    class(unit_cell_t), intent(out) :: unit_cell
    real(real_kind), intent(out), allocatable :: scatter_coefficients(:, :, :) ! Fourier coefficients (2,n_scatter_coeffs,n_scatter_types)
    real(real_kind), intent(out), allocatable :: atom_b_factor(:)
    real(real_kind), intent(out), allocatable :: atom_occupancy(:)
    integer, intent(out), allocatable :: atom_scatter_type(:)
    logical, intent(out), allocatable :: atom_is_not_bulk(:)
    integer, intent(out), allocatable :: atom_atomic_number(:)
    
    integer :: n_atom
    integer :: n_hkl
    integer :: n_scatter_coefficients
    integer :: n_scatter_types
    integer :: i, j, i_, j_
    real(real_kind) :: a, b, c, alpha, beta, gamma
    character(len = 1024) :: line
    
    integer :: unit_id = 19
    
    ASSERT(.not. allocated(hkl))
    ASSERT(.not. allocated(Fobs))
    ASSERT(.not. allocated(sigma_Fobs))
    ASSERT(.not. allocated(work_flag))
    ASSERT(.not. allocated(scatter_coefficients))
    ASSERT(.not. allocated(atom_b_factor))
    ASSERT(.not. allocated(atom_occupancy))
    ASSERT(.not. allocated(atom_scatter_type))
    ASSERT(.not. allocated(atom_is_not_bulk))
    ASSERT(.not. allocated(atom_atomic_number))
    
    open(unit = unit_id, FILE = filename, action = 'READ')
    read (unit_id, "(A)") line
    read (unit_id, *) n_atom
    read (unit_id, "(A)") line
    
    allocate(atom_b_factor(n_atom))
    allocate(atom_occupancy(n_atom))
    allocate(atom_scatter_type(n_atom))
    allocate(atom_is_not_bulk(n_atom))
    allocate(atom_atomic_number(n_atom))
    
    do i = 1, n_atom
      read (unit_id, *) i_, atom_b_factor(i), atom_occupancy(i), atom_scatter_type(i), atom_is_not_bulk(i), atom_atomic_number(i)
      ASSERT(i_ == i)
    end do
    read (unit_id, "(A)") line
!    ASSERT(line == "! n_scatter_coefficients, n_scatter_types")
    read (unit_id, *) n_scatter_coefficients, n_scatter_types
    read (unit_id, "(A)") line
!    ASSERT(line == "! i, j, scatter_coefficients(1:2, i, j)")
    allocate(scatter_coefficients(2, n_scatter_coefficients, n_scatter_types))
    do i = 1, n_scatter_coefficients
      do j = 1, n_scatter_types
        read (unit_id, *) i_, j_, scatter_coefficients(1:2, i, j)
        ASSERT(i_ == i)
        ASSERT(j_ == j)
      end do
    end do
    read (unit_id, "(A)") line
!    ASSERT(line == "! a, b, c, alpha, beta, gamma")
    read (unit_id, *) a, b, c, alpha, beta, gamma
    call unit_cell%init(a, b, c, alpha, beta, gamma)
    read (unit_id, "(A)") line
!    ASSERT(line == "! n_hkl")
    read (unit_id, *) n_hkl
    allocate(hkl(3, n_hkl))
    allocate(Fobs(n_hkl))
    allocate(sigma_Fobs(n_hkl))
    allocate(work_flag(n_hkl))
    read (unit_id, "(A)") line
!    ASSERT( line == "! i, h(i), k(i), l(i), Fobs(i), sigma_Fobs(i), is_work(i)")
    do i = 1, n_hkl
      read (unit_id, *) i_, hkl(:, i), Fobs(i), sigma_Fobs(i), work_flag(i)
      ASSERT(i_ == i)
    end do
    close(unit_id)
  
  end subroutine load

end module xray_debug_dump_module
