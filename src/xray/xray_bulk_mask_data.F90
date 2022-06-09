module xray_bulk_mask_data_module
  
  use xray_unit_cell_module, only: unit_cell_t
  
  implicit none
  
  public
  
  ! Arrays to hold the bulk solvent mask
  complex(8), dimension(:), allocatable :: mask_bs_grid_t_c
  ! Bulk solvent mask parameters
  double precision, dimension(3)   :: reciprocal_norms ! Norms of reciprocal vectors hkl=100; 010; 001
  double precision, dimension(3)   :: mask_grid_steps
  type(unit_cell_t), save :: unit_cell ! Copy of unit cell
  
  double precision, dimension(:), allocatable :: mask_cutoffs
  complex(8), dimension(:), allocatable :: f_mask
  
  
  ! hkl_indexing_bs_mask:     (H, K, L) set represented as a 1D array index of
  !                               FFT'd bulk solvent mask
  integer, dimension(:), allocatable :: hkl_indexing_bs_mask
  
  ! Convenient numerical constants
  double precision, parameter :: pi = 3.14159265359, zero = 0.0, d_tolerance = 1.e-10
  double precision :: k_sol = 0.35, b_sol = 46.0 ! FIXME: parameters are unused
  
  ! atom_types:    Type index for each atom, referring not to atom types in the
  !                standard MD topology but atom types for the SSF calculation
  ! mask_bs_grid:      Array to hold the bulk solvent mask
  ! grid_neighbors:    Array of relative positions for any grid neighbors
  !                    in the bulk
  integer, dimension(:), allocatable :: atom_types, mask_bs_grid, grid_neighbors
  
  ! Size of the mask grid (number of grid points)
  integer, save :: grid_dim(3)
  integer, save :: grid_size ! Product of grid_dim
  
  ! mask_update_period   how often to update the bs mask (in steps)

  
  ! grid_neighbors_size:    Number of grid neighbors that are within the
  !                         shrunken mask cutoff
  integer grid_neighbors_size

end module xray_bulk_mask_data_module
