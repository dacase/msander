module xray_bulk_mask_impl_gpu_module
  
  use xray_bulk_mask_data_module, only : f_mask
  use xray_contracts_module
  use xray_pure_utils, only : real_kind
  use xray_unit_cell_module, only : unit_cell_t
  
  implicit none
  private
  
  public :: update_f_bulk
  public :: finalize
  public :: init
  
  interface
    subroutine pmemd_xray_bulk_mask_init_gpu(&
        & grid_dim, &
        & reciprocal_norms, &
        & n_atom, mask_cutoffs, &
        & n_grid_points, mask_bs_grid, &
        & a, b, c, alpha_deg, beta_deg, gamma_deg, &
        & r_shrink, &
        & n_hkl, f_mask, hkl &
      ) bind(C)
      !    int *grid_dim,
      !    double *reciprocal_norms,
      !    int n_atom, double *mask_cutoffs,
      !    int n_grid_points, int *mask_bs_grid,
      !    double a, double b, double c, double alpha_deg, double beta_deg, double gamma_deg,
      !    double r_shrink
      use iso_c_binding
      implicit none
      integer(c_int), intent(in) :: grid_dim(3)
      real(c_double), intent(in) :: reciprocal_norms(3)
      integer(c_int), value :: n_atom
      real(c_double), intent(in), target :: mask_cutoffs(n_atom)
      integer(c_int), value :: n_grid_points
      integer(c_int), intent(inout), target :: mask_bs_grid(n_grid_points)
      real(c_double), value :: a
      real(c_double), value :: b
      real(c_double), value :: c
      real(c_double), value :: alpha_deg
      real(c_double), value :: beta_deg
      real(c_double), value :: gamma_deg
      real(c_double), value :: r_shrink
      integer(c_int), value :: n_hkl
      complex(c_double_complex) :: f_mask(n_hkl)
      integer(c_int) :: hkl(3, n_hkl)
    
    end subroutine pmemd_xray_bulk_mask_init_gpu
  end interface
  
  interface
    subroutine pmemd_xray_bulk_mask_update_f_bulk(n_atom, frac_xyz) bind(C)
      use iso_c_binding
      implicit none
      integer(c_int), value :: n_atom
      real(c_double), intent(in) :: frac_xyz(3, n_atom)
    end subroutine pmemd_xray_bulk_mask_update_f_bulk
  end interface
  
  interface
    subroutine pmemd_xray_bulk_mask_finalize_gpu() bind(C)
      implicit none
    end subroutine pmemd_xray_bulk_mask_finalize_gpu
  end interface

contains
  
  subroutine init(resolution_high, hkl, unit_cell_, atm_atomicnumber, &
      & solvent_mask_adjustment, solvent_mask_probe_radius)
    use xray_bulk_mask_impl_cpu_module, only : cpu_init => init
    implicit none
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell_
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    real(real_kind), intent(in) :: solvent_mask_adjustment
    real(real_kind), intent(in) :: solvent_mask_probe_radius
    
    call cpu_init(resolution_high, hkl, unit_cell_, atm_atomicnumber, &
        & solvent_mask_adjustment, solvent_mask_probe_radius)
    call gpu_init(hkl, solvent_mask_probe_radius)
  end subroutine init
  
  subroutine finalize()
    use xray_bulk_mask_impl_cpu_module, only : cpu_finalize => finalize
    implicit none
    call gpu_finalize()
    call cpu_finalize()
  end subroutine finalize
  
  subroutine update_f_bulk(frac, Fuser)
    use xray_interface2_data_module, only : new_order
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    complex(real_kind), allocatable, intent(in) :: Fuser(:)
    
    if( allocated( Fuser ) ) then
       f_mask = Fuser( new_order )
    else
       call check_precondition(size(frac, 1) == 3)
       call pmemd_xray_bulk_mask_update_f_bulk(size(frac, 2), frac)
    end if
    
  end subroutine update_f_bulk
  
  subroutine gpu_init(hkl, solvent_mask_probe_radius)
    use xray_bulk_mask_data_module
    implicit none
    integer, intent(in) :: hkl(:, :)
    real(real_kind), intent(in) :: solvent_mask_probe_radius
    
    call pmemd_xray_bulk_mask_init_gpu(&
        grid_dim, &
        reciprocal_norms, &
        size(mask_cutoffs), mask_cutoffs, &
        size(mask_bs_grid), mask_bs_grid, &
        & unit_cell%a(), unit_cell%b(), unit_cell%c(), unit_cell%alpha_deg(), unit_cell%beta_deg(), unit_cell%gamma_deg(), &
        solvent_mask_probe_radius, size(f_mask), f_mask, hkl &
    )
  end subroutine gpu_init
  
  subroutine gpu_finalize()
    implicit none
    call pmemd_xray_bulk_mask_finalize_gpu()
  end subroutine gpu_finalize

end module xray_bulk_mask_impl_gpu_module
