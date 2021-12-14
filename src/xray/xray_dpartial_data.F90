module xray_dpartial_data_module
  
  use xray_pure_utils, only : real_kind
  
  implicit none
  public
  
  integer, pointer, save :: hkl(:, :)
  real(real_kind), pointer, save :: mSS4(:)
  complex(real_kind), pointer, save :: Fcalc(:)
  real(real_kind), pointer, save :: abs_Fcalc(:)
  real(real_kind), allocatable, save :: atom_b_factor(:)
  integer, allocatable, save :: atom_scatter_type(:)

end module xray_dpartial_data_module