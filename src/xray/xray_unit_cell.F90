#include "../include/assert.fh"

module xray_unit_cell_module
  
  use xray_pure_utils, only : real_kind
  use xray_contracts_module
  use xray_pure_utils, only : to_degrees, to_radians, cross => cross_product, inverse
  
  implicit none
  private
  
  integer, parameter, private :: rk = real_kind
  
  public :: unit_cell_t
  
  type unit_cell_t
    private
    real(real_kind) :: m_a
    real(real_kind) :: m_b
    real(real_kind) :: m_c
    real(real_kind) :: m_alpha_deg
    real(real_kind) :: m_beta_deg
    real(real_kind) :: m_gamma_deg
    
    real(real_kind) :: m_frac_to_orth(3, 3)
    real(real_kind) :: m_orth_to_frac(3, 3)
    logical :: m_is_orth = .FALSE.
  
  contains
    
    procedure :: get_volume
    procedure :: init
    procedure :: is_orth

    procedure :: as_array
    procedure :: a => get_a
    procedure :: b => get_b
    procedure :: c => get_c
    procedure :: alpha_deg => get_alpha_deg
    procedure :: beta_deg => get_beta_deg
    procedure :: gamma_deg => get_gamma_deg
    procedure :: alpha_rad => get_alpha_rad
    procedure :: beta_rad => get_beta_rad
    procedure :: gamma_rad => get_gamma_rad

    generic :: get_s => impl_get_s_3i, impl_get_s_i3
    generic :: get_s2 => impl_get_s2_3i, impl_get_s2_i3, impl_get_s2_i3n
    generic :: to_frac => impl_to_frac_3, impl_to_frac_3n
    generic :: to_orth => impl_to_orth_3, impl_to_orth_3n
    generic :: to_orth_derivative => impl_to_orth_derivative_3n
    generic :: to_squared_orth_norm => impl_frac_to_orth_norm_square_3r, impl_frac_to_orth_norm_square_r3

    procedure, private :: impl_get_s2_3i
    procedure, private :: impl_get_s2_i3
    procedure, private :: impl_get_s2_i3n
    procedure, private :: impl_get_s_3i
    procedure, private :: impl_get_s_i3
    procedure, private :: impl_to_frac_3
    procedure, private :: impl_to_frac_3n
    procedure, private :: impl_to_orth_3
    procedure, private :: impl_to_orth_3n
    procedure, private :: impl_to_orth_derivative_3n
    procedure, private :: impl_frac_to_orth_norm_square_3r
    procedure, private :: impl_frac_to_orth_norm_square_r3
  end type unit_cell_t

contains
  
  subroutine init(this, a, b, c, alpha_deg, beta_deg, gamma_deg)
    implicit none
    class(unit_cell_t) :: this
    real(real_kind) :: a
    real(real_kind) :: b
    real(real_kind) :: c
    real(real_kind) :: alpha_deg
    real(real_kind) :: beta_deg
    real(real_kind) :: gamma_deg
    
    real(real_kind) :: ucell(3, 3)
    
    this%m_a = a
    this%m_b = b
    this%m_c = c
    this%m_alpha_deg = alpha_deg
    this%m_beta_deg = beta_deg
    this%m_gamma_deg = gamma_deg
    
    this%m_is_orth = (alpha_deg == 90 .and. beta_deg  == 90 .and. gamma_deg == 90)
    
    if (this%is_orth()) then
      ucell(:, :) = 0
      ucell(1, 1) = a
      ucell(2, 2) = b
      ucell(3, 3) = c
    else
      ucell(:, 1) = [ a, 0.0_rk, 0.0_rk]
      ucell(:, 2) = [ b * cos(to_radians(gamma_deg)), b * sin(to_radians(gamma_deg)), 0.0_rk]
      ucell(1, 3) = c * cos(to_radians(beta_deg))
      ucell(2, 3) = (b * c * cos(to_radians(alpha_deg)) - ucell(1, 3) * ucell(1, 2)) / ucell(2, 2)
      ucell(3, 3) = sqrt(c ** 2 - ucell(1, 3) ** 2 - ucell(2, 3) ** 2)
    end if
    
    this%m_frac_to_orth = ucell
    this%m_orth_to_frac = inverse(this%m_frac_to_orth)
  
  end subroutine init

  pure function as_array(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result(6)
    result = [ &
      & this%a(),  this%b(), this%c(), &
      & this%alpha_deg(),  this%beta_deg(), this%gamma_deg() &
    & ]
  end function as_array

  pure function is_orth(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    logical :: result
    result = this%m_is_orth
  end function is_orth
  
  pure function get_volume(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = dot_product(this%m_frac_to_orth(:, 1), &
        & cross(this%m_frac_to_orth(:, 2), this%m_frac_to_orth(:, 3)))
  end function get_volume
  
  pure function get_a(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_a
  end function get_a
  
  pure function get_b(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_b
  end function get_b
  
  pure function get_c(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_c
  end function get_c
  
  pure function get_alpha_deg(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_alpha_deg
  end function get_alpha_deg
  
  pure function get_alpha_rad(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = to_radians(this%alpha_deg())
  end function get_alpha_rad
  
  pure function get_beta_deg(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_beta_deg
  end function get_beta_deg
  
  pure function get_beta_rad(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = to_radians(this%beta_deg())
  end function get_beta_rad
  
  pure function get_gamma_deg(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = this%m_gamma_deg
  end function get_gamma_deg
  
  pure function get_gamma_rad(this) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind) :: result
    result = to_radians(this%gamma_deg())
  end function get_gamma_rad
  
  pure function impl_get_s_3i(this, h, k, l) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    integer, intent(in) :: h
    integer, intent(in) :: k
    integer, intent(in) :: l
    real(real_kind) :: result(3)
    result = this%m_orth_to_frac(1, :) * h &
         & + this%m_orth_to_frac(2, :) * k &
         & + this%m_orth_to_frac(3, :) * l
  end function impl_get_s_3i
  
  pure function impl_get_s_i3(this, hkl) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    integer, intent(in) :: hkl(3)
    real(real_kind) :: result(3)
    result = this%impl_get_s_3i(hkl(1), hkl(2), hkl(3))
  end function impl_get_s_i3
  
  pure function impl_get_s2_3i(this, h, k, l) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    integer, intent(in) :: h
    integer, intent(in) :: k
    integer, intent(in) :: l
    real(real_kind) :: result
    result = sum(this%get_s(h, k, l) ** 2)
  end function impl_get_s2_3i
  
  pure function impl_get_s2_i3(this, hkl) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    integer, intent(in) :: hkl(3)
    real(real_kind) :: result
    result = this%impl_get_s2_3i(hkl(1), hkl(2), hkl(3))
  end function impl_get_s2_i3
  
  function impl_get_s2_i3n(this, hkl) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    integer, intent(in) :: hkl(:, :)
    real(real_kind) :: result(size(hkl, 2))
    
    integer :: i
    
    ASSERT(size(hkl, 1) == 3)
    
    do i = 1, size(hkl, 2)
      result(i) = this%impl_get_s2_i3(hkl(:, i))
    end do
    
  end function impl_get_s2_i3n
  
  
  pure function impl_to_orth_3(this, frac) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: frac(3)
    real(real_kind) :: result(3)
    
    result = matmul(this%m_frac_to_orth, frac)
  end function impl_to_orth_3
  
  function impl_to_orth_3n(this, frac) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: frac(:, :)
    real(real_kind) :: result(size(frac, 1), size(frac, 2))
    ASSERT(size(frac, 1) == 3)
    
    result = matmul(this%m_frac_to_orth, frac)
  end function impl_to_orth_3n
  
  function impl_to_orth_derivative_3n(this, frac_derivative) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: frac_derivative(:, :)
    real(real_kind) :: result(size(frac_derivative, 1), size(frac_derivative, 2))
    ASSERT(size(frac_derivative, 1) == 3)
    
    result = matmul(transpose(this%m_orth_to_frac), frac_derivative)
  end function impl_to_orth_derivative_3n
  
  pure function impl_to_frac_3(this, orth) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: orth(3)
    real(real_kind) :: result(3)
    
    result = matmul(this%m_orth_to_frac, orth)
  end function impl_to_frac_3
  
  function impl_to_frac_3n(this, orth) result(result)
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: orth(:, :)
    real(real_kind) :: result(size(orth, 1), size(orth, 2))
    ASSERT(size(orth, 1) == 3)
    
    result = matmul(this%m_orth_to_frac, orth)
  end function impl_to_frac_3n
  
  pure function impl_frac_to_orth_norm_square_r3(this, frac) result(result)
    ! Returns square of real space norm
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: frac(3)
    real(real_kind) :: result
    
    ! Note: Can be optimzed via pre-calculated metric tensor
    !    result = tensor(1)*fx*fx + tensor(2)*fy*fy + &
    !           & tensor(3)*fz*fz + tensor(4)*fx*fy + &
    !           & tensor(5)*fx*fz + tensor(6)*fy*fz
    result = norm2(this%to_orth(frac)) ** 2
  end function impl_frac_to_orth_norm_square_r3
  
  pure function impl_frac_to_orth_norm_square_3r(this, fx, fy, fz) result(result)
    ! Returns square of real space norm
    implicit none
    class(unit_cell_t), intent(in) :: this
    real(real_kind), intent(in) :: fx
    real(real_kind), intent(in) :: fy
    real(real_kind), intent(in) :: fz
    real(real_kind) :: result
    
    result = this%impl_frac_to_orth_norm_square_r3([fx, fy, fz])
  end function impl_frac_to_orth_norm_square_3r

end module xray_unit_cell_module
