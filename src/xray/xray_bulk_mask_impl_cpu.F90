! <compile=optimized>
#include "../include/assert.fh"
module xray_bulk_mask_impl_cpu_module

  use xray_pure_utils, only: real_kind
  
  use xray_bulk_mask_data_module
  use xray_contracts_module
  use xray_unit_cell_module
  

  implicit none
  private

  public :: update_f_bulk
  public :: finalize
  public :: init

  double precision, dimension(:), allocatable :: rf_mask, if_mask

contains

  !----------------------------------------------------------------------------
  ! count_reduce:
  !----------------------------------------------------------------------------
  function count_reduce (red_n, factor) result(result)
    integer :: red_n, factor, result

    result = red_n
    do while(mod(result, factor) == 0)
      result = result / factor
    end do

  end function count_reduce

  !----------------------------------------------------------------------------
  ! check_max_prime: find the maximum prime of a number n, 
  !                  given a putative guess max_prime
  !----------------------------------------------------------------------------
  function check_max_prime (max_prime, n) result(result)
    integer :: max_prime, n, factor, n_
    logical :: result

    factor = 3
    n_ = n
    result = .true.
    n_ = count_reduce(n_, 2)
    do while(n_ > 1)
      if (factor > max_prime) then
        result = .false.
      end if
      n_ = count_reduce(n_, factor)
      factor = factor + 2
    end do

  end function check_max_prime

  !----------------------------------------------------------------------------
  ! adjust_gridding:  adjust number of grid point for FFT productivity
  !                   in this module max prime factor is set to 5
  !----------------------------------------------------------------------------
  function adjust_gridding (n, max_prime) result(n_fixed)
    integer :: n, n_, max_prime, n_fixed, mandatory_factor

    if (max_prime < 2) then
      max_prime = 0
    end if
    mandatory_factor = 1
    n_ = (n / mandatory_factor) * mandatory_factor
    if (n_ < n) then
      n_ = n + mandatory_factor
    end if

    do while(check_max_prime(max_prime, n_) .eqv. .false.)
      n_ = n_ + mandatory_factor
    end do
    n_fixed = n_

  end function adjust_gridding

  !----------------------------------------------------------------------------
  ! mod_grid:  compute the three dimensional grid index of an absolute 
  !            position in the corresponding linear array.
  !
  ! Arguments:
  !   x:        
  !   nx:       
  !----------------------------------------------------------------------------
  pure function mod_grid(x, nx)
    integer, intent(in) :: x, nx
    integer :: mod_grid

    mod_grid = x
    if (x < 0) then
      mod_grid = mod_grid + nx
    else
      mod_grid = mod(mod_grid, nx)
    end if

  end function mod_grid

  subroutine calc_grid_neighbors(mask_r_shrink)
    implicit none
    real(real_kind), intent(in) :: mask_r_shrink

    ! Create list of neighboring grid points
    integer :: low(3), high(3), i, n0, n1, n2, p0, m0, p1, m1, p2, m2, alloc_
    double precision :: x, shrink_truncation_radius_sq, frac(3), dist_sq

    alloc_ = 1
    do i = 1, 3
      x = mask_r_shrink * reciprocal_norms(i) * grid_dim(i)
      low(i) = int(floor(-x))
      high(i) = int(ceiling(x))
      alloc_ = alloc_ * (high(i) - low(i) + 1)
    end do
    allocate(grid_neighbors(alloc_))
    grid_neighbors = 0
    n0 = grid_dim(1)
    n1 = grid_dim(2)
    n2 = grid_dim(3)
    shrink_truncation_radius_sq = mask_r_shrink * mask_r_shrink
    grid_neighbors_size = 0
    do p0 = low(1), high(1)
      m0 = p0
      if (m0 < 0) then
        m0 = m0 + n0
      else
        m0 = mod(m0 , n0)
      end if
      frac(1) = dble(p0) / n0
      do p1 = low(2), high(2)
        m1 = p1
        if (m1 < 0) then
          m1 = m1 + n1
        else
          m1 = mod(m1 , n1)
        end if
        frac(2) = dble(p1) / n1
        do p2 = low(3), high(3)
          frac(3) = dble(p2) / n2
          dist_sq = unit_cell%to_squared_orth_norm(frac)
          if (dist_sq < shrink_truncation_radius_sq) then
            m2 = p2
            if (m2 < 0) then
              m2 = m2 + n2
            else
              m2 = mod(m2, n2)
            end if
            grid_neighbors_size = grid_neighbors_size + 1
            grid_neighbors(grid_neighbors_size) = m2 + m1 * n2 + m0 * n2 * n1
          end if
        end do
      end do
    end do
  end subroutine calc_grid_neighbors

  !--------------------------------------------------------------------------------------------
  ! init: intialize the mask to one, 'solvent present here.'
  !--------------------------------------------------------------------------------------------
  subroutine init(resolution_high, hkl, unit_cell_, atm_atomicnumber, &
      & solvent_mask_adjustment, solvent_mask_probe_radius)
    use xray_pure_utils, only: cross => cross_product
    implicit none
    double precision, intent(in) :: resolution_high
    class(unit_cell_t), intent(in) :: unit_cell_
    integer, intent(in) :: hkl(:, :)
    integer, intent(in) :: atm_atomicnumber(:)
    real(real_kind), intent(in) :: solvent_mask_adjustment
    real(real_kind), intent(in) :: solvent_mask_probe_radius

    integer :: num_atoms
    integer :: i, atomic_number, na, nb, nc
    double precision :: temp_grid, grid_stepX, grid_stepY, grid_stepZ
    double precision :: a,b,c,alpha,beta,gamma
    double precision :: cosa, sina, cosb, sinb, cosg, sing, V
    double precision, dimension(3) :: va, vb, vc, vas, vbs, vcs, s
    
    num_atoms = size(atm_atomicnumber)
    unit_cell = unit_cell_

    allocate(f_mask(size(hkl, 2)))
    allocate(hkl_indexing_bs_mask(size(hkl, 2)))

    allocate(atom_types(num_atoms))
    allocate(mask_cutoffs(num_atoms))
    do i = 1, num_atoms
      atomic_number = atm_atomicnumber(i)
      if( atomic_number == 6 ) then
          atom_types(i) = 2
          mask_cutoffs(i) = 1.775 + solvent_mask_adjustment
      elseif( atomic_number == 17 ) then
          atom_types(i) = 0
          mask_cutoffs(i) = 1.75 + solvent_mask_adjustment
      elseif ( atomic_number == 7 ) then
          atom_types(i) = 3
          mask_cutoffs(i) = 1.5 + solvent_mask_adjustment
      elseif( atomic_number == 11 ) then
          atom_types(i) = 0
          mask_cutoffs(i) = 2.27 + solvent_mask_adjustment
      elseif ( atomic_number == 8 ) then
        atom_types(i) = 4
        mask_cutoffs(i) = 1.45 + solvent_mask_adjustment
      elseif ( atomic_number == 16 ) then
        atom_types(i) = 5
        mask_cutoffs(i) = 1.8 + solvent_mask_adjustment
      else
        atom_types(i) = 1
        mask_cutoffs(i) = 1.2 + solvent_mask_adjustment
      endif
    end do

    ! work for grids related to solvent mask:

    ! TODO: next section is duplicated in init_ml()

    a = unit_cell%a()
    b = unit_cell%b()
    c = unit_cell%c()
    alpha = unit_cell%alpha_rad()
    beta = unit_cell%beta_rad()
    gamma = unit_cell%gamma_rad()
    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)
    sing = sin(gamma)
    cosg = cos(gamma)
    V = unit_cell%get_volume()
    
    vas = unit_cell%get_s(1, 0, 0)  ! h=1, k=0, l=0
    vbs = unit_cell%get_s(0, 1, 0)  ! h=0, k=1, l=0
    vcs = unit_cell%get_s(0, 0, 1)  ! h=0, k=0, l=1

    temp_grid = resolution_high / 4.0
    na = adjust_gridding((int(a / temp_grid)/2)*2+1, 5)
    nb = adjust_gridding((int(b / temp_grid)/2)*2+1, 5)
    nc = adjust_gridding((int(c / temp_grid)/2)*2+1, 5)
    
    va = unit_cell%to_orth([1.0_real_kind, 0.0_real_kind, 0.0_real_kind])
    vb = unit_cell%to_orth([0.0_real_kind, 1.0_real_kind, 0.0_real_kind])
    vc = unit_cell%to_orth([0.0_real_kind, 0.0_real_kind, 1.0_real_kind])
    
    grid_stepX = va(1) / na
    grid_stepY = vb(2) / nb
    grid_stepZ = vc(3) / nc

    ! Metric matrix (1:6), reciprocal cell lengths (7:9), 
    ! cart-to-frac matrix (10:15), cell volume (16)

    reciprocal_norms = [norm2(vas), norm2(vbs), norm2(vcs)]

    mask_grid_steps = (/grid_stepX, grid_stepY, grid_stepZ/)
    grid_dim = [na, nb, nc]
    grid_size =  na * nb * nc
#if 0
    write(6, *) 'resolution', resolution
    write(6, *) 'reciprocal_norms', reciprocal_norms
    write(6, *) 'mask_grid_steps', mask_grid_steps
    write(6, *) 'grid_size', grid_dim, grid_size
#endif
    allocate(mask_bs_grid(grid_size))
    allocate(mask_bs_grid_t_c(grid_dim(1) * grid_dim(2) * (grid_dim(3)/2 + 1)))


    do i = 1, size(hkl, 2)
      s(:) = unit_cell%get_s(hkl(:, i))
      
      hkl_indexing_bs_mask(i) = h_as_ih( hkl(1,i), hkl(2,i), hkl(3,i), na, nb, nc)
    end do

    call calc_grid_neighbors(solvent_mask_probe_radius)
    return

  end subroutine init


  subroutine finalize()
    deallocate(atom_types)
    deallocate(f_mask)
    deallocate(grid_neighbors)
    deallocate(hkl_indexing_bs_mask )
    deallocate(mask_bs_grid)
    deallocate(mask_bs_grid_t_c)
    deallocate(mask_cutoffs)
  end subroutine finalize

  pure function wrap_frac(x) result(result)
    real(real_kind), intent(in) :: x
    real(real_kind) :: result
    result = x
    if (x < -0.5) then
      result = x - floor(x)
    else if (x >= 0.5) then
      result = x - ceiling(x)
    end if
  end function wrap_frac
  
  pure function wrap_index(i, n) result(result)
    integer, intent(in) :: i
    integer, intent(in) :: n
    integer :: result
    if (i <= 0) then
      result = i - ((i - n) / n) * n
    else if (i > n) then
      result = i - ((i - 1) / n) * n
    else
      result = i
    end if
  end function wrap_index
  
  pure subroutine calc_grid_bounds(x, dx, grid_size, low, high)
    real(real_kind), intent(in) :: x
    real(real_kind), intent(in) :: dx
    integer, intent(in) :: grid_size
    integer, intent(out) :: low
    integer, intent(out) :: high
    
    low = int(floor((x - dx) * grid_size))
    high = int(ceiling((x + dx) * grid_size))
  end subroutine calc_grid_bounds

  !----------------------------------------------------------------------------
  ! grid_bulk_solvent: colors the bulk solvent mask
  !                    grid according to the presence of non-solvent atoms
  !
  ! Arguments:
  !   n_atom:  the total number of atoms
  !   crd:     atomic coordinates
  !----------------------------------------------------------------------------
  subroutine grid_bulk_solvent(frac)
    
    implicit none
    double precision, intent(in) :: frac(:, :)
    integer :: tid
    double precision :: atomX, atomY, atomZ, dx, dy, dz, cutoff, cutoffsq, distsq
    integer :: x_low, x_high, y_low, y_high, z_low, z_high, i, j, k, index, mdi, mdj, mdk
    real(real_kind) :: frac_crd(3)

    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(mask_cutoffs))

    mask_bs_grid = 1
    
    do tid = 1, size(frac, 2)

      atomX = frac(1, tid)
      atomY = frac(2, tid)
      atomZ = frac(3, tid)
      cutoff = mask_cutoffs(tid)
      cutoffsq = cutoff * cutoff

      ! Grid box around the atom
      call calc_grid_bounds(atomX, cutoff * reciprocal_norms(1), grid_dim(1), x_low, x_high)
      call calc_grid_bounds(atomY, cutoff * reciprocal_norms(2), grid_dim(2), y_low, y_high)
      call calc_grid_bounds(atomZ, cutoff * reciprocal_norms(3), grid_dim(3), z_low, z_high)

      ! Grid point is 0 if inside sphere and 1 if outside
      do i = x_low, x_high
        dx = wrap_frac(atomX - dble(i) / grid_dim(1))
        mdi = wrap_index(i, grid_dim(1)) - 1
        do j = y_low, y_high
          dy = wrap_frac(atomY - dble(j) / grid_dim(2))
          mdj = wrap_index(j, grid_dim(2)) - 1
          do k = z_low, z_high
            dz = wrap_frac(atomZ - dble(k) / grid_dim(3))
            distsq = unit_cell%to_squared_orth_norm(dx, dy, dz)
            
            if (distsq < cutoffsq) then
              mdk = wrap_index(k, grid_dim(3)) - 1
              index = mdk + grid_dim(3) * (mdj + mdi * grid_dim(2)) + 1
              mask_bs_grid(index) = 0
            end if
          end do
        end do
      end do
    end do
  end subroutine

  !----------------------------------------------------------------------------
  ! shrink_bulk_solvent: shave off the outermost layer of masking, expanding 
  !     the domain of bulk solvent in the simulation.  The grid is generally 
  !     four times the resolution of the X-ray data, implying for a 2.0A 
  !     resolution structure a grid of 0.5A or less.  This is pretty fine.
  !----------------------------------------------------------------------------
  subroutine shrink_bulk_solvent()

    integer :: tid, a, i0, j0, k0, i, j, k
    logical :: not_turned
    integer, allocatable ::  mask_bs_grid_tmp(:)
    
    mask_bs_grid_tmp = mask_bs_grid
    
    do tid = 0, grid_size - 1
      if (mask_bs_grid_tmp(tid + 1) == 0) then
        not_turned = .True.
        a = 1
        i0 =         tid / (grid_dim(3) * grid_dim(2))
        j0 =     mod(tid, (grid_dim(3) * grid_dim(2))) / grid_dim(3)
        k0 = mod(mod(tid, (grid_dim(3) * grid_dim(2))) , grid_dim(3))
        do while ((not_turned .and. (a .le. grid_neighbors_size)) .eqv. .True.)
          i = grid_neighbors(a) / (grid_dim(3) * grid_dim(2)) + i0
          j = mod(grid_neighbors(a), (grid_dim(3) * grid_dim(2))) / &
              grid_dim(3) + j0
          k = mod(mod(grid_neighbors(a), (grid_dim(3) * grid_dim(2))), &
                  grid_dim(3)) + k0
          
          ! If there is a "1" grid point in a sphere of shrinkage
          ! radius of a "0" grid point then it's turned to "1" also
          if (mask_bs_grid_tmp(mod_grid(k, grid_dim(3)) + &
              mod_grid(j, grid_dim(2)) * grid_dim(3) + &
              mod_grid(i, grid_dim(1)) * grid_dim(2) * &
              grid_dim(3) + 1) == 1) then
            mask_bs_grid(tid + 1) = 1
            not_turned = .False.
          end if
          a = a + 1
        end do
      end if
    end do
  end subroutine shrink_bulk_solvent

  !----------------------------------------------------------------------------
  ! fft_bs_mask: Transforms the bulk solvent (BS) mask into fourier space.
  !----------------------------------------------------------------------------
  subroutine fft_bs_mask()
    use xray_fft_interface_module, only: dft_i2c_3d
    implicit none
    call dft_i2c_3d(grid_dim, mask_bs_grid, mask_bs_grid_t_c)
  end subroutine fft_bs_mask
  
  !----------------------------------------------------------------------------
  ! h_as_ih:  represent a set of (h, k, l) as a 1D array index of FFT'd bulk 
  !           solvent mask, given grid dimensions (na, nb, nc)
  !----------------------------------------------------------------------------
  function h_as_ih (h, k, l, na, nb, nc) result(ih)
    integer :: h, k, l, na, nb, nc, m, ihh, ihk, ihl, ih
    logical :: error

    error = .false.

    ihh = h
    ihk = k
    ihl = l

    m = (na - 1) / 2
    if (-m > ihh .or. ihh > m) then
      error = .true.
    elseif (ihh < 0) then
      ihh = ihh + na
    end if
    m = (nb - 1) / 2
    if (-m > ihk .or. ihk > m) then
      error = .true.
    elseif (ihk < 0) then
      ihk = ihk + nb
    end if
    m = nc / 2 + 1
    if (0 > ihl .or. h >= m) then
      error = .true.
    end if

    if (error) then
      ih = -1
    else
      ih = ihh * nb * m + ihk * m + ihl
    end if

  end function h_as_ih
  
  
  subroutine update_f_bulk(frac, Fuser)
    
    use xray_interface2_data_module, only : new_order
    implicit none
    real(real_kind), intent(in) :: frac(:, :)
    complex(real_kind), allocatable, intent(in) :: Fuser(:)
    integer :: i
    ASSERT(size(frac, 1) == 3)
    ASSERT(size(frac, 2) == size(atom_types))
    ASSERT(size(frac, 2) == size(mask_cutoffs))
  
    if( allocated(Fuser) ) then
       f_mask = Fuser( new_order )
    else
       call grid_bulk_solvent(frac)
       call shrink_bulk_solvent()
       call fft_bs_mask()

       do i = 1, size(f_mask)
         ! High resolution reflexes are weighted with zero bulk scaling factor `k_bulk`
         ! therefore it should be fine to set them to zero
         if (hkl_indexing_bs_mask(i) /= -1) then
           f_mask(i) = conjg(mask_bs_grid_t_c(hkl_indexing_bs_mask(i)+1)) &
                       * unit_cell%get_volume() / grid_size
         else
           f_mask(i) = 0
         end if
       end do
    end if
  
  end subroutine update_f_bulk

  
end module xray_bulk_mask_impl_cpu_module
