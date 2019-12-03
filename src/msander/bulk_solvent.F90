! <compile=optimized>
module bulk_solvent_mod

  implicit none

  ! Arrays to hold the bulk solvent mask
  complex(8), dimension(:), allocatable :: f_mask, f_mask_tmp, mask_bs_grid_t_c, &
                                           mask_bs_grid_t_c_tmp
  ! Bulk solvent mask parameters
  double precision, dimension(16)   :: mask_cell_params
  double precision, dimension(3)    :: mask_grid_steps

  double precision, dimension(:), allocatable :: k_mask, k_scale, mask_cutoffs, &
                                                 s_squared, b_vector_mask

  ! hkl_indexing_bs_mask:     (H, K, L) set represented as a 1D array index of 
  !                               FFT'd bulk solvent mask
  integer, dimension(:), allocatable :: hkl_indexing_bs_mask

  ! Convenient numerical constants
  double precision, parameter :: pi = 3.14159265359, zero = 0.0, k_sol = 0.35, &
                                 b_sol = 46.0, mask_r_shrink = 0.9, &
                                 mask_r_probe = 1.11, d_tolerance = 1.e-10

  ! atom_types:               Type index for each atom, referring not to atom types in the
  !                           standard MD topology but atom types for the SSF calculation
  ! mask_bs_grid:             Array to hold the bulk solvent mask
  ! mask_bs_grid_tmp:         Array used in shrinking the bulk solvent mask when building it
  ! grid_neighbors:           Array of relative positions for any grid neighbors in the bulk
  integer, dimension(:), allocatable :: atom_types, mask_bs_grid, &
                                        mask_bs_grid_tmp, grid_neighbors

  ! Size of the mask grid (number of grid points)
  integer, dimension(4) :: mask_grid_size

  ! grid_neighbors_size:    Number of grid neighbors that are within the shrunken mask cutoff
  integer grid_neighbors_size

contains

  !--------------------------------------------------------------------------------------------
  ! count_reduce:
  !--------------------------------------------------------------------------------------------
  function count_reduce (red_n, factor) result(result)
    integer :: red_n, factor, result

    result = red_n
    do while(mod(result, factor) == 0)
      result = result / factor
    end do

  end function count_reduce

  !--------------------------------------------------------------------------------------------
  ! check_max_prime: find the maximum prime of a number n, given a putative guess max_prime
  !--------------------------------------------------------------------------------------------
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

  !--------------------------------------------------------------------------------------------
  ! adjust_gridding:  adjust number of grid point for FFT productivity
  !                   in this module max prime factor is set to 5
  !--------------------------------------------------------------------------------------------
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

  !--------------------------------------------------------------------------------------------
  ! mod_grid:  compute the three dimensional grid index of an absolute position in the
  !            corresponding linear array.
  !
  ! Arguments:
  !   x:        
  !   nx:       
  !--------------------------------------------------------------------------------------------
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

  subroutine calc_grid_neighbors()

    ! Create list of neighboring grid points
    integer :: low(3), high(3), i, n0, n1, n2, p0, m0, p1, m1, p2, m2, alloc_
    double precision :: x, shrink_truncation_radius_sq, frac(3), dist_sq

    alloc_ = 1
    do i = 1, 3
      x = mask_r_shrink * mask_cell_params(i + 6) * mask_grid_size(i);
      low(i) = int(floor(-x))
      high(i) = int(ceiling(x))
      alloc_ = alloc_ * (high(i) - low(i))
    end do
    allocate(grid_neighbors(alloc_))
    grid_neighbors = 0
    n0 = mask_grid_size(1)
    n1 = mask_grid_size(2)
    n2 = mask_grid_size(3)
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
          dist_sq = mask_cell_params(1)*frac(1)*frac(1) + &
                    mask_cell_params(2)*frac(2)*frac(2) + &
                    mask_cell_params(3)*frac(3)*frac(3) + &
                    mask_cell_params(4)*frac(1)*frac(2) + &
                    mask_cell_params(5)*frac(1)*frac(3) + &
                    mask_cell_params(6)*frac(2)*frac(3)
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
  ! init_bulk_solvent: intialize the bulk solvent mask to one, 'solvent present here.'
  !--------------------------------------------------------------------------------------------
  subroutine init_bulk_solvent(n_atom, NRF,resolution)

    use xray_globals_module, only: unit_cell
    use memory_module, only: i100, ix
    use ml_mod, only: cross
    implicit none
    integer, intent(in) :: n_atom, NRF
    double precision, intent(in) :: resolution

    integer :: i, atomic_number, na, nb, nc
    double precision :: temp_grid, grid_stepX, grid_stepY, grid_stepZ
    double precision :: a,b,c,alpha,beta,gamma
    double precision :: cosa, sina, cosb, sinb, cosg, sing, V
    double precision, dimension(3) :: va, vb, vc, vas, vbs, vcs

    mask_bs_grid = 1
    mask_bs_grid_tmp = 1

    allocate(f_mask(NRF))
    allocate(f_mask_tmp(NRF))
    allocate(k_mask(NRF))
    allocate(k_scale(NRF))
    allocate(s_squared(NRF))
    allocate(hkl_indexing_bs_mask(NRF))

    allocate(atom_types(n_atom))
    allocate(mask_cutoffs(n_atom))
    do i = 1, n_atom

      atomic_number = ix(i100+i)

      if( atomic_number == 6 ) then
          atom_types(i) = 2
          mask_cutoffs(i) = 1.775 + mask_r_probe

      elseif( atomic_number == 17 ) then
          atom_types(i) = 0
          mask_cutoffs(i) = 1.75 + mask_r_probe
      
      elseif ( atomic_number == 7 ) then
          atom_types(i) = 3
          mask_cutoffs(i) = 1.5 + mask_r_probe

      elseif( atomic_number == 11 ) then
          atom_types(i) = 0
          mask_cutoffs(i) = 2.27 + mask_r_probe

      elseif ( atomic_number == 8 ) then
        atom_types(i) = 4
        mask_cutoffs(i) = 1.45 + mask_r_probe

      elseif ( atomic_number == 16 ) then
        atom_types(i) = 5
        mask_cutoffs(i) = 1.8 + mask_r_probe

      else
        atom_types(i) = 1
        mask_cutoffs(i) = 1.2 + mask_r_probe

      endif
    end do

    ! work for grids related to solvent mask:

    ! TODO: next section is duplicated in init_ml()
    temp_grid = resolution / 4.0
    na = adjust_gridding((int(a / temp_grid)/2)*2+1, 5)
    nb = adjust_gridding((int(b / temp_grid)/2)*2+1, 5)
    nc = adjust_gridding((int(c / temp_grid)/2)*2+1, 5)

    a = unit_cell(1)
    b = unit_cell(2)
    c = unit_cell(3)
    alpha = 3.1415926536d0*unit_cell(4)/180.d0
    beta = 3.1415926536d0*unit_cell(5)/180.d0
    gamma = 3.1415926536d0*unit_cell(6)/180.d0
    sina = sin(alpha)
    cosa = cos(alpha)
    sinb = sin(beta)
    cosb = cos(beta)
    sing = sin(gamma)
    cosg = cos(gamma)
    V = a*b*c*sqrt( 1.d0 - cosa**2 - cosb**2 -cosg**2 &
        + 2.d0*cosa*cosb*cosg )

    va = (/ a, zero, zero/)
    vb = (/cosg * b, sing * b, zero/)
    vc = (/cosb * c, (cosa - cosb * cosg)/sing * c, &
    sqrt(1.0 - cosa * cosa - cosb * cosb - cosg * cosg + &
         2.0 * cosa * cosb * cosg)/sing * c/)
    vas = cross(vb, vc)
    vbs = cross(vc, va)
    vcs = cross(va, vb)
    vas(1:3) = vas(1:3) / V
    vbs(1:3) = vbs(1:3) / V
    vcs(1:3) = vcs(1:3) / V

    grid_stepX = a / na
    grid_stepY = vb(2) / nb
    grid_stepZ = vc(3) / nc

    ! Metric matrix (1:6), reciprocal cell lengths (7:9), cart-to-frac matrix (10:15),
    ! cell volume (16)

    mask_cell_params = (/a*a, b*b, c*c, 2*a*b*cos(gamma), 2*a*c*cos(beta), &
                        2*b*c*cos(alpha), norm2(vas), norm2(vbs), norm2(vcs), &
                        1.0 / a, -cos(gamma) / (sin(gamma)* a), &
                        (cos(alpha) * cos(gamma) - cos(beta)) / (V * sin(gamma)) * b * c, &
                        1 / (sin(gamma) * b), &
                        (cos(beta) * cos(gamma) - cos(alpha)) / (V * sin(gamma)) * a * c, &
                        a * b * sin(gamma) / V, V/)

    mask_grid_steps = (/grid_stepX, grid_stepY, grid_stepZ/)
    mask_grid_size = (/na, nb, nc, na*nb*nc/)
    write(6, *) 'resolution', resolution
    write(6, *) 'mask_cell_params', mask_cell_params
    write(6, *) 'mask_grid_steps', mask_grid_steps
    write(6, *) 'mask_grid_size', mask_grid_size
    allocate(mask_bs_grid(mask_grid_size(4)))
    allocate(mask_bs_grid_tmp(mask_grid_size(4)))
    allocate(mask_bs_grid_t_c(mask_grid_size(1) * mask_grid_size(2) * &
                              (mask_grid_size(3)/2 + 1)))
    allocate(mask_bs_grid_t_c_tmp(mask_grid_size(1) * mask_grid_size(2) * &
                                  (mask_grid_size(3) /2 + 1)))


#if 0
    !  need to get hkl here
    do i = 1, NRF
      hkl_indexing_bs_mask(i) = h_as_ih(hkl(i, 1), hkl(i, 2), hkl(i, 3), na, nb, nc)
      if (hkl_indexing_bs_mask(i) == -1) then
        stop 'Miller indices indexing failed'
      end if
    end do
#endif

    call calc_grid_neighbors()

    ! TODO: where should this go?
    !k_mask(i) = k_sol * exp(b_sol * s_squared(i))

    return

  end subroutine init_bulk_solvent

  !--------------------------------------------------------------------------------------------
  ! grid_bulk_solvent: CPU equivalent of the GPU kernel kGrid, colors the bulk solvent mask
  !                    grid according to the presence of non-solvent atoms
  !
  ! Arguments:
  !   n_atom:  the total number of atoms (not used but to define the size of crd)
  !   crd:     atomic coordinates
  !--------------------------------------------------------------------------------------------
  subroutine grid_bulk_solvent(n_atom, crd)

    implicit none
    integer :: tid, n_atom
    double precision :: atomX, atomY, atomZ, dx, dy, dz, cutoff, cutoffsq, &
                        distsq, coas, cobs, cocs
    integer :: x_low, x_high, y_low, y_high, z_low, z_high, i, j, k, index
    double precision :: frac(3)
    double precision :: crd(3, n_atom)

    do tid = 1, n_atom

      ! Cartesian to fractional coordinates
      atomX = mask_cell_params(10) * crd(1, tid) + &
              mask_cell_params(11) * crd(2, tid) + &
              mask_cell_params(12) * crd(3, tid)
      atomY = mask_cell_params(13) * crd(2, tid) + mask_cell_params(14) * crd(3, tid)
      atomZ = mask_cell_params(15) * crd(3, tid)
      cutoff = mask_cutoffs(tid)
      cutoffsq = cutoff * cutoff

      ! Grid box around the atom
      coas = cutoff * mask_cell_params(7)
      x_low = int(floor((atomX - coas) * mask_grid_size(1)))
      x_high = int(ceiling((atomX + coas) * mask_grid_size(1)))
      cobs = cutoff * mask_cell_params(8)
      y_low = int(floor((atomY - cobs) * mask_grid_size(2)))
      y_high = int(ceiling((atomY + cobs) * mask_grid_size(2)))
      cocs = cutoff * mask_cell_params(9)
      z_low = int(floor((atomZ - cocs) * mask_grid_size(3)))
      z_high = int(ceiling((atomZ + cocs) * mask_grid_size(3)))

      ! Grid point is 0 if inside sphere and 1 if outside
      do i = x_low, x_high
        frac(1) = dble(i) / mask_grid_size(1);
        dx = atomX - frac(1)
        do j = y_low, y_high
          frac(2) = dble(j) / mask_grid_size(2);
          dy = atomY - frac(2);
          do k = z_low, z_high
            frac(3) = dble(k) / mask_grid_size(3);
            dz = atomZ - frac(3);
            distsq = mask_cell_params(1)*dx*dx + mask_cell_params(2)*dy*dy + &
                     mask_cell_params(3)*dz*dz + mask_cell_params(4)*dx*dy + &
                     mask_cell_params(5)*dx*dz + mask_cell_params(6)*dy*dz
            if (distsq < cutoffsq) then
              index = mod_grid(k,mask_grid_size(3)) + &
                      mod_grid(j,mask_grid_size(2)) * mask_grid_size(3) + &
                      mod_grid(i,mask_grid_size(1)) * mask_grid_size(2) * mask_grid_size(3) + 1
              mask_bs_grid(index) = 0
              mask_bs_grid_tmp(index) = 0
            end if
          end do
        end do
      end do
    end do
  end subroutine

  !--------------------------------------------------------------------------------------------
  ! shrink_bulk_solvent: shave off the outermost layer of masking, expanding the domain of
  !                      bulk solvent in the simulation.  The grid is generally four times the
  !                      resolution of the X-ray data, implying for a 2.0A resolution
  !                      structure a grid of 0.5A or less.  This is pretty fine.
  !--------------------------------------------------------------------------------------------
  subroutine shrink_bulk_solvent()

    integer :: tid, a, i0, j0, k0, i, j, k
    logical :: not_turned

    do tid = 0, mask_grid_size(4) - 1
      if (mask_bs_grid_tmp(tid + 1) == 0) then
        not_turned = .True.
        a = 1
        i0 =         tid / (mask_grid_size(3) * mask_grid_size(2))
        j0 =     mod(tid, (mask_grid_size(3) * mask_grid_size(2))) / mask_grid_size(3)
        k0 = mod(mod(tid, (mask_grid_size(3) * mask_grid_size(2))) , mask_grid_size(3))
        do while ((not_turned .and. (a .le. grid_neighbors_size)) .eqv. .True.)
          i = grid_neighbors(a) / (mask_grid_size(3) * mask_grid_size(2)) + i0
          j = mod(grid_neighbors(a), (mask_grid_size(3) * mask_grid_size(2))) / &
              mask_grid_size(3) + j0
          k = mod(mod(grid_neighbors(a), (mask_grid_size(3) * mask_grid_size(2))), &
                  mask_grid_size(3)) + k0
          
          ! If there is a "1" grid point in a sphere of shrinkage
          ! radius of a "0" grid point then it's turned to "1" also
          if (mask_bs_grid_tmp(mod_grid(k, mask_grid_size(3)) + &
              mod_grid(j, mask_grid_size(2)) * mask_grid_size(3) + &
              mod_grid(i, mask_grid_size(1)) * mask_grid_size(2) * &
              mask_grid_size(3) + 1) == 1) then
            mask_bs_grid(tid + 1) = 1
            not_turned = .False.
          end if
          a = a + 1
        end do
      end if
    end do
  end subroutine

  !--------------------------------------------------------------------------------------------
  ! fft_bs_mask: FIX ME, THIS DOES NOT WORK.  This is the critical step in getting the CPU
  !              code to operate.  Transforms the bulk solvent (BS) mask into fourier space.
  !--------------------------------------------------------------------------------------------
  subroutine fft_bs_mask()

    use iso_c_binding
    implicit none
    include 'fftw3.f03'
    double precision :: mask_bs_grid_3d(mask_grid_size(1) * mask_grid_size(2) * mask_grid_size(3))
    double precision :: mask_bs_grid_3d_fft(mask_grid_size(1) * mask_grid_size(2) * mask_grid_size(3))
    integer :: i, j, k

    ! CHECK
    ! write(mdout, '(a,3i8)') '| mask_grid_size = ', mask_grid_size(1), mask_grid_size(2), &
    !                         mask_grid_size(3)
    ! END CHECK

#if 0
    call dfftw_plan_dft_r2c_3d(plan_forward, mask_grid_size(3), mask_grid_size(2), &
                               mask_grid_size(1), mask_bs_grid_3d, mask_bs_grid_3d_fft, &
                               FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan_forward, mask_bs_grid_3d, mask_bs_grid_3d_fft)
    mask_bs_grid_t_c = reshape(mask_bs_grid_3d_fft, &
                               (/(mask_grid_size(3) / 2 + 1) * mask_grid_size(2) * &
                                 mask_grid_size(1)/))
    call dfftw_destroy_plan(plan_forward)
#endif
  end subroutine fft_bs_mask
    

end module bulk_solvent_mod
