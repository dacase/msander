! <compile=optimized>
module ml_mod

  use file_io_dat
  ! use mdin_xray_dat_mod
  ! use mdin_ctrl_dat_mod
  ! use mdin_ewald_dat_mod
  ! use pme_slab_fft_mod
  
  implicit none

  integer :: sfactors_unit = 31, sf_weight = 32
  character(len=50) :: sfactors_name = ''

#if 0
#  include "fftw3.f"
#endif

  ! This is written out with nine lines of five entries, but really just one big line full
  ! of 45 entries.  The data is immediately reshaped into the requisite 9 col x 5 row array.
  double precision, dimension(9, 5), save :: &
  scat_factors = reshape((/ &
                            0.493002,   0.322912,   0.140191,   0.040810,  10.510900, &
                           26.125700,   3.142360,  57.799700,   0.003038,   2.310000, &
                            1.020000,   1.588600,   0.865000,  20.843900,  10.207500, &
                            0.568700,  51.651200,   0.215600,  12.212600,   3.132200, &
                            2.012500,   1.166300,   0.005700,   9.893300,  28.997500, &
                            0.582600, -11.529000,   3.048500,   2.286800,   1.546300, &
                            0.867000,  13.277100,   5.701100,   0.323900,  32.908900, &
                            0.250800,   6.905300,   5.203400,   1.437900,   1.586300, &
                            1.467900,  22.215100,   0.253600,  56.172000,   0.866900 /), &
                         shape(scat_factors)) ! H, C, N, O, S scattering factors, it1992 scattering table

  ! Space group, currently read in from the structure factors file.
  ! Not really used, since we use the expanded-to-P1 structure factors.
  character(len = 16) :: space_group

  ! Arrays to hold calculated structure factors and the bulk solvent mask
  complex(8), dimension(:), allocatable :: f_calc, f_mask, f_mask_tmp, mask_bs_grid_t_c, &
                                           mask_bs_grid_t_c_tmp

  ! Square root of -1 to calculate complex exponents
  complex(8), parameter :: imaginary_i = cmplx(0., 1.)

  ! Parameters for the structure factor computation, many of them supplied by the user
  ! adp_energy_weight, adp_energy_weight_incr are not used currently.
  ! Potentially, they can be used if we want to restrain B-factors
  double precision :: pseudo_energy_weight, pseudo_energy_weight_incr, scale_xray, &
                      F_obs_sq_sum_work, F_obs_sq_sum_free, r_work_factor_denominator, &
                      r_free_factor_denominator, adp_energy_weight, adp_energy_weight_incr, &
                      cell_volume, xray_r_work, xray_r_free

  ! Unit cell dimensions
  double precision, dimension(6)    :: uc_dimensions

  ! Bulk solvent mask parameters
  double precision, dimension(16)   :: mask_cell_params
  double precision, dimension(3)    :: mask_grid_steps

  ! 7 x 7 transformation matrices to anisotropically scale structure factors from calculated
  ! to experimental
  double precision, dimension(7, 7) :: Ucryst, MUcryst_inv

  ! Arrays used in the computation of structure factors
  double precision, dimension(:), allocatable :: f_obs, f_obs_weight, f_obs_sigmas, &
                                                 f_calc_abs, s_squared, BFactors, k_mask, &
                                                 k_scale, alpha_array, &
                                                 beta_array, delta_array, frc_adp, &
                                                 frc_adp_a_priori, mask_cutoffs, &
                                                 b_vector_base
  double precision, dimension(:,:), allocatable ::  s, scat_factors_precalc

  ! Arrays used in the estimation of maximum likelihood parameters,
  ! size is equal to the number of resolution bins/zones
  double precision, dimension(:), allocatable :: A_in_zones, B_in_zones, &
                                                 C_in_zones, q_in_zones, &
                                                 alpha_beta_bj, alpha_beta_OmegaI, alpha_beta_wi, &
                                                 t_optimal, alpha_in_zones, beta_in_zones

  ! Structure factor force computation scratch arrays
  double precision, dimension(:,:), allocatable :: frc_sf

  ! Convenient numerical constants
  double precision, parameter :: pi = 3.14159265359, zero = 0.0, k_sol = 0.35, &
                                 b_sol = 46.0, mask_r_shrink = 0.9, &
                                 mask_r_probe = 1.11, d_tolerance = 1.e-10

  ! NRF:                    The number of reflections
  ! NRF_work, NRF_work_sq:  Number of work reflections, and square thereof
  ! NRF_free:               Number of free reflections
  ! NAT_for_mask:           Number of mask atoms
  ! BFactors_len:           Number of B factors (given in the structure factors file),
  !                         has to be equal to NAT_for_mask, otherwise warning is thrown
  ! call_est:               Tracking counter for pseudo-energy weight estimation
  ! N_steps:                The total number of simulation steps (N_steps is set directly to
  !                         nstlim). Used to potentially control the restraints weight
  ! starting_N_step:        Step number to start on. Similar purpose as for N_steps
  ! total_N_steps:          Step number to finish on. Similar purpose as for N_steps
  ! mask_update_frequency:  Solvent mask update freuency
  ! grid_neighbors_size:    Number of grid neighbors that are within the shrunken mask cutoff
  ! n_bins:                 Number of reflections per resolution bin
  integer :: NRF, NRF_work, NRF_work_sq, NRF_free, NAT_for_mask, BFactors_len, &
             call_est, N_steps, starting_N_step, total_N_steps, &
             mask_update_frequency, grid_neighbors_size, n_bins

  ! Threshold for atom types that will be included in the calculation.  These
  ! type indexes are specific to the SF calculation with hard-wired parameters,
  ! for the moment. Hydrogens are included into structure factors if 0, otherwise 1
  integer, parameter :: h_threshold = 0

  ! atom_types:               Type index for each atom, referring not to atom types in the
  !                           standard MD topology but atom types for the SSF calculation
  ! mask_bs_grid:             Array to hold the bulk solvent mask
  ! mask_bs_grid_tmp:         Array used in shrinking the bulk solvent mask when building it
  ! hkl_indexing_bs_mask:     (H, K, L) set represented as a 1D array index of FFT'd bulk solvent mask
  ! h_sq, k_sq, l_sq:         Squares of H, K, and L indices (still integers)
  ! hk, kl, hl:               Products of H, K, and L indices
  ! grid_neighbors:           Array of relative positions for any grid neighbors in the bulk
  !                           solvent mask
  ! bins_work_population:     Number of work reflections in each resolution zone
  ! bins_free_population:     Number of free reflections in each resolution zone
  ! bins_free_start_indices:  Stores the indices of first free reflection in each resolution zone
  ! reflection_bin:           HKL's number of the corresponding resolution zone
  integer, dimension(:), allocatable :: atom_types, mask_bs_grid, &
                                        mask_bs_grid_tmp, hkl_indexing_bs_mask, &
                                        h_sq, k_sq, l_sq, hk, kl, hl, grid_neighbors, &
                                        bins_work_population, bins_free_population, &
                                        bins_free_start_indices, reflection_bin

  ! hkl:                      N columns by 3 rows array containing HKL indices of reflections
  ! bins_work_reflections:    Used to re-sort work reflections into the original order based on
  !                           resolution bin and relative position in it
  ! bins_free_reflections:    Used to re-sort free reflections into the original order based on
  !                           resolution bin and relative position in it
  integer, dimension(:,:), allocatable :: hkl, bins_work_reflections, bins_free_reflections

  ! Size of the mask grid (number of grid points)
  integer, dimension(4) :: mask_grid_size

  ! FFTW encoding... con this go away?
  !integer*8 plan_forward(1)

contains

  !--------------------------------------------------------------------------------------------
  ! cross:   Compute the cross product of vectors a and b.  Return the cross product.
  !--------------------------------------------------------------------------------------------
  function cross(a, b) result(cross_product)

    double precision, dimension(3) :: cross_product
    double precision, dimension(3), intent(in) :: a, b

    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

  
  !--------------------------------------------------------------------------------------------
  ! get_scat_factor:  get the scattering factor associated with an atom
  !
  ! Arguments:
  !   atom_type:    the type of the atom (indexed according to SF-specific array, not the
  !                 original topology numerical condification of atom types)
  !   s_i:          the atom index for which to get a scattering factor
  !--------------------------------------------------------------------------------------------
  function get_scat_factor(atom_type, s_i) result(scat_factor)

    double precision :: scat_factor
    integer :: i
    integer, intent(in) :: atom_type, s_i

    scat_factor = 0.0
    if (atom_type > h_threshold) then
      do i = 1, 4
        scat_factor = scat_factor + &
                      (scat_factors(i, atom_type) * exp(scat_factors(i + 4, atom_type) * &
                       s_squared(s_i)))
      end do
      scat_factor = scat_factor + scat_factors(9, atom_type)
    endif

  end function get_scat_factor

  !--------------------------------------------------------------------------------------------
  ! i1_over_i0:   approximation of modified Bessel functions of the first kind:
  !               first order over zero order
  !--------------------------------------------------------------------------------------------
  function i1_over_i0 (x) result(result)

    double precision, dimension(7) :: p  = (/  1.00000000,  3.51562290,  3.08994240, &
                                               1.20672920,  0.26597320,  0.03607680, &
                                               0.00458130 /), &
                                      pp = (/  0.50000000,  0.87890594,  0.51498869, &
                                               0.15084934,  0.02658733,  0.00301532, &
                                               0.00032411 /)
    double precision, dimension(9) :: q  = (/  0.39894228,  0.01328592,  0.00225319, &
                                              -0.00157565,  0.00916281, -0.02057706, &
                                               0.02635537, -0.01647633,  0.00392377 /), &
                                      qq = (/  0.39894228, -0.03988024, -0.00362018, &
                                               0.00163801, -0.01031555,  0.02282967, &
                                              -0.02895312,  0.01787654, -0.00420059 /)
    double precision :: x, y, result, be1, be0, pow_y_i, abs_x
    integer :: i
    be1 = 0
    be0 = 0
    abs_x = x
    if (abs_x < 0) then
      abs_x = -abs_x
    end if
    if (abs_x < 3.75) then
      y = x / 3.75
      y = y * y
      pow_y_i = 1
      do i = 1, 7
        be0 = be0 + p(i) * pow_y_i
        be1 = be1 + x * pp(i) * pow_y_i
        pow_y_i = pow_y_i * y
      end do
    else
      y = 3.75 / abs_x
      pow_y_i = 1
      do i = 1, 9
        be0 = be0 + q(i) * pow_y_i
        be1 = be1 + qq(i) * pow_y_i
        pow_y_i = pow_y_i * y
      end do
    end if
    result = be1/be0;
    if (x < 0.0 .and. result > 0.0) then
      result = -result
    end if

  end function i1_over_i0

  !--------------------------------------------------------------------------------------------
  ! ln_of_i0:   approximation of logarithm of modified Bessel function of the first kind,
  !             zero order
  !--------------------------------------------------------------------------------------------
  function ln_of_i0 (x) result(bessel_lni0)

    double precision, dimension(NRF_work) :: x, bessel_lni0
    double precision :: y, abs_x
    integer :: i

    do i = 1, NRF_work
      abs_x = x(i)
      if (abs_x < 0) then
        abs_x = -abs_x
      end if
      if (abs_x/3.75 < 1.0) then
        y = x(i) / 3.75
        y = y * y
        bessel_lni0(i) = 1.0 + y*(3.5156229 + &
                                  y*(3.0899424 + &
                                     y*(1.2067492 + &
                                        y*(0.2659732 + &
                                           y*(0.0360768 + y*0.0045813)))))
        bessel_lni0(i) = log(bessel_lni0(i));
      else
        y = 3.75 / abs_x
        y = 0.39894228 + y*(0.01328592 + &
                            y*(0.00225319 + &
                               y*(-0.00157565 + &
                                  y*(0.00916281 + &
                                     y*(-0.02057706 + &
                                        y*(0.02635537 + &
                                           y*(-0.01647633 + y*0.00392377)))))))
        bessel_lni0(i) = log(y) + abs_x - 0.5 * log(abs_x);
      end if
    end do

  end function ln_of_i0

  !--------------------------------------------------------------------------------------------
  ! h_as_ih:  represent a set of (h, k, l) as a 1D array index of FFT'd bulk solvent mask,
  !           given grid dimensions (na, nb, nc)
  !--------------------------------------------------------------------------------------------
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
  ! init_sf: gateway to structure factor computations.  Fires off open_sf_file and orders up
  !          reading for lots of parameters in the calculation.
  !
  ! Arguments:
  !   crd:    
  !--------------------------------------------------------------------------------------------
  subroutine init_sf(n_atom, nstlim)
    ! xxx prmtop_dat_mod
    ! xxx state_info_mod

    implicit none

    character(4) :: name
    integer :: i, n_atom, nstlim
    
    N_steps = nstlim
    allocate(atom_types(n_atom))
    allocate(mask_cutoffs(n_atom))
    do i = 1, n_atom

      ! FIXME: need to get atom names visible to this routine
      ! name = atm_igraph(i)

      if (name(1:1) == 'C') then
        if (name == 'Cl-') then
          atom_types(i) = 0
          mask_cutoffs(i) = 1.75 + mask_r_probe
          write(mdout, *) "Found CL atom"
        else
          atom_types(i) = 2
          mask_cutoffs(i) = 1.775 + mask_r_probe
        endif
      elseif (name(1:1) == 'N') then
        if (name == 'Na+') then
          atom_types(i) = 0
          mask_cutoffs(i) = 2.27 + mask_r_probe
          write(mdout, *) "Found NA atom"
        else
          atom_types(i) = 3
          mask_cutoffs(i) = 1.5 + mask_r_probe
        endif
      elseif (name(1:1) == 'O') then
        atom_types(i) = 4
        mask_cutoffs(i) = 1.45 + mask_r_probe
      elseif (name(1:1) == 'S') then
        atom_types(i) = 5
        mask_cutoffs(i) = 1.8 + mask_r_probe
      else
        atom_types(i) = 1
        mask_cutoffs(i) = 1.2 + mask_r_probe
      endif
      if (name == 'OXT') then ! Restrain protein only
        NAT_for_mask = i
        write(mdout, *) "NAT_for_mask = ", NAT_for_mask
      endif
    end do

    allocate(frc_sf(3, NAT_for_mask))
    allocate(frc_adp(NAT_for_mask))
    allocate(frc_adp_a_priori(NAT_for_mask))
    call open_sf_file(sfactors_unit, sfactors_name)
    write(mdout, *) "Parsed structure factors file"
    if (BFactors_len == NAT_for_mask) then
      write(mdout, *) "All good, Bfactors are read"
    else
      write(mdout, *) "Warning, not enough Bfactors"
    endif

    call_est = 1
    
    return

  end subroutine init_sf

  !--------------------------------------------------------------------------------------------
  ! finalize_ml_mod:  clean up the structure factors computation before termination of pmemd
  !--------------------------------------------------------------------------------------------
  subroutine finalize_ml_mod()

    implicit none

    complex(8), dimension(:), allocatable :: f_calc_tmp
    integer :: i, j, counter_sort, index_sort
        
    allocate(f_calc_tmp(NRF))
    
    counter_sort = 1
    do i = 1, n_bins
      do j = 1, bins_work_population(i)
        index_sort = bins_work_reflections(i, j)
        f_calc_tmp(index_sort) = f_calc(counter_sort)
        counter_sort = counter_sort + 1
      end do
    end do
    do i = 1, n_bins
      do j = 1, bins_free_population(i)
        index_sort = bins_free_reflections(i, j)
        f_calc_tmp(index_sort) = f_calc(counter_sort)
        counter_sort = counter_sort + 1
      end do
    end do

    open (unit=sf_weight, file="structure_factors.txt", status='replace', action='write')
    write (sf_weight, *) f_calc_tmp
    close (sf_weight)
    open (unit=sf_weight, file="final_weight.txt", status='replace', action='write')
    write (sf_weight, *) starting_N_step + N_steps, total_N_steps
    close (sf_weight)

  end subroutine finalize_ml_mod

  !--------------------------------------------------------------------------------------------
  ! file_line_count:  calculate the number of structure factors in the input file, stored in
  !                   line_num
  !--------------------------------------------------------------------------------------------
  subroutine file_line_count(input_unit, file_name, line_num)

    implicit none

    character(len = *) :: file_name
    character(len = 1023) :: line
    integer(kind = 4) :: ierror, input_status, input_unit, line_num
    
    write(mdout, *) "Opening structure factors file"

    open(unit = input_unit, file = file_name, status = 'old', iostat = input_status)
    if (input_status /= 0) then
      write ( *, '(a)') 'Fatal error! Cannot read structure factors file.'
      stop
    end if
    line_num = 0
    do
      read(input_unit, '(a)', iostat = input_status) line
      if (input_status /= 0) then
        ierror = line_num
        exit
      end if
      line_num = line_num + 1
    end do
    close(unit = input_unit)
    line_num = line_num - 6 ! gap for space group, dimensions, k-constants,
                            ! mask update frequency, BFactors_len, and Bfactors

    return

  end subroutine file_line_count

  !--------------------------------------------------------------------------------------------
  ! square:  compute the square of a number
  !--------------------------------------------------------------------------------------------
  pure function square(x)
    double precision, intent(in) :: x
    double precision :: square
    square = x * x
  end function

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

  !--------------------------------------------------------------------------------------------
  ! inverse: Matrix inverse, based on Doolittle LU factorization for Ax=b by Alex G.,
  !          December 2009.  The original matrix a(n,n) will be destroyed.
  !
  ! Arguments:
  !   a:       (n x n) array of coefficients for square matrix A
  !   n:       dimension of A
  !   c:       (n x n) inverse of matrix A
  !--------------------------------------------------------------------------------------------
  subroutine inverse(a,c,n)

    implicit none
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 allows such operations on matrices
    L = 0.0
    U = 0.0
    b = 0.0

    ! step 1: forward elimination
    do k = 1, n-1
      do i = k+1, n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j = k+1, n
          a(i,j) = a(i,j)-coeff*a(k,j)
        end do
      end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i = 1, n
      L(i,i) = 1.0
    end do

    ! U matrix is the upper triangular part of A
    do j = 1, n
      do i = 1, j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k = 1, n
      b(k) = 1.0
      d(1) = b(1)

      ! Step 3a: Solve Ld=b using the forward substitution
      do i = 2, n
        d(i) = b(i)
        do j = 1, i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do

      ! Step 3b: Solve Ux=d using the back substitution
      x(n) = d(n) / U(n,n)
      do i = n-1, 1, -1
        x(i) = d(i)
        do j = n, i+1, -1
          x(i) = x(i) - U(i,j)*x(j)
        end do
        x(i) = x(i) / u(i,i)
      end do

      ! Step 3c: fill the solutions x(n) into column k of C
      do i = 1, n
        c(i,k) = x(i)
      end do
      b(k) = 0.0
    end do
  end subroutine inverse

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
  ! bubble_sort:  go bubble sort!  Sort a list of real numbers a.
  !--------------------------------------------------------------------------------------------
  subroutine bubble_sort(a)
    double precision, dimension(:) :: a
    double precision :: temp
    integer :: i, j
    logical :: swapped

    do j = size(a)-1, 1, -1
      swapped = .FALSE.
      do i = 1, j
        if (a(i) > a(i+1)) then
          temp = a(i)
          a(i) = a(i+1)
          a(i+1) = temp
          swapped = .TRUE.
        end if
      end do
      if (.NOT. swapped) exit
    end do
  end subroutine bubble_sort

  !--------------------------------------------------------------------------------------------
  ! open_sf_file: open and read, unpack parameters, B-factors, and structure factors from a
  !               file.
  !--------------------------------------------------------------------------------------------
  subroutine open_sf_file(sfactors_unit, sfactors_name)

    implicit none

    character(*) :: sfactors_name
    double precision :: a, b, c, alpha, beta, gamma, V, grid_stepX, grid_stepY, grid_stepZ, &
                        temp_grid, resolution, d, d_star, reflections_per_bin, fo_fo
    double precision, dimension(3) :: va, vb, vc, vas, vbs, vcs
    double precision, dimension(:), allocatable :: f_obs_tmp, f_obs_sigma_tmp, d_star_sq, &
                                                   d_star_sq_sorted, bin_limits
    integer, dimension(:), allocatable :: r_free_flag_array, counter_w, counter_f
    integer (kind = 4) :: file_status
    integer :: d_i, i, j, k, na, nb, nc
    integer :: r_free_counter, r_free_flag, counter_sort, index_sort
    integer :: sfactors_unit
    integer, dimension(:,:), allocatable :: hkl_tmp
    
    call file_line_count(sfactors_unit, sfactors_name, NRF)
    write(mdout, '(A,I6,A)') 'Counted a total of', NRF, 'file lines'

    open(unit = sfactors_unit, file = sfactors_name, status = 'old', iostat = file_status)

    if (file_status /= 0) then
      write ( *, '(a)') 'Fatal error with structure factors file.'
      sfactors_unit = -1
      stop
    else
      allocate(f_calc(NRF))
      allocate(f_calc_abs(NRF))
      allocate(f_mask(NRF))
      allocate(f_mask_tmp(NRF))
      allocate(f_obs(NRF))
      allocate(f_obs_tmp(NRF))
      allocate(f_obs_weight(NRF))
      allocate(f_obs_sigmas(NRF))
      allocate(f_obs_sigma_tmp(NRF))
      allocate(r_free_flag_array(NRF))
      allocate(d_star_sq(NRF))
      allocate(reflection_bin(NRF))
      allocate(hkl(NRF, 3))
      allocate(hkl_indexing_bs_mask(NRF))
      allocate(hkl_tmp(NRF, 3))
      allocate(h_sq(NRF))
      allocate(k_sq(NRF))
      allocate(l_sq(NRF))
      allocate(hk(NRF))
      allocate(hl(NRF))
      allocate(kl(NRF))
      allocate(k_mask(NRF))
      allocate(k_scale(NRF))
      allocate(delta_array(NRF))
      allocate(alpha_array(NRF))
      alpha_array = 1.0
      allocate(beta_array(NRF))
      beta_array = 1.0
      allocate(s(3, NRF))
      allocate(s_squared(NRF))
      allocate(scat_factors_precalc(5, NRF))
      do i = 1, 7
        do j = 1, 7
          Ucryst(i, j) = zero
        end do
      end do
      F_obs_sq_sum_work = zero
      F_obs_sq_sum_free = zero
      r_work_factor_denominator = zero
      r_free_factor_denominator = zero
      r_free_counter = 0
      NRF_work = 0
      read(sfactors_unit, '(a)') space_group
      write(mdout, *) '| Original space group is ',space_group
      read(sfactors_unit, *) a, b, c, alpha, beta, gamma
      uc_dimensions = (/ a, b, c, alpha, beta, gamma /)
      read(sfactors_unit, *) pseudo_energy_weight, pseudo_energy_weight_incr
      write(mdout, '(a,f9.4,a,f9.4)') '| Initial weights = ', pseudo_energy_weight, ' ', &
                                      pseudo_energy_weight_incr
      read(sfactors_unit, *) mask_update_frequency
      write(mdout, *) '| Mask update frequency = ', mask_update_frequency
      alpha = alpha * pi / 180
      beta = beta * pi / 180
      gamma = gamma * pi / 180
      V = a * b * c * &
      sqrt(1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma) + &
           2.0 * cos(alpha) * cos(beta) * cos(gamma))
      va = (/ a, zero, zero/)
      vb = (/cos(gamma) * b, sin(gamma) * b, zero/)
      vc = (/cos(beta) * c, (cos(alpha) - cos(beta) * cos(gamma))/sin(gamma) * c, &
      sqrt(1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma) + &
           2.0 * cos(alpha) * cos(beta) * cos(gamma))/sin(gamma) * c/)
      vas = cross(vb, vc)
      vbs = cross(vc, va)
      vcs = cross(va, vb)
      vas(1:3) = vas(1:3) / V
      vbs(1:3) = vbs(1:3) / V
      vcs(1:3) = vcs(1:3) / V
      cell_volume = V
      read(sfactors_unit, *) BFactors_len
      allocate(BFactors(BFactors_len))
      read(sfactors_unit, *) BFactors

      resolution = 50.0
      do i = 1, NRF
        read(sfactors_unit, *) hkl_tmp(i,:), f_obs_tmp(i), f_obs_sigma_tmp(i), r_free_flag
        d_star =  (square(hkl_tmp(i, 1) * norm2(vas)) + square(hkl_tmp(i, 2) * norm2(vbs)) + &
                   square(hkl_tmp(i, 3) * norm2(vcs)) + &
                   2 * hkl_tmp(i, 2) * hkl_tmp(i, 3) * dot_product(vbs, vcs) + &
                   2 * hkl_tmp(i, 1) * hkl_tmp(i, 3) * dot_product(vas, vcs) + &
                   2 * hkl_tmp(i, 1) * hkl_tmp(i, 2) * dot_product(vbs, vas))
        d = sqrt(1.0 / d_star)
        if (d < resolution) then
          resolution = d
!         write(mdout, *) 'highest resolution', d, i
        end if
        if (r_free_flag == 0) then
          NRF_work = NRF_work + 1
          k = NRF_work
          r_work_factor_denominator = r_work_factor_denominator + f_obs_tmp(i)
          F_obs_sq_sum_work = F_obs_sq_sum_work + f_obs_tmp(i) * f_obs_tmp(i) !* &
          ! (1.0/f_obs_sigma_tmp(i))**2
        else
          k = NRF - r_free_counter
          r_free_counter = r_free_counter + 1
          r_free_factor_denominator = r_free_factor_denominator + f_obs_tmp(i)
          F_obs_sq_sum_free = F_obs_sq_sum_free + f_obs_tmp(i) * f_obs_tmp(i) !* &
          ! (1.0/f_obs_sigma_tmp(i))**2
        endif
        r_free_flag_array(i) = r_free_flag
        d_star_sq(i) = d_star
        f_obs_weight(k) = 1.0 ! (1.0/f_obs_sigma_tmp(i))**2
      enddo
      pseudo_energy_weight = pseudo_energy_weight / sqrt(maxval(f_obs_weight))
      pseudo_energy_weight_incr = pseudo_energy_weight_incr / sqrt(maxval(f_obs_weight))
      scale_xray = -2.0
      write(mdout, *) 'scaling coefficient = ', scale_xray
      NRF_free = NRF - NRF_work
      if (NRF_free /= r_free_counter) then
        write(mdout, *) 'STOP!!!!!!'
      endif
      write(mdout, *) 'NRFs', NRF_work, NRF_free, NRF

      ! Sort reflections for binning
      allocate(d_star_sq_sorted(NRF_free))
      r_free_counter = 0
      do i = 1, NRF
        if (r_free_flag_array(i) == 1) then
          r_free_counter = r_free_counter + 1
          d_star_sq_sorted(r_free_counter) = d_star_sq(i)
        end if
      end do
      call bubble_sort(d_star_sq_sorted)
      reflections_per_bin = 140.0
      if (reflections_per_bin > NRF_free) then
        reflections_per_bin = 1.0 * NRF_free
      end if
      n_bins = max(1, int(NRF_free / reflections_per_bin + 0.5))
      reflections_per_bin = 1.0 * NRF_free / n_bins
      write(mdout, *) "adjusted reflections per bin", reflections_per_bin
      allocate(bin_limits(n_bins + 1))
      bin_limits(1) = max(zero, d_star_sq_sorted(1) * (1 - d_tolerance))
      do i = 2, n_bins
        d_i = int((i-1) * reflections_per_bin + 0.5) + 1
        bin_limits(i) = d_star_sq_sorted(d_i) * (1 - d_tolerance)
        if (d_i == NRF_free) then
          write(mdout, *) "check your reflections file"
          stop
        end if
      end do
      bin_limits(n_bins + 1) = d_star_sq_sorted(NRF_free) * (1 + d_tolerance)
      write(mdout, *) "resolution bins"
      do i = 1, n_bins + 1
        write(mdout, *) 1.0/sqrt(bin_limits(i))
      end do
      reflection_bin = n_bins
      allocate(A_in_zones(n_bins))
      allocate(B_in_zones(n_bins))
      allocate(q_in_zones(n_bins))
      allocate(C_in_zones(n_bins))
      allocate(t_optimal(n_bins))
      allocate(alpha_beta_OmegaI(n_bins))
      allocate(alpha_beta_wi(n_bins))
      allocate(alpha_in_zones(n_bins))
      allocate(beta_in_zones(n_bins))
      allocate(alpha_beta_bj(NRF))
      A_in_zones = 0.0
      B_in_zones = 0.0
      q_in_zones = 0.0
      C_in_zones = 0.0
      t_optimal = 0.0
      alpha_beta_OmegaI = 0.0
      alpha_beta_wi = 0.0
      alpha_in_zones = 0.0
      beta_in_zones = 1.0
      alpha_beta_bj = 0.0

      allocate(bins_work_population(n_bins))
      allocate(bins_free_population(n_bins))
      bins_work_population = 0
      bins_free_population = 0
      do i = 1, NRF
         
        ! Reverse order is to compare with cctbx output
        do j = n_bins, 2, -1
          if (d_star_sq(i) < bin_limits(j)) then
            reflection_bin(i) = j - 1
          end if
        end do
        if (r_free_flag_array(i) == 0) then
          bins_work_population(reflection_bin(i)) = bins_work_population(reflection_bin(i)) + 1
        else
          bins_free_population(reflection_bin(i)) = bins_free_population(reflection_bin(i)) + 1
        end if
      end do
      allocate(bins_work_reflections(n_bins, maxval(bins_work_population)))
      allocate(bins_free_reflections(n_bins, maxval(bins_free_population)))
      allocate(counter_w(n_bins))
      allocate(counter_f(n_bins))
      bins_work_reflections = 0
      bins_free_reflections = 0
      counter_w = 1
      counter_f = 1
      do i = 1, NRF
        if (r_free_flag_array(i) == 0) then
          if (counter_w(reflection_bin(i)) > bins_work_population(reflection_bin(i))) then
            write(mdout, *) reflection_bin(i), "something weird during resorting work happened"
          end if
          bins_work_reflections(reflection_bin(i), counter_w(reflection_bin(i))) = i
          counter_w(reflection_bin(i)) = counter_w(reflection_bin(i)) + 1
        else
          if (counter_f(reflection_bin(i)) > bins_free_population(reflection_bin(i))) then
            write(mdout, *) reflection_bin(i), "something weird during resorting free happened"
          end if
          bins_free_reflections(reflection_bin(i), counter_f(reflection_bin(i))) = i
          counter_f(reflection_bin(i)) = counter_f(reflection_bin(i)) + 1
        end if
      end do
      counter_sort = 0
      do i = 1, n_bins
        do j = 1, bins_work_population(i)
          counter_sort = counter_sort + 1
          index_sort = bins_work_reflections(i, j)
          f_obs_sigmas(counter_sort) = f_obs_sigma_tmp(index_sort)
          f_obs(counter_sort) = f_obs_tmp(index_sort)
          hkl(counter_sort,:) = hkl_tmp(index_sort,:)
        end do
      end do
      allocate(bins_free_start_indices(n_bins))
      do i = 1, n_bins
        bins_free_start_indices(i) = counter_sort + 1
        write(mdout, *) i, bins_free_start_indices(i), bins_free_population(i)
        do j = 1, bins_free_population(i)
          counter_sort = counter_sort + 1
          index_sort = bins_free_reflections(i, j)
          f_obs_sigmas(counter_sort) = f_obs_sigma_tmp(index_sort)
          f_obs(counter_sort) = f_obs_tmp(index_sort)
          hkl(counter_sort,:) = hkl_tmp(index_sort,:)
        end do
      end do
      if (counter_sort .ne. NRF) then

        ! Might be an error in the code?
        write(mdout, *) "Binning went wrong. Check the reflections file"
        stop
      end if

      do i = 1, n_bins
        do j = bins_free_start_indices(i) + bins_free_population(i) - 1, &
               bins_free_start_indices(i), -1
          fo_fo = f_obs(j) * f_obs(j) ! / 2
          B_in_zones(i) = B_in_zones(i) + 2 * fo_fo ! / coef_norm
          q_in_zones(i) = q_in_zones(i) + 2 * fo_fo * fo_fo ! / coef_norm
        end do
        B_in_zones(i) = B_in_zones(i) / (2 * bins_free_population(i))
        q_in_zones(i) = q_in_zones(i) / (2 * bins_free_population(i))
      end do

      ! Re-sort due to binning ended
      temp_grid = resolution / 4.0
#if 0
      na = adjust_gridding((int(a / temp_grid)/2)*2+1, 5)
      nb = adjust_gridding((int(b / temp_grid)/2)*2+1, 5)
      nc = adjust_gridding((int(c / temp_grid)/2)*2+1, 5)
#endif

      ! FIXME: need to make nfft1,nfft2,nfft3 visible here:
      ! na = nfft1
      ! nb = nfft2
      ! nc = nfft3

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
      write(mdout, *) 'resolution', resolution
      write(mdout, *) 'mask_cell_params', mask_cell_params
      write(mdout, *) 'mask_grid_steps', mask_grid_steps
      write(mdout, *) 'mask_grid_size', mask_grid_size
      allocate(mask_bs_grid(mask_grid_size(4)))
      allocate(mask_bs_grid_tmp(mask_grid_size(4)))
      allocate(mask_bs_grid_t_c(mask_grid_size(1) * mask_grid_size(2) * &
                                (mask_grid_size(3)/2 + 1)))
      allocate(mask_bs_grid_t_c_tmp(mask_grid_size(1) * mask_grid_size(2) * &
                                    (mask_grid_size(3) /2 + 1)))
      allocate(b_vector_base(NRF_work))
      do i = 1, NRF
        s(:, i) = hkl(i, 1) * vas(:) + hkl(i, 2) * vbs(:) + hkl(i, 3) * vcs(:)
        s_squared(i) = -0.25 * (s(1, i) ** 2 + s(2, i) ** 2 + s(3, i) ** 2)
        s(:, i) = 2 * pi * s(:, i)
        do j = 1, 5
          scat_factors_precalc(j, i) = get_scat_factor(j, i)
        end do
        hkl_indexing_bs_mask(i) = h_as_ih(hkl(i, 1), hkl(i, 2), hkl(i, 3), na, nb, nc)
        if (hkl_indexing_bs_mask(i) == -1) then
          stop 'Miller indices indexing failed'
        end if
        k_mask(i) = k_sol * exp(b_sol * s_squared(i))
        h_sq(i) = hkl(i, 1) * hkl(i, 1)
        k_sq(i) = hkl(i, 2) * hkl(i, 2)
        l_sq(i) = hkl(i, 3) * hkl(i, 3)
        hk(i) = hkl(i, 1) * hkl(i, 2)
        hl(i) = hkl(i, 1) * hkl(i, 3)
        kl(i) = hkl(i, 2) * hkl(i, 3)
      end do
      NRF_work_sq = NRF_work * NRF_work
      Ucryst(1, 1) = 1.0 / NRF_work
      Ucryst(1, 2) = sum(1.0 * h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(1, 3) = sum(1.0 * k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(1, 4) = sum(1.0 * l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(1, 5) = sum(1.0 *   hk(1:NRF_work) / NRF_work_sq)
      Ucryst(1, 6) = sum(1.0 *   hl(1:NRF_work) / NRF_work_sq)
      Ucryst(1, 7) = sum(1.0 *   kl(1:NRF_work) / NRF_work_sq)

      ! In case if one needs only anisotropic scaling
!     Ucryst(1, 1) = 1.0
!     Ucryst(1, 2) = 0.0
!     Ucryst(1, 3) = 0.0
!     Ucryst(1, 4) = 0.0
!     Ucryst(1, 5) = 0.0
!     Ucryst(1, 6) = 0.0
!     Ucryst(1, 7) = 0.0

      Ucryst(2, 1) = Ucryst(1, 2)
      Ucryst(3, 1) = Ucryst(1, 3)
      Ucryst(4, 1) = Ucryst(1, 4)
      Ucryst(5, 1) = Ucryst(1, 5)
      Ucryst(6, 1) = Ucryst(1, 6)
      Ucryst(7, 1) = Ucryst(1, 7)

      Ucryst(2, 2) = sum(1.0 * h_sq(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(2, 3) = sum(1.0 * k_sq(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(2, 4) = sum(1.0 * l_sq(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(2, 5) = sum(1.0 *   hk(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(2, 6) = sum(1.0 *   hl(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(2, 7) = sum(1.0 *   kl(1:NRF_work)*h_sq(1:NRF_work) / NRF_work_sq)

      Ucryst(3, 2) = sum(1.0 * h_sq(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(3, 3) = sum(1.0 * k_sq(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(3, 4) = sum(1.0 * l_sq(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(3, 5) = sum(1.0 *   hk(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(3, 6) = sum(1.0 *   hl(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(3, 7) = sum(1.0 *   kl(1:NRF_work)*k_sq(1:NRF_work) / NRF_work_sq)

      Ucryst(4, 2) = sum(1.0 * h_sq(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(4, 3) = sum(1.0 * k_sq(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(4, 4) = sum(1.0 * l_sq(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(4, 5) = sum(1.0 *   hk(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(4, 6) = sum(1.0 *   hl(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)
      Ucryst(4, 7) = sum(1.0 *   kl(1:NRF_work)*l_sq(1:NRF_work) / NRF_work_sq)

      Ucryst(5, 2) = sum(1.0 * h_sq(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)
      Ucryst(5, 3) = sum(1.0 * k_sq(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)
      Ucryst(5, 4) = sum(1.0 * l_sq(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)
      Ucryst(5, 5) = sum(1.0 *   hk(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)
      Ucryst(5, 6) = sum(1.0 *   hl(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)
      Ucryst(5, 7) = sum(1.0 *   kl(1:NRF_work)*hk(1:NRF_work) / NRF_work_sq)

      Ucryst(6, 2) = sum(1.0 * h_sq(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)
      Ucryst(6, 3) = sum(1.0 * k_sq(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)
      Ucryst(6, 4) = sum(1.0 * l_sq(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)
      Ucryst(6, 5) = sum(1.0 *   hk(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)
      Ucryst(6, 6) = sum(1.0 *   hl(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)
      Ucryst(6, 7) = sum(1.0 *   kl(1:NRF_work)*hl(1:NRF_work) / NRF_work_sq)

      Ucryst(7, 2) = sum(1.0 * h_sq(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)
      Ucryst(7, 3) = sum(1.0 * k_sq(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)
      Ucryst(7, 4) = sum(1.0 * l_sq(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)
      Ucryst(7, 5) = sum(1.0 *   hk(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)
      Ucryst(7, 6) = sum(1.0 *   hl(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)
      Ucryst(7, 7) = sum(1.0 *   kl(1:NRF_work)*kl(1:NRF_work) / NRF_work_sq)

      call inverse(Ucryst, MUcryst_inv, 7)
      call calc_grid_neighbors()
    endif

    close(unit = sfactors_unit)

    return
  end subroutine open_sf_file

  !--------------------------------------------------------------------------------------------
  ! A_B_C_D_omega:  start alpha and beta estimation according to
  !                 https://doi.org/10.1107/S010876739500688X
  !
  ! Arguments:
  !   zone:       number of resolution bin
  !   f_calc_:    complex array of computed structure factors
  !--------------------------------------------------------------------------------------------
  subroutine A_B_C_D_omega(zone, f_calc_)

    integer :: zone, i, index_start, index_end
    complex(8) :: f_calc_(NRF)
    double precision :: coef_norm, fm_fm, fo_fm, D, p, q

    coef_norm = 2 * bins_free_population(zone)
    index_start = bins_free_start_indices(zone)
    index_end = bins_free_start_indices(zone) + bins_free_population(zone) - 1
    A_in_zones(zone) = 0.0
    C_in_zones(zone) = 0.0
    D = 0.0
    p = 0.0
    do i = index_end, index_start, -1  ! reverse order for compliance with cctbx
      fo_fm = f_obs(i) * abs(f_calc_(i))
      fm_fm = abs(f_calc_(i)) * abs(f_calc_(i))
      A_in_zones(zone) = A_in_zones(zone) + 2 * fm_fm ! / coef_norm
      C_in_zones(zone) = C_in_zones(zone) + 2 * fo_fm ! / coef_norm
      D = D + 2 * fo_fm * fo_fm ! / coef_norm
      p = p + 2 * fm_fm * fm_fm ! / coef_norm
      alpha_beta_bj(i) = fo_fm
    end do
    A_in_zones(zone) = A_in_zones(zone) / coef_norm
    C_in_zones(zone) = C_in_zones(zone) / coef_norm
    D = D / coef_norm
    p = p / coef_norm
    alpha_beta_OmegaI(zone) = 0.0

    if (((p > A_in_zones(zone) * A_in_zones(zone)) .and. &
         (q_in_zones(zone) > B_in_zones(zone) * B_in_zones(zone))) .or. &
        ((p < A_in_zones(zone) * A_in_zones(zone)) .and. &
         (q_in_zones(zone) < B_in_zones(zone) * B_in_zones(zone)))) then
      alpha_beta_OmegaI(zone) = (D - A_in_zones(zone) * B_in_zones(zone)) / &
                                (sqrt(abs(p - A_in_zones(zone) * A_in_zones(zone))) * &
                                 sqrt(abs(q_in_zones(zone) - &
                                          B_in_zones(zone) * B_in_zones(zone))))
    end if
    alpha_beta_wi(zone) = A_in_zones(zone) * B_in_zones(zone) - &
                          C_in_zones(zone) * C_in_zones(zone)

  end subroutine A_B_C_D_omega

  !--------------------------------------------------------------------------------------------
  ! blamm:  \Lambda(t) in the original paper
  !--------------------------------------------------------------------------------------------
  function blamm(t, zone) result(result)

    double precision :: t
    double precision :: result
    integer :: zone, i, index_start, index_end

    index_start = bins_free_start_indices(zone)
    index_end = bins_free_start_indices(zone) + bins_free_population(zone) - 1
    result = 0.0
    do i = index_start, index_end, 1
      result = result + alpha_beta_bj(i) * i1_over_i0(2.0 * t * alpha_beta_bj(i))
    end do
    result = result / bins_free_population(zone)
  end function blamm

  !--------------------------------------------------------------------------------------------
  ! t_optimal_function:  G(t) in the original paper
  !--------------------------------------------------------------------------------------------
  function t_optimal_function(t, zone) result(result)

    double precision :: t
    double precision :: result
    integer :: zone

    result = sqrt(1. + 4. * A_in_zones(zone) * B_in_zones(zone) * t * t) - 1. - &
             2. * t * blamm(t, zone)
  end function t_optimal_function

  !--------------------------------------------------------------------------------------------
  ! solvm:  binary search for a root of t_optimal_function(t) in a general case
  !--------------------------------------------------------------------------------------------
  subroutine solvm(zone)

    integer :: zone, n1, n2, nst1, nst2
    double precision :: eps, fgtopt, tl, fgtl, tr, fgtr

    nst1 = 10 * 5 + 1
    nst2 = 20 * 5 + 1
    eps = 0.0001/10.
    t_optimal(zone) = 0.
    fgtopt = 0.
    tl = C_in_zones(zone) / alpha_beta_wi(zone)
    fgtl = t_optimal_function(tl, zone)
    do n1 = 1, nst1
      tr = tl
      fgtr = fgtl
      tl = tr * 0.5
      fgtl = t_optimal_function(tl, zone)
      if (fgtl == 0.0) then
        t_optimal(zone) = tl
        fgtopt = 0.0
        return
      else if (fgtl < 0.0) then
        t_optimal(zone) = tr
        fgtopt = fgtr
        do n2 = 1, nst2
          if ((tr-tl) < eps * t_optimal(zone)) then
            return
          end if
          t_optimal(zone) = (tl * fgtr - tr * fgtl)/(fgtr-fgtl)
          fgtopt = t_optimal_function(t_optimal(zone), zone)
          if (fgtopt > 0.0) then
            tr = t_optimal(zone)
            fgtr = fgtopt
          else
            tl = t_optimal(zone)
            fgtl = fgtopt
          end if
        end do
      end if
    end do

  end subroutine solvm

  !--------------------------------------------------------------------------------------------
  ! estimate_t_optimal:  find root of G(t)
  !--------------------------------------------------------------------------------------------
  subroutine estimate_t_optimal(zone, f_calc_)

    integer :: zone
    complex(8) :: f_calc_(NRF)

    call A_B_C_D_omega(zone, f_calc_)
    if (alpha_beta_OmegaI(zone) <= 0.0) then
      t_optimal(zone) = 0.0
    else if (alpha_beta_wi(zone)/(A_in_zones(zone)*B_in_zones(zone)) <= 3.0E-07) then
      t_optimal(zone) = 1.0e+10
    else
      call solvm(zone)
    end if

  end subroutine estimate_t_optimal

  !--------------------------------------------------------------------------------------------
  ! smooth: smoothes t_optimal values over resuluion bins in maximum likelihood (ML) estimation
  !--------------------------------------------------------------------------------------------
  subroutine smooth(x)

    integer :: i
    double precision :: t_opt_1, t_opt_2, t_opt_3, x(n_bins)

    if (n_bins > 1) then
      t_opt_1 = x(1)
      t_opt_2 = x(2)
      do i = 2, n_bins - 1
        t_opt_3 = x(i + 1)
        x(i) = (t_opt_1+t_opt_2+t_opt_3)/3.0;
        t_opt_1 = t_opt_2;
        t_opt_2 = t_opt_3;
      end do
    end if

  end subroutine smooth

  !--------------------------------------------------------------------------------------------
  ! alpha_beta_in_zones: used in maximum likelihood (ML) estimation
  !--------------------------------------------------------------------------------------------
  subroutine alpha_beta_in_zones()

    integer :: i
    double precision :: Ai, Bi, tt, ww, hbeta

    do i = 1, n_bins
      Ai = A_in_zones(i)
      Bi = B_in_zones(i)
      if (t_optimal(i) == 0.0) then
        alpha_in_zones(i) = 0.0
        beta_in_zones(i)  = Bi
      else if (t_optimal(i) >= 1.0e+10) then
        alpha_in_zones(i) = sqrt(Ai/Bi)
        beta_in_zones(i) = 1.e-10
      else
        tt = 2.0 * t_optimal(i)
        ww = sqrt(1.0 + Ai * Bi * tt * tt)
        hbeta = Bi / (ww + 1.0)
        alpha_in_zones(i) = sqrt(hbeta * (ww-1.0) / Ai)
        beta_in_zones(i) = 2.0 * hbeta
      end if
    end do

  end subroutine alpha_beta_in_zones

  !--------------------------------------------------------------------------------------------
  ! alpha_beta_all:  expand alpha and beta ML parameters on all reflections
  !--------------------------------------------------------------------------------------------
  subroutine alpha_beta_all

    integer :: i

    do i = 1, NRF
      alpha_array(i) = alpha_in_zones(reflection_bin(i))
      beta_array(i)  = beta_in_zones(reflection_bin(i))
    end do

  end subroutine alpha_beta_all

  !--------------------------------------------------------------------------------------------
  ! estimate_alpha_beta:  estimate alpha and beta parameters in resolution bins,
  !                       as described in https://doi.org/10.1107/S010876739500688X,
  !                       Appendix A:
  !                       eq. 29     - estimate_t_optimal() estimates root of function G(t)
  !                       eq. 30, 31 - alpha_beta_in_zones() calculates alpha and beta from t
  !--------------------------------------------------------------------------------------------
  subroutine estimate_alpha_beta(f_calc_)

    integer :: i
    complex(8) :: f_calc_(NRF)

    do i = 1, n_bins
      call estimate_t_optimal(i, f_calc_)
    end do
    call smooth(t_optimal)
    call alpha_beta_in_zones()
    call smooth(alpha_in_zones)
    call smooth(beta_in_zones)
    call alpha_beta_all()

  end subroutine estimate_alpha_beta
  ! End alpha/beta estimation

  !--------------------------------------------------------------------------------------------
  ! init_bulk_solvent: intialize the bulk solvent mask to one, 'solvent present here.'
  !--------------------------------------------------------------------------------------------
  subroutine init_bulk_solvent()
    mask_bs_grid = 1
    mask_bs_grid_tmp = 1
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

    integer :: tid, n_atom
    double precision :: atomX, atomY, atomZ, dx, dy, dz, cutoff, cutoffsq, &
                        distsq, coas, cobs, cocs
    integer :: x_low, x_high, y_low, y_high, z_low, z_high, i, j, k, index
    double precision :: frac(3)
    double precision :: crd(3, n_atom)

    do tid = 1, NAT_for_mask

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

    double precision :: mask_bs_grid_3d(mask_grid_size(1) * mask_grid_size(2) * mask_grid_size(3))
    double precision :: mask_bs_grid_3d_fft(mask_grid_size(1) * mask_grid_size(2) * mask_grid_size(3))
    integer :: i, j, k

    ! CHECK
    write(mdout, '(a,3i8)') '| mask_grid_size = ', mask_grid_size(1), mask_grid_size(2), &
                            mask_grid_size(3)
    ! END CHECK

#if 0
    ! FIXME: need to work on the fft integration
    do i = 1, mask_grid_size(4)
      mask_bs_grid_3d(i) = mask_bs_grid(i)
    end do
    call slab_fft3drc_forward(mask_bs_grid_3d, mask_bs_grid_3d_fft, mask_grid_size(1), &
                              mask_grid_size(2), mask_grid_size(3))
    do i = 1, mask_grid_size(4)
      mask_bs_grid_t_c(i) = mask_bs_grid_3d_fft(i)
    end do

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
    
  !--------------------------------------------------------------------------------------------
  ! estimate_ml_parameters: currently this is only implemented on the CPU.  It is not all that
  !                         much time to pull down the array of structure factors and compute
  !                         these things on the CPU, but a GPU implementation of this
  !                         algorithm, looping over all reflections multiple times, would be
  !                         helpful.
  !
  ! Arguments:
  !   f_calc:    the initial computed structure factors (these will be modified with results
  !              from this routine)
  !   exay:      xray restraint energy
  !   nstep:     the number of steps, passed down all the way from runmd
  !--------------------------------------------------------------------------------------------
  subroutine estimate_ml_parameters(f_calc, f_obs, exray, nstep, NRF)
        
    !  xxx state_info_mod

    integer :: nstep, i, NRF
    double precision :: mean_delta
    complex(8)       :: f_calc(NRF)
    double precision :: f_obs(NRF)
    double precision :: k_scale(NRF)
    double precision :: r_work_factor_numerator, r_free_factor_numerator
    double precision :: exray
    
    if (mod(nstep, mask_update_frequency) == 0 .or. nstep == 0) then
      call estimate_alpha_beta(f_calc)
      delta_array = alpha_array / beta_array
      mean_delta = sum(delta_array) / NRF
    endif

    f_calc_abs = abs(f_calc)
    r_work_factor_numerator = sum(abs(f_calc_abs(1:NRF_work) - f_obs(1:NRF_work)))
    r_free_factor_numerator = sum(abs(f_calc_abs(NRF_work + 1:NRF) - f_obs(NRF_work + 1:NRF)))
    if (mod(nstep, ntpr) == 0 .or. nstep == 0) then
      xray_r_work = r_work_factor_numerator / r_work_factor_denominator
      xray_r_free = r_free_factor_numerator / r_free_factor_denominator
    end if
    pseudo_energy_weight = pseudo_energy_weight + pseudo_energy_weight_incr

    ! ML target function
    exray = exray + pseudo_energy_weight * &
                        sum(alpha_array(1:NRF_work) * delta_array(1:NRF_work) * &
                            f_calc_abs(1:NRF_work) * f_calc_abs(1:NRF_work) - &
                            ln_of_i0(2 * delta_array(1:NRF_work) * f_calc_abs(1:NRF_work) * &
                            f_obs(1:NRF_work)))

  end subroutine estimate_ml_parameters

#if 0
  !--------------------------------------------------------------------------------------------
  ! estimate_weight: currently not in use.
  !                  Handles the X-ray weight ramping if selected by the user
  !
  ! Arguments:
  !   n_atom:
  !   current_step:
  !--------------------------------------------------------------------------------------------
  subroutine estimate_weight(n_atom, current_step)

    ! xxx state_info_mod

    implicit none

    integer :: n_atom, current_step

    open (unit=sf_weight, file="final_weight.txt", status='old', action='read')
    read (sf_weight, *) starting_N_step, total_N_steps
    write(mdout, *) 'starting and final steps = ', starting_N_step, total_N_steps
    close (sf_weight)

    if (call_est == 1) then
      if (pseudo_energy_weight_incr > 0.0) then
        pseudo_energy_weight_incr = pseudo_energy_weight / total_N_steps
        pseudo_energy_weight = (1.d0 * starting_N_step) / total_N_steps + &
                               pseudo_energy_weight_incr
        write(mdout, *) 'starting weight and step = ', pseudo_energy_weight, &
                        pseudo_energy_weight_incr
      end if
    end if
    write(mdout, *) 'starting and final weights = ', scale_xray * pseudo_energy_weight, &
                    scale_xray * pseudo_energy_weight_incr * total_N_steps
    call_est = call_est + 1

  end subroutine
#endif

  !--------------------------------------------------------------------------------------------
  ! get_sf_force: CPU routine to encapsulate the subroutines for computing forces based on
  !               structure factor restraints.
  !
  ! Arguments:
  !   n_atom:   the number of atoms in the simulation (passed down from runmd)
  !   crd:      atomic coordinates for the whole system
  !   frc:      forces on all atoms in the system
  !   exray:    potential energies record, serves PME only
  !   nstep:    the step number (passed down from runmd)
  !--------------------------------------------------------------------------------------------
  subroutine get_sf_force(n_atom, crd, frc, exray, nstep)

    ! xxx state_info_mod

    implicit none

    double precision :: crd(3, n_atom)
    double precision :: frc(3, n_atom)

    double precision :: term, b(7), Uaniso(7)
    complex(8)       :: exp_of_rn_x_s_dot_product
    integer :: i, j, j_k, n_atom, nstep
    double precision :: exray

    call init_bulk_solvent()
    call grid_bulk_solvent(n_atom, crd)
    call shrink_bulk_solvent()
    call fft_bs_mask()
    
    do i = 1, NRF
      f_calc(i) = 0
      do j_k = 1, NAT_for_mask
        f_calc(i) = f_calc(i) + exp(BFactors(j_k) * s_squared(i)) * &
                    scat_factors_precalc(atom_types(j_k), i) * &
                    exp(imaginary_i * dot_product(s(:, i), crd(:, j_k)))
      end do
      f_mask(i) = conjg(mask_bs_grid_t_c(hkl_indexing_bs_mask(i) + 1)) * &
                        mask_cell_params(16) / mask_grid_size(4)
    end do
    f_calc = f_calc + k_mask * f_mask
    f_calc_abs = abs(f_calc)

    b_vector_base = log(f_obs(1:NRF_work) / f_calc_abs(1:NRF_work)) / NRF_work_sq
    b(1) = sum(b_vector_base)
    b(2) = sum(b_vector_base * h_sq(1:NRF_work))
    b(3) = sum(b_vector_base * k_sq(1:NRF_work))
    b(4) = sum(b_vector_base * l_sq(1:NRF_work))
    b(5) = sum(b_vector_base * hk(1:NRF_work))
    b(6) = sum(b_vector_base * hl(1:NRF_work))
    b(7) = sum(b_vector_base * kl(1:NRF_work))

    Uaniso = matmul(MUcryst_inv, b)

    k_scale = exp(Uaniso(1) + Uaniso(2)*h_sq + Uaniso(3)*k_sq + Uaniso(4)*l_sq + &
                              Uaniso(5)*hk   + Uaniso(6)*hl   + Uaniso(7)*kl)
    f_calc = f_calc * k_scale

    call estimate_ml_parameters(f_calc, f_obs, exray, nstep, NRF)

    do j_k = 1, NAT_for_mask
      frc_sf(:, j_k) = zero
      do i = 1, NRF_work
        exp_of_rn_x_s_dot_product = exp(BFactors(j_k) * s_squared(i)) * &
                                    scat_factors_precalc(atom_types(j_k), i) * &
                                    exp(imaginary_i * dot_product(s(:, i), crd(:, j_k)))

        ! CHECK (uncomment code below when this placeholder is no longer needed)
        !term = 0.0
        ! END CHECK
        
        term = k_scale(i) * (aimag(f_calc(i)) * real(exp_of_rn_x_s_dot_product) - &
                              real(f_calc(i)) * aimag(exp_of_rn_x_s_dot_product)) * &
               f_obs_weight(i) * &
               (delta_array(i) * (alpha_array(i) - &
                                  f_obs(i) * i1_over_i0(2 * delta_array(i) * &
                                                        f_calc_abs(i) * f_obs(i)) / &
                                  f_calc_abs(i)))
        frc_sf(:, j_k) = frc_sf(:, j_k) + s(:, i) * term
      end do
      frc(:, j_k) = frc(:, j_k) + scale_xray * pseudo_energy_weight * frc_sf(:, j_k)
    end do

    return

  end subroutine get_sf_force

end module ml_mod
