! Created by  on 12/29/19.

module k_scale_mod
    implicit none
    double precision, dimension(:), allocatable :: k_mask_bin_orig, k_iso, k_iso_test, k_aniso, k_aniso_test, &
            k_iso_exp, k_iso_exp_test, &
            s_squared_min_bin, s_squared_max_bin, s_squared_mean_bin
    double precision :: k_overall

    ! This is written out with nine lines of five entries, but really just one big line full
    ! of 45 entries.  The data is immediately reshaped into the requisite 9 col x 5 row array.
    double precision, dimension(9, 5), save :: &
            scat_factors = reshape((/ &
                    0.493002, 0.322912, 0.140191, 0.040810, 10.510900, &
                            26.125700, 3.142360, 57.799700, 0.003038, 2.310000, &
                            1.020000, 1.588600, 0.865000, 20.843900, 10.207500, &
                            0.568700, 51.651200, 0.215600, 12.212600, 3.132200, &
                            2.012500, 1.166300, 0.005700, 9.893300, 28.997500, &
                            0.582600, -11.529000, 3.048500, 2.286800, 1.546300, &
                            0.867000, 13.277100, 5.701100, 0.323900, 32.908900, &
                            0.250800, 6.905300, 5.203400, 1.437900, 1.586300, &
                            1.467900, 22.215100, 0.253600, 56.172000, 0.866900 /), &
                    shape(scat_factors)) ! H, C, N, O, S scattering factors, it1992 scattering table

    ! Space group, currently read in from the structure factors file.
    ! Not really used, since we use the expanded-to-P1 structure factors.
    character(len = 16) :: space_group

    ! Arrays to hold calculated structure factors and the bulk solvent mask
    complex(8), dimension(:), allocatable :: f_calc, f_mask, f_calc_tmp, f_mask_tmp, mask_bs_grid_t_c, &
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
    double precision, dimension(6) :: uc_dimensions

    ! Bulk solvent mask parameters
    double precision, dimension(16) :: mask_cell_params
    double precision, dimension(3) :: mask_grid_steps

    ! 7 x 7 transformation matrices to anisotropically scale structure factors from calculated
    ! to experimental
    double precision, dimension(7, 7) :: Ucryst, MUcryst_inv

    ! Arrays used in the computation of structure factors
    double precision, dimension(:), allocatable :: f_obs, f_obs_weight, f_obs_sigmas, &
            f_calc_abs, s_squared, s_squared_for_scaling, BFactors, k_mask, &
            k_scale, alpha_array, &
            beta_array, delta_array, frc_adp, &
            frc_adp_a_priori, mask_cutoffs, &
            b_vector_base
    double precision, dimension(:, :), allocatable :: s, scat_factors_precalc

    ! Arrays used in the estimation of maximum likelihood parameters,
    ! size is equal to the number of resolution bins/zones
    double precision, dimension(:), allocatable :: A_in_zones, B_in_zones, &
            C_in_zones, q_in_zones, &
            alpha_beta_bj, alpha_beta_OmegaI, alpha_beta_wi, &
            t_optimal, alpha_in_zones, beta_in_zones

    ! Structure factor force computation scratch arrays
    double precision, dimension(:, :), allocatable :: frc_sf

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
            mask_update_frequency, grid_neighbors_size, n_bins, &
            hi_res_shell_n

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
            bins_free_start_indices, reflection_bin, &
            scale_k1_indices

    ! hkl:                      N columns by 3 rows array containing HKL indices of reflections
    ! bins_work_reflections:    Used to re-sort work reflections into the original order based on
    !                           resolution bin and relative position in it
    ! bins_free_reflections:    Used to re-sort free reflections into the original order based on
    !                           resolution bin and relative position in it
    integer, dimension(:, :), allocatable :: hkl, bins_work_reflections, bins_free_reflections

    ! Size of the mask grid (number of grid points)
    integer, dimension(4) :: mask_grid_size

    ! FFTW encoding
    integer*8 plan_forward(1)

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

        double precision, dimension(7) :: p = (/  1.00000000, 3.51562290, 3.08994240, &
                1.20672920, 0.26597320, 0.03607680, &
                0.00458130 /), &
                pp = (/  0.50000000, 0.87890594, 0.51498869, &
                        0.15084934, 0.02658733, 0.00301532, &
                        0.00032411 /)
        double precision, dimension(9) :: q = (/  0.39894228, 0.01328592, 0.00225319, &
                -0.00157565, 0.00916281, -0.02057706, &
                0.02635537, -0.01647633, 0.00392377 /), &
                qq = (/  0.39894228, -0.03988024, -0.00362018, &
                        0.00163801, -0.01031555, 0.02282967, &
                        -0.02895312, 0.01787654, -0.00420059 /)
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
        result = be1 / be0;
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
            if (abs_x / 3.75 < 1.0) then
                y = x(i) / 3.75
                y = y * y
                bessel_lni0(i) = 1.0 + y * (3.5156229 + &
                        y * (3.0899424 + &
                                y * (1.2067492 + &
                                        y * (0.2659732 + &
                                                y * (0.0360768 + y * 0.0045813)))))
                bessel_lni0(i) = log(bessel_lni0(i));
            else
                y = 3.75 / abs_x
                y = 0.39894228 + y * (0.01328592 + &
                        y * (0.00225319 + &
                                y * (-0.00157565 + &
                                        y * (0.00916281 + &
                                                y * (-0.02057706 + &
                                                        y * (0.02635537 + &
                                                                y * (-0.01647633 + y * 0.00392377)))))))
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
    ! file_line_count:  calculate the number of structure factors in the input file, stored in
    !                   line_num
    !--------------------------------------------------------------------------------------------
    subroutine file_line_count(input_unit, file_name, line_num)

        implicit none

        character(len = *) :: file_name
        character(len = 1023) :: line
        integer(kind = 4) :: ierror, input_status, input_unit, line_num

        write(*, *) "Opening structure factors file"

        open(unit = input_unit, file = file_name, status = 'old', iostat = input_status)
        if (input_status /= 0) then
            write (*, '(a)') 'Fatal error! Cannot read structure factors file.'
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
    subroutine inverse(a, c, n)

        implicit none
        integer n
        double precision a(n, n), c(n, n)
        double precision L(n, n), U(n, n), b(n), d(n), x(n)
        double precision coeff
        integer i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 allows such operations on matrices
        L = 0.0d0
        U = 0.0d0
        b = 0.0d0

        ! step 1: forward elimination
        do k = 1, n - 1
            do i = k + 1, n
                coeff = a(i, k) / a(k, k)
                L(i, k) = coeff
                do j = k + 1, n
                    a(i, j) = a(i, j) - coeff * a(k, j)
                end do
            end do
        end do

        ! Step 2: prepare L and U matrices
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i = 1, n
            L(i, i) = 1.0d0
        end do

        ! U matrix is the upper triangular part of A
        do j = 1, n
            do i = 1, j
                U(i, j) = a(i, j)
            end do
        end do

        ! Step 3: compute columns of the inverse matrix C
        do k = 1, n
            b(k) = 1.0d0
            d(1) = b(1)

            ! Step 3a: Solve Ld=b using the forward substitution
            do i = 2, n
                d(i) = b(i)
                do j = 1, i - 1
                    d(i) = d(i) - L(i, j) * d(j)
                end do
            end do

            ! Step 3b: Solve Ux=d using the back substitution
            x(n) = d(n) / U(n, n)
            do i = n - 1, 1, -1
                x(i) = d(i)
                do j = n, i + 1, -1
                    x(i) = x(i) - U(i, j) * x(j)
                end do
                x(i) = x(i) / u(i, i)
            end do

            ! Step 3c: fill the solutions x(n) into column k of C
            do i = 1, n
                c(i, k) = x(i)
            end do
            b(k) = 0.0d0
        end do
    end subroutine inverse

    !--------------------------------------------------------------------------------------------
    ! bubble_sort:  go bubble sort!  Sort a list of real numbers a.
    !--------------------------------------------------------------------------------------------
    subroutine bubble_sort(a)
        double precision, dimension(:) :: a
        double precision :: temp
        integer :: i, j
        logical :: swapped

        do j = size(a) - 1, 1, -1
            swapped = .FALSE.
            do i = 1, j
                if (a(i) > a(i + 1)) then
                    temp = a(i)
                    a(i) = a(i + 1)
                    a(i + 1) = temp
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
                temp_grid, resolution, d, d_star, reflections_per_bin, fo_fo, low_res
        double precision, dimension(3) :: va, vb, vc, vas, vbs, vcs
        double precision, dimension(:), allocatable :: f_obs_tmp, f_obs_sigma_tmp, d_star_sq, &
                d_star_sq_sorted, bin_limits
        integer, dimension(:), allocatable :: r_free_flag_array, counter_w, counter_f, reflection_bin_tmp
        integer (kind = 4) :: file_status, file_status1, file_status2
        integer :: d_i, i, j, k, na, nb, nc
        integer :: r_free_counter, r_free_flag, counter_sort, index_sort
        integer :: sfactors_unit
        integer :: index_start, index_end, hi_res_shell_n_counter
        integer, dimension(:, :), allocatable :: hkl_tmp

        call file_line_count(sfactors_unit, sfactors_name, NRF)
        write(*, '(A,I6,A)') 'Counted a total of', NRF, 'file lines'

        open(unit = sfactors_unit, file = sfactors_name, status = 'old', iostat = file_status)
        open(unit = 77, file = 'f_calc', status = 'old', iostat = file_status1)
        open(unit = 78, file = 'f_mask', status = 'old', iostat = file_status2)

        if (file_status /= 0) then
            write (*, '(a)') 'Fatal error with structure factors file.'
            sfactors_unit = -1
            stop
        else
            allocate(f_calc(NRF))
            allocate(f_calc_tmp(NRF))
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
            allocate(reflection_bin_tmp(NRF))
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
            call initialize_scales()
            allocate(delta_array(NRF))
            allocate(alpha_array(NRF))
            alpha_array = 1.0
            allocate(beta_array(NRF))
            beta_array = 1.0
            allocate(s(3, NRF))
            allocate(s_squared(NRF))
            allocate(s_squared_for_scaling(NRF))
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
            write(*, *) '| Original space group is ', space_group
            read(sfactors_unit, *) a, b, c, alpha, beta, gamma
            uc_dimensions = (/ a, b, c, alpha, beta, gamma /)
            read(sfactors_unit, *) pseudo_energy_weight, pseudo_energy_weight_incr
            write(*, '(a,f9.4,a,f9.4)') '| Initial weights = ', pseudo_energy_weight, ' ', &
                    pseudo_energy_weight_incr
            read(sfactors_unit, *) mask_update_frequency
            write(*, *) '| Mask update frequency = ', mask_update_frequency
            alpha = alpha * pi / 180
            beta = beta * pi / 180
            gamma = gamma * pi / 180
            V = a * b * c * &
                    sqrt(1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma) + &
                            2.0 * cos(alpha) * cos(beta) * cos(gamma))
            va = (/ a, zero, zero/)
            vb = (/cos(gamma) * b, sin(gamma) * b, zero/)
            vc = (/cos(beta) * c, (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma) * c, &
                    sqrt(1.0 - cos(alpha) * cos(alpha) - cos(beta) * cos(beta) - cos(gamma) * cos(gamma) + &
                            2.0 * cos(alpha) * cos(beta) * cos(gamma)) / sin(gamma) * c/)
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
            low_res = 0.0
            do i = 1, NRF
                read(sfactors_unit, *) hkl_tmp(i, :), f_obs_tmp(i), f_obs_sigma_tmp(i), r_free_flag
                read(77, *) f_calc_tmp(i)
                read(78, *) f_mask_tmp(i)
                d_star = (square(hkl_tmp(i, 1) * norm2(vas)) + square(hkl_tmp(i, 2) * norm2(vbs)) + &
                        square(hkl_tmp(i, 3) * norm2(vcs)) + &
                        2 * hkl_tmp(i, 2) * hkl_tmp(i, 3) * dot_product(vbs, vcs) + &
                        2 * hkl_tmp(i, 1) * hkl_tmp(i, 3) * dot_product(vas, vcs) + &
                        2 * hkl_tmp(i, 1) * hkl_tmp(i, 2) * dot_product(vbs, vas))
                d = sqrt(1.0 / d_star)
                if (d > low_res) then
                    low_res = d
                end if

                if (d < resolution) then
                    resolution = d
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
                f_obs_weight(k) = 1.0d0 ! (1.0/f_obs_sigma_tmp(i))**2
            enddo
            write(*, *) low_res, resolution
            !pseudo_energy_weight = pseudo_energy_weight / sqrt(maxval(f_obs_weight))
            !pseudo_energy_weight_incr = pseudo_energy_weight_incr / sqrt(maxval(f_obs_weight))
            scale_xray = -2.0d0
            NRF_free = NRF - NRF_work
            if (NRF_free /= r_free_counter) then
                write(*, *) 'STOP!!!!!!'
            endif
            write(*, *) 'NRFs', NRF_work, NRF_free, NRF

            ! Sort reflections for binning
            allocate(d_star_sq_sorted(NRF_free))
            r_free_counter = 0
            hi_res_shell_n = 0
            do i = 1, NRF
                if (r_free_flag_array(i) == 1) then
                    r_free_counter = r_free_counter + 1
                    d_star_sq_sorted(r_free_counter) = d_star_sq(i)
                else if (1./sqrt(d_star_sq(i)) < min(4.0d0, resolution + 1.0d0)) then
                    hi_res_shell_n = hi_res_shell_n + 1
                end if
            end do
            write(*, *) 'hi_res_shell_n', hi_res_shell_n
            allocate(scale_k1_indices(hi_res_shell_n))
            scale_k1_indices = 0
            call bubble_sort(d_star_sq_sorted)
            reflections_per_bin = 140.0d0
            if (reflections_per_bin > NRF_free) then
                reflections_per_bin = 1.0d0 * NRF_free
            end if
            n_bins = max(1, int(NRF_free / reflections_per_bin + 0.5))
            write(*, *) 'n_bins', n_bins
            reflections_per_bin = 1.0 * NRF_free / n_bins
            write(*, *) "adjusted reflections per bin", reflections_per_bin
            allocate(bin_limits(n_bins + 1))
            bin_limits(1) = max(zero, d_star_sq_sorted(1) * (1 - d_tolerance))
            bin_limits(1) = 1 / (low_res*low_res) * (1 - d_tolerance) ! to include leftover work_reflections
            do i = 2, n_bins
                d_i = int((i - 1) * reflections_per_bin + 0.5) + 1
                bin_limits(i) = d_star_sq_sorted(d_i) * (1 - d_tolerance)
                if (d_i == NRF_free) then
                    write(*, *) "check your reflections file"
                    stop
                end if
            end do
            bin_limits(n_bins + 1) = d_star_sq_sorted(NRF_free) * (1 + d_tolerance)
            bin_limits(n_bins + 1) = 1/(resolution*resolution) * (1 + d_tolerance) ! leftover free
            write(*, *) "resolution bins"
            do i = 1, n_bins + 1
                write(*, *) 1.0d0 / sqrt(bin_limits(i))
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
            allocate(s_squared_min_bin(n_bins))
            allocate(s_squared_max_bin(n_bins))
            allocate(s_squared_mean_bin(n_bins))
            allocate(alpha_beta_bj(NRF))
            A_in_zones = 0.0d0
            B_in_zones = 0.0d0
            q_in_zones = 0.0d0
            C_in_zones = 0.0d0
            t_optimal = 0.0d0
            alpha_beta_OmegaI = 0.0d0
            alpha_beta_wi = 0.0d0
            alpha_in_zones = 0.0d0
            beta_in_zones = 1.0d0
            alpha_beta_bj = 0.0d0

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
            !write(*, *) 'reflection bin', maxval(reflection_bin)
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
                        write(*, *) reflection_bin(i), "something weird during resorting work happened"
                    end if
                    bins_work_reflections(reflection_bin(i), counter_w(reflection_bin(i))) = i
                    counter_w(reflection_bin(i)) = counter_w(reflection_bin(i)) + 1
                else
                    if (counter_f(reflection_bin(i)) > bins_free_population(reflection_bin(i))) then
                        write(*, *) reflection_bin(i), "something weird during resorting free happened"
                    end if
                    bins_free_reflections(reflection_bin(i), counter_f(reflection_bin(i))) = i
                    counter_f(reflection_bin(i)) = counter_f(reflection_bin(i)) + 1
                end if
            end do
            counter_sort = 0
            hi_res_shell_n_counter = 0
            do i = 1, n_bins
                do j = 1, bins_work_population(i)
                    counter_sort = counter_sort + 1
                    index_sort = bins_work_reflections(i, j)
                    f_obs_sigmas(counter_sort) = f_obs_sigma_tmp(index_sort)
                    f_obs(counter_sort) = f_obs_tmp(index_sort)
                    f_calc(counter_sort) = f_calc_tmp(index_sort)
                    f_mask(counter_sort) = f_mask_tmp(index_sort)
                    if (1./sqrt(d_star_sq(index_sort)) < min(4.0d0, resolution + 1.0d0)) then
                        hi_res_shell_n_counter = hi_res_shell_n_counter + 1
                        scale_k1_indices(hi_res_shell_n_counter) = counter_sort
                    end if
                    hkl(counter_sort, :) = hkl_tmp(index_sort, :)
                    reflection_bin_tmp(counter_sort) = reflection_bin(index_sort)
                end do
            end do
            if (minval(scale_k1_indices) == 0) then
                write(*, *) 'weird hi res shell indexing'
            end if
            allocate(k_mask_bin_orig(n_bins))
            allocate(bins_free_start_indices(n_bins))
            do i = 1, n_bins
                bins_free_start_indices(i) = counter_sort + 1
                !write(*, *) i, bins_free_start_indices(i), bins_free_population(i)
                do j = 1, bins_free_population(i)
                    counter_sort = counter_sort + 1
                    index_sort = bins_free_reflections(i, j)
                    f_obs_sigmas(counter_sort) = f_obs_sigma_tmp(index_sort)
                    f_obs(counter_sort) = f_obs_tmp(index_sort)
                    f_calc(counter_sort) = f_calc_tmp(index_sort)
                    f_mask(counter_sort) = f_mask_tmp(index_sort)
                    hkl(counter_sort, :) = hkl_tmp(index_sort, :)
                    reflection_bin_tmp(counter_sort) = reflection_bin(index_sort)
                end do
            end do
            reflection_bin = reflection_bin_tmp
            if (counter_sort /= NRF) then

                ! Might be an error in the code?
                write(*, *) "Binning went wrong. Check the reflections file"
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
            temp_grid = resolution / 4.0d0
            na = adjust_gridding((int(a / temp_grid)/2)*2+1, 5)
            nb = adjust_gridding((int(b / temp_grid)/2)*2+1, 5)
            nc = adjust_gridding((int(c / temp_grid)/2)*2+1, 5)
            allocate(b_vector_base(NRF_work))
            do i = 1, NRF
                s(:, i) = hkl(i, 1) * vas(:) + hkl(i, 2) * vbs(:) + hkl(i, 3) * vcs(:)
                s_squared_for_scaling(i) = s(1, i) ** 2 + s(2, i) ** 2 + s(3, i) ** 2
                s_squared_for_scaling(i) = s_squared_for_scaling(i) /4d0
                s_squared(i) = -1.0d0 * s_squared_for_scaling(i)

                s(:, i) = 2 * pi * s(:, i)
                do j = 1, 5
                    scat_factors_precalc(j, i) = get_scat_factor(j, i)
                end do
                hkl_indexing_bs_mask(i) = h_as_ih(hkl(i, 1), hkl(i, 2), hkl(i, 3), na, nb, nc)
                if (hkl_indexing_bs_mask(i) == -1) then
                    stop 'Miller indices indexing failed'
                end if
                ! k_mask(i) = k_sol * exp(b_sol * s_squared(i))
                h_sq(i) = hkl(i, 1) * hkl(i, 1)
                k_sq(i) = hkl(i, 2) * hkl(i, 2)
                l_sq(i) = hkl(i, 3) * hkl(i, 3)
                hk(i) = hkl(i, 1) * hkl(i, 2)
                hl(i) = hkl(i, 1) * hkl(i, 3)
                kl(i) = hkl(i, 2) * hkl(i, 3)
            end do
            index_end = 0
            do i = 1, n_bins
                index_start = index_end + 1
                index_end = index_start + bins_work_population(i) - 1
                s_squared_min_bin(i) = minval(s_squared_for_scaling(index_start:index_end))
                s_squared_max_bin(i) = maxval(s_squared_for_scaling(index_start:index_end))
                s_squared_mean_bin(i) = sum(s_squared_for_scaling(index_start:index_end))/bins_work_population(i)
            end do
            do i = 1, n_bins
                index_start = bins_free_start_indices(i)
                index_end = bins_free_start_indices(i) + bins_free_population(i) - 1
                s_squared_min_bin(i) = min(s_squared_min_bin(i), minval(s_squared_for_scaling(index_start:index_end)))
                s_squared_max_bin(i) = max(s_squared_max_bin(i), maxval(s_squared_for_scaling(index_start:index_end)))
                s_squared_mean_bin(i) = (s_squared_mean_bin(i) + &
                                         sum(s_squared_for_scaling(index_start:index_end))/bins_free_population(i))/2
            end do
            NRF_work_sq = NRF_work * NRF_work
!            Ucryst(1, 1) = 1.0 / NRF_work
!            Ucryst(1, 2) = sum(1.0 * h_sq(1:NRF_work) / NRF_work_sq)
!            Ucryst(1, 3) = sum(1.0 * k_sq(1:NRF_work) / NRF_work_sq)
!            Ucryst(1, 4) = sum(1.0 * l_sq(1:NRF_work) / NRF_work_sq)
!            Ucryst(1, 5) = sum(1.0 * hk(1:NRF_work) / NRF_work_sq)
!            Ucryst(1, 6) = sum(1.0 * hl(1:NRF_work) / NRF_work_sq)
!            Ucryst(1, 7) = sum(1.0 * kl(1:NRF_work) / NRF_work_sq)

            ! In case if one needs only anisotropic scaling
            Ucryst(1, 1) = 1.0d0
            Ucryst(1, 2) = 0.0d0
            Ucryst(1, 3) = 0.0d0
            Ucryst(1, 4) = 0.0d0
            Ucryst(1, 5) = 0.0d0
            Ucryst(1, 6) = 0.0d0
            Ucryst(1, 7) = 0.0d0

            Ucryst(2, 1) = Ucryst(1, 2)
            Ucryst(3, 1) = Ucryst(1, 3)
            Ucryst(4, 1) = Ucryst(1, 4)
            Ucryst(5, 1) = Ucryst(1, 5)
            Ucryst(6, 1) = Ucryst(1, 6)
            Ucryst(7, 1) = Ucryst(1, 7)

            Ucryst(2, 2) = sum(1.0d0 * h_sq(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(2, 3) = sum(1.0d0 * k_sq(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(2, 4) = sum(1.0d0 * l_sq(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(2, 5) = sum(1.0d0 * hk(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(2, 6) = sum(1.0d0 * hl(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(2, 7) = sum(1.0d0 * kl(1:NRF_work) * h_sq(1:NRF_work) / NRF_work_sq)

            Ucryst(3, 2) = sum(1.0d0 * h_sq(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(3, 3) = sum(1.0d0 * k_sq(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(3, 4) = sum(1.0d0 * l_sq(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(3, 5) = sum(1.0d0 * hk(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(3, 6) = sum(1.0d0 * hl(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(3, 7) = sum(1.0d0 * kl(1:NRF_work) * k_sq(1:NRF_work) / NRF_work_sq)

            Ucryst(4, 2) = sum(1.0d0 * h_sq(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(4, 3) = sum(1.0d0 * k_sq(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(4, 4) = sum(1.0d0 * l_sq(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(4, 5) = sum(1.0d0 * hk(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(4, 6) = sum(1.0d0 * hl(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)
            Ucryst(4, 7) = sum(1.0d0 * kl(1:NRF_work) * l_sq(1:NRF_work) / NRF_work_sq)

            Ucryst(5, 2) = sum(1.0d0 * h_sq(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)
            Ucryst(5, 3) = sum(1.0d0 * k_sq(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)
            Ucryst(5, 4) = sum(1.0d0 * l_sq(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)
            Ucryst(5, 5) = sum(1.0d0 * hk(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)
            Ucryst(5, 6) = sum(1.0d0 * hl(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)
            Ucryst(5, 7) = sum(1.0d0 * kl(1:NRF_work) * hk(1:NRF_work) / NRF_work_sq)

            Ucryst(6, 2) = sum(1.0d0 * h_sq(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)
            Ucryst(6, 3) = sum(1.0d0 * k_sq(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)
            Ucryst(6, 4) = sum(1.0d0 * l_sq(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)
            Ucryst(6, 5) = sum(1.0d0 * hk(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)
            Ucryst(6, 6) = sum(1.0d0 * hl(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)
            Ucryst(6, 7) = sum(1.0d0 * kl(1:NRF_work) * hl(1:NRF_work) / NRF_work_sq)

            Ucryst(7, 2) = sum(1.0d0 * h_sq(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)
            Ucryst(7, 3) = sum(1.0d0 * k_sq(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)
            Ucryst(7, 4) = sum(1.0d0 * l_sq(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)
            Ucryst(7, 5) = sum(1.0d0 * hk(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)
            Ucryst(7, 6) = sum(1.0d0 * hl(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)
            Ucryst(7, 7) = sum(1.0d0 * kl(1:NRF_work) * kl(1:NRF_work) / NRF_work_sq)

            call inverse(Ucryst, MUcryst_inv, 7)
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
        A_in_zones(zone) = 0.0d0
        C_in_zones(zone) = 0.0d0
        D = 0.0d0
        p = 0.0d0
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
        alpha_beta_OmegaI(zone) = 0.0d0

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
        result = 0.0d0
        do i = index_start, index_end, 1
            result = result + alpha_beta_bj(i) * i1_over_i0(2.0d0 * t * alpha_beta_bj(i))
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
        eps = 0.0001d0 / 10.
        t_optimal(zone) = 0.d0
        fgtopt = 0.d0
        tl = C_in_zones(zone) / alpha_beta_wi(zone)
        fgtl = t_optimal_function(tl, zone)
        do n1 = 1, nst1
            tr = tl
            fgtr = fgtl
            tl = tr * 0.5d0
            fgtl = t_optimal_function(tl, zone)
            if (fgtl == 0.0) then
                t_optimal(zone) = tl
                fgtopt = 0.0d0
                return
            else if (fgtl < 0.0) then
                t_optimal(zone) = tr
                fgtopt = fgtr
                do n2 = 1, nst2
                    if ((tr - tl) < eps * t_optimal(zone)) then
                        return
                    end if
                    t_optimal(zone) = (tl * fgtr - tr * fgtl) / (fgtr - fgtl)
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
            t_optimal(zone) = 0.0d0
        else if (alpha_beta_wi(zone) / (A_in_zones(zone) * B_in_zones(zone)) <= 3.0E-07) then
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
                x(i) = (t_opt_1 + t_opt_2 + t_opt_3) / 3.0d0;
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
                alpha_in_zones(i) = 0.00d0
                beta_in_zones(i) = Bi
            else if (t_optimal(i) >= 1.0e+10) then
                alpha_in_zones(i) = sqrt(Ai / Bi)
                beta_in_zones(i) = 1.e-10
            else
                tt = 2.0d0 * t_optimal(i)
                ww = sqrt(1.0d0 + Ai * Bi * tt * tt)
                hbeta = Bi / (ww + 1.0d0)
                alpha_in_zones(i) = sqrt(hbeta * (ww - 1.0d0) / Ai)
                beta_in_zones(i) = 2.0d0 * hbeta
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
            beta_array(i) = beta_in_zones(reflection_bin(i))
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
    ! estimate_ml_parameters: currently this is only implemented on the CPU.  It is not all that
    !                         much time to pull down the array of structure factors and compute
    !                         these things on the CPU, but a GPU implementation of this
    !                         algorithm, looping over all reflections multiple times, would be
    !                         helpful.
    !
    ! Arguments:
    !   f_calc:    the initial computed structure factors (these will be modified with results
    !              from this routine)
    !   pot_ene:   potential energies (this is geared towards PME only)
    !   nstep:     the number of steps, passed down all the way from runmd
    !--------------------------------------------------------------------------------------------
    subroutine estimate_ml_parameters(f_calc)

        integer :: i, j, index_sort, counter_sort
        double precision :: mean_delta, restraint
        complex(8) :: f_calc(NRF)
        double precision :: r_work_factor_numerator, r_free_factor_numerator

        call estimate_alpha_beta(f_calc)
        
        !do i = 1, NRF
        !    write(*, *) i, hkl(i, :), alpha_array(i), beta_array(i)
        !end do
        
        delta_array = alpha_array / beta_array
        mean_delta = sum(delta_array) / NRF

        f_calc_abs = abs(f_calc)
        r_work_factor_numerator = sum(abs(f_calc_abs(1:NRF_work) - f_obs(1:NRF_work)))
        r_free_factor_numerator = sum(abs(f_calc_abs(NRF_work + 1:NRF) - f_obs(NRF_work + 1:NRF)))
            xray_r_work = r_work_factor_numerator / r_work_factor_denominator
            xray_r_free = r_free_factor_numerator / r_free_factor_denominator
        restraint = sum(f_obs(1:NRF_work)**2 / beta_array(1:NRF_work) - &
                        log(2*f_obs(1:NRF_work)/beta_array(1:NRF_work)) +  alpha_array(1:NRF_work) * delta_array(1:NRF_work) * &
                        f_calc_abs(1:NRF_work) * f_calc_abs(1:NRF_work) - &
                        ln_of_i0(2 * delta_array(1:NRF_work) * f_calc_abs(1:NRF_work) * &
                                f_obs(1:NRF_work)))
        write(*, *) 'xray restraint: ', restraint

    end subroutine estimate_ml_parameters

    !--------------------------------------------------------------------------------------------
    ! get_sf_force: CPU routine to encapsulate the subroutines for computing forces based on
    !               structure factor restraints.
    !
    ! Arguments:
    !   n_atom:   the number of atoms in the simulation (passed down from runmd)
    !   crd:      atomic coordinates for the whole system
    !   frc:      forces on all atoms in the system
    !   pot_ene:  potential energies record, serves PME only
    !   nstep:    the step number (passed down from runmd)
    !--------------------------------------------------------------------------------------------
    subroutine get_sf_force()
        implicit none

        double precision :: term, current_r_work, r, k_overall_
        integer :: i, j, j_k, cycle

        ! write(*, *) sum(k_iso), sum(k_aniso), sum(k_mask), sum(f_obs), sum(f_calc), sum(f_mask)
        r = 1.0d0
        cycle = 0
        current_r_work = r_factor_w(f_calc)
        write(*, '(a,f8.6)') 'starting r_factor = ', current_r_work
        do while (r - current_r_work > 1.e-4 .and. cycle < 20)
            write(*, '(a,i3)') 'CYCLE ', cycle
            r = current_r_work
            if (cycle == 0) then
                call fit_k_iso_exp(current_r_work, sqrt(s_squared_for_scaling), f_obs, abs(f_calc))
                call k_mask_grid_search(current_r_work)
                call fit_k_iso_exp(current_r_work, sqrt(s_squared_for_scaling), f_obs, &
                abs(k_iso * k_aniso*(f_calc + k_mask*f_mask)))
            else
                call bulk_solvent_scaling(current_r_work)
            end if
            call anisotropic_scaling(current_r_work)
            cycle = cycle + 1
        end do

        f_calc = k_iso * k_iso_exp * k_aniso * (f_calc + k_mask * f_mask)
        write(*, '(a, f10.8)') '  xray rwork: ', &
        r_factor_w_selection(scale_selection(f_calc, 1,NRF_work)*f_calc, 1, NRF_work)
        write(*, '(a, f10.8)') '  xray rfree: ', &
        r_factor_w_selection(scale_selection(f_calc, 1 + NRF_work, NRF)*f_calc, 1 + NRF_work, NRF)

        call estimate_ml_parameters(f_calc)

        return

    end subroutine get_sf_force

    subroutine anisotropic_scaling(r_start)
        double precision :: r_start, r, b(7), Uaniso(7), u_star(6)
        f_calc_tmp = k_iso * k_iso_exp * (f_calc + k_mask * f_mask)
        f_calc_abs = abs(f_calc_tmp)
        write(*, '(a,f8.6)') 'k_aniso, starting r_factor = ', r_start

        b_vector_base = log(f_obs(1:NRF_work) / f_calc_abs(1:NRF_work)) / NRF_work_sq
        b(1) = sum(b_vector_base)
        b(2) = sum(b_vector_base * h_sq(1:NRF_work))
        b(3) = sum(b_vector_base * k_sq(1:NRF_work))
        b(4) = sum(b_vector_base * l_sq(1:NRF_work))
        b(5) = sum(b_vector_base * hk(1:NRF_work))
        b(6) = sum(b_vector_base * hl(1:NRF_work))
        b(7) = sum(b_vector_base * kl(1:NRF_work))

        Uaniso = matmul(MUcryst_inv, b) ! in Phenix it is u_star  multiplied by (-2 *pi *pi), b(5:7) are doubled
        u_star = Uaniso(2:7) / (-2.0d0 * pi * pi)
        u_star(4:6) = u_star(4:6) / 2.0d0
        write(*, '(a,f14.12,a,f14.12,a,f14.12,a,f14.12,a,f14.12,a,f14.12)') "u_star(11,22,33,12,13,23):", &
        u_star(1), ' ', u_star(2), ' ', u_star(3), ' ', u_star(4), ' ', u_star(5), ' ', u_star(6)
        

        k_aniso_test = exp(Uaniso(2) * h_sq + Uaniso(3) * k_sq + Uaniso(4) * l_sq + &
                           Uaniso(5) * hk + Uaniso(6) * hl + Uaniso(7) * kl) !&
                           ! + Uaniso(1))
        r = r_factor_w(k_aniso_test * f_calc_tmp)
        if (r < r_start) then
            r_start = r
            k_aniso = k_aniso_test
        end if
        write(*, '(a,f8.6)') 'k_aniso, final r_factor = ', r_start
    end subroutine

    subroutine fit_k_iso_exp(r_start, x, y, z)
        double precision, dimension(2) :: a
        double precision :: r_start, p, q, r, s, d, v, den, u
        double precision, dimension(NRF) :: x, y, z
        integer :: i

        a = 0d0
        p = 0d0
        q = 0d0
        r = 0d0
        s = 0d0
        write(*, '(a,f8.6)') 'k_iso_exp, starting r_factor = ', r_start

        ! fit over all work reflections
        do i = 1, NRF_work
            if (z(i) /= 0) then
                if (y(i) > 0) then
                    d = log(y(i) / z(i))
                    v = x(i) * x(i)
                    p = p + d
                    q = q + v
                    r = r + v * v
                    s = s + v * d
                end if
            end if
        end do

        if (r /= 0) then
            den = NRF_work - q * q / r
            if (den /= 0) then
                u = (p - s * q / r) / den
                a(2) = (u * q - s) / r
                a(1) = exp(u)
            end if
        end if

        !write(*, *) a

        if (a(2) > -100) then
            k_iso_exp_test = a(1) * exp(s_squared * a(2))
            r = special_r_factor(k_iso_exp_test * k_iso * k_aniso * (f_calc + k_mask * f_mask))
            if (r < r_start) then
                k_iso_exp = k_iso_exp_test
                r_start = r_factor_w_scale(k_iso_exp_test * k_iso * k_aniso * (f_calc + k_mask * f_mask), 1.0d0)
            end if
        end if
        write(*, '(a,f9.6,a,f9.6,a,f8.6,a,f8.6)') 'k_iso_exp: k_a = ', a(1), ' k_b = ', a(2),&
            ' r_factor in hi_res_shell = ', r, ' final r_factor = ', r_start

    end subroutine fit_k_iso_exp

    subroutine moving_average(x, result2)
        double precision :: x(n_bins), x_(n_bins+2), result1(n_bins+2), result2(n_bins)
        logical :: selection(n_bins+2)
        integer :: cycle, i
        x_(2:n_bins+1) = x
        x_(1) = x(1)
        x_(n_bins+2) = x(n_bins)
        
        do cycle = 1, 5
            result1 = x_
            selection = .False.
            do i = 1, n_bins
                if (i > 1 .and. i < n_bins + 2) then
                    if((result1(i-1)<result1(i) .and. result1(i+1)<result1(i)) .or. &
                       (result1(i-1)>result1(i) .and. result1(i+1)>result1(i))) then
                        selection(i) = .True.
                    end if
                end if
            end do
            do i = 1, n_bins
                if(i> 1 .and. i<n_bins+2 .and. selection(i)) then
                    result1(i) = (x_(i-1)+x_(i)+x_(i+1))/3.
                end if
            end do
            x_ = result1
        end do
        result2 = result1(2:n_bins+1)
    end subroutine

    subroutine smooth_k_mask(x)
        double precision, dimension(n_bins) :: x, result_, result1
        double precision :: r, d
        integer :: i
        call moving_average(x, result1)
        result_ = 0
        do i = 1, n_bins
            r = result1(i)
            d = 1/sqrt(s_squared_max_bin(i))/2
            if (r==0 .and. d<3) then
                exit
            end if
            result_(i)=r
        end do
        x = result_
    end subroutine

    subroutine k_mask_grid_search(r_start)
        double precision :: r_start, r, k_mask_best, k_overall_best, r_best, k_mask_per_bin_test, k_overall_tmp, &
                num, denum, k_mask_test(NRF), k_overall_
        double precision, dimension(14) :: k_mask_trial_range
        double precision, dimension(n_bins) :: k_mask_bin, k_mask_bin_
        double precision, dimension(NRF_work) :: tmp_scale
        integer :: i, j, sampling, index_start, index_end
        sampling = 14
        do i = 1, sampling
            k_mask_trial_range(i) = (i - 1) * 0.05d0
        end do
        k_mask_bin = 0d0
        tmp_scale = k_aniso(1:NRF_work) * k_iso(1:NRF_work) * k_iso_exp(1:NRF_work)
        index_end = 0
        write(*, '(a,f8.6)') 'k_mask_grid_search, starting r_factor = ', r_start
        do i = 1, n_bins
            !write(*, *) sum(f_mask(index_start:index_end))
            index_start = index_end + 1
            index_end = bins_work_population(i) + index_start - 1
            k_mask_best = 0.0d0
            k_overall_best = 1.0d0
            r_best = r_factor_w_selection_and_modified_fobs(f_calc, tmp_scale, index_start,index_end)
            !write(*,*) i, r_best
            do j = 1, sampling
                k_mask_per_bin_test = k_mask_trial_range(j)
                f_calc_tmp(index_start:index_end) = f_calc(index_start:index_end) + &
                                                    k_mask_per_bin_test * f_mask(index_start:index_end)
                k_overall_ = scale_selection_with_modified_fobs(f_calc_tmp, tmp_scale, index_start, index_end)
                r = r_factor_w_selection_and_modified_fobs_and_scale(f_calc_tmp, tmp_scale, k_overall_, index_start, index_end)
                if (r < r_best) then
                    k_mask_best = k_mask_per_bin_test
                    k_overall_best = k_overall_
                    r_best = r
                end if
            end do
            !write(*,*) i, r_best, k_mask_best, k_overall_best
            k_mask_bin(i) = k_mask_best
        end do
        k_mask_bin_orig = k_mask_bin
        !write(*, '(a)') 'bulk_solvent_scaling, k_mask in resolution bins'
        !do j = 1, n_bins
        !    write(*, '(i20,a,f14.12)') j, ' ', k_mask_bin(j)
        !end do 
        call smooth_k_mask(k_mask_bin)
        write(*, '(a)') 'k_mask_grid_search, smoothed k_mask in resolution bins'
        do i = 1, n_bins
            write(*, '(i20,a,f8.4)') i, ' ', k_mask_bin(i)
        end do 
        call populate_k_mask_linear_interpolation(k_mask_bin, k_mask)
        call bin_k_isotropic_as_scale_k1(r_start, k_iso, k_mask)
        r_start = r_factor_w_scale(k_iso * k_aniso * k_iso_exp * (f_calc + k_mask * f_mask), 1.0d0)
        ! TO-DO: Gauss fit of k_iso in bins based on mid sqrt(s_squared)
        ! if(n_bins>2) then
        ! end if
        write(*, '(a,f8.6)') 'k_mask_grid_search, final r_factor = ', r_start

    end subroutine k_mask_grid_search

    subroutine bin_k_isotropic_as_scale_k1(r_start, k_iso_in, k_mask_in)
        double precision :: r_start, r, num, denum, sc, k_iso_in(NRF), k_mask_in(NRF)
        integer :: i, j, index_start, index_end
        k_iso_test = 1d0
        index_end = 0
        f_calc_tmp = k_iso_exp * k_aniso * (f_calc + k_mask_in * f_mask)
        do i = 1, n_bins
            index_start = index_end + 1
            index_end = bins_work_population(i) + index_start - 1
            k_iso_test(index_start:index_end) = scale_selection(f_calc_tmp, index_start, index_end)
        end do
        r = r_factor_w(f_calc_tmp * k_iso_test)
        if (r < r_start) then
            index_end = 0
            r_start = r
            do i = 1, n_bins
                index_start = index_end + 1
                index_end = bins_work_population(i) + index_start - 1
                j = index_end
                sc = k_iso_test(index_start)
                k_iso_in(index_start:index_end) = k_iso_test(index_start:index_end)
                index_start = bins_free_start_indices(i)
                index_end = bins_free_start_indices(i) + bins_free_population(i) - 1
                k_iso_in(index_start:index_end) = sc
                index_end = j
            end do
        end if

    end subroutine bin_k_isotropic_as_scale_k1

    function linear_interpolation(x1, x2, y1, y2) result (k)
        double precision :: x1, x2, y1, y2
        double precision, dimension(2) :: k
        ! write(*, *) x1, x2, y1, y2
        k(1) = 0
        if(x1 /=x2) then
            k(1)=(y2-y1)/(x2-x1)
        end if
        k(2) = y1-k(1)*x1
    end function linear_interpolation

    subroutine populate_k_mask_linear_interpolation(k_mask_bin, k_mask_in)
        double precision, dimension(n_bins) :: k_mask_bin
        double precision, dimension(NRF) :: k_mask_tmp, k_mask_in
        integer :: i, index_start, index_end, index_start_, index_end_
        double precision, dimension(2) :: k
        double precision :: r0, r1, r2
        !do i = 1, NRF
        !    k_mask(i) = k_mask_bin(reflection_bin(i))
        !end do
        index_end = 0
        do i = 1, n_bins
            if (i == n_bins) then
                ! double check that! cctbx code is weird. for some reason last argument was n_bins - 2
                k = linear_interpolation(s_squared_min_bin(n_bins), s_squared_max_bin(n_bins), &
                        k_mask_bin(n_bins), k_mask_bin(n_bins-1)) !! as in cctbx. is it a mistake for the last index?
            else
                k = linear_interpolation(s_squared_min_bin(i), s_squared_max_bin(i), k_mask_bin(i), k_mask_bin(i+1))
            end if
            index_start = index_end + 1
            index_end = bins_work_population(i) + index_start - 1
            k_mask_tmp(index_start:index_end) = max(0.0d0, k(1) * s_squared_for_scaling(index_start:index_end) + k(2))
            
            r0 = 1.0d0 ! r_factor_w_selection(k_iso * k_iso_exp * k_aniso * (f_calc + k_mask * f_mask), index_start, index_end)
            r1 = r_factor_w_selection(k_iso * k_iso_exp * k_aniso * (f_calc + k_mask_tmp * f_mask), index_start, index_end)
            r2 = r_factor_w_selection(k_iso * k_iso_exp * k_aniso * (f_calc + k_mask_bin(i) * f_mask), &
            index_start, index_end)
            ! write(*, *) k_mask_bin(i), k, r1, r2, r0
            if (r0 > min(r1, r2)) then
                if (r1 < r2) then
                    k_mask_in(index_start:index_end) = k_mask_tmp(index_start:index_end)
                    index_start_ = bins_free_start_indices(i)
                    index_end_ = bins_free_start_indices(i) + bins_free_population(i) - 1
                    k_mask_in(index_start_:index_end_) = k(1) * s_squared_for_scaling(index_start_:index_end_) + k(2)
                else
                    k_mask_in(index_start:index_end) = k_mask_bin(i)
                    index_start_ = bins_free_start_indices(i)
                    index_end_ = bins_free_start_indices(i) + bins_free_population(i) - 1
                    k_mask_in(index_start_:index_end_) = k_mask_bin(i)
                end if
            end if
        end do

    end subroutine populate_k_mask_linear_interpolation
    
    function solve_cubic_equation(a, b, c) result(solution)
        double precision :: a, b, c, solution(3), Disc, theta, S, T, p, q, arg, eps
        solution = 0d0
        eps = d_tolerance * 10
        p = ( 3*b - square(a) ) / 3
        q = ( 2*square(a)*a - 9*a*b + 27*c ) / 27
        Disc = square(p/3)*p/3 + square(q/2)
        if (abs(p)<eps .and. abs(q)<eps .and. abs(Disc)<eps) then
            solution(1:3) = - sign(abs(c)**(1./3.), c)
        else if (Disc < 0) then
            ! theta = acos(-q/2/(sqrt(abs((p/3)**3))))
            arg = (q/p)**2*27/(4*p)
            if (abs(1-abs(arg))<1.e-9) then
                arg = 1
            end if
            if (q>0) then
                theta = acos(-arg)
            else
                theta = acos(arg)
            end if
            solution(1) = 2 *(sqrt(abs(p/3)) * cos(theta/3)) - a/3
            solution(2) = - 2 *(sqrt(abs(p/3)) * cos((theta+2*pi)/3)) - a/3
            solution(3) = - 2 *(sqrt(abs(p/3)) * cos((theta-2*pi)/3)) - a/3
        else
            S = sign(abs(-q/2+sqrt(Disc))**(1./3.), -q/2+sqrt(Disc))
            T = sign(abs(-q/2-sqrt(Disc))**(1./3.), -q/2-sqrt(Disc))
            solution(1:3) = S + T - a/3
            if (Disc <= eps) then
                solution(2:3) = -(S+T)/2 - a/3
            end if
        end if
    end function
    
    subroutine bulk_solvent_scaling(r_start)
        double precision :: r_start, min_f_mask, a2, b2, c2, y2, a3, b3, c3, d3, y3, p, q, r, t, I, v, w, u, den, &
        k_m_test(3), k_best, r_best, shift, k_mask_per_bin_test, k_overall_, inc, k_best_test, k_mask_bin(n_bins), &
        tmp_scale(NRF_work), a, b, c, scale_k1, upper_limit, k_mask_test(NRF)
        integer :: l, j, index_start, index_end
        write(*, '(a,f8.6)') &
        'bulk_solvent_scaling (simulteneous k_iso and k_mask analytical determination), starting r_factor = ', r_start
        shift = 0.05d0
        index_end=0
        scale_k1 = estimate_scale_k1(k_iso_exp * k_aniso * (f_calc + k_mask * f_mask))
        !write(*, *) 'k1', scale_k1
        tmp_scale = k_aniso(1:NRF_work) * scale_k1 * k_iso_exp(1:NRF_work)
        
        do l = 1, n_bins
            index_start = index_end + 1
            index_end = bins_work_population(l) + index_start - 1
            min_f_mask = minval(tmp_scale(index_start:index_end) * abs(f_mask(index_start:index_end)))
            k_best = k_mask_bin_orig(l)
            f_calc_tmp(index_start:index_end) = tmp_scale(index_start:index_end) * (f_calc(index_start:index_end) + &
                                                    k_best * f_mask(index_start:index_end))
            k_overall_ = scale_selection(f_calc_tmp, index_start, index_end)
            r_best = r_factor_w_selection(f_calc_tmp, index_start, index_end)
            !write(*, *) 'previous round', k_best, r_best
            a2 = 0.0d0
            b2 = 0.0d0
            c2 = 0.0d0
            y2 = 0.0d0
            a3 = 0.0d0
            b3 = 0.0d0
            c3 = 0.0d0
            d3 = 0.0d0
            y3 = 0.0d0
            if (min_f_mask > 1.0e-9) then
                do j = index_start, index_end
                    p = real(tmp_scale(j) * f_calc(j))
                    r = aimag(tmp_scale(j) * f_calc(j))
                    q = real(tmp_scale(j) * f_mask(j))
                    t = aimag(tmp_scale(j) * f_mask(j))
                    I = f_obs(j) * f_obs(j)
                    v = p*q + r*t
                    w = q*q + t*t
                    u = p*p + r*r
                    a2 = a2 + u*I
                    b2 = b2 + 2.*v*I
                    c2 = c2 + w*I
                    y2 = y2 + I*I
                    a3 = a3 + u*v
                    b3 = b3 + 2.*v*v+u*w
                    c3 = c3 + 3.*v*w
                    d3 = d3 + w*w
                    y3 = y3 + v*I
                end do
                den = d3*y2-c2*c2
                !! ASSERT(den != 0.0);
                !! coefficients of x**3 + ax**2 + bx + c = 0
                a = (c3*y2-c2*b2-c2*y3)/den
                b = (b3*y2-c2*a2-y3*b2)/den
                c = (a3*y2-y3*a2)/den
                k_m_test = solve_cubic_equation(a, b, c)
                !write(*, *) k_m_test
                do j = 1, 3
                    if (k_m_test(j) >= 0) then
                        r = r_factor_w_selection(tmp_scale * (f_calc + k_m_test(j) * f_mask), index_start, index_end)
                        if (r < r_best) then
                            r_best = r
                            k_best = k_m_test(j)
                        end if
                    end if
                end do
                !write(*, *) 'select the best', k_best, r_best
            end if
            if (k_best > 1) then
                k_best = 1.0d0
            end if
            r_best = r_factor_w_selection(tmp_scale * (f_calc + k_best * f_mask), index_start, index_end)
            upper_limit = k_best+shift+1.0e-3
            ! fine-tune by grid search around
            k_mask_per_bin_test = max(0.0d0,k_best-shift)
            !write(*, *) k_best, r_best, upper_limit, k_mask_per_bin_test
            do while (k_mask_per_bin_test <= upper_limit)
                f_calc_tmp(index_start:index_end) = tmp_scale(index_start:index_end) * (f_calc(index_start:index_end) + &
                                                    k_mask_per_bin_test * f_mask(index_start:index_end))
                k_overall_ = scale_selection(f_calc_tmp, index_start, index_end)
                r = r_factor_w_selection_and_scale(f_calc_tmp, index_start, index_end, k_overall_)
                !write(*,*) k_mask_per_bin_test, r
                if (r < r_best) then
                    ! write(*,*) k_mask_per_bin_test, r
                    k_best = k_mask_per_bin_test
                    r_best = r
                end if
                k_mask_per_bin_test = k_mask_per_bin_test + 0.01d0
            end do
            !write(*, *) 'best', k_best, r_best
            k_mask_bin(l) = k_best
        end do
        !write(*, '(a)') 'bulk_solvent_scaling, k_mask in resolution bins'
        !do j = 1, n_bins
        !    write(*, '(i20,a,f14.12)') j, ' ', k_mask_bin(j)
        !end do 
        call smooth_k_mask(k_mask_bin)
        write(*, '(a)') 'bulk_solvent_scaling, smoothed k_mask in resolution bins'
        do j = 1, n_bins
            write(*, '(i20,a,f8.4)') j, ' ', k_mask_bin(j)
        end do 
        call populate_k_mask_linear_interpolation(k_mask_bin, k_mask_test)
        call bin_k_isotropic_as_scale_k1(r_start, k_iso_test, k_mask_test)
        r = r_factor_w_scale(k_iso_test * k_aniso * k_iso_exp * (f_calc + k_mask_test * f_mask), 1.0d0)
        if (r - r_start<0.5d0/100) then
            k_iso = k_iso_test
            k_mask = k_mask_test
            r_start = r
            k_mask_bin_orig = k_mask_bin
        end if
        write(*, '(a,f8.6)') 'bulk_solvent_scaling, final r_factor = ', r_start
        
    end subroutine bulk_solvent_scaling

    function estimate_scale_k1(f_m) result (scale_k1)
        double precision :: scale_k1, cutoff, high_res, num, denum
        integer :: i, min_reflections
        complex(8), dimension(NRF) :: f_m
        ! overall scaling on min(best 1A shell, 500 reflection)
        min_reflections=500
        scale_k1 = 1d0
        !write(*, *) 'hi_res_shell_n', sum(f_m(scale_k1_indices(1:hi_res_shell_n)))
        if (hi_res_shell_n > min_reflections) then
            denum = sum(abs(f_m(scale_k1_indices(1:hi_res_shell_n))) * abs(f_m(scale_k1_indices(1:hi_res_shell_n))))
            num = sum(f_obs(scale_k1_indices(1:hi_res_shell_n)) * abs(f_m(scale_k1_indices(1:hi_res_shell_n))))
            if (denum == 0) then
                ! make test for zero f_obs in sf.F90
                scale_k1 = 0d0
            else
                scale_k1 = num / denum
            end if
        end if
        !write(*, *) 'hi_res_shell_k1', scale_k1

    end function estimate_scale_k1

    subroutine initialize_scales()
        integer :: i

        k_overall = 1d0
        k_mask = 0d0
        allocate(k_iso(NRF))
        allocate(k_iso_test(NRF))
        allocate(k_iso_exp(NRF))
        allocate(k_iso_exp_test(NRF))
        k_iso = 1d0
        k_iso_test = 1d0
        k_iso_exp = 1d0
        k_iso_exp_test = 1d0
        allocate(k_aniso(NRF))
        allocate(k_aniso_test(NRF))
        k_aniso = 1d0
        k_aniso_test = 1d0

    end subroutine initialize_scales

    function r_factor_w(f_m) result(r)
        double precision :: r, sc
        complex(8), dimension(NRF) :: f_m
        sc = scale(f_m)
        r = r_factor_w_scale(f_m, sc)

    end function r_factor_w
    
    function r_factor_w_selection(f_m, index_start, index_end) result(r)
        double precision :: r, num, denum, sc
        complex(8), dimension(NRF) :: f_m
        integer :: index_start, index_end
        sc = scale_selection(f_m, index_start, index_end)
        num = sum(abs(f_obs(index_start:index_end) - sc * abs(f_m(index_start:index_end))))
        denum = sum(f_obs(index_start:index_end))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 9999999
        else
            r = num / denum
        end if
    end function r_factor_w_selection
    
    function r_factor_w_selection_and_modified_fobs(f_m, modifier, index_start, index_end) result(r)
        double precision :: r, num, denum, sc, modifier(NRF)
        complex(8), dimension(NRF) :: f_m
        integer :: index_start, index_end
        sc = scale_selection_with_modified_fobs(f_m, modifier, index_start, index_end)
        num = sum(abs(f_obs(index_start:index_end)/modifier(index_start:index_end) - sc * abs(f_m(index_start:index_end))))
        denum = sum(f_obs(index_start:index_end)/modifier(index_start:index_end))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 9999999
        else
            r = num / denum
        end if
    end function

    function r_factor_w_selection_and_modified_fobs_and_scale(f_m, modifier, sc, index_start, index_end) result(r)
        double precision :: r, num, denum, sc, modifier(NRF)
        complex(8), dimension(NRF) :: f_m
        integer :: index_start, index_end
        num = sum(abs(f_obs(index_start:index_end)/modifier(index_start:index_end) - sc * abs(f_m(index_start:index_end))))
        denum = sum(f_obs(index_start:index_end)/modifier(index_start:index_end))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 9999999
        else
            r = num / denum
        end if
    end function
        
    function r_factor_w_selection_and_scale(f_m, index_start, index_end, sc) result(r)
        double precision :: r, num, denum, sc
        complex(8), dimension(NRF) :: f_m
        integer :: index_start, index_end
        num = sum(abs(f_obs(index_start:index_end) - sc * abs(f_m(index_start:index_end))))
        denum = sum(f_obs(index_start:index_end))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 9999999
        else
            r = num / denum
        end if
    end function
    
    function r_factor_w_scale(f_m, sc) result(r)
        double precision :: r, num, sc
        complex(8), dimension(NRF) :: f_m
        integer :: i
        num = sum(abs(f_obs(1:NRF_work) - sc * abs(f_m(1:NRF_work))))
        if (r_work_factor_denominator == 0) then
            ! make test for zero f_obs in sf.F90
            r = 9999999
        else
            r = num / r_work_factor_denominator
        end if

    end function r_factor_w_scale
    
    function scale(f_m) result(r)
        double precision :: r
        complex(8), dimension(NRF) :: f_m
        r = scale_selection(f_m, 1, NRF)
    end function
    
    function scale_selection(f_m, index_start, index_end) result(r)
        double precision :: r, num, denum
        complex(8), dimension(NRF) :: f_m
        integer :: i, index_start, index_end
        num = sum(f_obs(index_start: index_end) * abs(f_m(index_start: index_end)))
        denum = sum(abs(f_m(index_start: index_end)) * abs(f_m(index_start: index_end)))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 0d0
        else
            r = num / denum
        end if
    end function
    
    function scale_selection_with_modified_fobs(f_m, modifier, index_start, index_end) result(r)
        double precision :: r, num, denum, modifier(NRF)
        complex(8), dimension(NRF) :: f_m
        integer :: i, index_start, index_end
        num = sum(f_obs(index_start: index_end)/modifier(index_start: index_end) * abs(f_m(index_start: index_end)))
        denum = sum(abs(f_m(index_start: index_end)) * abs(f_m(index_start: index_end)))
        if (denum == 0) then
            ! make test for zero f_obs in sf.F90
            r = 0d0
        else
            r = num / denum
        end if
    end function
    
    function special_r_factor(f_m) result(r_best)
        double precision :: r, r_best, scale, scale_best, scale_, offset_factor, step_factor, step, num, denum
        complex(8), dimension(NRF) :: f_m
        offset_factor = 3.0
        step_factor = 20.0
        num = sum(f_obs(1:NRF_work) * abs(f_m(1:NRF_work)))
        denum =  sum(abs(f_m(1:NRF_work)) * abs(f_m(1:NRF_work)))
        if (denum > 0)  then
            scale = num/denum
        else
            scale = 0.0d0
        end if
        scale_ = scale - scale/offset_factor
        r_best = r_factor_w_scale(f_m, scale_)
        step = scale / step_factor
        do while (scale_ <= scale + scale/offset_factor)
            r = r_factor_w_scale(f_m, scale_)
            if (r < r_best) then
                r_best = r
                scale_best = scale_
            end if
            scale_ = scale_ + step
        end do
    end function

end module k_scale_mod
