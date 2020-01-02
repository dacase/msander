! <compile=optimized>
#include "../include/assert.fh"

module ml_mod

  use file_io_dat
  implicit none

  ! 7 x 7 transformation matrices to anisotropically scale structure factors 
  !       from calculated to experimental
  double precision, dimension(7, 7) :: Ucryst, MUcryst_inv

  ! Arrays used in the computation of structure factors and ML parameters
  double precision, dimension(:), allocatable :: alpha_array, beta_array, &
                                                 delta_array

  ! Arrays used in the estimation of maximum likelihood parameters,
  ! size is equal to the number of resolution bins/zones
  double precision, dimension(:), allocatable :: A_in_zones, B_in_zones, &
                       C_in_zones, q_in_zones, &
                       alpha_beta_bj, alpha_beta_OmegaI, alpha_beta_wi, &
                       t_optimal, alpha_in_zones, beta_in_zones, &
                       b_vector_base

  ! Convenient numerical constants
  double precision, parameter :: pi = 3.14159265359, zero = 0.0, d_tolerance = 1.e-10

  ! h_sq, k_sq, l_sq:         Squares of H, K, and L indices (still integers)
  ! hk, kl, hl:               Products of H, K, and L indices
  integer, dimension(:), allocatable :: h_sq, k_sq, l_sq, hk, kl, hl

  ! NRF:                    The number of reflections
  ! NRF_work, NRF_work_sq:  Number of work reflections, and square thereof
  ! NRF_free:               Number of free reflections
  ! call_est:               Tracking counter for pseudo-energy weight estimation
  ! N_steps:                The total number of simulation steps
  !                         Used to potentially control the restraints weight
  ! starting_N_step:        Step number to start on.
  ! total_N_steps:          Step number to finish on.
  ! n_bins:                 Number of reflections per resolution bin
  integer :: NRF,NRF_work, NRF_work_sq, NRF_free, &
             call_est, N_steps, starting_N_step, total_N_steps, n_bins

  ! bins_work_population:     Number of work reflections in each resolution zone
  ! bins_free_population:     Number of free reflections in each resolution zone
  ! bins_free_start_indices:  Stores the indices of first free reflection 
  !                           in each resolution zone
  ! reflection_bin:           HKL's number of the corresponding resolution zone
  integer, dimension(:), allocatable :: bins_work_population, &
           bins_free_population, bins_free_start_indices, reflection_bin

  ! hkl:                      3 by N carray containing HKL indices
  ! bins_work_reflections:    Used to re-sort work reflections into the 
  !                           original order based on resolution bin and 
  !                           relative position in it
  ! bins_free_reflections:    Used to re-sort free reflections
  integer, dimension(:,:), allocatable :: bins_work_reflections, &
                                          bins_free_reflections

#ifdef MPI
#     include "parallel.h"
#else
      integer :: mytaskid = 0
#endif
contains

  !----------------------------------------------------------------------------
  ! cross:   Compute the cross product of vectors a and b.  
  !----------------------------------------------------------------------------
  function cross(a, b) result(cross_product)

    double precision, dimension(3) :: cross_product
    double precision, dimension(3), intent(in) :: a, b

    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)

  end function cross

  !----------------------------------------------------------------------------
  ! i1_over_i0:   approximation of modified Bessel functions of the first kind:
  !               first order over zero order
  !----------------------------------------------------------------------------
  function i1_over_i0 (x) result(result)

    double precision, dimension(7) :: &
          p  = (/  1.00000000,  3.51562290,  3.08994240, &
                   1.20672920,  0.26597320,  0.03607680,  0.00458130 /), &
          pp = (/  0.50000000,  0.87890594,  0.51498869, &
                   0.15084934,  0.02658733,  0.00301532,  0.00032411 /)
    double precision, dimension(9) :: &
          q  = (/  0.39894228,  0.01328592,  0.00225319, &
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

  !----------------------------------------------------------------------------
  ! ln_of_i0:   approximation of logarithm of modified Bessel function of 
  !             the first kind, zero order
  !----------------------------------------------------------------------------
  function ln_of_i0 (x) result(bessel_lni0)

    double precision, intent(in) :: x
    double precision :: bessel_lni0
    double precision :: y, abs_x
    integer :: i

    abs_x = abs(x)
    if (abs_x/3.75 < 1.0) then
      y = x / 3.75
      y = y * y
      bessel_lni0 = 1.0 + y*(3.5156229 + &
                                y*(3.0899424 + &
                                   y*(1.2067492 + &
                                      y*(0.2659732 + &
                                         y*(0.0360768 + y*0.0045813)))))
      bessel_lni0 = log(bessel_lni0);
    else
      y = 3.75 / abs_x
      y = 0.39894228 + y*(0.01328592 + &
                          y*(0.00225319 + &
                             y*(-0.00157565 + &
                                y*(0.00916281 + &
                                   y*(-0.02057706 + &
                                      y*(0.02635537 + &
                                         y*(-0.01647633 + y*0.00392377)))))))
      bessel_lni0 = log(y) + abs_x - 0.5 * log(abs_x);
    end if

  end function ln_of_i0

  !----------------------------------------------------------------------------
  ! init_ml: gateway to the maximum-likelihood target function
  !          this routine also intializes anisotropic scaling data
  !----------------------------------------------------------------------------
  subroutine init_ml(target, nstlim, num_hkl, hkl_tmp, &
                     f_obs_tmp, f_obs_sigma_tmp, test_flag, d_star_sq_out, &
                     resolution)

    use xray_globals_module, only: unit_cell
    implicit none

    character(len=4), intent(in) :: target
    integer, intent(in) :: nstlim, num_hkl
    integer, dimension(3,num_hkl), intent(inout) :: hkl_tmp
    double precision, dimension(num_hkl), intent(inout) :: f_obs_tmp, &
            f_obs_sigma_tmp
    integer, dimension(num_hkl), intent(inout) :: test_flag
    double precision, dimension(num_hkl), intent(out) :: d_star_sq_out
    double precision, intent(out) :: resolution

    double precision :: a, b, c, alpha, beta, gamma, V, &
                        d, d_star, reflections_per_bin, fo_fo
    double precision :: cosa, sina, cosb, sinb, cosg, sing
    double precision, dimension(3) :: va, vb, vc, vas, vbs, vcs
    double precision, dimension(:), allocatable :: d_star_sq, &
                                                   d_star_sq_sorted, bin_limits
    integer, dimension(:), allocatable :: counter_w, counter_f
    integer (kind = 4) :: file_status
    integer :: d_i, i, j, k, na, nb, nc, mdout = 6
    integer :: r_free_counter, r_free_flag, counter_sort, index_sort
    integer:: hkl(3,num_hkl)
    double precision :: f_obs(num_hkl), fobs_weight(num_hkl), f_obs_sigmas(num_hkl)

    N_steps = nstlim
    NRF = num_hkl

    allocate(d_star_sq(NRF))
    allocate(reflection_bin(NRF))
    allocate(delta_array(NRF))
    allocate(alpha_array(NRF))
    alpha_array(:) = 1.0
    allocate(beta_array(NRF))
    beta_array(:) = 1.0

    r_free_counter = 0
    NRF_work = 0

    ! Following block of work might be done somewhere else.
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

    ! start of block to separate work and free reflections, and sort into bins

    resolution = 50.0
    do i = 1, NRF
      d_star =  (square(hkl_tmp(1,i) * norm2(vas)) + square(hkl_tmp(2,i) * norm2(vbs)) + &
                 square(hkl_tmp(3,i) * norm2(vcs)) + &
                 2 * hkl_tmp(2,i) * hkl_tmp(3,i) * dot_product(vbs, vcs) + &
                 2 * hkl_tmp(1,i) * hkl_tmp(3,i) * dot_product(vas, vcs) + &
                 2 * hkl_tmp(1,i) * hkl_tmp(2,i) * dot_product(vbs, vas))
      d = sqrt(1.0 / d_star)
      resolution = min( d, resolution )
      if (test_flag(i) == 1) then
        NRF_work = NRF_work + 1
      else
        r_free_counter = r_free_counter + 1
      endif
      d_star_sq(i) = d_star
      ! f_obs_weight(i) = 1.0 ! (1.0/f_obs_sigma_tmp(i))**2
    enddo

    NRF_free = NRF - NRF_work
    REQUIRE( NRF_free == r_free_counter )
    if (mytaskid == 0 ) &
      write(6,'(a,3i8)') '| number of reflections: ', NRF_work, NRF_free, NRF

    if (target(1:2) == 'ml') then

    ! Sort reflections for binning
    allocate(b_vector_base(NRF_work))
    allocate(d_star_sq_sorted(NRF_free))
    r_free_counter = 0
    do i = 1, NRF
      if (test_flag(i) == 0) then
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

    allocate(bin_limits(n_bins + 1))
    bin_limits(1) = max(zero, d_star_sq_sorted(1) * (1 - d_tolerance))
    do i = 2, n_bins
      d_i = int((i-1) * reflections_per_bin + 0.5) + 1
      bin_limits(i) = d_star_sq_sorted(d_i) * (1 - d_tolerance)
      REQUIRE (d_i /= NRF_free)
    end do
    bin_limits(n_bins + 1) = d_star_sq_sorted(NRF_free) * (1 + d_tolerance)
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
      if (test_flag(i) == 1) then
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
      if (test_flag(i) == 1) then
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
        hkl(:,counter_sort) = hkl_tmp(:,index_sort)
        d_star_sq_out(counter_sort) = d_star_sq(index_sort)
      end do
    end do
    allocate(bins_free_start_indices(n_bins))
    do i = 1, n_bins
      bins_free_start_indices(i) = counter_sort + 1
      ! write(mdout, *) i, bins_free_start_indices(i), bins_free_population(i)
      do j = 1, bins_free_population(i)
        counter_sort = counter_sort + 1
        index_sort = bins_free_reflections(i, j)
        f_obs_sigmas(counter_sort) = f_obs_sigma_tmp(index_sort)
        f_obs(counter_sort) = f_obs_tmp(index_sort)
        hkl(:,counter_sort) = hkl_tmp(:,index_sort)
        d_star_sq_out(counter_sort) = d_star_sq(index_sort)
      end do
    end do

    ! need to put fobs into f_obs_tmp to send back to calling 
    !   program; same for f_obs_sigmas, f_obs_weight(?), hkl and test_flag

    f_obs_tmp(:) = f_obs(:)
    f_obs_sigma_tmp(:) = f_obs_sigmas(:)
    ! f_obs_weight_tmp(:) = f_obs_weight(:)   !unused for now
    hkl_tmp(:,:) = hkl(:,:)
    test_flag(1:NRF_work) = 1
    test_flag(NRF_work+1:NRF) = 0

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

    endif  ! only done for target=='ml'

    ! Setting up anisotropic scaling parameters:

    allocate(h_sq(NRF))
    allocate(k_sq(NRF))
    allocate(l_sq(NRF))
    allocate(hk(NRF))
    allocate(hl(NRF))
    allocate(kl(NRF))

    do i = 1, NRF
      h_sq(i) = hkl(1,i) * hkl(1,i)
      k_sq(i) = hkl(2,i) * hkl(2,i)
      l_sq(i) = hkl(3,i) * hkl(3,i)
      hk(i) = hkl(1,i) * hkl(2,i)
      hl(i) = hkl(1,i) * hkl(3,i)
      kl(i) = hkl(2,i) * hkl(3,i)
    end do
    NRF_work_sq = NRF_work * NRF_work
    Ucryst(1, 1) = 1.0 / NRF_work
    Ucryst(1, 2) = sum(1.0 * h_sq(1:NRF_work) / NRF_work_sq)
    Ucryst(1, 3) = sum(1.0 * k_sq(1:NRF_work) / NRF_work_sq)
    Ucryst(1, 4) = sum(1.0 * l_sq(1:NRF_work) / NRF_work_sq)
    Ucryst(1, 5) = sum(1.0 *   hk(1:NRF_work) / NRF_work_sq)
    Ucryst(1, 6) = sum(1.0 *   hl(1:NRF_work) / NRF_work_sq)
    Ucryst(1, 7) = sum(1.0 *   kl(1:NRF_work) / NRF_work_sq)

    ! In case if one needs only isotropic scaling
    ! Ucryst(1, 1) = 1.0
    ! Ucryst(1, 2) = 0.0
    ! Ucryst(1, 3) = 0.0
    ! Ucryst(1, 4) = 0.0
    ! Ucryst(1, 5) = 0.0
    ! Ucryst(1, 6) = 0.0
    ! Ucryst(1, 7) = 0.0

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

    call_est = 1

    return

  end subroutine init_ml

  !----------------------------------------------------------------------------
  ! square:  compute the square of a number
  !----------------------------------------------------------------------------
  pure function square(x)
    double precision, intent(in) :: x
    double precision :: square
    square = x * x
  end function

  !----------------------------------------------------------------------------
  ! inverse: Matrix inverse, based on Doolittle LU factorization for Ax=b 
  !   by Alex G., December 2009.  The original matrix a(n,n) will be destroyed.
  !
  ! Arguments:
  !   a:       (n x n) array of coefficients for square matrix A
  !   n:       dimension of A
  !   c:       (n x n) inverse of matrix A
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! bubble_sort:  go bubble sort!  Sort a list of real numbers a.
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! A_B_C_D_omega:  start alpha and beta estimation according to
  !                 https://doi.org/10.1107/S010876739500688X
  !
  ! Arguments:
  !   zone:       number of resolution bin
  !   f_calc_:    complex array of computed structure factors
  !----------------------------------------------------------------------------
  subroutine A_B_C_D_omega(zone, f_calc_,f_obs)

    
    integer, intent(in) :: zone
    complex(8), intent(in) :: f_calc_(NRF)
    double precision, intent(in) :: f_obs(NRF)
    integer :: i, index_start, index_end
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

  !----------------------------------------------------------------------------
  ! blamm:  \Lambda(t) in the original paper
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! t_optimal_function:  G(t) in the original paper
  !----------------------------------------------------------------------------
  function t_optimal_function(t, zone) result(result)

    double precision :: t
    double precision :: result
    integer :: zone

    result = sqrt(1. + 4. * A_in_zones(zone) * B_in_zones(zone) * t * t) - 1. - &
             2. * t * blamm(t, zone)
  end function t_optimal_function

  !----------------------------------------------------------------------------
  ! solvm:  binary search for a root of t_optimal_function(t) in a general case
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! estimate_t_optimal:  find root of G(t)
  !----------------------------------------------------------------------------
  subroutine estimate_t_optimal(zone, f_calc_, f_obs)

    integer, intent(in) :: zone
    complex(8), intent(in) :: f_calc_(NRF)
    double precision, intent(in) :: f_obs(NRF)

    call A_B_C_D_omega(zone, f_calc_,f_obs)
    if (alpha_beta_OmegaI(zone) <= 0.0) then
      t_optimal(zone) = 0.0
    else if (alpha_beta_wi(zone)/(A_in_zones(zone)*B_in_zones(zone)) <= 3.0E-07) then
      t_optimal(zone) = 1.0e+10
    else
      call solvm(zone)
    end if

  end subroutine estimate_t_optimal

  !----------------------------------------------------------------------------
  ! smooth: smoothes t_optimal values over resuluion bins in maximum 
  !         likelihood (ML) estimation
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! alpha_beta_in_zones: used in maximum likelihood (ML) estimation
  !----------------------------------------------------------------------------
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

  !----------------------------------------------------------------------------
  ! alpha_beta_all:  expand alpha and beta ML parameters on all reflections
  !----------------------------------------------------------------------------
  subroutine alpha_beta_all

    integer :: i

    do i = 1, NRF
      alpha_array(i) = alpha_in_zones(reflection_bin(i))
      beta_array(i)  = beta_in_zones(reflection_bin(i))
    end do

  end subroutine alpha_beta_all

  !----------------------------------------------------------------------------
  ! estimate_alpha_beta:  estimate alpha and beta in resolution bins,
  !      as described in https://doi.org/10.1107/S010876739500688X, Appendix A:
  !      eq. 29     - estimate_t_optimal() estimates root of function G(t)
  !      eq. 30, 31 - alpha_beta_in_zones() calculates alpha and beta from t
  !----------------------------------------------------------------------------
  subroutine estimate_alpha_beta(f_calc_, f_obs)

    complex(8), intent(in) :: f_calc_(NRF)
    double precision, intent(in) :: f_obs(NRF)
    integer :: i

    do i = 1, n_bins
      call estimate_t_optimal(i, f_calc_, f_obs)
    end do
    call smooth(t_optimal)
    call alpha_beta_in_zones()
    call smooth(alpha_in_zones)
    call smooth(beta_in_zones)
    call alpha_beta_all()

  end subroutine estimate_alpha_beta

#if 0
  !----------------------------------------------------------------------------
  ! estimate_ml_parameters: currently this is only implemented on the CPU.
  !                         Also, get the xray restraint energy
  !
  ! Arguments:
  !   f_calc:    computed structure factors 
  !   exay:      xray restraint energy
  !   nstep:     the number of steps, passed down all the way from runmd
  !----------------------------------------------------------------------------
  subroutine estimate_ml_parameters(f_calc, f_obs, exray, nstep)
        
    use xray_globals_module, only: xray_weight
    implicit none
    complex(8), intent(in)   :: f_calc(NRF)
    double precision, intent(in) :: f_obs(NRF)
    double precision, intent(out) :: exray
    integer, intent(in) :: nstep
    double precision :: f_calc_abs(NRF)
    integer :: i
    
    if (mod(nstep, ml_update_frequency) == 0 .or. nstep == 0) then
      call estimate_alpha_beta(f_calc, f_obs)
      delta_array = alpha_array / beta_array
    endif

    f_calc_abs = abs(f_calc)


    ! ML target function
    exray = xray_weight * &
            sum(alpha_array(1:NRF_work) * delta_array(1:NRF_work) * &
            f_calc_abs(1:NRF_work) * f_calc_abs(1:NRF_work) - &
            ln_of_i0(2 * delta_array(1:NRF_work) * f_calc_abs(1:NRF_work) * &
            f_obs(1:NRF_work)))

  end subroutine estimate_ml_parameters
#endif

end module ml_mod
