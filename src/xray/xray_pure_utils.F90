! This module contains functions/subrotines without sideffects
! to facilitate unit testing

#include "../include/assert.fh"
module xray_pure_utils
  ! MUST NOT depend on any other module
  
  use xray_contracts_module
  
  implicit none

  double precision, private :: r_dummy ! FIXME: enforce same kind as in xray_global_module
                                       !        without injecting explicit module dependency
  integer, parameter :: real_kind = kind(r_dummy)
  real(8), parameter :: PI = 3.1415926535897932384626433832795d0 ! FIXME: import from global constants module
  real(4), private, parameter :: PI_real_4 = 3.1415926535897932384626433832795 ! FIXME: import from global constants module

  interface cross_product
    module procedure cross_product_real_4
    module procedure cross_product_real_8
  end interface cross_product
  
  interface is_sorted
    module procedure is_sorted_real_4
    module procedure is_sorted_real_8
    module procedure is_sorted_integer_4
    module procedure is_sorted_integer_8
  end interface is_sorted

  interface to_radians
    module procedure to_radians_4
    module procedure to_radians_8
  end interface to_radians

  interface to_degrees
    module procedure to_degrees_4
    module procedure to_degrees_8
  end interface to_degrees
  
contains
  
  !--------------------------------------------------------------------------------------------
  ! ln_of_i0:   approximation of logarithm of modified Bessel function of the first kind,
  !             zero order
  !--------------------------------------------------------------------------------------------
  elemental real(real_kind)  function ln_of_i0 (x) result(bessel_lni0)
    implicit none
    real(real_kind), intent(in):: x
    real(real_kind) :: y, abs_x
    
    abs_x = abs(x)
    
    if (abs_x/3.75d0 < 1.0d0) then
      y = x / 3.75d0
      y = y * y
      bessel_lni0 = 1.0d0 + y*(3.5156229d0 + &
          y*(3.0899424d0 + &
              y*(1.2067492d0 + &
                  y*(0.2659732d0 + &
                      y*(0.0360768d0 + y*0.0045813d0)))))
      bessel_lni0 = log(bessel_lni0);
    else
      y = 3.75d0 / abs_x
      y = 0.39894228d0 + y*(0.01328592d0 + &
          y*(0.00225319d0 + &
              y*(-0.00157565d0 + &
                  y*(0.00916281d0 + &
                      y*(-0.02057706d0 + &
                          y*(0.02635537d0 + &
                              y*(-0.01647633d0 + y*0.00392377d0)))))))
      bessel_lni0 = log(y) + abs_x - 0.5d0 * log(abs_x);
    end if
  
  
  end function ln_of_i0
  
  !--------------------------------------------------------------------------------------------
  ! i1_over_i0:   approximation of modified Bessel functions of the first kind:
  !               first order over zero order
  !--------------------------------------------------------------------------------------------
  elemental real(real_kind) function i1_over_i0 (x) result(result)
    implicit none
    real(real_kind), dimension(7), parameter :: p = (/  1.00000000d0, 3.51562290d0, 3.08994240d0, &
        1.20672920d0, 0.26597320d0, 0.03607680d0, &
        0.00458130d0 /), &
        pp = (/  0.50000000d0, 0.87890594d0, 0.51498869d0, &
            0.15084934d0, 0.02658733d0, 0.00301532d0, &
            0.00032411d0 /)
    real(real_kind), dimension(9), parameter :: q = (/  0.39894228d0, 0.01328592d0, 0.00225319d0, &
        -0.00157565d0, 0.00916281d0, -0.02057706d0, &
        0.02635537d0, -0.01647633d0, 0.00392377d0 /), &
        qq = (/  0.39894228d0, -0.03988024d0, -0.00362018d0, &
            0.00163801d0, -0.01031555d0, 0.02282967d0, &
            -0.02895312d0, 0.01787654d0, -0.00420059d0 /)
    real(real_kind), intent(in) :: x
    real(real_kind) :: y, be1, be0, pow_y_i, abs_x
    integer :: i
    be1 = 0d0
    be0 = 0d0
    abs_x = x
    if (abs_x < 0d0) then
      abs_x = -abs_x
    end if
    if (abs_x < 3.75d0) then
      y = x / 3.75d0
      y = y * y
      pow_y_i = 1d0
      do i = 1, 7
        be0 = be0 + p(i) * pow_y_i
        be1 = be1 + x * pp(i) * pow_y_i
        pow_y_i = pow_y_i * y
      end do
    else
      y = 3.75d0 / abs_x
      pow_y_i = 1d0
      do i = 1, 9
        be0 = be0 + q(i) * pow_y_i
        be1 = be1 + qq(i) * pow_y_i
        pow_y_i = pow_y_i * y
      end do
    end if
    result = be1 / be0;
    if (x < 0.0d0 .and. result > 0.0d0) then
      result = -result
    end if
  
  end function i1_over_i0
  
  !------------------------------------------------------------------
  ! Partition `indices` according to `move_to_start` inidicator array
  !------------------------------------------------------------------
  subroutine index_partition(move_to_start, indices)
    logical, dimension(:), intent(in) :: move_to_start
    integer, dimension(:), intent(inout) :: indices
    integer :: j, first, temp_i
    
    first = size(move_to_start) + 1
    
    do j = 1, size(move_to_start)
      if (.not. move_to_start(indices(j))) then
        first = j
        exit
      end if
    end do
    
    do j = first + 1, size(move_to_start)
      if (move_to_start(indices(j))) then
        temp_i = indices(first)
        indices(first) = indices(j)
        indices(j) = temp_i
        first = first + 1
      end if
    end do
  
  end subroutine index_partition
  
  pure function sorted(a) result(result)
    real(real_kind), intent(in) :: a(:)
    real(real_kind) :: result(size(a))
    
    integer :: indices(size(a))
    integer :: i
    
    do i = 1, size(indices)
      indices(i) = i
    end do
    
    call index_sort(a, indices)
    result = a(indices)
  
  end function sorted
  
  !--------------------------------------------------------------
  ! sort real-valued array `a` indiced by `indices`
  !--------------------------------------------------------------
  pure subroutine index_sort(a, indices)
    real(real_kind), intent(in) :: a(:)
    integer, intent(inout) :: indices(:)
    
    call index_hybrid_quick_sort(a, indices)
    
  end subroutine index_sort
  
  
  pure subroutine index_insertion_sort(a, indices)
    real(real_kind), intent(in) :: a(:)
    real(real_kind) :: key
    integer, intent(inout) :: indices(:)
    integer :: temp_i
    integer :: i, j
    
    do i = 1, size(indices)
      key = a(indices(i))
      temp_i = indices(i)
      j = i - 1
      do while (j > 0)
        if (a(indices(j)) <= key) then
          exit
        end if
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      indices(j + 1) = temp_i
    end do
  end subroutine index_insertion_sort
  
  pure recursive subroutine index_quick_sort(a, indices)
    real(real_kind), intent(in) :: a(:)
    integer, intent(inout) :: indices(:)
    real(real_kind) :: pivot
    integer :: i, j, tmp

    pivot = a( indices(size(indices) / 2) )
    i = 1
    j = size(indices)
    do
      do while (a(indices(i)) < pivot)
        i = i + 1
      end do
      do while (pivot < a(indices(j)))
        j = j - 1
      end do
      if (i >= j) exit
      tmp = indices(i);  indices(i) = indices(j);  indices(j) = tmp
      i = i + 1
      j = j - 1
    end do
    if (1 < i - 1) then
      call index_quick_sort(a, indices(:i - 1))
    end if
    if (j + 1 < size(indices)) then
      call index_quick_sort(a, indices(j + 1:))
    end if
  end subroutine index_quick_sort
  
  pure recursive subroutine index_hybrid_quick_sort(a, indices)
    real(real_kind), intent(in) :: a(:)
    integer, intent(inout) :: indices(:)
    real(real_kind) :: pivot
    integer :: i, j, tmp
    
    if ( size(indices) < 16 ) then
      call index_insertion_sort(a, indices)
      return
    end if
    
    pivot = a( indices(size(indices) / 2) )
    i = 1
    j = size(indices)
    do
      do while (a(indices(i)) < pivot)
        i = i + 1
      end do
      do while (pivot < a(indices(j)))
        j = j - 1
      end do
      if (i >= j) exit
      tmp = indices(i);  indices(i) = indices(j);  indices(j) = tmp
      i = i + 1
      j = j - 1
    end do
    if (1 < i - 1) then
      call index_hybrid_quick_sort(a, indices(:i - 1))
    end if
    if (j + 1 < size(indices)) then
      call index_hybrid_quick_sort(a, indices(j + 1:))
    end if
  end subroutine index_hybrid_quick_sort
  
  ! Given a logical-selection array, return an array of indices into the
  ! selection's target array, and the number of selected indices.
  pure function pack_index(selection) result(result)
    implicit none
    logical, intent(in) :: selection(:)
    integer, allocatable :: result(:)
    integer :: i, n
    allocate(result(count(selection)))
    n = 0
    do i = 1, size(selection)
      if (selection(i)) then
        n = n + 1
        result(n) = i
      end if
    end do
  end function pack_index
  
  !--------------------------------------------------------------------------------------------
  ! t_optimal_function:  G(t) in the original paper
  !--------------------------------------------------------------------------------------------
  pure function t_optimal_function(t, A, B, fo_fm) result(result)
    
    double precision, intent(in) :: t, A, B
    real(real_kind), intent(in) :: fo_fm(:)
    double precision :: result
    
    result = sqrt(1.0d0 + 4.0d0 * A * B * t ** 2) - 1.0d0 - 2.0d0 * t * blamm(t, fo_fm)
  end function t_optimal_function
  
  
  !--------------------------------------------------------------------------------------------
  ! blamm:  \Lambda(t) in the original paper
  !--------------------------------------------------------------------------------------------
  pure function blamm(t, alpha_beta_bj) result(result)
    implicit none
    real(real_kind), intent(in) :: t
    real(real_kind), intent(in) :: alpha_beta_bj(:)
    real(real_kind) :: result
    
    result = sum(alpha_beta_bj * i1_over_i0(2 * t * alpha_beta_bj)) / size(alpha_beta_bj)
    
  end function blamm
  
  
  !--------------------------------------------------------------------------------------------
  ! solvm:  binary search for a root of t_optimal_function(t) in a general case
  !--------------------------------------------------------------------------------------------
  pure function solvm(A, B, C, alpha_beta_wi, fo_fm) result(t_optimal)
    implicit none
    real(real_kind), intent(in) :: A, B, C, alpha_beta_wi
    real(real_kind), intent(in) :: fo_fm(:)
    real(real_kind) :: t_optimal

    integer, parameter :: nst1 = 10 * 5 + 1  ! TODO: add docstring
    integer, parameter :: nst2 = 20 * 5 + 1  ! TODO: add docstring

    real(real_kind) :: eps, fgtopt, tl, fgtl, tr, fgtr
    integer :: n1, n2

    eps = 0.00001d0
    t_optimal = 0.0d0
    fgtopt = 0.0d0
    tl = C / alpha_beta_wi
    fgtl = t_optimal_function(tl, A, B, fo_fm)
    do n1 = 1, nst1
      tr = tl
      fgtr = fgtl
      tl = tr * 0.5d0
      fgtl = t_optimal_function(tl, A, B, fo_fm)
      if (fgtl == 0.0d0) then
        t_optimal = tl
        fgtopt = 0.0d0
        return
      else if (fgtl < 0.0d0) then
        t_optimal = tr
        fgtopt = fgtr
        do n2 = 1, nst2
          if ((tr-tl) < eps * t_optimal) then
            return
          end if
          t_optimal = (tl * fgtr - tr * fgtl)/(fgtr-fgtl)
          fgtopt = t_optimal_function(t_optimal, A, B, fo_fm)
          if (fgtopt > 0.0d0) then
            tr = t_optimal
            fgtr = fgtopt
          else
            tl = t_optimal
            fgtl = fgtopt
          end if
        end do
      end if
    end do
  
  end function solvm
  
  !--------------------------------------------------------------------------------------------
  ! estimate_t_optimal:  find root of G(t) estimation according to
  !                      https://doi.org/10.1107/S010876739500688X
  !
  !--------------------------------------------------------------------------------------------
  pure function estimate_t_optimal(abs_Fobs, abs_Fcalc, B, q, A, p) result(result)
    implicit none
    
    real(real_kind), intent(in) :: abs_Fobs(:)               !< Experimental structure moduli factor
    real(real_kind), intent(in) :: abs_Fcalc(size(abs_Fobs)) !< Model structure moduli factor
    real(real_kind), intent(in) :: B                         !< Sum of abs_Fobs ** 2 (B in eq. 18)
    real(real_kind), intent(in) :: q                         !< Sum of abs_Fobs ** 2
    real(real_kind), intent(in) :: A                         !< Sum of abs_Fcalc ** 2 (A in eq. 18)
    real(real_kind), intent(in) :: p                         !< Sum of abs_Fcalc ** 2
    real(real_kind) :: result
  
    real(real_kind) :: C, D, r
    real(real_kind) :: OmegaI, wi
    
    real(real_kind) :: fo_fm(size(abs_Fobs))
    
    fo_fm = abs_Fobs * abs_Fcalc
  
    C = sum(fo_fm) / size(abs_Fobs)
    D = sum(fo_fm ** 2) / size(abs_Fobs)
  
    r = (p - A ** 2) * (q - B ** 2)
  
    if (r > 0) then
      OmegaI = (D - A * B) / sqrt(r)
    else
      OmegaI = 0
    end if

    wi = A * B - C ** 2
    
    if (OmegaI <= 0.0d0) then
      result = 0
    else if (wi / (A * B) <= 3.0E-07) then
      result = 1.0e+10
    else
      result = solvm(A, B, C, wi, fo_fm)
    end if
  
  end function estimate_t_optimal
  
  
  !--------------------------------------------------------------------------------------------
  ! calc_bin_alpha_beta: used in maximum likelihood (ML) estimation
  !--------------------------------------------------------------------------------------------
  pure subroutine calc_bin_alpha_beta(t, A, B, alpha, beta)
    
    real(real_kind), intent(in) :: t
    real(real_kind), intent(in) :: A
    real(real_kind), intent(in) :: B
    real(real_kind), intent(out) :: alpha
    real(real_kind), intent(out) :: beta
    
    double precision :: tt, ww, hbeta
    
    if (t == 0.0d0) then
      alpha = 0.0d0
      beta = B
    else if (t >= 1.0e+10) then
      alpha = sqrt(A / B)
      beta = 1.e-10
    else
      tt = 2.0d0 * t
      ww = sqrt(1.0d0 + A * B * tt ** 2)
      hbeta = B / (ww + 1.0d0)
      alpha = sqrt(hbeta * (ww - 1.0d0) / A)
      beta = 2.0d0 * hbeta
    end if
  
  end subroutine calc_bin_alpha_beta
  
  subroutine create_equisized_bins(n_elements, n_bins, bin_start, bin_count)
    
    implicit none
    integer, intent(in) :: n_elements
    integer, intent(in) :: n_bins
    integer, intent(out) :: bin_start(n_bins)
    integer, intent(out) :: bin_count(n_bins)
    
    real(real_kind) :: n_elements_per_bin
    integer :: i
    
    n_elements_per_bin = n_elements * 1.0 / n_bins
    
    bin_start(1) = 1
    do i = 1, n_bins - 1
      bin_start(i + 1) = int(n_elements_per_bin * i + 0.5) + 1
      bin_count(i) = int(n_elements_per_bin * i + 0.5) - int(n_elements_per_bin * (i - 1) + 0.5)
    end do
    bin_count(n_bins) = n_elements - bin_start(n_bins) + 1
  
  end subroutine create_equisized_bins
  

  subroutine create_logspace_resolution_bins(resolution, &
                                             n_reflections_in_worst_resolution_bin, &
                                             min_bin_size, &
                                             bin_start, bin_size, &
                                             n_bins)
    ! Create resolution bins on a logarithmic scale.
    !
    ! Note: implementation is slightly different from cctbx (log_binning)
    !
    ! See Urzhumtsev et al. (2009) Acta Crystallogr D Biol Crystallogr. 65:1283-91.
    !     https://doi.org/10.1107/S0907444909039638
    !
    
    implicit none
    real(real_kind), intent(in) :: resolution(:)
    integer, intent(in) :: n_reflections_in_worst_resolution_bin
    integer, intent(in) :: min_bin_size
    integer, intent(out) :: bin_start(:)
    integer, intent(out) :: bin_size(size(bin_start))
    integer, intent(out) :: n_bins

    ! Preconditions:
    ASSERT(is_sorted(resolution))
    ASSERT(all(resolution > 0))
    ASSERT(size(bin_start) > 0)
    ASSERT(min_bin_size > 0)
    
    if (size(resolution) <= n_reflections_in_worst_resolution_bin * 2 &
        .or. size(bin_start) == 1 &
      ) then
      bin_start(1) = 1
      bin_size(1) = size(resolution)
      n_bins = 1
      return
    end if
    
    ! Partition all resolutions but
    call create_equiwide_bins( &
        log(resolution(:size(resolution) - n_reflections_in_worst_resolution_bin)), &
        min_bin_size, &
        log(resolution(size(resolution) - n_reflections_in_worst_resolution_bin)) - log(resolution(size(resolution) - 2 * n_reflections_in_worst_resolution_bin)), &
        bin_start(:size(bin_start) - 1), &
        bin_size(:size(bin_start) - 1), &
        n_bins &
    )
    
    ! Fix-ups for last bin with worst resolution
    n_bins = n_bins + 1
    bin_start(n_bins) = size(resolution) - n_reflections_in_worst_resolution_bin + 1
    bin_size(n_bins) = n_reflections_in_worst_resolution_bin
  
  end subroutine create_logspace_resolution_bins
  
  
  subroutine create_equiwide_bins(values, min_bin_size, min_step, bin_start, bin_size, n_bins)
    ! Partition `values` into bins of equal width
    
    implicit none
    real(real_kind), intent(in) :: values(:)           !< Sorted values
    integer, intent(in) :: min_bin_size                !< Minimum bin size
    real(real_kind), intent(in) :: min_step
    integer, intent(out) :: bin_start(:)
    integer, intent(out) :: bin_size(size(bin_start))
    integer, intent(out) :: n_bins
    
    real(real_kind) :: step
    real(real_kind) :: right_bound
    
    integer :: i, bin_end

    ! Preconditions:
    ASSERT(is_sorted(values))
    ASSERT(size(values) > 0)
    ASSERT(size(bin_start) > 0)
    ASSERT(size(bin_start) == size(bin_size))
    ASSERT(min_bin_size > 0)
    ASSERT(min_step >= 0)
    
    step = max((values(size(values)) - values(1)) / size(bin_start), min_step)
    i = size(values)
    n_bins = 0
    
    ! Find bins in reverse order
    do while (i >= min_bin_size .and. n_bins < size(bin_start))
      n_bins = n_bins + 1
      bin_end = i
      right_bound = values(i)
      i = i - min_bin_size
      do while (i > 1)
        if (.not. (values(i) > right_bound - step)) then
          exit
        end if
        i = i - 1
      end do
      bin_start(n_bins) = i + 1
      bin_size(n_bins) = bin_end - bin_start(n_bins) + 1
    end do
    
    if (bin_start(n_bins) > 1) then
      bin_size(n_bins) = bin_size(n_bins) + (bin_start(n_bins) - 1)
      bin_start(n_bins) = 1
    end if

    if (n_bins > 2) then
      ! Unite two last bins if last one is smaller than quarter of previous
      if (bin_size(n_bins - 1) / 4 > bin_size(n_bins)) then
        bin_size(n_bins - 1) = bin_size(n_bins - 1)  + bin_size(n_bins)
        n_bins = n_bins - 1
      end if
    end if

    bin_start(n_bins) = 1

    ! Fix-up: reverse order to return bins in increasing order
    bin_start(1:n_bins) = bin_start(n_bins:1:-1)
    bin_size(1:n_bins) = bin_size(n_bins:1:-1)
  
  end subroutine create_equiwide_bins
  
  
  subroutine assign_resolution_bin_indices(resolution, max_bin_resolution, resolution_bin_index)
    real(real_kind), intent(in) :: resolution(:)
    real(real_kind), intent(in) :: max_bin_resolution(:)
    integer, intent(out) :: resolution_bin_index(size(resolution))
    
    integer :: i, j
    
    ! Preconditions:
    ASSERT(is_sorted(resolution))
    ASSERT(is_sorted(max_bin_resolution))
    
    i = 1 ! bin index
    do j = 1, size(resolution)
      if (resolution(j) > max_bin_resolution(i)) then
        i = min(i + 1, size(max_bin_resolution))
      end if
      resolution_bin_index(j) = i
    end do
  
  end subroutine assign_resolution_bin_indices
  
  
  subroutine exponential_fit_1d_analytical(y, z, x, a, b)
    ! Fit
    !     y = a * z * exp(-b * x)
    ! via (a, b)
    real(real_kind), intent(in) :: y(:)
    real(real_kind), intent(in) :: z(size(y))
    real(real_kind), intent(in) :: x(size(y))
    
    real(real_kind), intent(out) :: a  !< Multiplier
    real(real_kind), intent(out) :: b  !< Exponential decay
    
    real(real_kind) :: r, p, q, s, den, u
    real(real_kind) :: y_over_z
    real(real_kind) :: log_yz
    integer :: i, n

    ! Preconditions:
    ASSERT(count(y > 0) > 0)
    ASSERT(count(z > 0) > 0)
    
    p = 0
    q = 0
    r = 0
    s = 0
    n = 0 ! Number of valid points
    
    do i = 1, size(y)
      y_over_z = y(i) / z(i)
      if (y_over_z > 0) then
        log_yz = log(y_over_z)
        p = p + log_yz
        q = q + x(i)
        r = r + x(i) ** 2
        s = s + x(i) * log_yz
        n = n + 1
      end if
    end do
    
    if (r /= 0) then
      den = n - q ** 2 / r
      if (den /= 0) then
        u = (p - s * q / r) / den
        b = (u * q - s) / r
        a = exp(u)
        return
      end if
    end if
    
    ! Fitting failed, return dummy exponet
    ! Should return nans instead?
    b = 0
    a = 1
  
  end subroutine exponential_fit_1d_analytical
  
  
  !---------------------------------------------------------------------------------------
  ! Calculate k_overall that minimizes least square residual between k*|Fcalc| and |Fobs|
  !
  ! Ref: eq. (5) from https://doi.org/10.1107/S0907444913000462
  !---------------------------------------------------------------------------------------
  function calc_k_overall(abs_Fobs, abs_Fcalc) result(result)
    implicit none
    real(real_kind), intent(in) :: abs_Fobs(:) !< Magnitude of experimnal structure factors
    real(real_kind), intent(in) :: abs_Fcalc(size(abs_Fobs)) !< Magnitude of model structure factors
    real(real_kind):: result
    real(real_kind):: denom
    
    ! Precondition
    ASSERT(size(abs_Fobs) == size(abs_Fcalc))
    ASSERT(all(abs_Fobs >= 0))
    ASSERT(all(abs_Fcalc >= 0))
    
    denom = sum(abs_Fcalc ** 2)
    if (denom > 0) then
      result = sum(abs_Fobs * abs_Fcalc) / denom
    else
      result = 1
    end if
  
  end function calc_k_overall
  
  
  !---------------------------------------------------------------------------------------
  ! Calculate k_overall that minimizes least square residual between k*Fcalc and Fobs
  !
  !---------------------------------------------------------------------------------------
  function calc_k_overallc(Fobs, Fcalc) result(result)
    implicit none
    complex(real_kind), intent(in) :: Fobs(:) !< experimnal structure factors
    complex(real_kind), intent(in) :: Fcalc(:) !< model structure factors
    real(real_kind):: result
    real(real_kind):: denom
    
    ! Precondition
    ASSERT(size(Fobs) == size(Fcalc))
    
    denom = sum(Fcalc(:) * conjg(Fcalc(:)))
    if (denom > 0) then
       result= sum( real(Fobs(:)*conjg(Fcalc(:))) ) / denom
    else
      result = 1
    end if
  
  end function calc_k_overallc
  
  
  function calc_unscaled_r_factor(abs_Fobs, abs_Fcalc) result (result)
    real(real_kind), intent(in) :: abs_Fobs(:)            !< Magnitude of expreimental structure factors
    real(real_kind), intent(in) :: abs_Fcalc(size(abs_Fobs)) !< Magnitude of model structure factors
    real(real_kind) :: result
    real(real_kind), parameter :: bad_r_factor = 999
    real(real_kind) :: denum

    ! Precondition
    ASSERT(size(abs_Fobs) == size(abs_Fcalc))
    ASSERT(all(abs_Fobs >= 0))
    ASSERT(all(abs_Fcalc >= 0))
    
    denum = sum(abs_Fobs)
    
    if (denum == 0) then
      result = bad_r_factor
    else
      result = sum(abs(abs_Fobs - abs_Fcalc)) / denum
    end if
  end function calc_unscaled_r_factor
  
  
  function calc_r_factor(abs_Fobs, abs_Fcalc) result (result)
    real(real_kind), intent(in) :: abs_Fobs(:)            !< Magnitude of expreimental structure factors
    real(real_kind), intent(in) :: abs_Fcalc(size(abs_Fobs)) !< Magnitude of model structure factors
    real(real_kind) :: result
    real(real_kind) :: scale

    ! Precondition
    ASSERT(size(abs_Fobs) == size(abs_Fcalc))
    ASSERT(all(abs_Fobs >= 0))
    ASSERT(all(abs_Fcalc >= 0))
    
    scale = calc_k_overall(abs_Fobs, abs_Fcalc)
    result = calc_unscaled_r_factor(abs_Fobs, abs_Fcalc * scale)
  end function calc_r_factor
  
  
  function linspace(start, stop, num, endpoint) result(result)
    implicit none
    
    real(real_kind), intent(in) :: start, stop
    integer, intent(in) :: num
    logical, optional, intent(in) :: endpoint
    real(real_kind) :: result(num)
    
    real(real_kind) :: step
    logical :: use_endpoint
    integer ::i

    ! Precondition
    ASSERT(num > 0)
    
    use_endpoint = .FALSE.
    
    if (present(endpoint)) then
      use_endpoint = endpoint
    end if
    
    if (use_endpoint) then
      step = (stop - start) / (num - 1)
    else
      step = (stop - start) / num
    end if
    
    do i = 1, num
      result(i) = start + step * (i - 1)
    end do
  
  end function linspace
  
  
  function calc_k_bulk_cubic(absFobs, Fprot, Fbulk) result(result)
    ! Returns -1.0 in place of imaginary roots
    implicit none
    
    real(real_kind), intent(in) :: absFobs(:)
    complex(real_kind), intent(in) :: Fprot(size(absFobs))
    complex(real_kind), intent(in) :: Fbulk(size(absFobs))
    real(real_kind) :: result(3)
    integer :: j
    real(real_kind) :: p, r, q, t, I, u, v, w
    real(real_kind) :: a2, b2, c2, y2, a3, b3, c3, d3, y3, den
    real(real_kind) :: a, b, c
    
    ! Precondition
    ASSERT(size(absFobs) == size(Fprot))
    ASSERT(size(absFobs) == size(Fbulk))
    
    a2 = 0.0d0
    b2 = 0.0d0
    c2 = 0.0d0
    y2 = 0.0d0
    a3 = 0.0d0
    b3 = 0.0d0
    c3 = 0.0d0
    d3 = 0.0d0
    y3 = 0.0d0
    
    do j = 1, size(absFobs)
      p = real(Fprot(j))
      r = aimag(Fprot(j))
      q = real(Fbulk(j))
      t = aimag(Fbulk(j))
      I = absFobs(j) ** 2
      v = p * q + r * t
      w = q * q + t * t
      u = p * p + r * r
      a2 = a2 + u * I
      b2 = b2 + 2 * v * I
      c2 = c2 + w * I
      y2 = y2 + I * I
      a3 = a3 + u * v
      b3 = b3 + 2 * v * v + u * w
      c3 = c3 + 3 * v * w
      d3 = d3 + w * w
      y3 = y3 + v * I
    end do
    den = d3 * y2 - c2 * c2
    !! ASSERT(den != 0.0);
    !! coefficients of x**3 + ax**2 + bx + c = 0
    a = (c3 * y2 - c2 * b2 - c2 * y3) / den
    b = (b3 * y2 - c2 * a2 - y3 * b2) / den
    c = (a3 * y2 - y3 * a2) / den
    result = solve_cubic_equation_real(a, b, c, sentinel=-1.0_real_kind)
  end function calc_k_bulk_cubic
  
  
  pure function solve_cubic_equation_real(a, b, c, sentinel) result(solution)
    ! Analytical solution of x**3 + a x**2 + bx + c = 0.
    ! Returns `sentinel` in place of imaginary roots (default=0)
    implicit none
    real(real_kind), intent(in):: a, b, c
    real(real_kind), optional, intent(in) :: sentinel
    real(real_kind) :: sentinel_value
    real(real_kind) :: solution(3)
    real(real_kind) :: Disc, theta, S, T, p, q, arg, eps
    real(real_kind), parameter :: d_tolerance = 1.e-10
    real(real_kind), parameter :: pi = 3.141592653589793238_real_kind
    real(real_kind), parameter :: one_third = 1.0_real_kind / 3
    
    if (present(sentinel)) then
      sentinel_value = sentinel
    else
      sentinel_value = 0
    end if
    
    eps = d_tolerance * 10
    p = (3 * b - a ** 2) / 3
    q = (2 * a ** 3 - 9 * a * b + 27 * c) / 27
    Disc = (p / 3) ** 3 + (q / 2) ** 2
    if (abs(p)<eps .and. abs(q)<eps .and. abs(Disc)<eps) then
      solution(1:3) = - sign(abs(c) ** one_third, c)
    else if (Disc < 0) then
      arg = (q / p) ** 2 * 27 / (4 * p)
      arg = sqrt(abs(arg))
      if (abs(1 - abs(arg)) < 1.e-9) then
        arg = 1
      end if
      theta = acos(-sign(arg, q))
      solution(1) = 2 * (sqrt(abs(p / 3)) * cos(theta / 3)) - a / 3
      solution(2) = 2 * (sqrt(abs(p / 3)) * cos((theta + 2 * pi) / 3)) - a / 3
      solution(3) = 2 * (sqrt(abs(p / 3)) * cos((theta - 2 * pi) / 3)) - a / 3
    else
      S = sign(abs(-q / 2 + sqrt(Disc)) ** one_third, -q / 2 + sqrt(Disc))
      T = sign(abs(-q / 2 - sqrt(Disc)) ** one_third, -q / 2 - sqrt(Disc))
      solution(1) = S + T - a / 3
      if (Disc <= eps) then
        solution(2:3) = -(S + T) / 2 - a / 3
      else
        solution(2:3) = sentinel
      end if
    end if
  end function solve_cubic_equation_real
  
  
  pure subroutine linear_interpolation(x1, x2, y1, y2, k, b)
    ! Interploate y = k*x + b
    real(real_kind), intent(in) :: x1, x2, y1, y2
    real(real_kind), intent(out) :: k, b
    k = 0.0d0
    if (x1 /= x2) then
      k = (y2 - y1) / (x2 - x1)
    end if
    b = y1 - k * x1
  end subroutine linear_interpolation
  
  pure function smooth_with_moving_average(x) result(result)
    implicit none
    real(real_kind), intent(in):: x(:)
    real(real_kind) :: result(size(x))
    real(real_kind) :: x_(size(x) + 2), y(size(x) + 2)
    integer :: cycle, i, n_bins
    n_bins = size(x)
    x_(2:n_bins + 1) = x
    x_(1) = x(1)
    x_(n_bins + 2) = x(size(x))
    do cycle = 1, 5
      y = x_
      do i = 2, n_bins + 1
        if ((y(i - 1) < y(i) .and. y(i + 1) < y(i)) .or. &
            (y(i - 1) > y(i) .and. y(i + 1) > y(i))) then
          x_(i) = (y(i - 1) + y(i) + y(i + 1)) / 3.0d0
        end if
      end do
    end do
    result = x_(2:n_bins + 1)
  end function smooth_with_moving_average
 
  
  pure function smooth_k_bulk(k_bulk_in_bins, d_in_bins) result(result)
    implicit none
    real(real_kind), intent(in) :: k_bulk_in_bins(:)
    real(real_kind), intent(in) :: d_in_bins(size(k_bulk_in_bins))
    real(real_kind) :: tmp(size(k_bulk_in_bins))
    real(real_kind) :: result(size(k_bulk_in_bins))
    real(real_kind) :: k, d
    integer :: i
    
    tmp = smooth_with_moving_average(k_bulk_in_bins)
    
    result = 0
    do i = size(k_bulk_in_bins), 1, -1
      k = tmp(i)
      d = d_in_bins(i)
      ! TODO: Check suspicious floating point exact comparison: "r==0"
      !       Probably should be replaced with "r < eps"
      if (k == 0 .and. d < 3) then
        exit
      end if
      result(i) = k
    end do
  end function smooth_k_bulk
  
  
  !--------------------------------------------------------------------------------------------
  ! inverse: Matrix inverse, based on Doolittle LU factorization for Ax=b by Alex G.,
  !          December 2009.  The original matrix a(n,n) will be destroyed.
  !
  ! Arguments:
  !   a:       (n x n) array of coefficients for square matrix A
  ! Returns:
  !   c:       (n x n) inverse of matrix A
  !--------------------------------------------------------------------------------------------
  function inverse(z) result(c)
    
    implicit none
    real(real_kind), intent(in) :: z(:, :)
    real(real_kind), dimension(size(z, 1), size(z, 1)) :: c
    
    real(real_kind), dimension(size(z, 1), size(z, 1)) ::  L, U, a
    real(real_kind), dimension(size(z, 1)) :: b, d, x
    real(real_kind) coeff
    integer i, j, k, n
    
    ! assert size(a,1) == size(a,2)
    n = size(z, 1)
    a = z
    
    
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
  end function inverse

  ! cctbx analogue of symmetric matrix inversion
  pure function inverse_symmetric(z) result(c)

    implicit none
    real(real_kind), intent(in) :: z(:, :)
    real(real_kind), dimension(size(z, 1), size(z, 1)) ::  c
    real(real_kind), dimension(size(z, 1)) :: eigenvalues, diagonal_elements
    real(real_kind), dimension(size(z, 1) * (size(z, 1) + 1) / 2) :: a
    real(real_kind), dimension(size(z, 1) * size(z, 1)) :: eigenvectors
    real(real_kind), parameter :: epsilon_ = 1e-9
    integer i, j, k, n
    integer jac, ik, il, ilq, ilr, im, imq, imr, ind, iq, km, l, ll, lm, lq, m, mm, mq
    real(real_kind) anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y, denominator

    ! initialize lower triangle matrix and other necessary values
    n = size(z, 1)
    anorm = 0
    iq = 0
    k = 1
    eigenvectors = 0
    do i = 1, n
      do j = 1, i
        a(k) = z(i, j)
        k = k + 1
        if (i /= j) then
          anorm = anorm + z(i, j) ** 2
        else
          eigenvectors(i + (i - 1) * n) = 1
        end if
      end do
    end do

    anorm = sqrt(2 * anorm)
    anrmx = epsilon_ * anorm / n;
    if (anrmx < epsilon_) then
      anrmx = epsilon_
    end if

    if (anorm > 0) then
      thr = anorm
      do while (thr > anrmx)
        thr = thr / n
        ind = 1
        do while (ind > 0)
          ind = 0
          l = 0
          do while (l /= n - 1)
            lq = l * (l + 1) / 2
            ll = l + lq
            m = l + 1
            ilq = n * l
            do while (m /= n)
              mq = m * (m + 1) / 2
              lm = l + mq
              if (a(lm + 1) * a(lm + 1)>thr * thr) then
                ind = 1
                mm = m + mq
                x = 0.5_real_kind * (a(ll + 1) - a(mm + 1))
                denominator = sqrt(a(lm + 1) * a(lm + 1) + x * x)
                ! ASSERT(denominator /= 0);
                y = -a(lm + 1) / denominator
                if (x<0) then
                  y = -y
                end if
                sinx = y / sqrt(2 * (1 + (sqrt(1 - y * y))))
                sinx2 = sinx * sinx
                cosx = sqrt(1 - sinx2)
                cosx2 = cosx * cosx
                sincs = sinx * cosx
                imq = n * m
                do i = 0, n - 1
                  iq = i * (i + 1) / 2
                  if (i /= l .and. i /= m) then
                    if (i<m) then
                      im = i + mq
                    else
                      im = m + iq
                    end if
                    if (i<l) then
                      il = i + lq
                    else
                      il = l + iq
                    end if
                    x = a(il + 1) * cosx - a(im + 1) * sinx
                    a(im + 1) = a(il + 1) * sinx + a(im + 1) * cosx
                    a(il + 1) = x
                  end if
                  ilr = ilq + i
                  imr = imq + i
                  x = eigenvectors(ilr + 1) * cosx - eigenvectors(imr + 1) * sinx
                  eigenvectors(1 + imr) = eigenvectors(1 + ilr) * sinx + eigenvectors(1 + imr) * cosx
                  eigenvectors(1 + ilr) = x
                end do
                x = 2 * a(lm + 1) * sincs
                y = a(ll + 1) * cosx2 + a(mm + 1) * sinx2 - x
                x = a(ll + 1) * sinx2 + a(mm + 1) * cosx2 + x
                a(lm + 1) = (a(ll + 1) - a(mm + 1)) * sincs + a(lm + 1) * (cosx2 - sinx2)
                a(ll + 1) = y
                a(mm + 1) = x
              end if
              m = m + 1
            end do
            l = l + 1
          end do
        end do
      end do
    end if

    ! Sort eigenvalues & eigenvectors
    k = 0
    do i = 0, n - 2
      im = i
      km = k
      x = a(k + 1)
      l = 0
      do j = 0, n - 1
        if (j > i .and. a(l + 1) > x) then
          im = j
          km = l
          x = a(l + 1)
        end if
        l = l + j + 2
      end do
      if (im /= i) then
        a(km + 1) = a(k + 1)
        a(k + 1) = x
        l = n * i
        m = n * im
        do j = 0, n - 1
          x = eigenvectors(1 + l)
          eigenvectors(1 + l) = eigenvectors(1 + m);
          eigenvectors(1 + m) = x
          l = l + 1
          m = m + 1
        end do
      end if
      k = k + i + 2
    end do

    k = 0
    do j = 0, n - 1
      eigenvalues(j + 1) = a(k + 1);
      k = k + j + 2
    end do

    ! generalized_inverse_as_packed_u preparation
    do i = 0, n - 1
      x = eigenvalues(i + 1)
      if (abs(x) < anrmx) then
        diagonal_elements(i + 1) = 0
      else
        diagonal_elements(i + 1) = 1 / x;
      end if
    end do

    ! transpose_multiply_diagonal_multiply_as_packed_u analogue
    ik = 0
    do i = 0, n - 1
      x = eigenvectors(i + 1) * diagonal_elements(1)
      do k = i, n - 1
        a(ik + 1) = x * eigenvectors(k + 1)
        ik = ik + 1
      end do
    end do
    jac = n
    do j = 1, n - 1
      ik = 0
      do i = 0, n - 1
        x = eigenvectors(jac + i + 1) * diagonal_elements(j + 1)
        do k = i, n - 1
          a(ik + 1) = a(ik + 1) + x * eigenvectors(jac + k + 1)
          ik = ik + 1
        end do
      end do
      jac = jac + n
    end do

    ! populate the whole matrix
    k = 1
    do i = 1, n
      do j = i, n
        c(i, j) = a(k)
        c(j, i) = a(k)
        k = k + 1
      end do
    end do

  end function inverse_symmetric
  
  pure subroutine smooth_resolution_bins(values)
    implicit none
    integer :: i, num_bins
    real(real_kind) :: t_opt_1, t_opt_2, t_opt_3
    real(real_kind), dimension(:), intent(inout) :: values(:)
    
    num_bins = size(values)
    
    if (num_bins > 1) then
      t_opt_1 = values(1)
      t_opt_2 = values(2)
      do i = 2, num_bins - 1
        t_opt_3 = values(i + 1)
        values(i) = (t_opt_1+t_opt_2+t_opt_3)/3.0d0;
        t_opt_1 = t_opt_2;
        t_opt_2 = t_opt_3;
      end do
    end if
  
  end subroutine smooth_resolution_bins
  
  
  !------------------------------------------------------------------
  ! Calculate resolution for given miller indices and lattice vectors
  !------------------------------------------------------------------
  pure function calc_resolution(hkl, recip) result(result)
    implicit none
    double precision, intent(in) ::  recip(3, 3)          ! Reciprocal lattice vectors (3,3)
    integer, dimension(:,:), intent(in) :: hkl            ! Array of miller indices (3,N)
    real(real_kind), dimension(size(hkl, 2)) :: result    ! Output resolution array (N)
    integer :: i
    
    real(real_kind), dimension(3) :: vas, vbs, vcs
    real(real_kind) :: d_star
    
    vas = recip(:,1)
    vbs = recip(:,2)
    vcs = recip(:,3)
    
    do i = 1, size(result)
      d_star = ((hkl(1, i) * norm2(vas)) ** 2 + (hkl(2, i) * norm2(vbs)) ** 2 + &
          (hkl(3, i) * norm2(vcs))**2 + &
          2d0 * hkl(2, i) * hkl(3, i) * dot_product(vbs, vcs) + &
          2d0 * hkl(1, i) * hkl(3, i) * dot_product(vas, vcs) + &
          2d0 * hkl(1, i) * hkl(2, i) * dot_product(vbs, vas))
      result(i) = sqrt(1.0d0 / d_star)
    end do
  
  end function
  
  !----------------------------------------------------------------------------
  ! cross:   Compute the cross product of vectors a and b.
  !----------------------------------------------------------------------------
  pure function cross_product_real_4(a, b) result(result)
    
    real(4), intent(in) :: a(3), b(3)
    real(4) :: result(3)

    result(1) = a(2) * b(3) - a(3) * b(2)
    result(2) = a(3) * b(1) - a(1) * b(3)
    result(3) = a(1) * b(2) - a(2) * b(1)
  
  end function cross_product_real_4
  
  pure function cross_product_real_8(a, b) result(result)
    
    real(8), intent(in) :: a(3), b(3)
    real(8) :: result(3)
    
    result(1) = a(2) * b(3) - a(3) * b(2)
    result(2) = a(3) * b(1) - a(1) * b(3)
    result(3) = a(1) * b(2) - a(2) * b(1)
  
  end function cross_product_real_8
  
  subroutine set_start_size_from_bin_index(bin_index, n_bins, bin_start, bin_size)
    implicit none
    integer, intent(in) :: bin_index(:)
    integer, intent(in) :: n_bins
    integer, intent(out) :: bin_start(n_bins)
    integer, intent(out) :: bin_size(n_bins)
    integer :: last, i
    
    ! Preconditions:
    ASSERT(size(bin_index) > 0)
    ASSERT(is_sorted(bin_index))
    ASSERT(bin_index(1) == 1)
    ASSERT(bin_index(size(bin_index)) <= n_bins)
    
    last = 0
    do i = 1, size(bin_index)
      do while (bin_index(i) /= last)
        last = last + 1
        bin_start(last) = i
      end do
    end do
    
    bin_size(last) = size(bin_index) - bin_start(last) + 1
    bin_size(:last - 1) = bin_start(2:last) - bin_start(:last - 1)
    
    if (last /= n_bins) then
      bin_start(last + 1:n_bins) = size(bin_index)
      bin_size(last + 1:n_bins) = 0
    end if
  
  end subroutine set_start_size_from_bin_index
  
  pure function is_sorted_real_4(array) result(result)
    real(4), intent(in) :: array(:)
    logical :: result
    integer :: i
    
    do i = 1, size(array) - 1
      if (array(i) > array(i + 1)) then
        result = .FALSE.
        return
      end if
    end do
    result = .TRUE.
  end function is_sorted_real_4
  
  pure function is_sorted_real_8(array) result(result)
    real(8), intent(in) :: array(:)
    logical :: result
    integer :: i
    
    do i = 1, size(array) - 1
      if (array(i) > array(i + 1)) then
        result = .FALSE.
        return
      end if
    end do
    result = .TRUE.
  end function is_sorted_real_8
  
  pure function is_sorted_integer_4(array) result(result)
    integer(4), intent(in) :: array(:)
    logical :: result
    integer :: i
    
    do i = 1, size(array) - 1
      if (array(i) > array(i + 1)) then
        result = .FALSE.
        return
      end if
    end do
    result = .TRUE.
  end function is_sorted_integer_4
  
  pure function is_sorted_integer_8(array) result(result)
    integer(8), intent(in) :: array(:)
    logical :: result
    integer :: i
    
    do i = 1, size(array) - 1
      if (array(i) > array(i + 1)) then
        result = .FALSE.
        return
      end if
    end do
    result = .TRUE.
  end function is_sorted_integer_8
  
  elemental function to_radians_4(degrees) result(result)
    implicit none
    real(4), intent(in) :: degrees
    real(4), parameter :: radian_over_degree = PI_real_4 / 180
    real(4) :: result
    result = degrees * radian_over_degree
  end function to_radians_4
  
  elemental function to_radians_8(degrees) result(result)
    implicit none
    real(8), intent(in) :: degrees
    real(8), parameter :: radian_over_degree = PI / 180
    real(8) :: result
    result = degrees * radian_over_degree
  end function to_radians_8
  
  elemental function to_degrees_4(radians) result(result)
    implicit none
    real(4), intent(in) :: radians
    real(4), parameter :: degree_over_radian = 180 / PI_real_4
    real(4) :: result
    result = radians * degree_over_radian
  end function to_degrees_4
  
  elemental function to_degrees_8(radians) result(result)
    implicit none
    real(8), intent(in) :: radians
    real(8), parameter :: degree_over_radian = 180 / PI
    real(8) :: result
    result = radians * degree_over_radian
  end function to_degrees_8
  
end module xray_pure_utils
