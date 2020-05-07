!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the
! kinetic energy, which can be used to rescale the velocities.  The
! procedure can be applied to all atoms or to smaller groups.  If it is
! applied to intersecting groups in sequence, the kinetic energy that is
! given as an input (kk) has to be up-to-date with respect to the previous
! rescalings.
!
! When applied to the entire system, and when performing standard
! molecular dynamics (fixed c.o.m. (center of mass)) the degrees of
! freedom of the c.o.m. have to be discarded in the calculation of
! ndeg, and the c.o.m. momentum HAS TO BE SET TO ZERO.  When applied to
! subgroups, one can chose to: (a) calculate the subgroup kinetic energy
! in the usual reference frame, and count the c.o.m. in ndeg (b) calculate
! the subgroup kinetic energy with respect to its c.o.m. motion, discard
! the c.o.m. in ndeg and apply the rescale factor with respect to the
! subgroup c.o.m. velocity.  They should be almost equivalent.  If the
! subgroups are expected to move one respect to the other, the choice (b)
! should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces
! to an istantaneous randomization of the kinetic energy, as described in
! paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator
! against discretization errors.  The easiest recipe is: htilde = h +
! conint where h is the total energy (kinetic + potential) and conint is
! a quantity accumulated along the trajectory as minus the sum of all the
! increments of kinetic energy due to the thermostat.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/assert.fh"

module resamplekin_mod

contains
  function resamplekin(kk, sigma, ndeg, taut)
    implicit none
    real(kind = 8) :: resamplekin
    real(kind = 8), intent(in) :: kk    ! present value of the kinetic energy 
                                        ! of the atoms to be thermalized 
                                        ! (in arbitrary units)
    real(kind = 8), intent(in) :: sigma ! target average value of the kinetic 
                                        ! energy (ndeg k_b T/2)  
                                        ! (in the same units as kk)
    integer(kind = 4), intent(in) :: ndeg  ! number of degrees of freedom 
                                           ! of the atoms to be thermalized
    real(kind = 8), intent(in) :: taut  ! relaxation time of the thermostat, 
                                ! in units of 'how often this routine is called'
    real(kind = 8) :: factor, rr
    if(taut>0.1) then
      factor = exp(-1.0 / taut)
    else
      factor = 0.0
    end if
    rr = gasdev()
    resamplekin = kk + (1.0 - factor) * (sigma * (sumnoises(ndeg - 1) + rr**2) / ndeg - kk) &
        + 2.0 * rr * sqrt(kk * sigma / ndeg * (1.0 - factor) * factor)
  
  contains
    
    double precision function sumnoises(nn)
      implicit none
      integer, intent(in) :: nn
      ! returns the sum of n independent gaussian noises squared
      ! (i.e. equivalent to summing the square of the return values 
      !  of nn calls to gasdev)
      if(nn==0) then
        sumnoises = 0.0
      else if(nn==1) then
        sumnoises = gasdev()**2
      else if(modulo(nn, 2)==0) then
        sumnoises = 2.0 * gamdev(nn / 2)
      else
        sumnoises = 2.0 * gamdev((nn - 1) / 2) + gasdev()**2
      end if
    end function sumnoises
  
  end function resamplekin
  
  double precision function gamdev(ia)
    ! gamma-distributed random number
    
    implicit none
    integer, intent(in) :: ia
    integer j
    real(kind = 8) am, e, s, v1, v2, x, y

    ASSERT( ia > 0 )
    if( ia<6 ) then
      x = 1.d0
      do j = 1, ia
        x = x * ran1()
      end do
      x = -log(x)
    else
      e = 9999.d0
      do while( ran1() < e ) 
        v1 = 2.d0 * ran1() - 1.d0
        v2 = 2.d0 * ran1() - 1.d0
        if(v1**2 + v2**2 > 1.d0) cycle
        y = v2 / v1
        am = ia - 1
        s = sqrt(2.d0 * am + 1.d0)
        x = s * y + am
        if(x <= 0.d0) cycle
        e = (1.d0 + y**2) * exp(am * log(x / am) - s * y)
      end do
    endif
    gamdev = x
  end
  
  double precision function gasdev()
    use random, only : gauss
    call gauss(0.0d0, 1.0d0, gasdev)
  end function gasdev
  
  double precision function ran1()
    use random, only : amrand
    call amrand(ran1)
  end function ran1

end module resamplekin_mod
