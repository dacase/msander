subroutine nose_hoover_init_LES ( mass, v, f )
! initialize Nose-Hoover chain thermostats
  use les_data, only: cnum
  use constants, only: pi, hbar, kB
  use nose_hoover_module, only: nose_hoover_module_init,  &
                                Thermostat_init,  &
                                Thermostat_link,  &
                                Thermostat_type
  use nose_hoover_module, only: file_nhc, nchain, thermo, nthermo, tau
  
  implicit none
# include "../include/memory.h"
# include "../include/md.h"
# include "parallel.h"
  _REAL_, intent(in) :: mass( natom )
  _REAL_ :: v( 3, natom )
  _REAL_ :: f( 3, natom )
  _REAL_ :: kT, beta
  integer :: idim, iatom

  if ( worldrank.eq.0 ) then
     open( file_nhc, file = "NHC.dat" )
  endif

  kT = kB * temp0
  beta = 1.d0 / kT

  call nose_hoover_module_init
  nthermo = 3 * natom

  allocate(thermo(3,natom))

  do iatom = 1, natom
     if (cnum(iatom).eq.0) then
        tau = hbar * beta
     else
        tau = hbar*beta/sqrt(dble(nbead))
     endif
     do idim  = 1, 3
        call Thermostat_init(nchain,thermo(idim,iatom),1,0,kT,tau,.false.,.false.)
        call Thermostat_link(thermo(idim,iatom),mass(iatom), &
              v(idim,iatom),f(idim,iatom))
     enddo
  enddo

end subroutine nose_hoover_init_LES


