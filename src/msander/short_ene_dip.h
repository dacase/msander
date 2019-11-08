
!------------------------------------------------------------------------------
! short_ene_dip: calculate the direct Ewald component of the potentials 
!                with polarizabilities.
!
! Arguments:
!   See the descriptions in short_ene above for most formal arguments.  The 
!   new additions are:
!   mpoltype:  a code fro the type of polarizability to employ
!   pol:
!   pol2:
!   dipole:    
!------------------------------------------------------------------------------
subroutine short_ene_dip(i, xk, yk, zk, ipairs, numtot, numvdw, ewaldcof, &
                         eedtbdns, eed_cub, eed_lin, charge, dipole, ntypes, &
                         iac, ico, cn1, cn2, cn6, filter_cut, eelt, epol, &
                         evdw, ehb, frc, field, pol, pol2, dir_vir, ee_type, &
                         eedmeth, mpoltype, dxdr, eedvir, cn3, cn4, cn5)

#ifdef LES
  use les_data, only: lestmp, nlesty, lfac, lesfac, lestyp
#endif
  use nblist, only: bckptr,imagcrds,tranvec
  use constants, only: zero, one, two, three, four, five, &
                        six, twelve, third, half
  use pol_gauss

  implicit none
#include "../include/md.h"
  _REAL_ xk, yk, zk
  integer i, numvdw, numtot
  integer ipairs(*), ee_type, eedmeth, mpoltype
  _REAL_ ewaldcof, eed_cub(4,*), eed_lin(2,*), charge(*), dipole(3,*), &
         dir_vir(3,3), eedvir
  _REAL_ eedtbdns, filter_cut, dxdr
  integer ntypes,iac(*),ico(*)
  _REAL_ cn1(*), cn2(*), cn6(*), eelt, epol, evdw, ehb, frc(3,*), &
         field(3,*), pol(*)
  _REAL_ pol2(*)
  _REAL_ cn3(*), cn4(*), cn5(*)
  integer ic,j,m,n,ind,iaci
  _REAL_ del
  _REAL_ switch, d_switch_dx
  _REAL_ ee_vir_iso
  _REAL_ edx, edy, edz
  _REAL_ b0, b1, b2, b3, fac, dotir, dotjr, dotij, fact
  _REAL_ dphii_dx, dphii_dy, dphii_dz, dphij_dx, dphij_dy, dphij_dz
  _REAL_ dphii_dx_cor, dphii_dy_cor, dphii_dz_cor
  _REAL_ dphij_dx_cor, dphij_dy_cor, dphij_dz_cor
  _REAL_ term, term0, term1, termi, termj, cgj
  _REAL_ filter_cut2, xx
  _REAL_ xktran(3,18)
  integer, parameter :: mask27 = 2**27 - 1
  _REAL_ mr4, f4
  _REAL_ delx, dely, delz, delr, delr2, f10, r10, cgi, delr2inv, r6, f6, &
         f12, df, dx, x, dfx, dfy, dfz, dumx, dumy, dumz
  integer itran
  _REAL_ au3, exp_au3, lambda3, lambda5, lambda7, b1_o
  _REAL_ au, a2u2, a3u3, a4u4, exp_au, v, v3, v4, termi_o, termj_o
   
  fac = two * ewaldcof * ewaldcof
  ee_vir_iso = zero
  del = one / eedtbdns
  dumx = zero
  dumy = zero
  dumz = zero
  edx = zero
  edy = zero
  edz = zero
  f4 = zero
  filter_cut2 = filter_cut * filter_cut
  cgi = charge(i)
  iaci = ntypes * (iac(i) - 1)
  do m = 1, 18
    xktran(1,m) = tranvec(1,m) - xk
    xktran(2,m) = tranvec(2,m) - yk
    xktran(3,m) = tranvec(3,m) - zk
  end do
#ifdef LES
  lestmp=nlesty*(lestyp(i)-1)
#endif
  do m = 1, numvdw
    n = ipairs(m)
    itran = ishft(n, -27)
    n = iand(n, mask27)
    j = bckptr(n)
    delx = imagcrds(1,n) + xktran(1,itran)
    dely = imagcrds(2,n) + xktran(2,itran)
    delz = imagcrds(3,n) + xktran(3,itran)
    delr2 = delx*delx + dely*dely + delz*delz
    if (delr2 < filter_cut2) then
      delr = sqrt(delr2)
      delr2inv = one/delr2
      x = dxdr * delr
      cgj = charge(j)
      if (eedmeth == 1) then

        ! Cubic spline on switch
        ind = int(eedtbdns*x) + 1
        dx = x - (ind-one)*del
        switch = eed_cub(1,ind) + dx*(eed_cub(2,ind) + &
                                      dx*(eed_cub(3,ind) + &
                                          dx*eed_cub(4,ind)*third)*half)
        d_switch_dx = eed_cub(2,ind) + dx*(eed_cub(3,ind) + &
                                           dx*eed_cub(4,ind)*half)
      else if (eedmeth == 2) then

        ! Linear lookup on switch, deriv
        xx = eedtbdns*x + 1
        ind = int(xx)
        dx = xx - ind
        switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
        d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      else if (eedmeth == 3) then
         
        ! Direct function call:
        call get_ee_func(x, switch, d_switch_dx, ee_type)
      else if ( eedmeth == 4 ) then
          
        ! Use un-modified Coulomb interaction, no switch
        switch = one
        d_switch_dx = zero
      else
        write(6,*) 'bad eedmeth in ew_short_dip: ',eedmeth
        call mexit(6, 1)
      end if
      ! End switch over eedmeth electrostatics options
 
      ! Calculate Thole damping function
      if (pol(i) == 0 .or. pol(j) == 0) then
        lambda3 = one
        lambda5 = one
        lambda7 = one
      else if (mpoltype == 1) then
        lambda3 = one
        lambda5 = one
        lambda7 = one
      else if (mpoltype == 2) then
        au3 = delr*delr2 / (pol2(i)*pol2(j))
        exp_au3 = exp(-au3)
        lambda3 = one - exp_au3
        lambda5 = one - (one+au3)*exp_au3
        lambda7 = one - (one + au3 + three/five*au3*au3)*exp_au3
      else if (mpoltype == 3) then
        au = delr / (pol2(i)*pol2(j))
        exp_au = exp(-au)
        a2u2 = au*au
        a3u3 = a2u2*au
        a4u4 = a3u3*au
        lambda3=one-(a2u2/two+au+one)*exp_au
        lambda5=lambda3-a3u3/six*exp_au
        lambda7=lambda5-a4u4/six/five*exp_au
      else if (mpoltype == 4) then
        v = delr / (pol2(i)*pol2(j))
        if ( v < one ) then
          v3 = v*v*v
          v4 = v3*v
          lambda3 = four*v3 - three*v4
          lambda5 = v4
          lambda7 = v4 / five
        else
          lambda3 = one
          lambda5 = one
          lambda7 = one
        end if
      end if
      ! End switch over polarization styles

      ! Tom Darden Got the idea for B_l from Walter Smith's CCP5 article 1982
      ! Ewald for point multipoles
      ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
      ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)
      b0 = switch*delr*delr2inv
      fact = d_switch_dx*dxdr
      b1 = (b0 - fact)*delr2inv
      fact = fac*fact
      b2 = (three*b1 - fact)*delr2inv
      fact = fac*fact
      b3 = (five*b2 - fact)*delr2inv
      b1_o = b1
      b1 = b1 * lambda3
      b2 = b2 * lambda5
      b3 = b3 * lambda7
         
      ! B1 = (B0 - d_switch_dx*dxdr)*delr2inv
      ! B2 = (Three*B1 - fac*ewaldcof*d_switch_dx)*delr2inv
      ! B3 = (Five*B2 - fac*fac*ewaldcof*d_switch_dx)*delr2inv
       
      ! epol = dip_i dot grad_i of phii
      ! phii is direct sum electrostatic potential at i due to j
      ! so phii = cgj*B0 + dipj dot gradj of B0 = cgj*B0 - dotjr*B1
      ! dphii_dx etc are derivatives with respect to r_i
      ! phij is direct sum electrostatic potential at j due to i
      ! dphij_dx etc are derivatives with respect to r_j
      dotjr = dipole(1,j)*delx + dipole(2,j)*dely + dipole(3,j)*delz
      dotir = dipole(1,i)*delx + dipole(2,i)*dely + dipole(3,i)*delz
      dotij = dipole(1,i)*dipole(1,j) + dipole(2,i)*dipole(2,j) + &
              dipole(3,i)*dipole(3,j)
         
      ! gradi phii = cgj*rij*B1 + dipj*B1 - dotjr*rij*B2
      ! so epol = -cgi*dotjr*B1 + (cgj*B1 - dotjr*B2)*dotir + dotij*B1
      eelt = eelt + cgi*cgj*b0
      term = cgj*dotir - cgi*dotjr + dotij
      epol = epol + term*b1 - dotir*dotjr*b2
      term0 = cgi*cgj*b0 + term*b1 - dotir*dotjr*b2
         
      ! so ene = ene + term0; dfx = dterm0_dx etc
      ! grad term0 = term1*rij + B1*grad_i term - grad_i dotir*dotjr*B2
      ! grad_i term = -cgj*dip_i + cgi*dip_j
      ! grad_i dotir = -dip_i; similar for dotjr
      ! grad_i term0 = term1*rij + (-cgj*B1+dotjr*B2)*dip_i +
      !                (cgi*B1+dotir*B2)*dip_j
      term1 = cgi*cgj*b1_o + term*b2 - dotir*dotjr*b3
      termi = cgi*b1+dotir*b2
      termj = cgj*b1-dotjr*b2
      termi_o = cgi*b1_o+dotir*b2
      termj_o = cgj*b1_o-dotjr*b2
      dfx = term1*delx + termi*dipole(1,j) - termj*dipole(1,i)
      dfy = term1*dely + termi*dipole(2,j) - termj*dipole(2,i)
      dfz = term1*delz + termi*dipole(3,j) - termj*dipole(3,i)
      ic = ico(iaci+iac(j))
      r6 = delr2inv * delr2inv * delr2inv
#ifdef LES
      lfac = lesfac(lestmp+lestyp(j))
      f6 = cn2(ic) * r6 * lfac
      f12 = cn1(ic) * (r6*r6) * lfac
#else
      if (vdwmodel == 0) then
        f6 = cn2(ic)*r6
        f12 = cn1(ic)*(r6*r6)
      else
        f6 = cn5(ic)*r6
        f12 = cn4(ic)/exp(cn3(ic)/sqrt(delr2inv))
      endif
#endif
      if (lj1264 == 1) then
        mr4 = delr2inv*delr2inv
        f4 = cn6(ic)*mr4
        evdw = evdw + f12 - f6 - f4
      else
        evdw = evdw + f12 - f6
      end if
         
      ! Force related quantities
      df = (twelve*f12 - six*f6)*delr2inv
      dfx = dfx + df*delx
      dfy = dfy + df*dely
      dfz = dfz + df*delz

      ! Inserting Gaussian correction to energy
      if (ipolg == 1) then
        call Calc_Dip_Gauss_Correc(i, j, delx, dely, delz, delr2inv, delr, &
                                   delr2, eelt, epol, dfx, dfy, dfz, dotjr, &
                                   dotir, dotij, dipole(1,j), dipole(2,j), &
                                   dipole(3,j), dipole(1,i), dipole(2,i), &
                                   dipole(3,i), dphii_dx_cor, dphii_dy_cor, &
                                   dphii_dz_cor, dphij_dx_cor, dphij_dy_cor, &
                                   dphij_dz_cor, eed_cub)
      end if

      frc(1,j) = frc(1,j) + dfx
      frc(2,j) = frc(2,j) + dfy
      frc(3,j) = frc(3,j) + dfz
      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz
         
      ! Field related quantities
      dphii_dx = termj_o*delx + b1*dipole(1,j)
      dphii_dy = termj_o*dely + b1*dipole(2,j)
      dphii_dz = termj_o*delz + b1*dipole(3,j)
      dphij_dx = -termi_o*delx + b1*dipole(1,i)
      dphij_dy = -termi_o*dely + b1*dipole(2,i)
      dphij_dz = -termi_o*delz + b1*dipole(3,i)

      ! Gaussian Correction to field
      if (ipolg == 1) then
        dphii_dx = dphii_dx + dphii_dx_cor
        dphii_dy = dphii_dy + dphii_dy_cor
        dphii_dz = dphii_dz + dphii_dz_cor
        dphij_dx = dphij_dx + dphij_dx_cor
        dphij_dy = dphij_dy + dphij_dy_cor
        dphij_dz = dphij_dz + dphij_dz_cor
      end if
      edx = edx + dphii_dx
      edy = edy + dphii_dy
      edz = edz + dphii_dz
      field(1,j) = field(1,j) - dphij_dx
      field(2,j) = field(2,j) - dphij_dy
      field(3,j) = field(3,j) - dphij_dz
    end if
    ! End branch for range check
  end do
  ! End loop over interacting particles
   
  do m = numvdw+1, numtot
    n = ipairs(m)
    itran = ishft(n, -27)
    n = iand(n, mask27)
    j = bckptr(n)
    delx = imagcrds(1,n) + xktran(1,itran)
    dely = imagcrds(2,n) + xktran(2,itran)
    delz = imagcrds(3,n) + xktran(3,itran)
    delr2 = delx*delx + dely*dely+delz*delz
    if (delr2 < filter_cut2) then
      delr = sqrt(delr2)
      delr2inv = one / delr2
      x = dxdr * delr
      cgj = charge(j)
      if (eedmeth == 1) then

        ! Cubic spline on switch
        ind = int(eedtbdns*x) + 1
        dx = x - (ind-one)*del
        switch = eed_cub(1,ind) + dx*(eed_cub(2,ind) + &
                                      dx*(eed_cub(3,ind) + &
                                          dx*eed_cub(4,ind)*third)*half)
        d_switch_dx = eed_cub(2,ind) + dx*(eed_cub(3,ind) + &
                                           dx*eed_cub(4,ind)*half)
      else if (eedmeth == 2) then
            
        ! Linear lookup on switch, deriv
        xx = eedtbdns*x + 1
        ind = int(xx)
        dx = xx - ind
        switch = (one - dx)*eed_lin(1,ind) + dx*eed_lin(1,ind+1)
        d_switch_dx = (one - dx)*eed_lin(2,ind) + dx*eed_lin(2,ind+1)
      else if ( eedmeth == 3 )then
          
        ! Direct function call:
        call get_ee_func(x, switch, d_switch_dx, ee_type)
      else if (eedmeth == 4) then

        ! Use un-modified Coulomb interaction, no switch
        switch = one
        d_switch_dx = zero
      else
        write(6,*) 'bad eedmeth in ew_short_dip: ', eedmeth
        call mexit(6, 1)
      end if  ! ( eedmeth == 1 )
      if (pol(i) == 0 .or. pol(j) == 0) then
        lambda3 = one
        lambda5 = one
        lambda7 = one
      else if (mpoltype == 1) then
        lambda3 = one
        lambda5 = one
        lambda7 = one
      else if (mpoltype == 2) then
        au3 = delr*delr2 / (pol2(i)*pol2(j))
        exp_au3 = exp(-au3)
        lambda3 = one - exp_au3
        lambda5 = one - (one+au3)*exp_au3
        lambda7 = one - (one + au3 + three/five*au3*au3)*exp_au3
      else if ( mpoltype == 3 ) then
        au = delr / (pol2(i)*pol2(j))
        exp_au = exp(-au)
        a2u2 = au*au
        a3u3 = a2u2*au
        a4u4 = a3u3*au
        lambda3 = one - (a2u2/two + au + one)*exp_au
        lambda5 = lambda3 - a3u3/six*exp_au
        lambda7 = lambda5 - a4u4/six/five*exp_au
      else if ( mpoltype == 4 ) then
        v = delr / (pol2(i)*pol2(j))
        if ( v < one ) then
          v3 = v*v*v
          v4 = v3*v
          lambda3 = four*v3 - three*v4
          lambda5 = v4
          lambda7 = v4 / five
        else
          lambda3 = one
          lambda5 = one
          lambda7 = one
        end if
      end if
      b0 = switch * delr * delr2inv
      fact = d_switch_dx * dxdr
      b1 = (b0 - fact) * delr2inv
      fact = fac * fact
      b2 = (three*b1 - fact)*delr2inv
      fact = fac * fact
      b3 = (five*b2 - fact)*delr2inv
      b1_o = b1
      b1 = b1 * lambda3
      b2 = b2 * lambda5
      b3 = b3 * lambda7

      ! B0 = switch*delr*delr2inv
      ! B1 = (B0 - d_switch_dx*dxdr)*delr2inv
      ! B2 = (Three*B1 - fac*ewaldcof*d_switch_dx)*delr2inv
      ! B3 = (Five*B2 - fac*fac*ewaldcof*d_switch_dx)*delr2inv
      dotjr = dipole(1,j)*delx + dipole(2,j)*dely + dipole(3,j)*delz
      dotir = dipole(1,i)*delx + dipole(2,i)*dely + dipole(3,i)*delz
      dotij = dipole(1,i)*dipole(1,j) + dipole(2,i)*dipole(2,j) + &
              dipole(3,i)*dipole(3,j)
      eelt = eelt + cgi*cgj*b0
      term = cgj*dotir - cgi*dotjr + dotij
      epol = epol + term*b1 - dotir*dotjr*b2
      term0 = cgi*cgj*b0 + term*b1 - dotir*dotjr*b2
      term1 = cgi*cgj*b1_o + term*b2 - dotir*dotjr*b3
      termi = cgi*b1 + dotir*b2
      termj = cgj*b1 - dotjr*b2
      termi_o = cgi*b1_o + dotir*b2
      termj_o = cgj*b1_o - dotjr*b2
      dfx = (term1)*delx + termi*dipole(1,j) - termj*dipole(1,i)
      dfy = (term1)*dely + termi*dipole(2,j) - termj*dipole(2,i)
      dfz = (term1)*delz + termi*dipole(3,j) - termj*dipole(3,i)
      df = zero
      ehb = zero
      f10 = zero
      r10 = zero
      ! Force related quantities
      dfx = dfx + df*delx
      dfy = dfy + df*dely
      dfz = dfz + df*delz
      frc(1,j) = frc(1,j) + dfx
      frc(2,j) = frc(2,j) + dfy
      frc(3,j) = frc(3,j) + dfz
      dumx = dumx + dfx
      dumy = dumy + dfy
      dumz = dumz + dfz
        
      ! Field related quantities
      dphii_dx = termj_o*delx + b1*dipole(1,j)
      dphii_dy = termj_o*dely + b1*dipole(2,j)
      dphii_dz = termj_o*delz + b1*dipole(3,j)
      dphij_dx = -termi_o*delx + b1*dipole(1,i)
      dphij_dy = -termi_o*dely + b1*dipole(2,i)
      dphij_dz = -termi_o*delz + b1*dipole(3,i)
      edx = edx + dphii_dx
      edy = edy + dphii_dy
      edz = edz + dphii_dz
      field(1,j) = field(1,j) - dphij_dx
      field(2,j) = field(2,j) - dphij_dy
      field(3,j) = field(3,j) - dphij_dz
    end if  ! ( delr2 < filter_cut2 )
  end do  !  m = numvdw+1,numtot
  frc(1,i) = frc(1,i) - dumx
  frc(2,i) = frc(2,i) - dumy
  frc(3,i) = frc(3,i) - dumz
  field(1,i) = field(1,i) - edx
  field(2,i) = field(2,i) - edy
  field(3,i) = field(3,i) - edz
  return

end subroutine short_ene_dip 

