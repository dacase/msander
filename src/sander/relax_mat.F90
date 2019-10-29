! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

module relax_mat

public noeread, noecalc
private caldis, calrate, dinten, drates, indexn, remarc

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Get distance-dependent terms for relaxation matrix analysis
#ifdef NMODE
subroutine caldis(x,ddep,dddep,newf,amass)
#else
subroutine caldis(x,ddep,dddep)
#endif
   
   
   !  Subroutine CALculate DIStances:
   
   !     Calculates distance information bewteen all hydrogen atoms:
   !           nath = number of H atoms
   !           x    = array of atom coordinates
   !           ddep  = distance dependent part of rate matrix element
   !           dddep = derivative of ddep with respect to proton m,i
   
   !      Treats methyl rotors as three equivalent protons, and
   !      aromatic D and E as equivalent.  The distances are
   !      calculated assuming fast motion for the methyl group
   !      and slow (compared to correlation time) motion for aromatic flip,
   !      which is the same as 1/6 averaging.
   
   !      Variables "i" and "j" are indexed in the nath scheme; "ii"
   !      and "jj" are the corresponding entries in the natmet scheme.
   use constants, only : zero
   implicit none
#  include "nmr.h"

!Passed in
   _REAL_ :: x(*),ddep(ma,ma),dddep(3,ma,ma)
#ifdef NMODE
   _REAL_ :: amass(*)
   logical newf
#endif

!Local
   _REAL_ :: dddepp(3)
   _REAL_ :: dmet(3,3),rmet(3,3,3)
   _REAL_ :: sum2, sumeff, dprod, denom, term, diffij, work
   integer :: i, jj, ii, j, m2ij, kl, k, l, klmn, m, n, nn, lk, lkm

   ! --- constants:
   
   _REAL_, parameter :: x5o81=(5.0d0/81.0d0)
   _REAL_, parameter :: x2o81=(2.0d0/81.0d0)
   _REAL_, parameter :: x6o81=(6.0d0/81.0d0)
   _REAL_, parameter :: x5o9=(5.0d0/9.0d0)
   _REAL_, parameter :: x2o9=(2.0d0/9.0d0)
   _REAL_, parameter :: x6o9=(6.0d0/9.0d0)
   _REAL_, parameter :: x5o18=(5.0d0/18.0d0)
   _REAL_, parameter :: x1o9=(1.0d0/9.0d0)
   _REAL_, parameter :: x1o3=(1.0d0/3.0d0)
   _REAL_, parameter :: x1o18=(1.0d0/18.0d0)
   _REAL_, parameter :: x5o36=(5.0d0/36.0d0)
   ! --- subscript arrays to collapse mutliple do-loops into one:
   integer, dimension(9), parameter :: ksub = (/1,2,3,1,2,3,1,2,3/)
   integer, dimension(9), parameter :: lsub = (/1,1,1,2,2,2,3,3,3/)
   integer, dimension(81), parameter ::  k1sub = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &
                                                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, &
                                                   3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3/)
   integer, dimension(81), parameter ::  l1sub = (/1,1,1,1,1,1,1,1,1, &
                                                   2,2,2,2,2,2,2,2,2, &
                                                   3,3,3,3,3,3,3,3,3, &
                                                   1,1,1,1,1,1,1,1,1, &
                                                   2,2,2,2,2,2,2,2,2, &
                                                   3,3,3,3,3,3,3,3,3, &
                                                   1,1,1,1,1,1,1,1,1, &
                                                   2,2,2,2,2,2,2,2,2, &
                                                   3,3,3,3,3,3,3,3,3/)
   integer, dimension(81), parameter ::  m1sub = (/1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3, &
                                                   1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3, &
                                                   1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3/)
   integer, dimension(81), parameter ::  n1sub = (/1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3, &
                                                   1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3, &
                                                   1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3, &
                                                   1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3, &
                                                   1,2,3,1,2,3,1,2,3/)
   integer, dimension(6), parameter :: l2sub = (/1,1,1,2,2,2/)
   integer, dimension(6), parameter :: k2sub = (/1,2,3,1,2,3/)
   integer, dimension(4), parameter :: k3sub = (/1,1,2,2/)
   integer, dimension(4), parameter :: l3sub = (/1,2,1,2/)
   integer, dimension(18), parameter :: l4sub = (/1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2/)
   integer, dimension(18), parameter :: k4sub = (/1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3/)
   integer, dimension(18), parameter ::  m4sub = (/1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3/)
   
   do i=1,nath
      do jj=1,natmet
         dddep(1,i,jj) = zero
         dddep(2,i,jj) = zero
         dddep(3,i,jj) = zero
      end do
   end do
   ii = 0
   do i=1,nath
      
      !                **first 2 methyl protons of a methyl group
      !                   and first aromatic D or E skipped
      
      if(m2(i) == 1.or.m2(i) == 4) cycle
      ii = ii + 1
      
      jj = 0
      do j=1,nath
         
         !                  **first 2 methyl protons of a methyl group
         !                    and first aromatic D or E skipped
         
         if(m2(j) == 1.or.m2(j) == 4) cycle
         jj = jj + 1
         
         m2ij=m2(i)*m2(j)
         if(m2ij == 9) then !             **methyl - methyl:
            
            if(ii == jj) then !           **intramethyl:

               ddep(ii,jj) = 8.13d-3  ! this is 1/(4*1.77**6)
               
            else !                         **intermethyl:
               
               do kl=1,9
                  k = ksub(kl)
                  l = lsub(kl)
                  rmet(k,l,1)=x(3*(i-k)+1)-x(3*(j-l)+1)
                  sum2=rmet(k,l,1)**2
                  rmet(k,l,2)=x(3*(i-k)+2)-x(3*(j-l)+2)
                  sum2=sum2+rmet(k,l,2)**2
                  rmet(k,l,3)=x(3*(i-k)+3)-x(3*(j-l)+3)
                  sum2=sum2+rmet(k,l,3)**2
                  dmet(k,l)=sum2
               end do
               sumeff=0.0d0
               do klmn=1,81
                  k = k1sub(klmn)
                  l = l1sub(klmn)
                  m = m1sub(klmn)
                  n = n1sub(klmn)
                  dprod=(rmet(k,l,1)*rmet(m,n,1)) &
                        +(rmet(k,l,2)*rmet(m,n,2)) &
                        +(rmet(k,l,3)*rmet(m,n,3))
                  denom = 1.0d0/sqrt((dmet(k,l)*dmet(m,n))**5)
                  term=(3.0d0*dprod**2-dmet(k,l)*dmet(m,n))*denom
                  sumeff=sumeff + term
                  do nn=1,3
                     dddep(nn,i-k+1,jj) = dddep(nn,i-k+1,jj) &
                           - x5o81*term*rmet(k,l,nn)/dmet(k,l) &
                           - (x2o81*dmet(m,n)*rmet(k,l,nn) - &
                           x6o81*dprod*rmet(m,n,nn))*denom
                  end do
               end do
               ddep(ii,jj)=sumeff/162.0d0
            end if
         
         else if (m2ij == 6) then !     **methyl to ordinary:
         
            if(m2(i) == 3) then !       ** i is the methyl:
               do k=1,3
                  rmet(k,1,1)=x(3*(i-k)+1)-x(3*(j-1)+1)
                  sum2=rmet(k,1,1)**2
                  rmet(k,1,2)=x(3*(i-k)+2)-x(3*(j-1)+2)
                  sum2=sum2+rmet(k,1,2)**2
                  rmet(k,1,3)=x(3*(i-k)+3)-x(3*(j-1)+3)
                  sum2=sum2+rmet(k,1,3)**2
                  dmet(k,1)=sum2
               end do
               sumeff=0.0d0
               do kl=1,9
                  k = ksub(kl)
                  l = lsub(kl)
                  dprod=rmet(k,1,1)*rmet(l,1,1) &
                        +rmet(k,1,2)*rmet(l,1,2) &
                        +rmet(k,1,3)*rmet(l,1,3)
                  denom = 1.0d0/sqrt((dmet(k,1)*dmet(l,1))**5)
                  term=(3.0d0*dprod**2-dmet(k,1)*dmet(l,1))*denom
                  sumeff = sumeff + term
                  do n=1,3
                     dddep(n,i-k+1,jj) = dddep(n,i-k+1,jj) &
                           - x5o9*term*rmet(k,1,n)/dmet(k,1) &
                           - (x2o9*dmet(l,1)*rmet(k,1,n) - &
                           x6o9*dprod*rmet(l,1,n))*denom
                  end do
               end do
               ddep(ii,jj)=sumeff/18.0d0
         
            else !                     ** j is the methyl:
               do k=1,3
                  rmet(k,1,1)=x(3*(i-1)+1) - x(3*(j-k)+1)
                  sum2=rmet(k,1,1)**2
                  rmet(k,1,2)=x(3*(i-1)+2) - x(3*(j-k)+2)
                  sum2=sum2+rmet(k,1,2)**2
                  rmet(k,1,3)=x(3*(i-1)+3) - x(3*(j-k)+3)
                  sum2=sum2+rmet(k,1,3)**2
                  dmet(k,1)=sum2
               end do
               sumeff=0.0d0
               do kl=1,9
                  k = ksub(kl)
                  l = lsub(kl)
                  dprod=rmet(k,1,1)*rmet(l,1,1) &
                        +rmet(k,1,2)*rmet(l,1,2) &
                        +rmet(k,1,3)*rmet(l,1,3)
                  denom = 1.0d0/sqrt((dmet(k,1)*dmet(l,1))**5)
                  term=(3.0d0*dprod**2-dmet(k,1)*dmet(l,1))*denom
                  sumeff = sumeff + term
                  do n=1,3
                     dddep(n,i,jj) = dddep(n,i,jj) &
                        - x5o18*term*(rmet(k,1,n)/dmet(k,1) &
                        + rmet(l,1,n)/dmet(l,1)) &
                        - ( dmet(l,1)*rmet(k,1,n) &
                        +dmet(k,1)*rmet(l,1,n) &
                        -3.0d0*dprod*(rmet(k,1,n)+rmet(l,1,n)))*denom*x1o9
                  end do
               end do
               ddep(ii,jj)=sumeff/18.0d0
      
            end if
   
         else if (m2ij == 15) then !        *** aromatic to methyl:
   
            if(m2(i) == 3) then !           **i is methyl, j is aromatic:
               do lk=1,6
                  l = l2sub(lk)
                  k = k2sub(lk)
                  rmet(k,l,1)=x(3*(i-k)+1)-x(3*(j-l)+1)
                  sum2=rmet(k,l,1)**2
                  rmet(k,l,2)=x(3*(i-k)+2)-x(3*(j-l)+2)
                  sum2=sum2+rmet(k,l,2)**2
                  rmet(k,l,3)=x(3*(i-k)+3)-x(3*(j-l)+3)
                  sum2=sum2+rmet(k,l,3)**2
                  dmet(k,l)=sum2
               end do
               sumeff=0.0d0
               do lkm=1,18
                  l = l4sub(lkm)
                  k = k4sub(lkm)
                  m = m4sub(lkm)
                  dprod=(rmet(k,l,1)*rmet(m,l,1)) &
                        +(rmet(k,l,2)*rmet(m,l,2)) &
                        +(rmet(k,l,3)*rmet(m,l,3))
                  denom = 1.0d0/sqrt((dmet(k,l)*dmet(m,l))**5)
                  term=(3.0d0*dprod**2-dmet(k,l)*dmet(m,l))*denom
                  sumeff = sumeff + term
                  do n=1,3
                     dddep(n,i-k+1,jj) = dddep(n,i-k+1,jj) &
                           - x5o18*term*rmet(k,l,n)/dmet(k,l) &
                           - (x1o9*dmet(m,l)*rmet(k,l,n) - &
                           x1o3*dprod*rmet(m,l,n))*denom
                  end do
               end do
              ddep(ii,jj)=sumeff/36.0d0
   
            else !                          **i is aromatic, j is methyl:
   
               do lk=1,6
                  l = l2sub(lk)
                  k = k2sub(lk)
                  rmet(k,l,1)=x(3*(i-l)+1)-x(3*(j-k)+1)
                  sum2=rmet(k,l,1)**2
                  rmet(k,l,2)=x(3*(i-l)+2)-x(3*(j-k)+2)
                  sum2=sum2+rmet(k,l,2)**2
                  rmet(k,l,3)=x(3*(i-l)+3)-x(3*(j-k)+3)
                  sum2=sum2+rmet(k,l,3)**2
                  dmet(k,l)=sum2
               end do
               sumeff=0.0d0
               do lkm=1,18
                  l = l4sub(lkm)
                  k = k4sub(lkm)
                  m = m4sub(lkm)
                  dprod=(rmet(k,l,1)*rmet(m,l,1)) &
                        +(rmet(k,l,2)*rmet(m,l,2)) &
                        +(rmet(k,l,3)*rmet(m,l,3))
                  denom = 1.d0/sqrt((dmet(k,l)*dmet(m,l))**5)
                  term=(3.0d0*dprod**2-dmet(k,l)*dmet(m,l))*denom
                  sumeff = sumeff + term
                  do n=1,3
                     dddep(n,i-l+1,jj) = dddep(n,i-l+1,jj) &
                           - x5o36*term*(rmet(k,l,n)/dmet(k,l) &
                           + rmet(m,l,n)/dmet(m,l)) &
                           - ( dmet(m,l)*rmet(k,l,n) &
                           +dmet(k,l)*rmet(m,l,n) &
                           -3.0d0*dprod*(rmet(k,l,n)+rmet(m,l,n)))*denom*x1o18
                   end do
               end do
               ddep(ii,jj)=sumeff/36.0d0

            end if

         else if (m2ij == 25) then !           **aromatic to aromatic:

            if(ii /= jj) then !                **inter
               sumeff=0.0d0
               do kl=1,4
                  k = k3sub(kl)
                  l = l3sub(kl)
                  sum2=0.0d0
                  do n=1,3
                     diffij = x(3*(i-k)+n)-x(3*(j-l)+n)
                     sum2=sum2 + diffij**2
                     dddepp(n) = -6.0d0*diffij
                  end do
                  sum2 = 1.0d0/sum2
                  term=sum2**3
                  sumeff = sumeff + term
                  dddep(1,i-k+1,jj) = dddep(1,i-k+1,jj) + &
                        0.25d0*term*dddepp(1)*sum2
                  dddep(2,i-k+1,jj) = dddep(2,i-k+1,jj) + &
                        0.25d0*term*dddepp(2)*sum2
                  dddep(3,i-k+1,jj) = dddep(3,i-k+1,jj) + &
                        0.25d0*term*dddepp(3)*sum2
               end do
               ddep(ii,jj)=0.25d0*sumeff
      
            else !                                **intra:

               ddep(ii,jj) = 1.582d-4  !  this is 1/4.30**6

            end if
   
         else if(m2ij == 10) then !              **aromatic to ordinary:

            if(m2(i) /= 2) then !                **i is aromatic:
               sumeff=0.0d0
               do k=1,2
                  sum2=0.0d0
                  do n=1,3
                     diffij = x(3*(i-k)+n)-x(3*(j-1)+n)
                     sum2=sum2 + diffij**2
                     dddep(n,i-k+1,jj) = -6.0d0*diffij
                  end do
                  sum2 = 1.0d0/sum2
                  term = sum2**3
                  sumeff = sumeff + term
                  dddep(1,i-k+1,jj) = 0.5d0*term*dddep(1,i-k+1,jj)*sum2
                  dddep(2,i-k+1,jj) = 0.5d0*term*dddep(2,i-k+1,jj)*sum2
                  dddep(3,i-k+1,jj) = 0.5d0*term*dddep(3,i-k+1,jj)*sum2
               end do
               ddep(ii,jj) = 0.5d0*sumeff
      
            else !                                **j is aromatic:

               sumeff=0.0d0
               do k=1,2
                  sum2=0.0d0
                  do n=1,3
                     diffij = x(3*(i-1)+n)-x(3*(j-k)+n)
                     sum2=sum2 + diffij**2
                     dddepp(n) =  - 6.0d0*diffij
                  end do
                  sum2 = 1.0d0/sum2
                  term = sum2**3
                  sumeff = sumeff + term
                  dddep(1,i,jj) = dddep(1,i,jj)+0.5d0*term*dddepp(1)*sum2
                  dddep(2,i,jj) = dddep(2,i,jj)+0.5d0*term*dddepp(2)*sum2
                  dddep(3,i,jj) = dddep(3,i,jj)+0.5d0*term*dddepp(3)*sum2
               end do
               ddep(ii,jj) = 0.5d0*sumeff
            end if
   
         else if (m2(i)==2.and.m2(j)==2) then !  *** both i and j are ordinary:
   
            if (i == j) then !                   **intra:

               ddep(ii,jj) = 0.0d0
            else if (i < j) then !               **inter,first time:

               sum2=0.0d0
               do n=1,3
                  diffij = x(3*(i-1)+n) - x(3*(j-1)+n)
                  sum2 = sum2 + diffij**2
                  dddep(n,i,jj) = -6.0d0*diffij
               end do
               work = 1.0d0/sum2**4
               dddep(1,i,jj) = dddep(1,i,jj)*work
               dddep(2,i,jj) = dddep(2,i,jj)*work
               dddep(3,i,jj) = dddep(3,i,jj)*work
               ddep(ii,jj) = work*sum2
#ifdef NMODE
      
               !  --- get correction factor for this spin pair:
               
               call corf(x,i,j,ihyp(i),ihyp(j),newf,amass)
               newf = .false.
               ddep(ii,jj) = ddep(ii,jj)*gamma_nmr
               dddep(1,i,jj) = gamma_nmr*dddep(1,i,jj)+ &
                     ddep(ii,jj)*dgamma(1)
               dddep(2,i,jj) = gamma_nmr*dddep(2,i,jj)+ &
                     ddep(ii,jj)*dgamma(2)
               dddep(3,i,jj) = gamma_nmr*dddep(3,i,jj)+ &
                     ddep(ii,jj)*dgamma(3)
               do n=1,iscale
                  dratg(ii,jj,n) = ddep(ii,jj)*dgamma(n+3)/gamma_nmr
               end do
#endif
            else !                               **inter, use previous:

               dddep(1,i,jj) = -dddep(1,j,ii)
               dddep(2,i,jj) = -dddep(2,j,ii)
               dddep(3,i,jj) = -dddep(3,j,ii)
               ddep(ii,jj) = ddep(jj,ii)
#ifdef NMODE
               do n=1,iscale
                  dratg(ii,jj,n) = dratg(jj,ii,n)
               end do
#endif
            end if
   
         else
            write(6,*) 'caldis error:',i,j,m2(i),m2(j)
            call mexit(6, 1)
         end if

      end do
   end do

   return
end subroutine caldis

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine calrate here]
subroutine calrate(ddep,rate,trp)
   
   
   !   Subroutine CALculate RATE matrix.
   
   !   Correlation times are calculated according to:
   !          tau(i,j) = tau(i)*tau(j)/(tau(i)+tau(j))
   !   That is they assume that the motion of the whole molecule, as well
   !   as internal motion is isotropic.
   
   !   The spectral densities default to a resonant frequency of 500
   !   MHz (see OMEGA in subroutine mdread).
   !   Units of time for correlation times are nsec;
   !   units for the rate matrix elements are sec**-1.
   
   !-----modified P F Yip 8/9/89, by dac 10/25/89, Ross Walker 08/02/05
   use constants, only : TWOPI 
   implicit none
#  include "nmr.h"
   
   _REAL_ :: gsh, gyr, hbar, tenm16, con0, con1, con2
   parameter(gyr=2.6751965d4,hbar=1.0545919d0,gsh=gyr*hbar*gyr)
   parameter(tenm16=1.d-16,con0=gsh*gsh*tenm16)
   parameter(con1=3.d0*con0,con2=6*con0)
   
   _REAL_ :: rate(ma,ma),trp(ma,ma)
   _REAL_ :: spd(lst,3)
   _REAL_ :: ddep(ma,ma)

!Local
   _REAL_ :: omega2, taui, tauj, tauw, tsws, ts2ws
   _REAL_ :: popisq, popjsq, popij, s0, s1, s2, ootc3, tc3
   integer :: ij, i, j, k, ii, jj

   !  --- square of Larmor frequency in (radians/nsec)**2 :
   
   omega2 = (TWOPI*omega*1.d-3)**2
   
   !---- the first block ( the 80 do loop ) calculates the spec densities
   
   ij = 0
   do i=1,nath
      if (m2(i) == 1 .or. m2(i) == 4) cycle
      taui = tau(i)
      do j=1,i
         if (m2(j) == 1 .or. m2(j) == 4) cycle
         ij=ij+1
         tauj = tau(j)
         tauw = (taui*tauj)/(taui + tauj)
         tsws=tauw*tauw*omega2
         ts2ws=4.d0*tsws
         spd(ij,1) = con0*tauw
         spd(ij,2) = con1*tauw/(1.d0 + tsws)
         spd(ij,3) = con2*tauw/(1.d0 + ts2ws)
      end do
   end do
   
   
   ! Calculate transition probabilites from distances and spectral densities.
   
   ! Note: populations are actually sqrt(pop(i)) and the spectral densities
   ! are calculated to generate the symmetrized rate matrix:
   
   !                                -1
   !                     R_prime = P   R  P
   
   ! where P is a diagonal matrix containing the square roots of the
   ! populations of each group of equivalent nuclei.
   
   !      The fundamental equations are from Macura and Ernst, Mol. Phys.
   !       41, 95-117 (1980), see esp. section 4.  The symmetrization
   !       procedure has been derived by many people; an explicit
   !       presentation that corresponds to the way we do things (i.e.,
   !       with jump models for methyl groups) is given by Olejniczak,
   !       J. Magn. Res. 81, 392-394 (1989).
   
   !       The rate matrix R ( nonsymmetric )  is:
   
   !         rate(i,i)=selfrel(i,i)+{SUM(j not i )otherrel(i,j)}
   !         rate(i,j)=crossrel(i,j) = population(i)*(W2(i,j)-W0(i,j))
   !       [Note that population(i) is the same as popisq and also
   !       ( pop(i)**2 ).]
   
   !        After the transformation above,
   !              rate(i,j)=pop(i)*pop(j)*(W2(i,j)-W0(i,j))
   !        This is computed in the (else loop of if(ii.eq.jj))
   
   !        The transformed rate matrix is now symmetric and can be
   !        diagonalized.  Later, in subroutine remarc, when the
   !        intensities are calcualated, they are multiplied by
   !        population factors that effect the inverse transformation
   !        back to the original magnetization Bloch equations.
   
   
   k=1
   ii = 0
   do i=1,nath
      if (m2(i) == 1 .or. m2(i) == 4) cycle
      ii = ii + 1
      if (m2(i) == 3 ) then
         popisq = 3.0d0
      else if (m2(i) == 5) then
         popisq = 2.0d0
      else
         popisq = 1.0d0
      end if
      jj = 0
      do j=1,i
         if(m2(j) == 1 .or. m2(j) == 4) cycle
         jj = jj + 1
         if(m2(j) == 3) then
            popjsq = 3.0d0
         else if (m2(j) == 5) then
            popjsq = 2.0d0
         else
            popjsq = 1.0d0
         end if
         popij = pop(i)*pop(j)
         
         s0 = spd(k,1)
         s1 = spd(k,2)
         s2 = spd(k,3)
         if (ii == jj) then
            
            !       For the equivalent protons in the group self-relaxtion must
            !       be accounted for: (W1 + W2)
            
            if (popisq == 1.0d0) then
               rate(ii,ii) = 0.0d0
            else if (popisq == 2.0d0) then
               if (iroesy > 0) then
                  rate(ii,ii)=0.5d0*(9.d0*s0+3.d0*s1+2.*s2)*ddep(ii,ii)
               else
                  rate(ii,ii)=2.0d0*(0.50d0*s1+s2)*ddep(ii,ii)
               end if
            else if (popisq == 3.0d0) then
               if (iroesy > 0) then
                  
                  ! --- just use the fast methyl motion limit here:
                  
                  rate(ii,ii) = 0.5d0*(13.d0*s0 + 3.d0*s1 + &
                        3.d0*s2)*ddep(ii,ii)
               else
                  
                  !      For methyl relaxation,
                  !  --- use Kalk & Berendsen formula [JMR 24, 343 (1978)], Eq. 11:
                  
                  taui = 0.5d0* tau(i)
                  ootc3 = 1./taui + 1./taumet
                  tc3 = 1./ootc3
                  rate(ii,ii)=2.0d0*5.372d0*(0.25d0*( &
                        taui/(1.+omega2*taui*taui) &
                        + 4.d0*taui/(1.+omega2*4.d0*taui*taui) ) &
                        + 0.75d0*( tc3/(1.d0+omega2*tc3*tc3) &
                        + 4.d0*tc3/(1.d0+omega2*4.d0*tc3*tc3) ) )
#if 0
                  write(6,*) 'for methyl: ',ii,rate(ii,ii),taui,tc3
                  term1 = 2.0d0*5.372d0*0.25d0* &
                        taui/(1.+omega2*taui*taui)
                  term2 = 2.0d0*5.372d0*0.25d0* &
                        4.d0*taui/(1.+omega2*4.d0*taui*taui)
                  term3 = 2.0d0*5.372d0* &
                        0.75d0*tc3/(1.d0+omega2*tc3*tc3)
                  term4 = 2.0d0*5.372d0* &
                        0.75d0*4.d0*tc3/(1.d0+omega2*4.d0*tc3*tc3)
                  write(6,*) '            ', term1,term2,term3,term4
#endif
               end if
            else
               write(6,*) 'bad value for popisq:', popisq,i,ii
               call mexit(6, 1)
            end if  ! (popisq == 1.0d0)
            trp(ii,ii) = 0.0d0
         else
            
            !         Calculate the cross relaxation terms in rate: (W2 - W0)
            
            if (iroesy > 0) then
               rate(ii,jj) = ddep(ii,jj)*(s2 + 4.0d0*s0)*popij*0.5d0
            else
               rate(ii,jj) = ddep(ii,jj)*(s2 - s0)*popij
            end if
            
            !       Now the off-diagonal terms which contribute to the diagonal rate
            !       elements: (W0 + 2W1 + W2). These are summed below.
            !     Note that trp is not symmetric, therefore trp2 is needed;
            !     trp is the lower diag and trp2 is the upper diag.
            
            if (iroesy > 0) then
               trp(ii,jj) = &
                     ddep(ii,jj)*(3.d0*s1 + s2 + 5.d0*s0)*popjsq*0.5d0
               trp(jj,ii) = &
                     ddep(ii,jj)*(3.d0*s1 + s2 + 5.d0*s0)*popisq*0.5d0
            else
               trp(ii,jj) = ddep(ii,jj)*(s1 + s2 + s0)*popjsq
               trp(jj,ii) = ddep(ii,jj)*(s1 + s2 + s0)*popisq
            end if
#ifdef DEBUG_NMR
            !           write(6,*) 'trp:',ii,jj,trp(ii,jj),trp(jj,ii)
#endif
            
         end if  ! (ii == jj)
         k = k+1
      end do
   end do
   
   !     The diagonal self-relaxation and all the cross-peak terms
   !     of the transition probabilities have already been assigned
   !     to the rate matrix.  Here the summation over i.ne.j is
   !     performed to generate the diagonal terms of the rate matrix.
   
   do ii=1,natmet
      do jj=1,natmet
         rate(ii,ii)=rate(ii,ii)+trp(ii,jj)
      end do
#ifdef DEBUG_NMR
      write(6,*) 'rate:',ii,ii,rate(ii,ii)
#endif
   end do
   
   ! --- filling up the whole rate matrix:
   
   do ii=1,natmet
      do jj=ii+1,natmet
         rate(ii,jj)=rate(jj,ii)
#ifdef DEBUG_NMR
         write(6,*) 'rate:',ii,jj,rate(ii,jj)
#endif
      end do
   end do
   
   return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine amatg here]
subroutine amatg(vecs,eig,taum,amat,iamat)
   
   !   --- compute the intensity matrix at Gaussian points needed for
   !         the integral form of the exact derivative
   use constants, only : zero 
   implicit none
#  include "nmr.h"
   _REAL_    vecs(ma,ma), eig(ma), amat(ma,ma,5), s(5)
   _REAL_    et(ma,5)
   integer   iamat(ma)
   
   integer   i, j, k
   integer   ig, igmax
   _REAL_    taum, taus
   _REAL_    v
   save s,igmax
   data s/0.0469100770d0,0.2307653449d0,0.5d0,0.7692346551d0, &
         0.9530899230d0/
   data igmax /5/
   
   do ig=1,igmax
      taus = s(ig)*taum
      do k=1,natmet
         et(k,ig) = exp(-taus*eig(k))
      end do
   end do
   do i=1,natmet
      if (iamat(i) == 0) cycle
      do j=1,natmet
         amat(i,j,1) = zero
         amat(i,j,2) = zero
         amat(i,j,3) = zero
         amat(i,j,4) = zero
         amat(i,j,5) = zero
         do k=1,natmet
            v = vecs(i,k)*vecs(j,k)
            amat(i,j,1) = amat(i,j,1) + v*et(k,1)
            amat(i,j,2) = amat(i,j,2) + v*et(k,2)
            amat(i,j,3) = amat(i,j,3) + v*et(k,3)
            amat(i,j,4) = amat(i,j,4) + v*et(k,4)
            amat(i,j,5) = amat(i,j,5) + v*et(k,5)
         end do
      end do
   end do
   return
end subroutine amatg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine noecalc here]
subroutine noecalc(x,f,xx,ix)
   
   !  Subroutine NOEsy intensity CALCulation:
   
   !  --- computes relaxation matrix approximation to noe
   !      intensities and their derivatives.  Places energy
   !      into enoe and the negative of the gradient is used
   !      to update f.
   
   !   Inputs:
   
   !         x        --  coordinates
   !         f        --  forces
   !         xx,ix    -- main storage arrays
   
   !   Outputs:
   
   !         enoe     -- contains the penaly energy for noe-related restraints
   !         f        -- is updated with new derivatives
   
   !   Overview:
   
   !      This is mostly a "driver" routine.  It reads input for the
   !      submolecules, in namelist "&noeexp", then calls indexn, caldis
   !      and remarc to do the brunt of the work.  The computation
   !      and printing of various statistics is also handled by this routine.
   
   !      There is a major additional section, indicated by #ifdef NMODE,
   !      that allows vibrational modes to be read in and used to compute
   !      motional correction factors for the NOEs.  This is still
   !      experimental code, and is *not* currently described in the
   !      users manual.
   
   !-----------------------------------------------------------------------------
   use file_io_dat
   implicit none
#ifdef NMODE
   integer:: line, n
   _REAL_ :: consq, dev, dpen, earg, edev, penalty
#endif
   integer:: i, im, imet, imix, indxs, ip, ix, k, ksub, l, nuse
   _REAL_ :: f, prob, r, rms, tauc, ten, x, xhyd, xmet, xrfac6, xrfact, xx, z
   logical newf
   
#  include "nmr.h"
#  include "../include/memory.h"
#  include "../include/md.h"
#  include "def_time.h"
#ifdef NMODE
   
   real freq4, xw(3*matom)
   _REAL_ kt
   parameter (kt = 0.6d0)
   parameter (consq = 2.39805d-3)
   !                  ----CONSQ = hc/2kT in cm, with T=300K
   !                       (use for quantum, Bose statistics)
#endif
#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra.h"
   
   _REAL_ frespa
   dimension xx(*),ix(*)
   dimension xmet(isubr),imet(isubi)
   equivalence (nath,imet(1))
   equivalence (tau(1),xmet(1))
   dimension x(*),f(*)
   dimension xhyd(3*ma)
   
   ten = 10.d0
   
   ! --- set up printer output:
   
   1000 format(1x,'Wnoesy = ', f8.3)
   1010 format(1x,79('-')/4x,'sub mix    Proton 1',8x,'Proton 2',6x, &
         '  Exp.      Calc.     Penalty   ID'/1x,79('-'))
   if (master .and. iprint /= 0) then
      rewind (81)
      write(81,1000) wnoesy
      write(81,1010)
      write(6,1000) wnoesy
   end if
   
   ! --- zero out the forces arising from the "extra" degrees of freedom:
   
   do im=1,iscale
      f(3*natom+im) = 0.0d0
   end do
   
   ! --- use respa impulses for this force every noeskp steps:
   
   enoe = 0.0d0
   if (mod(irespa,noeskp) /= 0) return
   frespa = noeskp
   
   ! --- read in relaxation matrix parameters, desired mixing times, and
   !      experimental intensities
   
   call timer_start(TIME_NOECALC1)
   ntot = 0
   ntota = 0
   ntotb = 0
#ifdef MPI
   ksub = mytaskid + 1 - numtasks
#else
   ksub = 0
#endif
   enoe = 0.0d0
   do i=1,3*natom+iscale
      xx(l110 -1 + i) = 0.0d0
   end do
   
   140 continue
   
   ! --- grab submolecule info from where it was stored in noeread:
   
#ifdef MPI
   ksub = ksub + numtasks
#else
   ksub = ksub + 1
#endif
   if (ksub > maxsub) then
      call timer_stop(TIME_NOECALC1)
      goto 280
   end if
   indxs = l105 + (ksub-1)*isubr - 1
   do i=1,isubr
      xmet(i) = xx(indxs+i)
   end do
   indxs = i65 + (ksub-1)*isubi - 1
   do i=1,isubi
      imet(i) = ix(indxs+i)
   end do
   
   
   ! --- initialize coordinates of the sub-molecule:
   
   if (ihet == 0) then
      l = 0
      do i=1,nath
         k = 3*(ihyp(i)-1)
         xhyd(l+1)=x(k+1)
         xhyd(l+2)=x(k+2)
         xhyd(l+3)=x(k+3)
         l=l+3
      end do
   end if
   
   ! --- overall rotational correlation time in h2o is taurot; here
   !      we hard-wire in the assumption that this is increased by
   !      23% in d2o.
   
   tauc = taurot*2.d0
   if (id2o == 1) tauc = 1.23d0*tauc
   
   ! --- in calrate, the relative tumbling time for each proton pair
   !      is calculated as taui*tauj/(taui+tauj).  This gives a constant
   !      in the current code, but is included for (future?) models in
   !      which the effective correlation time is not the same for all
   !      atoms.  Here we set up the tau(i) array as twice the rotational
   !      correlation time.
   
   do im=1,nath
      tau(im) = tauc
   end do
   
   !----calculate  penalty function and its gradient:  caldis computes the
   !       distance-dependent parts of the rate matrix, and their derivatives,
   !       and remarc completes the work (q.v.).  Subroutine remhet is
   !       called for vibrational calculations of Lipari-Szabo order
   !       parameters for heteronuclear relaxation.  Logical variable "newf"
   !       is true if the frequencies are changed since the last time
   !       order parameters were evaluated.
   
   call timer_stop(TIME_NOECALC1)
   if (ksub == 1) then
      newf = .true.
   else
      newf = .false.
   end if
   if (ihet == 1) then
#ifdef NMODE
      call remhet(xx(l110),xx(lcrd),ksub,newf,xx(lmass))
#endif
   else
      call timer_start(TIME_CALDIS)
#ifdef NMODE
      call caldis(xhyd,xx(l115),xx(l120),newf,xx(lmass))
#else
      call caldis(xhyd,xx(l115),xx(l120))
#endif
      call timer_stop(TIME_CALDIS)
      call remarc(xx(l115),xx(l120),xx(l110),ksub, &
            xx(l125),xx(l130),xx(l135),xx(l140),xx(l145))
   end if
   
   !----redo calc. for another submolecule
   
   call timer_start(TIME_NOECALC1)
   goto 140
   
   !---if all sub mol. are done, compute and print some statistics,
   !     update the forces and return:
   
   280 continue
   
   ! --- get linear correlation coefficient between cacl. and obs.
   
   1050 format(/50x,'Total NOESY penalty: ',f7.2/)
   1060 format(/21x,'   #   Pearson r  rms error   R1      Rr      ', &
         /1x,75('-'))
   call timer_start(TIME_NOECALC2)
   if (master .and. iprint /= 0) then
      write(81,1050) enoe
      if (ntot < 3) goto 350
      write(81,1060)
      call pearsn(exper,calc,ntot,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            'Full  Correlation: = ',nuse,r,rms,xrfact,xrfac6
      call pearsn(expera,calca,ntota,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            'Intra Correlation: = ',nuse,r,rms,xrfact,xrfac6
      call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            'Inter Correlation: = ',nuse,r,rms,xrfact,xrfac6
      ntotb= 0
      do ip=1,ntot
         if (exper(ip) < 0.1) cycle
         ntotb = ntotb + 1
         experb(ntotb) = exper(ip)
         calcb(ntotb) = calc(ip)
      end do
      call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            '>0.1  Correlation: = ',nuse,r,rms,xrfact,xrfac6
      
      ntotb= 0
      do ip=1,ntot
         if (exper(ip) < 0.05 .or. exper(ip) > 0.1) cycle
         ntotb = ntotb + 1
         experb(ntotb) = exper(ip)
         calcb(ntotb) = calc(ip)
      end do
      call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            '0.05-0.1 Correl  : = ',nuse,r,rms,xrfact,xrfac6
      
      ntotb= 0
      do ip=1,ntot
         if (exper(ip) > 0.05 .or. exper(ip) < 0.025) cycle
         ntotb = ntotb + 1
         experb(ntotb) = exper(ip)
         calcb(ntotb) = calc(ip)
      end do
      call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            '0.025-0.05 Correl: = ',nuse,r,rms,xrfact,xrfac6
      
      ntotb= 0
      do ip=1,ntot
         if (exper(ip) > 0.025) cycle
         ntotb = ntotb + 1
         experb(ntotb) = exper(ip)
         calcb(ntotb) = calc(ip)
      end do
      call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
      write(81,'(a21,i5,4f10.5)') &
            '<0.025 Correlation:= ',nuse,r,rms,xrfact,xrfac6
      
      !  --statistics for each mixing time:
      
      do imix=1,nummt
         ntotb= 0
         do ip=1,ntot
            if (ipmix(ip) /= imix) cycle
            ntotb = ntotb + 1
            experb(ntotb) = exper(ip)
            calcb(ntotb) = calc(ip)
         end do
         call pearsn(experb,calcb,ntotb,r,prob,z,rms,xrfact,xrfac6,nuse)
         write(81,'(i10,a11,i5,4f10.5)') &
               imix,' Correl: = ',nuse,r,rms,xrfact,xrfac6
      end do
      350 continue
   end if
   call timer_stop(TIME_NOECALC2)

   ! --- update full derivatives

#ifdef NMODE
   do i=1,3*natom
      f(i) = 0.0
   end do
   do i=3*natom+1,3*natom+iscale
      f(i) = f(i) + frespa*xx(l110 -1 + i)
   end do
   
   1080 format(' ',78('-'))
   if (iprint /= 0) then
      write(6,1080)
      write(6,*) 'Heteronuclear normal mode analysis:'
      write(6,*)
      if( iscale > 1) then
         write(6,*) 'frequencies:'
         write(6,'(6f12.6)') (x(3*natom+i),i=1,iscale-1)
      end if
      write(6,1080)
      rewind (81)
      390 read(81,'(a80)',end=400) line
      write(6,'(a80)') line
      goto 390
      400 rewind (81)
      write(6,1080)
   end if
   
   !  --- now set up penalties to keep the effective frequencies close
   !         to their "true" values:  right now assumes first six
   !         frequencies are zero.
   
   !      set the penalty function to:
   
   !      penalty = (xdev/freq_true)(freq_ref - freq_true)**2
   !               + exp(-6*freq_ref)
   
   !      where freq_exp is the "true" frequency, and freq_ref is the
   !          current value of the "refined" frequency.  The first term
   !          is a usual quadratic penalty (with weight equal to the
   !          inverse of freq_true,) and the second term should serve to
   !          keep all frequencies positive.
   
   !     xdev = 0.05  is the default
   
   edev = 0.0
   do n=1,iscale-1
      dev = x(3*natom+n) - freq(n)
      earg = min(-6.0*x(3*natom+n),ten)
      earg = max(-ten-ten,earg)
      penalty = xdev*dev**2/freq(n) + exp(earg)
      enoe = enoe + penalty
      edev = edev + penalty
      dpen = 2.0*xdev*dev/freq(n) - 6.*exp(earg)
      f(3*natom+n) = f(3*natom+n) - dpen
   end do
   if (iprint /= 0) then
      write(6,'(a23,f12.5)') 'global scaling factor: ', &
            x(3*natom+iscale)
      write(6,'(a29,e14.6)') 'frequency deviation penalty: ',edev
   end if
#else
   do i=1,3*natom+iscale
      f(i) = f(i) + frespa*xx(l110 -1 + i)
   end do
#endif
   
   return
   
#ifdef DEBUG_NMR
   1030 format(1x,79('=')/'Data for submolecule',i4,': id2o = ', &
         i1,', oscale = ',e12.5,', taumet = ',e12.5,/ &
         36x,'omega = ',e12.5, 'taurot =',e12.5/)
#endif
   
end subroutine noecalc 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine noeread here]
subroutine noeread(xx,ix,ih)
   
   use file_io_dat
   implicit none
#ifdef NMODE
   integer:: im, ivec, ivform, k, l, lmax, lmin, n3, nf
   _REAL_ :: consq, omecut, omescl, xkappa
#endif
   integer:: i, ifind, iin, imet, imix, imol, indxs, ix, j, ksub, numpks
   _REAL_ :: xmet, xx
#  include "nmr.h"
#  include "../include/memory.h"
#  include "../include/md.h"
#  include "parallel.h"
   
   dimension xx(*),ix(*)
   character(len=4) ih(*)
   dimension xmet(isubr),imet(isubi)
   equivalence (nath,imet(1))
   equivalence (tau(1),xmet(1))
   logical hetero
   character(len=80) line
   namelist /noeexp/ npeak,emix,ihp,jhp,aexp,awt,arange,id2o,oscale, &
         taumet,omega,invwt1,invwt2,hetero,iroesy,taurot,peakid
#ifdef NMODE
   logical hsfull
   integer g98_vecs,jaguar_vecs
   character(len=4) star
   character(len=2) spacer
   real freq4
   _REAL_ kt
   namelist /nmode/ nvect,hsfull,bose,xdev,xkappa,ivform,omegax, &
         iusev,jaguar_vecs,g98_vecs,vtemp,per_mode,nmsnap
#endif

   
   maxsub = 0
   imol = 0
   ksub = 0
   id2o = 0
   iroesy = 0
   hetero = .false.
   taurot = 1.0d0
   oscale = 1.0d0
   invwt1 = 1.0d0
   invwt2 = 1.0d0
   if (iredir(4) /= 0) then
      call amopen(35,redir(4)(1:iredir(4)),'O','F','R')
      iin = 35
      write(6,1020) redir(4)(1:iredir(4))
   else
      return
   end if
   
   !  --- read and echo title from noeexp file:
   
   write(6,*) 'Here are comments from the NOEsy input file:'
   42 read(iin,'(a)',end=140) line
   if (line(1:1) == '#') then
      write(6,*) line
      goto 42
   end if
   backspace (iin)
   write(6,*)
   
#ifdef NMODE
   
   !    --- first, read in namelist nmode, after setting defaults:
   
   nvect = 99999
   hsfull = .true.
   bose = .true.
   xkappa = 1.0
   omegax = 1000.
   xdev = 0.05
   ivform = 1
   jaguar_vecs = 0
   g98_vecs = 0
   per_mode = .false.
   nmsnap = 0
   vtemp = 300.
   
   call nmlsrc('nmode',iin,ifind)
   if (ifind == 0) then
      write(6,*) 'Unable to fine namelist "nmode"'
      call mexit(6,1)
   end if
   read(iin,nml=nmode,err=460)
   write(6,*) 'read nmode: ',nvect,hsfull,bose,xdev,xkappa,ivform, &
         vtemp,omegax
   if (hsfull) then
      ihsful = 1
   else
      ihsful = 0
   end if
   kt = vtemp*0.002
   !                  ----CONSQ = hc/2kT in cm
   !                       (use for quantum, Bose statistics)
   consq = 0.71942/vtemp
   xx(lcrd-1 +3*natom+iscale) = xkappa
   
   !   --- read in the normal mode frequencies and vectors:
   
   nf = 60
   
   if( jaguar_vecs > 0 ) then
      call amopen(nf,vecs,'O','F','R')
      do lmin=1,nvect,7
         lmax = min(nvect, lmin+6)
         read(nf, '(13x,7f9.2)') (freq(l), l=lmin,lmax)
         write(6, '(13x,7f9.2)') (freq(l), l=lmin,lmax)
         do k=1,3*natom
            read(nf, '(13x,7f9.5)') (vect(k,l), l=lmin,lmax)
         end do
         read(nf, '(a2)' ) spacer
      end do
   else if( g98_vecs > 0 ) then
      call amopen(nf,vecs,'O','F','R')
      do lmin=1,nvect,5
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         lmax = min(nvect, lmin+4)
         read(nf, '(23x,5f10.4)') (freq(l), l=lmin,lmax)
         write(6, '(23x,5f10.4)') (freq(l), l=lmin,lmax)
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         read(nf, '(a2)' ) spacer
         do k=1,3*natom
            read(nf, '(23x,5f10.5)') (vect(k,l), l=lmin,lmax)
         end do
      end do
   else
      if (ivform == 0) then
         call amopen(nf,vecs,'O','U','R')
         read(nf,err=470) n3
      else
         call amopen(nf,vecs,'O','F','R')
         read(nf,'(a40)',err=470) title
         read(nf,'(i5)',err=470) n3
      end if
      if(n3 /= natom*3) then
         write(6,*) 'number of atoms wrong in vecs: ',n3,natom*3
         call mexit(6, 1)
      end if
      
      !  --- do not need coords, put temporarily into vect(,1):
      
      if (ivform == 0) then
         read(nf,err=470)
      else
         read(nf,'(7f11.5)',err=470) (vect(j,1),j=1,n3)
      end if
      do l = 1,nvect
         if (ivform == 0) then
            read(nf,end=70,err=470) ivec,freq4
            freq(l) = freq4
            read(nf,err=470) (vect(k,l), k=1,n3)
         else
            read(nf,'(a4)',end=70,err=470) star
            read(nf,'(i5,f12.6)',err=470) ivec,freq(l)
            read(nf,'(7f11.5)',err=470) (vect(k,l),k=1,n3)
         end if
         if (freq(l) < 0.5) freq(l) = 10000.
      end do
      goto 80
      70 nvect = l-1
   end if
   80 write(6,*) 'Found',nvect,' eigenvectors in vecs file'
   close(nf)
   
   !  --- hard-wired: scale frequencies as rafael suggested:
   
   if (freqe /= 'dummy') then
      call amopen(nf,freqe,'O','F','R')
      read(nf,*) omecut,omescl
      do i=1,nvect
         if(freq(i) < omecut) freq(i) = omescl*freq(i)
      end do
   end if
   
   !  --- initially set scaling factors to "true" frequencies,
   
   do im=1,iscale-1
      xx(lcrd-1 +3*natom+im) = freq(im)
   end do
   
   !  --- kludge for now: zero out velocities of all real atoms:
   
   do i=1,3*natom
      xx(lvel-1+i) = 0.0
   end do
#endif
   
   140 numpks = 0
   
   !   --zero out the arange and npeak arrays:
   
   do i=1,mxtau
      npeak(i) = 0
      do j=1,mxp
         arange(i,j) = 0.0d0
      end do
   end do
   call nmlsrc('noeexp',iin,ifind)
   if (ifind == 0) goto 280
   read(iin,nml=noeexp,end=280,err=450)
   imol = imol + 1
   if (hetero) then
      ihet=1
   else
      ihet=0
   end if
   do imix=1,mxtau
      if(npeak(imix) < 0) goto 190
      if (npeak(imix) > mxp) then
         write(6,*) 'Npeak is too big.'
         call mexit(6, 1)
      end if
      
      do i=1,npeak(imix)
         numpks = numpks + 1
      end do
   end do
   
   ! ------ if we fall out of the do-loop, set the number of mixing times
   !         to MXTAU:
   
   nummt = mxtau
   goto 200
   
   ! ------come here when negative peak number is encountered;
   !         total number of mixing times is nummt
   
   190 nummt=imix-1
   200 write(6,1040) nummt, numpks
   if (nummt <= 0) goto 280
   
   ! --- done reading experimental info; get various information about
   !       protons in this submolecule into the appropriate arrays:
   
   if (.not.hetero) call indexn(ix,ih,iin)
   
   ! --- Store info on this submolecule into appropriate locations in the
   !       xx and ix arrays:
   
   ksub = ksub + 1
   if (ksub > mxsub) then
      write(6,*) 'Too many sub-molecules!'
      call mexit(6, 1)
   end if
   indxs = l105 + (ksub-1)*isubr - 1
   do i=1,isubr
      xx(indxs+i) = xmet(i)
   end do
   indxs = i65 + (ksub-1)*isubi - 1
   do i=1,isubi
      ix(indxs+i) = imet(i)
   end do
   
   !   ---cycle back for another sub-molecule:
   
   goto 140
   
   280 maxsub = ksub
   return
   
   ! namelist error on read:
   
   450 write(6,*) 'Namelist reports error in reading noeexp'
   write(6,*) '-- Subscript out of range implies dimensioning problem'
   write(6,*) '-- (see nmr.h)'
   call mexit(6,1)
#ifdef NMODE
   460 write(6,*) 'Namelist reports error in reading nmode'
   call mexit(6,1)
   470 write(6,*) 'Error reported in reading frequency file'
   call mexit(6,1)
#endif
   
   1020 format(' Noesy volumes will be read from file: ',a)
   1040 format(1x,'Read ', i5,' mixing times with ', i5, &
         ' total peaks.')
   
end subroutine noeread 

end module relax_mat
