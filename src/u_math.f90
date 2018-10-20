module u_math

  ! ModAEM 1.4pre1
  ! Copyright (c) 2001 WHPA Inc.
  !
  ! This program is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU General Public License
  ! as published by the Free Software Foundation; either version 2
  ! of the License, or (at your option) any later version.
  !
  ! This program is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ! GNU General Public License for more details.
  !
  ! You should have received a copy of the GNU General Public License
  ! along with this program; if not, write to the Free Software
  ! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
  !
  ! Contact the author by e-mail at: vic@wittmanhydro.com
  ! Or by regular mail at:
  ! WHPA, Inc
  ! 320 W 8th St
  ! Bloomington, IN 47401

! This module contains useful computational routines and special functions required
! by certain modules

use u_constants

implicit none

public :: rExpIntE1, &
          rPolynomial

  ! Constants
  real (kind=ModAEM_Double),public,parameter :: rGAMMA = .5772156649_ModAEM_Double
  ! The cutoff for exp(-rU) ~ 0.0 was adopted from Haitjema - GFLOW
  real (kind=ModAEM_Double),public,parameter :: rU_CUTOFF = 150.0_ModAEM_Real

contains

function rExpIntE1(rU) result (rE1)
  ! Returns the exponential integral of the first kind, E1(U)
  ! From Abramowicz and Stegun, 5.1.53 and 5.1.56
  real (kind=ModAEM_Real),intent(in) :: rU
  real (kind=ModAEM_Real) :: rE1
  ! Locals
  real (kind=ModAEM_Real) :: rU0

  ! Local computations are performed in double precision, return value is single
  ! Coefficient arrays for polynomial computations
  ! Polynomial for 0 <= x <= 1 [A&S 5.1.53]
  real (kind=ModAEM_Double ),dimension(0:5) :: &
        rA_SmallX = &
              (/  0.00107857_ModAEM_Double, & 
                 -0.00976004_ModAEM_Double, &
                  0.05519968_ModAEM_Double, &
                 -0.24991055_ModAEM_Double, &
                  0.99999193_ModAEM_Double, &
                 -0.57721566_ModAEM_Double /)
  ! Polynomial for the numerator of the rational expression for x > 1 [A&S 5.1.56]
  real (kind=ModAEM_Double ),dimension(0:4) :: &
        rA_LargeX_Numerator = &
              (/  1.0000000000_ModAEM_Double, &
                  8.5733287401_ModAEM_Double, &
                 18.0590169730_ModAEM_Double, &
                  8.6347608925_ModAEM_Double, &
                  0.2677737343_ModAEM_Double /)
  ! Polynomial for the denominator of the rational expression for x > 1 [A&S 5.1.56]
  real (kind=ModAEM_Double ),dimension(0:4) :: &
        rA_LargeX_Denominator = &
              (/  1.0000000000_ModAEM_Double, &
                  9.5733223454_ModAEM_Double, &
                 25.6329561486_ModAEM_Double, &
                 21.0996530827_ModAEM_Double, &
                  3.9584969228_ModAEM_Double /)
  
  ! Check the argument
  if ( rU <= rZERO ) then
    write(unit=*,fmt="(""rfExpIntE1: Error - argument is less than zero"")")
    rE1 = -HUGE(ModAEM_Real)
    return
  endif
  
  ! Compute the E1(U) function
  if ( rU == rZERO ) then
    rE1 = -HUGE(ModAEM_Real)
  else if ( rU < rONE ) then
    rE1 = rPolynomial(rU,rA_SmallX,5) - log(rU) 
  else
    ! The cutoff for exp(-rU) ~ 0.0 was adopted from Haitjema - GFLOW
    if ( rU < rU_CUTOFF ) then
      rU0 = rU
    else
      rU0 = rZERO
    end if
    rE1 = rPolynomial(rU,rA_LargeX_Numerator,4) /   &
               ( rPolynomial(rU,rA_LargeX_Denominator,4) * rU * exp(rU0) )
  endif

  return
end function rExpIntE1

function rPolynomial(rX,rA,iOrder) result (rY)
  ! Computes the iOrder order polynomial rA(0)x^iOrder + rA(1)x^(iOrder-1) + ...
  ! rA = vector of coefficients, dimension(0:iOrder)
  ! Allcomputations are in double precision
  ! Arguments
  real (kind=ModAEM_Double),intent(in) :: rX
  real (kind=ModAEM_Double),intent(in),dimension(0:) :: rA
  integer (kind=ModAEM_Integer),intent(in) :: iOrder
  real (kind=ModAEM_Double) :: rY
  ! Locals
  integer (kind=ModAEM_Integer) :: i

  ! Range checking
  if ( ubound(rA,1) < iOrder ) then
    rY = -HUGE(ModAEM_Double)
    return
  endif

  rY = rA(0)
  do i=1,iOrder
    rY = rY*rX + rA(i)
  end do

  return
end function rPolynomial

end module u_math

