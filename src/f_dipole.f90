module f_dipole

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


  !! module f_dipole (FDP)
  !! 
  !! Written by Victor A.Kelson
  !! 
  !! Revision History:
  !!   1.0.0   17 March 1999
  !!           First "source code release" version, adapted from
  !!           previous testing version.
  !!
  !! Module of data structures and functions for second order
  !! dipoles of complex strength.  The real part of the strength
  !! is the 'dipole' strength (jump in the streamfunction) and the
  !! imaginary part is the 'doublet' strength (jump in the 
  !! potential)
  !!
  !! This module encapsulates containers of dipole functions
  !! in layers.  The various data structures are organized for 
  !! efficient computation on vector or parallel machinery.
  !!
  !! Limitations:
  !!   At present, this module has no hooks for series expansion
  !!   block implementation.  No known computational bugs exist.
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   io         --  Universal ModAEM I/O functions and constants
  !!   i_dipole   --  Influence function module for dipoles
  !!

use u_constants
use u_io
use i_dipole
use i_linesink

implicit none

public

  type,public :: FDP_DIPOLE
    !! type FDP_DIPOLE
    !!
    !! Type that holds information for one dipole
    !!
    !! Members:
    !!   complex :: cZC,cZL
    !!     The center and directed length vector for the dipole
    !!   complex :: cRho(3)
    !!     The complex strength coefficients for the dipole
    !!     cRho(1) is the strength at the first end, cRho(2) is the
    !!     strength at the center and cRho(3) is the strength at
    !!     the second end.
    !! 
    complex (kind=ModAEM_Real) :: cZC
    complex (kind=ModAEM_Real) :: cZL
    complex (kind=ModAEM_Real),dimension(3) :: cRho
  end type FDP_DIPOLE

  type,public :: FDP_COLLECTION
    !! type FDP_COLLECTION
    !!
    !! Type that holds the dipole entries for a layer
    !!
    !! Members:
    !!   type (FDP_DIPOLE),pointer :: Dipoles(:)
    !!     FDP_DIPOLE structures that hold dipole information
    !!   integer :: iCount
    !!     Number of FDP_DIPOLE structures currently in use
    !! 
    type (FDP_DIPOLE),dimension(:),pointer :: Dipoles
    integer (kind=ModAEM_Integer) :: iCount
  end type FDP_COLLECTION

contains

function FDP_Create(iMax,io) result (fdp)
  !! function FDP_Create
  !!
  !! Creates a new FDP_COLLECTION object
  !!
  !! Calling Sequence:
  !!    fdp => FDP_Create(iMax)
  !!
  !! Arguments:
  !!    (in)      integer :: iMax
  !!                The maximum number of dipoles to be stored
  !!
  !! Return Value:
  !!   On success, fdp points to a new FDP_COLLECTION object
  !!   On failure (allocation error), fatal error
  !!
  ! [ ARGUMENTS ]
  integer (kind=ModAEM_Integer),intent(in) :: iMax
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  type ( FDP_COLLECTION ),pointer :: fdp
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat

  allocate (fdp,stat=iStat)
  if (IO_Assert( (iStat == 0), "FDP_Create: allocation failed",io )) return

  allocate (fdp%Dipoles(iMax),stat=iStat)
  if (IO_Assert( (iStat == 0), "FDP_Create: allocation failed",io )) return

  fdp%Dipoles = FDP_DIPOLE(cZERO,cmplx(rONE,rONE,ModAEM_Real),(/cZERO,cZERO,cZERO/))
  fdp%iCount = 0

  return
end function FDP_Create

subroutine FDP_New(fdp,cZ1,cZ2,cRho,iRV,io)
  !! subroutine FDP_New
  !!
  !! Makes a new dipole entry. On call, the geometry and strengths of 
  !! the dipole are provided. The internal dipole structures are then
  !! set up. Returns iRV=the index into the dipole table on success or
  !! iRV<0 on failure.
  !!
  !! Calling Sequence:
  !!    call FDP_New(fdp,cZ1,cZ2,cRho)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    complex :: cZ1,cZ2
  !!             Complex coordinates of the ends of the dipole
  !!   (in)    complex :: cRho(3)
  !!             The complex strengths at cZ1, at the center, and
  !!             at cZ2, respectively.
  !!   (out)   integer :: iRV
  !!             On success, the index in fdp%Dipoles used
  !! 
  !! Note: On failure, forces a fatal error
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  complex (kind=ModAEM_Real),intent(in) :: cZ1,cZ2
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cRho
  integer (kind=ModAEM_Integer),intent(out) :: iRV
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_New: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_New: FDP_Alloc has not been called",io )) return
    if (IO_Assert( (fdp%iCount < size(fdp%Dipoles)), "FDP_New: Space exhausted",io )) return
    if (IO_Assert( (size(cRho) == 3), "FDP_New: Illegal strength vector",io )) return
    !! Note: the next assertion ensures that the caller has checked the length appropriately
    if (IO_Assert( (abs(cZ2-cZ1) > rTINY ), "FDP_New: Dipole length is too short",io )) return
  endif

  ! Compute the center and directed length parameters
  fdp%iCount = fdp%iCount + 1
  dip => fdp%Dipoles(fdp%iCount)
  dip = FDP_DIPOLE(rHALF*(cZ2+cZ1),rHALF*(cZ2-cZ1),cRho)
  iRV = fdp%iCount

  return
end subroutine FDP_New

subroutine FDP_Update(fdp,iDP,cRho,io)
  !! subroutine FDP_Update
  !!
  !! Updates the strength coefficients for a dipole entry.
  !!
  !! Calling Sequence:
  !!    call FDP_Update(fdp,iDP,cRho,iRV)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to be used
  !!   (in)    integer :: iDP 
  !!             The index for the dipole
  !!   (in)    complex :: cRho(3)
  !!             The complex strengths at cZ1, at the center, and
  !!             at cZ2, respectively.
  !!
  !! Note: On failure, forces a fatal error
  !!
  ! [ ARGUMENTS ]
  type ( FDP_COLLECTION ),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in) :: iDP
  complex (kind=ModAEM_Real),dimension(3),intent(in) :: cRho(3)
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
    if (IO_Assert( (iDP <= size(fdp%Dipoles)), "FDP_Update: Space exhausted",io )) return
    if (IO_Assert( (size(cRho) == 3), "FDP_Update: Illegal strength vector",io )) return
  endif

  dip => fdp%Dipoles(iDP)
  dip%cRho = cRho

  return
end subroutine FDP_Update

subroutine FDP_GetInfluence_IDP(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cF,io)
  !! subroutine FDP_GetInfluence
  !!
  !! Retrieves arrays of influence functions for use in matrix generation
  !!
  !! Calling Sequence:
  !!    call FDP_GetInfluence(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cF)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ) :: fdp
  !!             The FDP_COLLECTION to be used
  !!   (in)    integer :: iWhich
  !!             The influence function to be computed;  iWhich values are
  !!                INFLUENCE_P   - Complex potential
  !!                INFLUENCE_Q   - Complex discharge
  !!                INFLUENCE_F   - Integrated flux
  !!                INFLUENCE_G   - Areal infiltration
  !!                INFLUENCE_Q   - Extraction rate
  !!                INFLUENCE_J   - Jump magnitude
  !!                INFLUENCE_D   - Difference in potential 
  !!                INFLUENCE_Z   - All zeroes
  !!   (in)    integer :: iDP1
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP
  !!             The number of consecutive dipoles to be computed
  !!   (in)    complex :: cPathZ(:)
  !!             Complex coordinates of the control point(s) to be used. For
  !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
  !!             cPathZ(1) is used; for iWhich=INFLUENCE_F, the influence 
  !!             function is computed along the path cPathZ(:)
  !!   (in)    complex :: cOrientation
  !!             Orientation normal vector for iWhich=INFLUENCE_W
  !!   (out)   complex :: cF(1:iNDP,3)
  !!             The returned influence functions.  Indexes 1:iNDP relate
  !!             to dipole indices iDP1:iDP1+iNDP-1, respectively.  Elements
  !!             cF(:,1), cF(:,2), cF(:,3) are the coefficients for the first
  !!             end, center and second end of the dipole, respectively.
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in) :: iWhich,iDP1,iNDP
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  complex (kind=ModAEM_Real),dimension(:,:,:),intent(out) :: cF
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDP2,i,j,iStat
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ1
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ2
  real (kind=ModAEM_Real),dimension(:),allocatable :: rX
  complex (kind=ModAEM_Real) :: cUnitNormal
  complex (kind=ModAEM_Real),dimension(3,1) :: cP,cR,cS,cG,cQ,cJ
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_GetInfluence_IDP: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_GetInfluence_IDP: FDP_Alloc has not been called",io )) return
    if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                    "FDP_GetInfluence_IDP: Bad index range",io ) ) return
    if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                    "FDP_GetInfluence_IDP: Bad index range",io ) ) return             
    if (IO_Assert( (size(cF,1)>=iNDP .and. size(cF,2)>=3 .and. size(cF,3)>=1), &
                    "FDP_GetInfluence_IDP: Invalid result vector",io ) ) return             
  endif

  ! It is assumed that the caller has eliminated the potential for a singularity
  ! by selecting appropriate control points
  iDP2 = iDP1+iNDP-1
  allocate ( cMapZ1(iDP1:iDP2),cMapZ2(iDP1:iDP2),rX(iDP1:iDP2),stat=iStat)
  if (IO_Assert( (iStat==0), "FDP_GetInfluence_IDP: Allocation failed",io)) return
  !
  select case ( iWhich )
    case ( INFLUENCE_P )      ! Potential
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        cP = cIDP_InfluenceP(cMapZ1(iDP1+i-1))
        cF(i,:,1) = cP(:,1)
      end do
    case ( INFLUENCE_D )      ! Potential difference
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      cMapZ2(:) = (cPathZ(2)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        cP = cIDP_InfluenceP(cMapZ1(iDP1+i-1)) - cIDP_InfluenceP(cMapZ2(iDP1+i-1))
        cF(i,:,1) = cP(:,1)
      end do
    case ( INFLUENCE_W )      ! Discharge
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        cUnitNormal = cOrientation/abs(cOrientation)
        cR = cIDP_InfluenceW(cMapZ1(iDP1+i-1))
        ! Rotate the discharge vector by the unit vector in the direction of caOrientation
        ! and reverse the real and imaginary parts for use in influence functions...
        cF(i,:,1) = (cUnitNormal * ( cR(:,1)/fdp%Dipoles(iDP1+i-1)%cZL ))
      end do
    case ( INFLUENCE_F )      ! Integrated Flow
      cF = cZERO
      do j=1,size(cPathZ)-1
        cMapZ1(:) = (cPathZ(j)-fdp%Dipoles(iDP1:iDP2)%cZC ) / fdp%Dipoles(iDP1:iDP2)%cZL
        cMapZ2(:) = (cPathZ(j+1)-fdp%Dipoles(iDP1:iDP2)%cZC ) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i=1,iNDP
          cS = cIDP_InfluenceF(cMapZ1(iDP1+i-1),cMapZ2(iDP1+i-1))
          cF(i,:,1) = cF(i,:,1) + cS(:,1)
        end do
      end do
    case ( INFLUENCE_G )      ! Recharge rate
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        cG = cIDP_InfluenceG(cMapZ1(iDP1+i-1))
        cF(i,:,1) = cG(:,1)
      end do
    case ( INFLUENCE_Q )      ! Extraction rate
      do i=1,iNDP
        cQ = cIDP_InfluenceQ()
        cF(i,:,1) = cQ(:,1)
      end do
    case ( INFLUENCE_J )      ! Jump
      rX(:) = real((cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL)
      do i=1,iNDP
        cJ = cIDP_InfluenceJ(rX(iDP1+i-1))
        cF(i,:,1) = cJ(:,1)
      end do
    case ( INFLUENCE_Z )
      cF = cZERO
  end select
  deallocate (cMapZ1,cMapZ2)

  return
end subroutine FDP_GetInfluence_IDP

subroutine FDP_GetInfluence_ILS(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cF,io)
  !! subroutine FDP_GetInfluence_ILS
  !!
  !! Retrieves arrays of influence functions for use in matrix generation,
  !! using the first-order linesink function, instead of the dipole functions.
  !!
  !! Calling Sequence:
  !!    call FDP_GetInfluence(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cF)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    integer :: iWhich
  !!             The influence function to be computed;  iWhich values are
  !!                INFLUENCE_P   - Complex potential
  !!                INFLUENCE_Q   - Complex discharge
  !!                INFLUENCE_F   - Integrated flux
  !!                INFLUENCE_G   - Areal infiltration
  !!                INFLUENCE_Q   - Extraction rate
  !!                INFLUENCE_D   - Difference in potential 
  !!                INFLUENCE_Z   - All zeroes
  !!   (in)    integer :: iDP1
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP
  !!             The number of consecutive dipoles to be computed
  !!   (in)    complex :: cPathZ(:)
  !!             Complex coordinates of the control point(s) to be used. For
  !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
  !!             cPathZ(1) is used; for iWhich=INFLUENCE_F, the influence 
  !!             function is computed along the path cPathZ(:)
  !!   (in)    complex :: cOrientation
  !!             Orientation normal vector for iWhich=INFLUENCE_W
  !!   (out)   complex :: cF(1:iNDP,1)
  !!             The returned influence functions.  Indexes 1:iNDP relate
  !!             to dipole indices iDP1:iDP1+iNDP-1, respectively.  The 
  !!             influence function is in terms of the sink density of the
  !!             linesink.
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in) :: iWhich,iDP1,iNDP
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  complex (kind=ModAEM_Real),dimension(:,:,:),intent(out) :: cF
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDP2,i,j,iStat
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ1
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ2
  complex (kind=ModAEM_Real) :: cUnitNormal
  complex (kind=ModAEM_Real),dimension(1,1) :: cP,cR,cS,cG,cQ
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_GetInfluence_ILS: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_GetInfluence_ILS: FDP_Alloc has not been called",io )) return
    if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                    "FDP_GetInfluence_ILS: Bad index range",io ) ) return             
    if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and.(iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                    "FDP_GetInfluence_ILS: Bad index range",io ) ) return                    
    if (IO_Assert( (size(cF,1)>=iNDP .and. size(cF,2)>=1 .and. size(cF,3)>=1), &
                    "FDP_GetInfluence_ILS: Invalid result vector",io ) ) return             
  endif

  ! Allocate the Map-Z arrays
  iDP2 = iDP1+iNDP-1
  allocate (cMapZ1(iDP1:iDP2),cMapZ2(iDP1:iDP2),stat=iStat)
  if (IO_Assert( (iStat==0), "FDP_GetInfluence_ILS: Allocation failed",io)) return

  ! It is assumed that the caller has eliminated the potential for a singularity
  ! by selecting appropriate control points
  select case ( iWhich )
    case ( INFLUENCE_P )
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        dip => fdp%Dipoles(iDP1+i-1)
        cP = cILS_InfluenceP(cMapZ1(iDP1+i-1),dip%cZL)
        cF(i,:,1) = cP(:,1)
      end do
    case ( INFLUENCE_D )
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      cMapZ2(:) = (cPathZ(2)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        dip => fdp%Dipoles(iDP1+i-1)
        cP = cILS_InfluenceP(cMapZ1(iDP1+i-1),dip%cZL) - cILS_InfluenceP(cMapZ2(iDP1+i-1),dip%cZL)
        cF(i,:,1) = cP(:,1)
      end do
    case ( INFLUENCE_W )
      cMapZ1(:) = (cPathZ(1)-fdp%Dipoles(iDP1:iDP2)%cZC) / fdp%Dipoles(iDP1:iDP2)%cZL
      do i=1,iNDP
        dip => fdp%Dipoles(iDP1+i-1)
        ! Rotate the discharge vector by the unit vector in the direction of cOrientation
        cUnitNormal = cOrientation/abs(cOrientation)
        cR = cILS_InfluenceW(cMapZ1(iDP1+i-1),dip%cZL)
        cF(i,:,1) = cUnitNormal * ( cR(:,1)/fdp%Dipoles(iDP1+i-1)%cZL )
      end do
    case ( INFLUENCE_F )
      cF = cZERO
      do j=1,size(cPathZ)-1
        cMapZ1(:) = (cPathZ(j)-fdp%Dipoles(iDP1:iDP2)%cZC ) / fdp%Dipoles(iDP1:iDP2)%cZL
        cMapZ2(:) = (cPathZ(j+1)-fdp%Dipoles(iDP1:iDP2)%cZC ) / fdp%Dipoles(iDP1:iDP2)%cZL
        do i=1,iNDP
          dip => fdp%Dipoles(iDP1+i-1)
          cS = cILS_InfluenceF(cMapZ1(iDP1+i-1),cMapZ2(iDP1+i-1),dip%cZL)
          cF(i,:,1) = cF(i,:,1) + cS(:,1)
        end do
      end do
    case ( INFLUENCE_G )
      do i=1,iNDP
        cG = cILS_InfluenceG()
        cF(i,:,1) = cG(:,1)
      end do
    case ( INFLUENCE_Q )
      do i=1,iNDP
        dip => fdp%Dipoles(i)
        cQ = cILS_InfluenceQ(dip%cZL)
        cF(i,:,1) = cQ(:,1)
      end do
    case ( INFLUENCE_Z )
      cF = cZERO
  end select

  deallocate (cMapZ1,cMapZ2)

  return
end subroutine FDP_GetInfluence_ILS

function cFDP_Potential(fdp,cZ,io,iDP1,iNDP) result(cOmega)
  !! complex function cFDP_Potential
  !!
  !! Computes the complex potential due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    cOmega = cFDP_Potential(fdp,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!   (in)    integer :: iDP1 [OPTIONAL]
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP [OPTIONAL]
  !!             The number of consecutive dipoles to be computed
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP1,iNDP
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cOmega
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd,iStat
  complex (kind=ModAEM_Real),dimension(3,1) :: cP
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iDP1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return             
    endif
    iStart = iDP1
    if ( present(iNDP) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return            
      endif
      iEnd = iStart+iNDP-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fdp%iCount
  end if

  ! Sum up the contribution of all dipoles
  allocate (cMapZ(iStart:iEnd),stat=iStat)
  if (IO_Assert( (iStat==0), "cFDP_Potential: Allocation failed",io)) return
  cOmega = cZERO
  cMapZ = (cZ - fdp%Dipoles%cZC) / fdp%Dipoles%cZL
  do i=iStart,iEnd
    dip => fdp%Dipoles(i)
    cP = cIDP_InfluenceP(cMapZ(i))
    cOmega = cOmega + sum(dip%cRho*cP(:,1))
  end do
  deallocate (cMapZ)

  return
end function cFDP_Potential

function cFDP_Discharge(fdp,cZ,io,iDP1,iNDP) result(cQ)
  !! complex function cFDP_Discharge
  !!
  !! Computes the complex discharge due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    cQ = cFDP_Discharge(fdp,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!   (in)    integer :: iDP1 [OPTIONAL]
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP [OPTIONAL]
  !!             The number of consecutive dipoles to be computed
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP1,iNDP
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cQ
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd,iStat
  complex (kind=ModAEM_Real),dimension(3,1) :: cP
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iDP1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return              
    endif
    iStart = iDP1
    if ( present(iNDP) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return            
      endif
      iEnd = iStart+iNDP-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fdp%iCount
  end if

  ! Sum up the contribution of all dipoles
  allocate (cMapZ(iStart:iEnd),stat=iStat)
  if (IO_Assert( (iStat==0), "cFDP_Potential: Allocation failed",io)) return
  cQ = cZERO
  cMapZ = (cZ - fdp%Dipoles%cZC) / fdp%Dipoles%cZL
  do i=iStart,iEnd
    dip => fdp%Dipoles(i)
    cP = cIDP_InfluenceW(cMapZ(i))
    cQ = cQ + conjg(sum(dip%cRho*cP(:,1)) / dip%cZL)
  end do
  deallocate (cMapZ)

  return
end function cFDP_Discharge

function rFDP_Flow(fdp,cPathZ,io,iDP1,iNDP) result(rFlow)
  !! real function cFDP_Flow
  !!
  !! Computes the integrated flow due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    rFlow = fFDP_Flow(fdp,cPathZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    complex :: cPathZ(:)
  !!             The path across which the flow is desired
  !!   (in)    integer :: iDP1 [OPTIONAL]
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP [OPTIONAL]
  !!             The number of consecutive dipoles to be computed
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP1,iNDP
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rFlow
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,j,iStart,iEnd,iStat
  complex (kind=ModAEM_Real),dimension(3,1) :: cS
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ1,cMapZ2
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip


  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iDP1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return             
    endif
    iStart = iDP1
    if ( present(iNDP) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                      "FDP_GetInfluence_ILS: Bad index range",io ) ) return            
      endif
      iEnd = iStart+iNDP-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fdp%iCount
  end if

  ! Sum up the contribution of all dipoles
  allocate (cMapZ1(iStart:iEnd),cMapZ2(iStart:iEnd),stat=iStat)
  if (IO_Assert( (iStat==0), "cFDP_Potential: Allocation failed",io)) return
  rFlow = rZERO
  do j=1,ubound(cPathZ,1)-1
    cMapZ1 = (cPathZ(j)-fdp%Dipoles%cZC ) / fdp%Dipoles%cZL
    cMapZ2 = (cPathZ(j+1)-fdp%Dipoles%cZC ) / fdp%Dipoles%cZL
    do i=iStart,iEnd
      dip => fdp%Dipoles(i)
      cS = cIDP_InfluenceF(cMapZ1(i),cMapZ2(i))
      rFlow = rFlow + sum(real(dip%cRho * cS(:,1)))
    end do
  end do
  deallocate (cMapZ1,cMapZ2)

  return
end function rFDP_Flow

function rFDP_Recharge(fdp,cZ,io,iDP1,iNDP) result(rGamma)
  !! complex function cFDP_Recharge
  !!
  !! Computes the recharge rate due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    rGamma = cFDP_Recharge(fdp,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    complex :: cZ
  !!             The point at which to determine the recharge
  !!   (in)    integer :: iDP1 [OPTIONAL]
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP [OPTIONAL]
  !!             The number of consecutive dipoles to be computed
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP1,iNDP
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: rGamma
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd,iStat
  complex (kind=ModAEM_Real),dimension(3,1) :: cG
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iDP1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                     "FDP_Recharge: Bad index range",io) ) return
    endif
    iStart = iDP1
    if ( present(iNDP) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                      "FDP_Recharge: Bad index range",io ) ) return            
      endif
      iEnd = iStart+iNDP-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fdp%iCount
  end if

  ! Sum up the contribution of all dipoles
  allocate (cMapZ(iStart:iEnd),stat=iStat)
  if (IO_Assert( (iStat==0), "FDP_Recharge: Allocation failed",io)) return
  rGamma = rZERO
  do i=iStart,iEnd
    dip => fdp%Dipoles(i)
    cG = cIDP_InfluenceG(cMapZ(i))
    rGamma = rGamma + sum(dip%cRho*cG(:,1))
  end do
  deallocate (cMapZ)

  return
end function rFDP_Recharge

function rFDP_Extraction(fdp,io,iDP1,iNDP) result(rQ)
  !! real function rFDP_Extraction
  !!
  !! Computes the recharge rate due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    rGamma = rFDP_Extraction(fdp,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    integer :: iDP1 [OPTIONAL]
  !!             The index for the first dipole to be used
  !!   (in)    integer :: iNDP [OPTIONAL]
  !!             The number of consecutive dipoles to be computed
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP1,iNDP
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: rQ
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd,iStat
  complex (kind=ModAEM_Real),dimension(3,1) :: cQ
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cMapZ
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iDP1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iDP1>=lbound(fdp%Dipoles,1) .and. iDP1<=ubound(fdp%Dipoles,1)), &
                      "FDP_Recharge: Bad index range",io ) ) return
    endif
    iStart = iDP1
    if ( present(iNDP) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iDP1+iNDP-1)>=lbound(fdp%Dipoles,1) .and. (iDP1+iNDP-1)<=ubound(fdp%Dipoles,1)), &
                      "FDP_Recharge: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNDP-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fdp%iCount
  end if

  ! Sum up the contribution of all dipoles
  rQ = rZERO
  do i=iStart,iEnd
    dip => fdp%Dipoles(i)
    cQ = cIDP_InfluenceQ()
    rQ = rQ + sum(dip%cRho*cQ(:,1))
  end do

  return
end function rFDP_Extraction

function rFDP_PotentialJump(fdp,io,iDP,cZ) result(rJump)
  !! real function rFDP_Jump
  !!
  !! Computes the potential jump at the position cZ in dipole iDP of FDP_COLLECTION fdp
  !!
  !! Calling Sequence:
  !!    rJump = rFDP_Jump(fdp,iDP,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!   (in)    integer :: iDP [OPTIONAL]
  !!             The index for the dipole to be used
  !!   (in)    complex :: cZ
  !!             The location to be investigated
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in),optional :: iDP
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: rJump
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  complex (kind=ModAEM_Real),dimension(3,1) :: cJ
  complex (kind=ModAEM_Real) :: cMapZ
  complex (kind=ModAEM_Real) :: cA
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Jump: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Jump: FDP_Alloc has not been called",io )) return
    if (IO_Assert( (iDP > 0 .and. iDP <= fdp%iCount), "FDP_Jump: Illegal index",io )) return
  endif

  ! Sum up the contribution of all dipoles
  dip => fdp%Dipoles(iDP)
  cMapZ = (cZ - dip%cZC) / dip%cZL
  cJ = cIDP_InfluenceJ(real(cMapZ))
  rJump = real(sum(dip%cRho*cJ(:,1)))

  return
end function rFDP_PotentialJump

subroutine FDP_Report(fdp,io)
  !! subroutine FDP_Report
  !!
  !! Writes a report of all dipole information to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call FDP_Report(fdp)
  !! 
  !! Arguments:
  !!    (in)   type ( FDP_COLLECTION ),pointer :: fdp
  !!             The FDP_COLLECTION to use
  !!
  ! [ ARGUMENTS ]
  type (FDP_COLLECTION),pointer :: fdp
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type (FDP_DIPOLE),pointer :: dip

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fdp)), "FDP_Update: FDP_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp%Dipoles)), "FDP_Update: FDP_Alloc has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module FDP'//'Line-dipole function information')" )

  ! OK. Write the report...
  write ( unit=kIOOutputLU, &
          fmt="(""Number of dipole/doublet functions: "",i5/)" &
        ) fdp%iCount

  if ( fdp%iCount > 0 ) then
    write ( unit=kIOOutputLU, &
            fmt="("" Geometric parameters""/""Index"",t10,""ZC"",t40,""ZL"")" &
          )
    do i=1,fdp%iCount
      dip => fdp%Dipoles(i)
      write ( unit=kIOOutputLU, &
              fmt="(i5,t10,4(d12.5,3x))" &
            ) i,dip%cZc,dip%cZl
    end do

    write ( unit=kIOOutputLU, &
            fmt="(/""Strength parameters""/""Index"",t10,""Rho1"",t40," // &
                " ""Mu*"",t70,""Rho3"")" &
          )
    do i=1,fdp%iCount
      dip => fdp%Dipoles(i)
      write ( unit=kIOOutputLU, &
              fmt="(i5,t10,6(d12.5,3x))" & 
            ) i,dip%cRho
    end do
  end if

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine FDP_Report

end module f_dipole
