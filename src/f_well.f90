module f_well

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


  !! module f_well (fwl)
  !! 
  !! Written by Victor A.Kelson
  !! 
  !! Revision History:
  !!   1.0.0   22 March 1999
  !!           First "source code release" version, adapted from
  !!           previous testing version.
  !!
  !! Module of data structures and functions for 2-D, steady-state
  !! wells.  
  !!
  !! This module encapsulates containers of well functions
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
  !!   i_well     --  Influence function module for wells
  !!
use u_constants
use u_io
use i_well

implicit none

public

  type,public :: FWL_WELL
    !! type FWL_WELL
    !!
    !! PUBLIC type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: rDischarge
    !!     Discharge of the well; positive indicates a well withdrawal
    !!     and negative and injection rate.
    !! 
    complex (kind=ModAEM_Real) :: cZC
    real (kind=ModAEM_Real) :: rRadius
    real (kind=ModAEM_Real) :: rDischarge
  end type FWL_WELL

  type,public :: FWL_COLLECTION
    !! type FWL_COLLECTION
    !!
    !! PUBLIC type that holds the well entries for a layer
    !!
    !! Members:
    !!   type (fwl_WELL),pointer :: Wells(:)
    !!     FWL_WELL structures that hold well information
    !!   integer :: iCount
    !!     Number of FWL_WELL structures currently in use
    !! 
    type ( FWL_WELL ),dimension(:),pointer :: Wells
    integer (kind=ModAEM_Integer) :: iCount
  end type

contains

function FWL_Create(iMax,io) result (fwl)
  !! function FWL_Create
  !!
  !! Creates a new FWL_COLLECTION object
  !!
  !! Calling Sequence:
  !!    fwl => FWL_Create(iMax)
  !!
  !! Arguments:
  !!    (in)      integer :: iMax
  !!                The maximum number of wells to be stored
  !!
  ! [ ARGUMENTS ]
  integer (kind=ModAEM_Integer),intent(in) :: iMax
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  type ( FWL_COLLECTION ),pointer :: fwl
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (fwl,stat=iStat)
  if ( IO_Assert( (iStat == 0), "FWL_Create: allocation failed",io )) return

  allocate (fwl%Wells(iMax),stat=iStat)
  if ( IO_Assert( (iStat == 0), "FWL_Create: allocation failed",io )) return
  fwl%iCount = 0

  return
end function FWL_Create

subroutine FWL_Alloc(fwl,iMax,io) 
  !! subroutine FWL_Alloc
  !!
  !! Dimensions the internal buffers for iaMax wells in layer iL
  !!
  !! Calling Sequence:
  !!    call FWL_Alloc(fwl,iMax)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             The FWL_COLLECTION object to be used
  !!   (in)    integer :: iMax
  !!             The maximum number of wells in fwl
  !!
  !! Note: If allocation fails, causes a fatal error
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  integer (kind=ModAEM_Integer),intent(in) :: iMax
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Alloc: FWL_New has not been called",io )) return
    if (IO_Assert( (.not. associated(fwl%Wells)), "FWL_Alloc: Wells already allocated",io )) return
  endif

  ! Now, allocate space for the specified layer and initialize
  allocate (fwl%Wells(iMax),stat=iStat)
  if (IO_Assert( (iStat==0), "FWL_Alloc: Allocation failed",io)) return
  fwl%Wells = FWL_WELL(cZERO,rZERO,rZERO)
  fwl%iCount = 0

  return
end subroutine FWL_Alloc

subroutine FWL_New(fwl,cZc,rDischarge,rRadius,iRV,io)
  !! subroutine FWL_New
  !!
  !! Makes a new well entry. On call, the geometry and discharge of 
  !! the well are provided. The internal well structures are then
  !! set up. Returns iRV=the index into the well table on success or
  !! iRV<0 on failure.
  !!
  !! Calling Sequence:
  !!    call FWL_New(fwl,cZC,rDischarge,iRV)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    complex :: cZC
  !!             Complex coordinates of the center of the well
  !!   (in)    real :: rDischarge
  !!             The discharge of the well
  !!   (out)   integer :: iRV
  !!             On success, the index in fwl%Wells used
  !!
  !! Note: On failure, forces a fatal error
  !!
  ! [ ARGUMENTS ]
  type ( FWL_COLLECTION ),pointer :: fwl
  complex (kind=ModAEM_Real),intent(in) :: cZc
  real (kind=ModAEM_Real), intent(in) :: rRadius
  real (kind=ModAEM_Real), intent(in) :: rDischarge
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer),intent(out) :: iRV
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), &
                    "FWL_New: FWL_Create has not been called",io ) ) return
    if (IO_Assert( (associated(fwl%Wells)), &
                    "FWL_New: FWL_Alloc has not been called",io ) ) return
    if (IO_Assert( (fwl%iCount < size(fwl%Wells)), &
                    "FWL_New: Space exhausted",io ) ) return
  endif

  fwl%iCount = fwl%iCount + 1
  wel => fwl%Wells(fwl%iCount)
  wel = FWL_WELL(cZc,rRadius,rDischarge)
  iRV = fwl%iCount

  return
end subroutine FWL_New

subroutine FWL_Update(fwl,iWL,rDischarge,io) 
  !! subroutine FWL_Update
  !!
  !! Updates the discharge for a well entry.
  !!
  !! Calling Sequence:
  !!    call FWL_Update(fwl,iWL,rDischarge)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    integer :: iWL
  !!             The index for the well in fwl
  !!   (in)    real :: rDischarge
  !!             The discharge of the well
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  integer (kind=ModAEM_Integer),intent(in) :: iWL
  real (kind=ModAEM_Real),intent(in) :: rDischarge
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
    if (IO_Assert( (iWL <= size(fwl%Wells)), "FWL_Update: Space exhausted",io )) return
  endif

  wel => fwl%Wells(iWL)
  wel%rDischarge = rDischarge

  return
end subroutine FWL_Update

subroutine FWL_GetInfluence(fwl,iWhich,iWL1,iNWL,cPathZ,cOrientation,cF,io)
  !! subroutine FWL_GetInfluence
  !!
  !! Retrieves arrays of influence functions for use in matrix generation
  !!
  !! Calling Sequence:
  !!    call FWL_GetInfluence(fwl,iWhich,iWL1,iNWL,cPathZ,cOrientation,cF)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    integer :: iWhich
  !!             The influence function to be computed;  iWhich values are
  !!                INFLUENCE_P   - Complex potential
  !!                kInfluenceQ   - Complex discharge
  !!                kInfluenceF   - Integrated flux
  !!                kInfluenceG   - Areal infiltration
  !!                kInfluenceQ   - Extraction rate
  !!                kInfluenceD   - Difference in potential
  !!                kInfluenceZ   - All zeroes
  !!   (in)    integer :: iWL1
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL
  !!             The number of consecutive wells to be computed
  !!   (in)    complex :: cPathZ(:)
  !!             Complex coordinates of the control point(s) to be used. For
  !!             iWhich = (INFLUENCE_P, INFLUENCE_W and INFLUENCE_G) only
  !!             cPathZ(1) is used; for iWhich=INFLUENCE_F, the influence 
  !!             function is computed along the path cPathZ(:)
  !!   (in)    complex :: cOrientation
  !!             Orientation normal vector for iWhich=INFLUENCE_W
  !!   (out)   complex :: cF(1:iNWL,3)
  !!             The returned influence functions.  Indexes 1:iNWL relate
  !!             to well indices iWL1:iWL1+iNWL-1, respectively.  
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  integer (kind=ModAEM_Integer),intent(in) :: iWhich,iWL1,iNWL
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  complex (kind=ModAEM_Real),dimension(:,:,:),intent(out) :: cF
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iWL2,i,j
  real (kind=ModAEM_Real),dimension(1,1) :: rF
  complex (kind=ModAEM_Real),dimension(1,1) :: cW,cUnitNormal
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_GetInfluence: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_GetInfluence: FWL_Alloc has not been called",io )) return
    if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                    "FWL_GetInfluence: Bad index range",io ) ) return
    if  (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                    "FWL_GetInfluence: Bad index range",io ) ) return
    if (IO_Assert( (size(cF,1)>=iNWL .and. size(cF,2)>=1 ), "FWL_GetInfluence: Invalid result vector",io )) return
  endif

  select case ( iWhich )
    case ( INFLUENCE_P )
      do i=1,iNWL
        wel => fwl%Wells(i)
        cF(i,:,:) = cIWL_InfluenceP(cPathZ(1),wel%cZC)
      end do
    case ( INFLUENCE_D )
      do i=1,iNWL
        wel => fwl%Wells(i)
        cF(i,:,:) = cIWL_InfluenceP(cPathZ(1),wel%cZC) - cIWL_InfluenceP(cPathZ(2),wel%cZC)
      end do
    case ( INFLUENCE_W )
      do i=1,iNWL
        wel => fwl%Wells(i)
        cUnitNormal = cOrientation/abs(cOrientation)
        cW = cIWL_InfluenceW(cPathZ(1),wel%cZC)
        cF(i,:,:) = cUnitNormal * conjg( cW )
      end do
    case ( INFLUENCE_F )
      do i=1,iNWL
        wel => fwl%Wells(i)
        rF = rZERO
        do j=1,ubound(cPathZ,1)-1
          rF = rF + cIWL_InfluenceF(cPathZ(j),cPathZ(j+1),wel%cZC)
        end do
        cF(i,:,:) = -cmplx(rZERO,rF,ModAEM_Real)
      end do
    case ( INFLUENCE_G )
       cF(i,:,:) = rZERO
    case ( INFLUENCE_Q )
      do i=1,iNWL
        cF(i,:,:) = rONE
      end do
    case ( INFLUENCE_Z )
      cF = rZERO
  end select

  return
end subroutine FWL_GetInfluence

function cFWL_Potential(fwl,cZ,io,iWL1,iNWL) result(cOmega)
  !! complex function cFWL_Potential
  !!
  !! Computes the complex potential due to the current set of wells. 
  !!
  !! Calling Sequence:
  !!    cOmega = cFWL_Potential(fwl,cZ,iWL1,iNWL)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!   (in)    integer :: iWL1 [OPTIONAL]
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL [OPTIONAL]
  !!             The number of consecutive wells to be computed
  !! Note:
  !!   If iWL1 is not provided, all wells will be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iWL1,iNWL
  complex (kind=ModAEM_Real) :: cOmega
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd
  complex (kind=ModAEM_Real),dimension(1,1) :: cP
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iWL1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
    endif
    iStart = iWL1
    if ( present(iNWL) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNWL-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fwl%iCount
  end if

  ! Sum up the contributions of all wells
  cOmega = cZERO
  do i=iStart,iEnd
    wel => fwl%Wells(i)
    cP = cIWL_InfluenceP(cZ,wel%cZC)
    cOmega = cOmega + wel%rDischarge*cP(1,1)
  end do
  
  return
end function cFWL_Potential

function cFWL_Discharge(fwl,cZ,io,iWL1,iNWL) result(cQ)
  !! complex function cFWL_Discharge
  !!
  !! Computes the complex discharge due to the current set of wells. 
  !!
  !! Calling Sequence:
  !!    cQ = cFWL_Discharge(fwl,cZ,io,iWL1,iNWL)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!   (in)    integer :: iWL1 [OPTIONAL]
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL [OPTIONAL]
  !!             The number of consecutive wells to be computed
  !! Note:
  !!   If iWL1 is not provided, all wells will be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iWL1,iNWL
  complex (kind=ModAEM_Real) :: cQ
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd
  complex (kind=ModAEM_Real),dimension(1,1) :: cW
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iWL1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
    endif
    iStart = iWL1
    if ( present(iNWL) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNWL-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fwl%iCount
  end if

  ! Sum up the contributions of all wells
  cQ = cZERO
  do i=iStart,iEnd
    wel => fwl%Wells(i)
    cW = cIWL_InfluenceW(cZ,wel%cZC)
    cQ = cQ + conjg(wel%rDischarge * cW(1,1))
  end do

  return
end function cFWL_Discharge

function rFWL_Flow(fwl,cPathZ,io,iWL1,iNWL) result(rFlow)
  !! complex function cFWL_Discharge
  !!
  !! Computes the complex discharge due to the current set of wells. 
  !!
  !! Calling Sequence:
  !!    rFlow = fFWL_Flow(fwl,cPathZ,iWL1,iNWL)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    complex :: cPathZ(:)
  !!             The path across which the flow is desired
  !!   (in)    integer :: iWL1 [OPTIONAL]
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL [OPTIONAL]
  !!             The number of consecutive wells to be computed
  !! Note:
  !!   If iWL1 is not provided, all wells will be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iWL1,iNWL
  real (kind=ModAEM_Real) :: rFlow
  type ( IO_STATUS ),pointer :: io
  ! Locals
  integer (kind=ModAEM_Integer) :: i,j,iStart,iEnd
  real (kind=ModAEM_Real),dimension(1,1) :: rS
  complex (kind=ModAEM_real) :: cZC,cZL,cZMap
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iWL1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                     "FWL_GetInfluence_ILS: Bad index range",io ) ) return
    endif
    iStart = iWL1
    if ( present(iNWL) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNWL-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fwl%iCount
  end if

  rFlow = rZERO
  ! Sum up the contributions of all wells
  do i=iStart,iEnd
    wel => fwl%Wells(i)
    do j=1,ubound(cPathZ,1)-1
      ! SPECIAL CASE: Does the line segment pass within the well bore?
      cZC = rHALF * ( cPathZ(j) + cPathZ(j+1) )
      cZL = rHALF * ( cPathZ(j) - cPathZ(j+1) )
      cZMap = ( wel%cZC - cZC ) / cZL
      if ( wel%rRadius >= abs(aimag(cZMap))*abs(cZL) ) then
        rS = rZERO
      else
        rS = cIWL_InfluenceF(cPathZ(j),cPathZ(j+1),wel%cZC)
      endif
      rFlow = rFlow + wel%rDischarge*rS(1,1)
    end do
  end do

  return
end function rFWL_Flow

function rFWL_Recharge(fwl,cZ,io,iWL1,iNWL) result(rG)
  !! complex function rFWL_Recharge
  !!
  !! Computes the complex potential due to the current set of wells. 
  !!
  !! Calling Sequence:
  !!    rG = rFWL_Recharge(fwl,cZ,iWL1,iNWL)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!   (in)    integer :: iWL1 [OPTIONAL]
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL [OPTIONAL]
  !!             The number of consecutive wells to be computed
  !! Note:
  !!   If iWL1 is not provided, all wells will be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in),optional :: iWL1,iNWL
  real (kind=ModAEM_Real) :: rG
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd
  complex (kind=ModAEM_Real),dimension(1,1) :: cG
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iWL1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
    endif
    iStart = iWL1
    if ( present(iNWL) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNWL-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fwl%iCount
  end if

  ! Sum up the contributions of all wells
  rG = cZERO
  do i=iStart,iEnd
    wel => fwl%Wells(i)
    cG = cIWL_InfluenceG(cZ,wel%cZC)
    rG = rG + wel%rDischarge*cG(1,1)
  end do
  
  return
end function rFWL_Recharge

function rFWL_Extraction(fwl,io,iWL1,iNWL) result(rQ)
  !! complex function rFWL_Extraction
  !!
  !! Computes the complex potential due to the current set of wells. 
  !!
  !! Calling Sequence:
  !!    rQ = rFWL_Extraction(fwl,cZ,iWL1,iNWL)
  !!
  !! Arguments:
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!   (in)    integer :: iWL1 [OPTIONAL]
  !!             The index for the first well to be used
  !!   (in)    integer :: iNWL [OPTIONAL]
  !!             The number of consecutive wells to be computed
  !! Note:
  !!   If iWL1 is not provided, all wells will be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  integer (kind=ModAEM_Integer),intent(in),optional :: iWL1,iNWL
  real (kind=ModAEM_Real) :: rQ
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,iStart,iEnd
  complex (kind=ModAEM_Real),dimension(1,1) :: cQ
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  ! Set the range, using the optional arguments
  if ( present(iWL1) ) then 
    if ( io%lDebug ) then
      if (IO_Assert( (iWL1>=lbound(fwl%Wells,1) .and. iWL1<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
    endif
    iStart = iWL1
    if ( present(iNWL) ) then 
      if ( io%lDebug ) then
        if (IO_Assert( ((iWL1+iNWL-1)>=lbound(fwl%Wells,1) .and. (iWL1+iNWL-1)<=ubound(fwl%Wells,1)), &
                      "FWL_GetInfluence_ILS: Bad index range",io ) ) return
      endif
      iEnd = iStart+iNWL-1
    else
      iEnd = iStart
    end if
  else
    iStart = 1
    iEnd = fwl%iCount
  end if

  ! Sum up the contributions of all wells
  rQ = cZERO
  do i=iStart,iEnd
    wel => fwl%Wells(i)
    cQ = cIWL_InfluenceQ()
    rQ = rQ + wel%rDischarge*cQ(1,1)
  end do
  
  return
end function rFWL_Extraction

subroutine FWL_Report(fwl,io)
  !! subroutine FWL_Report
  !!
  !! Writes a report of all well information to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call FWL_Report(fwl)
  !!
  !! Arguments:
  !!    (in)   type ( FWL_COLLECTION ),pointer :: fwl
  !!             The FWL_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type (fwl_COLLECTION),pointer :: fwl
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type (fwl_WELL),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(fwl)), "FWL_Update: FWL_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl%Wells)), "FWL_Update: FWL_Alloc has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module FWL'//'Well function information')" )

  write (unit=kIOOutputLU,fmt="(""Number of well functions: "",i5/)") fwl%iCount

  if ( fwl%iCount > 0 ) then
    write (unit=kIOOutputLU,fmt="(""Geometric parameters""/""Index"",t10,""ZC"",t40,""Radius"")")
    do i=1,fwl%iCount
      wel => fwl%Wells(i)
      write(unit=kIOOutputLU,fmt="(i5,t10,3(d12.5,3x))") i,wel%cZc,wel%rRadius
    end do

    write (unit=kIOOutputLU,fmt="(""Strength parameters""/""Index"",t10,""Discharge"")")
    do i=1,fwl%iCount
      wel => fwl%Wells(i)
      write(unit=kIOOutputLU,fmt="(i5,t10,(d12.5,3x))") i,wel%rDischarge
    end do
  end if

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine FWL_Report

end module f_well
