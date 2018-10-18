module m_wl1

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

  !! module m_wl1
  !!
  !! Element module for 2-D head specified wells
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!
  
use u_constants
use u_io
use f_well
use u_matrix
use m_aqu

implicit none

public

  type :: WL1_WELL
    !! type WL1_WELL
    !!
    !! Type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: cZHead
    !!     The point at which the head is specified
    !!   real :: rHead
    !!     The given head
    !!   real :: rHeadCorrection
    !!     The correction on the head as computed by another module (for derived packages)
    !!   integer :: iID
    !!     Identification label for the well (for interaction with e.g. GUIs)
    !!   integer :: iFWLIndex
    !!     Index for the well entry in the FWL module
    !!   real :: rStrength
    !!     Computed strength of the element
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rRadius
    complex (kind=ModAEM_Real) ::cZHead
    real (kind=ModAEM_Real) :: rHead
    real (kind=ModAEM_Real) :: rHeadCorrection
    integer (kind=ModAEM_Integer) :: iID
    integer (kind=ModAEM_Integer) :: iFWLIndex
    real (kind=ModAEM_Real) :: rStrength
  end type WL1_WELL

  type :: WL1_COLLECTION
    !! type WL1_COLLECTION
    !!
    !! Type that holds head-specified wells in a layer
    !!
    !! Members:
    !!   type ( WL1_WELL ),dimension(:),pointer :: Wells
    !!     Array of WL1_WELL objects for the layer; dimensioned for the maximum
    !!     number of wells according to the input file (see WL1_Read)
    !!   integer :: iCount
    !!     The number of strings actually in use
    !! 
    type ( WL1_WELL),dimension(:),pointer :: Wells
    integer (kind=ModAEM_Integer) :: iCount
    integer (kind=ModAEM_Integer) :: iRegenerate
  end type WL1_COLLECTION

  ! Module flags for matrix generator routines
  integer (kind=ModAEM_Integer),private,parameter :: WL1Vertex=1

contains

function WL1_Create(io) result (wl1)
  !! function WL1_Create
  !!
  !! Creates a new WL1_COLLECTION object
  !!
  !! Calling Sequence:
  !!    WL1 => WL1_Create()
  !!
  !! Arguments:
  !!
  ! [ ARGUMENTS ]
  ! [ RETURN VALUE ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (wl1,stat=iStat)
  if (IO_Assert( (iStat == 0), "WL1_Create: allocation failed",io )) return
  nullify(wl1%Wells)
  wl1%iCount = 0
  wl1%iRegenerate = 1
 
  return
end function WL1_Create

subroutine WL1_Alloc(wl1,iNWL,io)
  !! Subroutine WL1_Alloc
  !! 
  !! Allocates Wells for the WL1_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call WL1_Alloc(wl1, iNStr)
  !!
  !! Arguments:
  !!    (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!              The WL1_COLLECTION object to be used
  !!    (in)    integer :: iNWL
  !!              The number of wells to make space for
  !!
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  integer (kind=ModAEM_Integer),intent(in) :: iNWL
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (wl1%Wells(iNWL),stat=iStat)
  if (IO_Assert( (iStat == 0), "WL1_Alloc: allocation failed",io )) return

end subroutine WL1_Alloc


subroutine WL1_Destroy(wl1,io)
  !! subroutine WL1_Destroy
  !!
  !! Frees memory allocated for wl1 Wells and the WL1 Collection object
  !!
  !! Calling Sequence:
  !!     call WL1_Destroy(wl1)
  !!
  !! Arguments:
  !!  type ( WL1_COLLECTION ),pointer :: wl1
  !!              Pointer to the WL1_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( wl1_COLLECTION ),pointer :: wl1
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat


  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Destroy: WL1_Create has not been called",io )) return
  endif

  if ( associated(wl1%Wells) ) then
    deallocate (wl1%Wells,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                     "WL1_Destroy: deallocation of Wells failed",io )) return
  end if
  deallocate (wl1,stat=iStat)
  if (IO_Assert( (iStat == 0), "WL1_Destroy: deallocation failed",io )) return
  
  return
end subroutine WL1_Destroy


subroutine WL1_New(wl1,Well,io) 
  !! function WL1_New
  !!
  !! Adds a new WL1_WELL object to the WL1_COLLECTION 'wl1'
  !!
  !! Calling Sequence:
  !!    call WL1_New(wl1,Well)
  !!
  !! Arguments:
  !!    (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!              The WL1_COLLECTION object to be used
  !!    (in)    type ( WL1_WELL), pointer :: Well
  !!              The new well object
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( WL1_WELL ) :: Well
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_New: WL1_Create has not been called",io )) return
  endif

  if (IO_Assert( (wl1%iCount < size(wl1%Wells)), &
                  "WL1_New: Space exhausted",io )) return

  wl1%iCount = wl1%iCount + 1
  wl1%Wells(wl1%iCount) = Well

  return
end subroutine WL1_New

function iWL1_GetInfo(wl1,iOption,io) result(iValue)
  !! function WL1_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iWL1_GetInfo(wl1,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             WL1_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_GetInfo: WL1_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FWL )
      iValue = wl1%iCount
    case ( SIZE_EQUATIONS )
      iValue = wl1%iCount
    case ( SIZE_UNKNOWNS )
      iValue = wl1%iCount
    case ( INFO_REGENERATE )
      iValue = wl1%iRegenerate
    case default
      iValue = wl1%iCount
  end select

  return
end function iWL1_GetInfo

subroutine WL1_SetupFunctions(wl1,fwl,io)
  !! subroutine WL1_SetupFunctions
  !!
  !! This routine sets up the functions in f_well and f_dipole for the line-sinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call WL1_Setup(wl1,fwl,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i, iWell, iWL
  real (kind=ModAEM_Real) :: rStrength,rDisch,rHead1,rHead2,rHead
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  complex (kind=ModAEM_Real),dimension(3) :: cCPResult
  type ( WL1_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Setup: WL1_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "WL1_Setup: Illegal FWL_COLLECTION object",io )) return
  endif

  do iWell = 1,wl1%iCount
    wel => wl1%Wells(iWell)

    call FWL_New(fwl,wel%cZ,rZERO,wel%rRadius,iWL,io)
    if ( io%lError ) return
    wel%iFWLIndex = iWL
  end do

  return
end subroutine WL1_SetupFunctions

subroutine WL1_SetupMatrix(wl1,aqu,mat,io)
  !! subroutine WL1_SetupMatrix
  !!
  !! This routine sets up the functions in f_well and f_dipole for the line-sinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call WL1_Setup(wl1,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    type ( AQU_COLLECTION ),pointer
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i, iWell, iWL
  real (kind=ModAEM_Real) :: rStrength,rDisch,rHead1,rHead2,rHead
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  complex (kind=ModAEM_Real),dimension(3) :: cCPResult
  type ( WL1_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Setup: WL1_Create has not been called",io )) return
    if (IO_Assert( (associated(aqu)), &
                    "WL1_Setup: Illegal AQU_COLLECTION object",io )) return
    if (IO_Assert( (associated(mat)), &
                    "WL1_Setup: Illegal MAT_MATRIX object",io )) return
  endif

  ! Build matrix generator entries for all segments
  do iWell = 1,wl1%iCount
    ! Set up the unknown variables
    wel => wl1%Wells(iWell) 
    call MAT_CreateVariable(mat,ELEM_WL1,iWell,0,WL1Vertex,io)
    if ( io%lError ) return

    call MAT_CreateEquation(mat,(/ wel%cZHead /),EQN_HEAD,ELEM_WL1, &
                            0,iWell,0, &
                            rAQU_HeadToPotential(aqu,wel%rHead,wel%czHead,io),cZERO,io)
    if ( io%lError ) return
  end do

  return
end subroutine WL1_SetupMatrix

subroutine WL1_Prepare(wl1,io) 
  !! subroutine WL1_Prepare
  !! 
  !! Prepares the module for a new iteration
  !!
  !! Do-nothing for m_wl1
  !!
  !! Calling Sequence:
  !!    call WL1_Setup(wl1,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( IO_STATUS ),pointer :: io

  return
end subroutine WL1_Prepare

subroutine WL1_ComputeCoefficients(wl1,fwl,cPathZ,iEqType,iElementType,iElementString, &
                                   iElementVertex,iElementFlag,cOrientation,rARow,io)
  !! subroutine WL1_ComputeCoefficients
  !!
  !! Computes a row of matrix coefficients (with no corrections) for the WL1 
  !! elements in layer iL.
  !!
  !! Calling Sequence:
  !!    call WL1_ComputeCoefficients(wl1,cPathZ,iEqType,cOrientation,rRow)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!   (in)    complex :: cPathZ(:)
  !!             The control point (or control path) to be used
  !!   (in)    integer :: iEqType
  !!             The equation type
  !!   (in)    integer :: iElementType
  !!             The element type that created this equation
  !!   (in)    integer :: iElementString
  !!             The element string corresponding to this equation
  !!   (in)    integer :: iElementVertex
  !!             The element vertex corresponding to this equation
  !!   (in)    integer :: iElementFlag
  !!             The element flag (if any) for this equation
  !!   (in)    complex :: cOrientation
  !!             Orientation unit vector (for discharge-based equations)
  !!   (out)   real :: rARow(:)
  !!             The output row of coefficients (to be concatenated with
  !!             row portions for other element modules.
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( FWL_COLLECTION ),pointer :: fwl
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  integer (kind=ModAEM_Integer), intent(in) :: iEqType
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real),dimension(:),intent(out) :: rARow
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat,iCol,iVtx,iWL1,iNWL,iWhich
  complex (kind=ModAEM_Real),dimension(:,:,:),allocatable :: cWLF
  type ( WL1_WELL ),pointer :: wel

  if ( io%lDebug ) then 
    if (IO_Assert( (associated(wl1)), &
                    "WL1_ComputeCoefficients: WL1_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "WL1_Setup: Illegal FWL_COLLECTION object",io )) return
  endif

  if (wl1%iCount > 0) then

    iCol = 0
    rARow = rZERO
    ! ASSUMES that WL1_Setup routine created consecutive well entries 
    wel => wl1%Wells(1)
    iWL1 = wel%iFWLIndex
    iNWL = wl1%iCount
    allocate (cWLF(1:iNWL,1,1),stat=iStat)
    if (IO_Assert( (iStat==0), "WL1_ComputeCoefficients: Allocation failed",io )) return

    ! Get the appropriate incluence functions for the boundary condition type
    select case ( iEqType )
      case (EQN_HEAD)
        iWhich = INFLUENCE_P
      case (EQN_FLOW)
        iWhich = INFLUENCE_F
      case (EQN_INHO)
        iWhich = INFLUENCE_P
      case (EQN_DISCHARGE)
        iWhich = INFLUENCE_W
      case (EQN_RECHARGE)
        iWhich = INFLUENCE_G
      case (EQN_CONTINUITY)
        iWhich = INFLUENCE_Q
      case (EQN_POTENTIALDIFF)
        iWhich = INFLUENCE_D
      case (EQN_TOTALFLOW)
        iWhich = INFLUENCE_Z
    end select
   
    call FWL_GetInfluence(fwl,iWhich,iWL1,iNWL,cPathZ,cOrientation,cWLF(1:iNWL,:,:),io)
    if ( io%lError ) return
    do iVtx=1,iNWL
      iCol = iCol+1
      rARow(iCol) = real( cWLF(iVtx,1,1) )
    end do
 
    deallocate(cWLF)

  endif

  return
end subroutine WL1_ComputeCoefficients

function rWL1_ComputeRHS(wl1,iElementType,iElementString,iElementVertex, &
                         iElementFlag,rSpecValue,rCheck,io) result(rRHS)
  !! function rWL1_ComputeRHS
  !!
  !! Computes the right-hand side value for the solution
  !!
  !! Calling Sequence:
  !!   rRHS = rWL1_ComputeRHS(wl1,rValue,iElementType,iElementString,iElementVertex, &
  !!                          iElementFlag,rSpecValue,rCheck)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             WL1_COLLECTION object to be used
  !!   (in)    integer :: iElementType
  !!             Element type (either ELAM_AQU or ELEM_IN0)
  !!   (in)    integer :: iElementString
  !!             Element string number
  !!   (in)    integer :: iElementVertex
  !!             Element vertex number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!   (in)    real :: rSpecValue
  !!             The new result value from the solution vector
  !!   (in)    real :: rCheck
  !!             The new result value from the solution vector
  !!
  !! Return Value:
  !!   real :: rRHS
  !!     The RHS value for the module
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real), intent(in) :: rSpecValue
  real (kind=ModAEM_Real), intent(in) :: rCheck
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rRHS

  ! For WL1, this is easy...
  rRHS = rSpecValue - rCheck

  return
end function rWL1_ComputeRHS

subroutine WL1_StoreResult(wl1,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
  !! subroutine WL1_StoreResult
  !!
  !! Stores the results of a solution for a single equation associated with
  !! the WL1 module.
  !!
  !! Calling Sequence:
  !!    WL1_StoreResult(wl1,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    real :: rValue
  !!             The new result value from the solution vector
  !!   (in)    integer :: iElementType
  !!             Element type (always ELEM_WL1)
  !!   (in)    integer :: iElementWell
  !!             Element well number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  real (kind=ModAEM_Real), intent(in) :: rValue
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type ( WL1_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_StoreResult: WL1_Create has not been called",io )) return
    if (IO_Assert( (iElementString>=1 .and. iElementString<=wl1%iCount), &
                    "WL1_StoreResult: Bad element well ID",io)) return
  endif

  wel => wl1%Wells(iElementString)
  wel%rStrength = wel%rStrength + rValue 

  return
end subroutine WL1_StoreResult

function rWL1_ComputeCheck(wl1,iElementWell,cCPZ,rPot,io) result(rCheck)
  !! function rWL1_ComputeCheck
  !!
  !! Returns the check value for the specified domain and vertex at the point cZ
  !!
  !! Calling Sequence:
  !!    rCheck = rWL1_ComputeCheck(in0,iElementString,iElementVertex,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             WL1_COLLECTION object to be used
  !!   (in)    integer :: iElementString
  !!             The domain number to be examined
  !!   (in)    integer :: iElementVertex
  !!             The vertex number to be examined
  !!   (in)    complex,dimension(:) :: cCPZ
  !!             Control point(s) to be used in coefficient calculations
  !!   (in)    real :: rPot
  !!             The potential at cCPZ(1)
  !!
  !! Return Value:
  !!   (out)   real :: rCheck
  !!             The check value
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  integer (kind=ModAEM_Integer),intent(in) :: iElementWell
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cCPZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rCheck
  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_ComputeCheck: WL1_Create has not been called",io )) return
  endif

  rCheck = rPot

  return
end function rWL1_ComputeCheck

subroutine WL1_Update(wl1,fwl,io)
  !! subroutine WL1_Update
  !!
  !! Updates the underlying function objects for the specified layer.
  !!
  !! Calling Sequence:
  !!    WL1_Update(wl1)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iWell
  type ( WL1_WELL ),pointer :: wel


  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Update: WL1_Create has not been called",io )) return
  endif

  do iWell = 1,wl1%iCount
    wel => wl1%Wells(iWell)
    call FWL_Update(fwl,wel%iFWLIndex,wel%rStrength,io)
    if ( io%lError ) return
  end do

  ! Make sure I never force a matrix regeneration
  wl1%iRegenerate = 0

  return
end subroutine WL1_Update

function lWL1_CheckPoint(wl1,cZ,rTol,io) result(lSing)
  !! logical function lWL1_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lWL1_CheckPoint(wl1,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             The WL1_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (wl1_COLLECTION),pointer :: wl1
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iWell
  type ( WL1_WELL ),pointer :: wel

  lSing = .false.
  do iWell = 1,wl1%iCount
    wel => wl1%Wells(iWell)
    if ( abs(real(wel%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(wel%cZ-cZ))<rTol ) then
      lSing = .true.
      exit
    endif
  end do
  
  
  return
end function lWL1_CheckPoint

function lWL1_CheckProximity(wl1,cZ,rProximity,iDir,iElementID,io) result(lFound)
  !! function lWL1_CheckProximity
  !!
  !! Checks to see if the point specified is within the radius of a well, then
  !! moves the point outside well perimeter if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    lFound = WL1_CheckProximity(wl1,cZ,rProximity,iDir,iElementID) result(lFound)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             WL1_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (in)    integer :: iDir
  !!             Flag for checking wells:
  !!               iDir =  0   -  Check without regard to the direction of travel
  !!               iDir = <0   -  Check for back tracing (ignore pumping wells)
  !!               iDir = >0   -  Check for formard tracing (ignore injection wells)
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the line-sink found (if lFound is .true.)
  !!   (out)   logical :: lFound
  !!             .true. if a well is found
  !!             .false. if a well is not found
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lFound
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  real (kind=ModAEM_Real) :: rDist,rMinDist
  type ( WL1_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Inquiry: WL1_Create has not been called",io )) return
  endif

  ! Search through the wells, according to the iDir
  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do i=1,wl1%iCount
    wel => wl1%Wells(i)
    if ( iDir < 0 .and. wel%rStrength >= rZERO ) then
      ! skip pumping wells
      cycle
    endif
    if ( iDir > 0 .and. wel%rStrength <= rZERO ) then
      ! skip injection wells
      cycle
    endif
    ! check the well
    rDist = abs(cZ - wel%cZ)
    if ( rDist < rProximity .or. rDist < wel%rRadius ) then
      if ( rDist < rMinDist ) then
        rMinDist = rDist
        iElementID = wel%iID
        lFound = .true.
      endif
    endif
  end do

  return
end function lWL1_CheckProximity

subroutine WL1_Read(wl1,io)
  !! subroutine WL1_Read
  !!
  !! Populates an WL1_COLLECTION using data from kIOInputLU
  !!
  !! Calling Sequence:
  !!    call WL1_Read(wl1)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer :: wl1
  !!             WL1_COLLECTION to be populated
  !!
  !! The format of the WL1 section of the input file appears as follows:
  !! WL1 NWells
  !!     x y r x y head id
  !!     ... Up to NWells
  !! END
  !!
  !! NOTE: It is assumed that the WL1 line was found by the caller

  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCAL DIRECTIVES ]
  type (DIRECTIVE),dimension(3),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD /)
  ! [ LOCALS ]
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  real (kind=ModAEM_Real) :: rRad, rHead
  complex (kind=ModAEM_Real) :: cZ, cZHead
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iID
  logical (kind=ModAEM_Integer) :: lFlag
  type ( WL1_WELL ), pointer :: wel

  call IO_MessageText("  Reading WL1 module input",io)

  if (IO_Assert( (associated(wl1)), "WL1_Read: WL1_Create has not been called",io )) return

  ! Process input   
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "WL1_Read: I/O Error",io)) return
        exit
      case (kOpFileEOF)
        ! EOF is unexpected for all ModWL1 "ifXXXRead" routines. 
        ! Report the condition, but proceed as if EOD was found.
        if (IO_Assert( .false., "WL1_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(wl1%Wells)), "WL1_Read: No space allocated",io)) return
        if (IO_Assert( (wl1%iCount<size(wl1%Wells)), "WL1_Read: Space exhausted",io)) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ,rRad,cZHead,rHead,iID
        if (IO_Assert( (iStat==0), "WL1_Read: I/O Error",io)) return
        wl1%iCount = wl1%iCount+1
        wel => wl1%Wells(wl1%iCount)
        wel%cZ = cZ
        wel%rRadius = rRad
        wel%cZHead = cZHead
        wel%rHead = rHead
        wel%rHeadCorrection = rZERO
        wel%iID = iID
        ! No FWL is declared here; see WL1_Setup
        wel%iFWLIndex = -1
        wel%rStrength = rZero
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        exit
      case (kOpDBG)
        ! Change the io%lDebug flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "WL1_Read: I/O Error",io )) return
        call IO_SetDebug(lFlag,io)
        if ( io%lError ) return
      case (kOpPCD)
        ! Change the IO_Proceed flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "WL1_Read: I/O Error",io )) return
        call IO_SetProceed(lFlag,io)
        if ( io%lError ) return
    end select
  end do  

  call IO_MessageText("  Leaving WL1 module",io)

end subroutine WL1_Read

subroutine WL1_Inquiry(wl1,iLU,io)
  !! subroutine WL1_Inquiry
  !!
  !! Writes an inquiry report for all line-sinks to iLU
  !!
  !! Calling Sequence:
  !!    call WL1_Inquiry(wl1,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( WL1_WELL ),pointer :: wel
  
  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Inquiry: WL1_Create has not been called",io )) return
  endif

  do i=1,wl1%iCount
    wel => wl1%Wells(i)
    write ( unit=iLU, &
            fmt="(""WL1"",2("","",i9),6("","",e14.6))" &
          ) wel%iID,wel%cZ,wel%rRadius,wel%cZHead,wel%rHead
  end do

  return
end subroutine WL1_Inquiry

subroutine WL1_Report(wl1,io)
  !! subroutine WL1_Report
  !!
  !! Writes a debugging report for all line-sinks to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call WL1_Report(wl1)
  !!
  !! Arguments:
  !!   (in)    type ( WL1_COLLECTION ),pointer
  !!             WL1_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL1_COLLECTION ),pointer :: wl1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type ( WL1_WELL ),pointer :: wel
  integer (kind=ModAEM_Integer) :: iWell
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN

  if ( io%lDebug ) then 
    if (IO_Assert( (associated(wl1)), &
                    "WL1_Report: WL1_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module WL1'//'Head-specified well element information')" )

  if ( .not. associated(wl1%Wells) ) then
    write ( unit=kIOOutputLU, &
            fmt="("" No Wells Allocated""/)" &
          )
  else
    write ( unit=kIOOutputLU, &
            fmt="("" Number of wells: "",i5,""  used: "",i5)" &
          ) ubound(wl1%Wells,1),wl1%iCount

    ! Write the detail information for all wells...
    write ( unit=kIOOutputLU, &
            fmt="("" Wells:"")" &
          )
    write ( unit=kIOOutputLU, &
            fmt= "(""    ID"",t15,""       X"",t30,""       Y"",t45, " // &
                 " ""  Radius"",t60,""  X(Head)"",t75,""  Y(Head)"",t90, " // &
                 " ""    Head"",t105,""Discharge"",t120,""FWL Index"")" &
          )
    do iWell=1,wl1%iCount
      wel => wl1%Wells(iWell)
      write( unit=kIOOutputLU, &           
             fmt="(i10,t15,d12.5,t30,d12.5,t45,d12.5,t60,d12.5,t75,d12.5,t90,d12.5,t105,d12.5,t120,i10)" &
           ) wel%iID,wel%cZ,wel%rRadius,wel%cZHead,wel%rHead,wel%rStrength,wel%iFWLIndex
    end do


    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iWL1_GetInfo(wl1,SIZE_FWL,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iWL1_GetInfo(wl1,SIZE_FPD,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iWL1_GetInfo(wl1,SIZE_FDP,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iWL1_GetInfo(wl1,SIZE_EQUATIONS,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iWL1_GetInfo(wl1,SIZE_UNKNOWNS,io)
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine WL1_Report

end module m_wl1

