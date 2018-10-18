module m_ls1

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

  !! module m_ls1
  !!
  !! Element module for 2-D head specified line-sinks
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! line-sink elements. Elements are defined as strings of points, with the
  !! discharge rate (per length) specified for each segment.
  !!
  !! Note: ModAEM uses the "traditional" line-sink function (Strack, 1989)
  !! for matrix generation, but uses strings of dipoles terminated by a well
  !! for computational performance, once a solution is achieved.

use u_constants
use u_io
use f_well
use f_dipole
use u_matrix
use m_aqu

implicit none

public

  type :: LS1_VERTEX
    !! type LS1_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rHead
    !!     The specified head at the vertex
    !!   real :: rHeadCorrection
    !!     Correction on the specified head computed by an outside module
    !!   real :: rDPStrength
    !!     The dipole strength at the vertex
    !!   real :: rLength
    !!     The segment length
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rHead
    real (kind=ModAEM_Real) :: rHeadCorrection
    real (kind=ModAEM_Real) :: rDPStrength
    real (kind=ModAEM_Real) :: rStrength
    real (kind=ModAEM_Real) :: rLength
    integer (kind=ModAEM_Integer) :: iFDPIndex
  end type LS1_VERTEX

  type :: LS1_STRING
    !! type LS1_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   type ( LS1_VERTEX ),dimension(:),pointer :: Vertices
    !!     A vector of LS1_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !! 
    type ( LS1_VERTEX ),dimension(:),pointer :: Vertices
    integer (kind=ModAEM_Integer) :: iNPts
    integer (kind=ModAEM_Integer) :: iFWLIndex
    integer (kind=ModAEM_Integer) :: iID
  end type LS1_STRING

  type :: LS1_COLLECTION
    !! type LS1_COLLECTION
    !!
    !! Type that holds information for all LS1 elements in a layer
    !!
    !! Members:
    !!   type ( LS1_STRING ),dimension(:),pointer :: Strings
    !!     A vector of LS1_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !! 
    type ( LS1_STRING ),dimension(:),pointer :: Strings
    integer (kind=ModAEM_Integer) :: iNStr
    integer (kind=ModAEM_Integer) :: iRegenerate
  end type LS1_COLLECTION

  ! This buffer holds LS1 data for all layers
  type ( LS1_COLLECTION ),dimension(:),allocatable,target,private :: ls1ers

  ! Module flags for matrix generator routines
  integer (kind=ModAEM_Integer),private,parameter :: LS1Vertex=1

contains

!**pd Modified to allocate only the collection object
function LS1_Create(io) result (ls1)
  !! function LS1_Create
  !!
  !! Creates a new LS1_COLLECTION object
  !!
  !! Calling Sequence:
  !!    ls1 => LS1_Create()
  !!
  !! Arguments:
  !!
  ! [ ARGUMENTS ]
  ! [ RETURN VALUE ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (ls1,stat=iStat)
  if (IO_Assert( (iStat == 0), "LS1_Create: allocation failed",io )) return
  nullify(ls1%Strings)
  ls1%iNStr = 0
  ls1%iRegenerate = 1
 
  return
end function LS1_Create

!**pd New LS1_Alloc subroutine
subroutine LS1_Alloc(ls1,iNStr,io)
  !! Subroutine LS1_Alloc
  !! 
  !! Allocates Strings for the LS1_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call LS1_Alloc(ls1, iNStr)
  !!
  !! Arguments:
  !!    (in)    type ( LS1_COLLECTION ),pointer :: ls0
  !!              The LS1_COLLECTION object to be used
  !!    (in)    integer :: iNStr
  !!              The number of strings to make space for
  !!
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  integer (kind=ModAEM_Integer),intent(in) :: iNStr
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (ls1%Strings(iNStr),stat=iStat)
  if (IO_Assert( (iStat == 0), "LS1_Alloc: allocation failed",io )) return

end subroutine LS1_Alloc


!**pd  New LS1_Destroy subroutine
subroutine LS1_Destroy(ls1,io)
  !! subroutine LS1_Destroy
  !!
  !! Frees memory allocated for ls1 Linesinks and strings of vertices 
  !! and the ls1 Collection object
  !!
  !! Calling Sequence:
  !!     call LS1_Destroy(ls1)
  !!
  !! Arguments:
  !!  type ( ls1_COLLECTION ),pointer :: ls1
  !!              Pointer to the ls1_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( ls1_COLLECTION ),pointer :: ls1
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  integer (kind=ModAEM_Integer):: i
  type (ls1_STRING), pointer :: str


  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Destroy: LS1_Create has not been called",io )) return
  endif

  ! First deallocate each string of vertices
  do i = 1, ls1%iNStr
    str => ls1%Strings(i)
    deallocate (str%Vertices,stat=iStat)
    if (IO_Assert( (iStat == 0), "LS1_Destroy: deallocation of Vertices failed",io )) return
  end do
  ! Then deallocate the strings table
  if ( associated(ls1%Strings) ) then
    deallocate (ls1%Strings,stat=iStat)
    if (IO_Assert( (iStat == 0), "LS1_Destroy: deallocation of Strings failed",io )) return
  endif
  ! Now deallocate the collection
  deallocate (ls1,stat=iStat)
  if (IO_Assert( (iStat == 0), "LS1_Destroy: deallocation failed",io )) return
  
  return
end subroutine LS1_Destroy


subroutine LS1_New(ls1,Vertices,iNPts,io) 
  !! function LS1_New
  !!
  !! Adds a new LS1_STRING object to the LS1_COLLECTION 'ls1'
  !!
  !! Calling Sequence:
  !!    call LS1_New(ls1,Vertices,iNPt)
  !!
  !! Arguments:
  !!    (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!              The LS1_COLLECTION object to be used
  !!    (in)    type ( LS1_VERTEX ) :: Vertices(:)
  !!              Vector that defines the points along the barrier
  !!    (in)    integer :: iNPt
  !!              The number of vertices in the string
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( LS1_VERTEX ),dimension(:) :: Vertices
  integer (kind=ModAEM_Integer),intent(in) :: iNPts
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  type ( LS1_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_New: LS1_Create has not been called",io )) return
  endif

  if (IO_Assert( (ls1%iNStr < size(ls1%Strings)), &
                  "LS1_New: Space exhausted",io )) return
  if (IO_Assert( (iNPts <= size(Vertices)), &
                  "LS1_New: Size of provided vertices is inconsistent",io )) return

  ls1%iNStr = ls1%iNStr + 1
  str => ls1%Strings(ls1%iNStr)
  allocate ( str%Vertices(iNPts),stat=iStat )
  if (IO_Assert( (iStat==0), "LS1_New: Allocation failed",io )) return
  str%Vertices = Vertices(1:iNPts)
  str%iNPts = iNPts

  return
end subroutine LS1_New

function iLS1_GetInfo(ls1,iOption,io) result(iValue)
  !! function LS1_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iLS1_GetInfo(ls1,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!             LS1_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue
  integer (kind=ModAEM_Integer) :: iStr
  type ( LS1_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_GetInfo: LS1_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FWL )
      iValue = ls1%iNStr
    case ( SIZE_FDP )
      do iStr=1,ls1%iNStr
        str => ls1%Strings(ls1%iNStr)
        iValue = iValue + str%iNPts-1
      end do
    case ( SIZE_EQUATIONS )
      do iStr=1,ls1%iNStr
        str => ls1%Strings(ls1%iNStr)
        iValue = iValue + str%iNPts-1
      end do
    case ( SIZE_UNKNOWNS )
      do iStr=1,ls1%iNStr
        str => ls1%Strings(ls1%iNStr)
        iValue = iValue + str%iNPts-1
      end do
    case ( INFO_REGENERATE )
      iValue = ls1%iRegenerate
    case default
      iValue = 0
  end select

  return
end function iLS1_GetInfo

subroutine LS1_SetupFunctions(ls1,fwl,fdp,io)
  !! subroutine LS1_Setup
  !!
  !! This routine sets up the functions in f_well and f_dipole for the line-sinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call LS1_Setup(ls1)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer
  !!             FDP_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx,i,iDP,iWL
  real (kind=ModAEM_Real) :: rStrength,rDisch,rHead1,rHead2,rHead
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  complex (kind=ModAEM_Real),dimension(3) :: cCPResult
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx,last_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Setup: LS1_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "LS0_Setup: Illegal FWL_COLLECTION object",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "LS1_Setup: Illegal FDP_COLLECTION object",io )) return
  endif

  do iStr = 1,ls1%iNStr
    str => ls1%Strings(iStr)
    ! Build dipoles for all segments
    do iVtx = 1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      this_vtx%rLength = abs(next_vtx%cZ - this_vtx%cZ)
      call FDP_New(fdp,this_vtx%cZ,next_vtx%cZ,(/cZERO,cZERO,cZERO/),this_vtx%iFDPIndex,io)
      if ( io%lError ) return
    end do

    ! Put a well at the end of the string
    last_vtx => str%Vertices(str%iNPts)
    call FWL_New(fwl,last_vtx%cZ,rZERO,rZERO,iWL,io)
    if ( io%lError ) return
    str%iFWLIndex = iWL
  end do

  return
end subroutine LS1_SetupFunctions

subroutine LS1_SetupMatrix(ls1,aqu,mat,io)
  !! subroutine LS1_SetupMatrix
  !!
  !! This routine sets up the matrix entries for the linesinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call LS1_Setup(ls1)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    type ( AQU_COLLECTION ),pointer
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx,i,iDP,iWL
  real (kind=ModAEM_Real) :: rStrength,rDisch,rHead1,rHead2,rHead
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  complex (kind=ModAEM_Real),dimension(3) :: cCPResult
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx,last_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Setup: LS1_Create has not been called",io )) return
    if (IO_Assert( (associated(aqu)), &
                    "LS1_Setup: Illegal AQU_COLLECTION object",io )) return
    if (IO_Assert( (associated(mat)), &
                    "LS1_Setup: Illegal MAT_MATRIX object",io )) return
  endif

  ! Build matrix generator entries for all segments
  do iStr = 1,ls1%iNStr
    ! Set up the unknown variables 
    ! Vertex entry -- all vertices 
    str => ls1%Strings(iStr)
    do iVtx = 1,str%iNPts-1
      call MAT_CreateVariable(mat,ELEM_LS1,iStr,iVtx,LS1Vertex,io)
      if ( io%lError ) return
    end do

    ! Set up control points and equations -- One equation per segment
    do iVtx = 1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      rHead1 = this_vtx%rHead
      rHead2 = next_vtx%rHead
      ! Compute control points for the segment.  There is one unknown
      ! per segment.
      call MAT_ComputeControlPoints(this_vtx%cZ,next_vtx%cZ,1,cCPResult,rZERO,io)
      if ( io%lError ) return
      ! Now, create the equation entry...
      rHead = rHead1 + (rHead2-rHead1)*(cCPResult(2)-this_vtx%cZ)/(next_vtx%rHead-this_vtx%cZ)
      call MAT_CreateEquation(mat,(/ cCPResult(2) /),EQN_HEAD,ELEM_LS1, &
                             iStr,iVtx,0,rAQU_HeadToPotential(aqu,rHead,cCPResult(2),io),cZERO,io)
      if ( io%lError ) return
    end do
  end do

  return
end subroutine LS1_SetupMatrix

subroutine LS1_Prepare(ls1,io) 
  !! subroutine LS1_Prepare
  !! 
  !! Prepares the module for a new iteration
  !!
  !! Do-nothing for m_ls1
  !!
  !! Calling Sequence:
  !!    call LS1_Setup(ls1,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( IO_STATUS ),pointer :: io

  return
end subroutine LS1_Prepare

subroutine LS1_ComputeCoefficients(ls1,fwl,fdp,cPathZ,iEqType,iElementType,iElementString, &
                                   iElementVertex,iElementFlag,cOrientation,rARow,io)
  !! subroutine LS1_ComputeCoefficients
  !!
  !! Computes a row of matrix coefficients (with no corrections) for the LS1 
  !! elements in layer iL.
  !!
  !! Calling Sequence:
  !!    call LS1_ComputeCoefficients(ls1,cPathZ,iEqType,cOrientation,rRow)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer
  !!             FDP_COLLECTION object to be used
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
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
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
  integer (kind=ModAEM_Integer) :: iStat,iCol,iStr,iVtx,iDP1,iNDP,iWhich
  complex (kind=ModAEM_Real),dimension(:,:,:),allocatable :: cDPF
  complex (kind=ModAEM_Real),dimension(1) :: cWLF1
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: first_vtx

  if ( io%lDebug ) then 
    if (IO_Assert( (associated(ls1)), &
                    "LS1_ComputeCoefficients: LS1_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "LS0_Setup: Illegal FWL_COLLECTION object",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "LS1_ComputeCoefficients: Illegal FDP_COLLECTION object",io )) return
  endif

  iCol = 0
  rARow = rZERO
  do iStr=1,ls1%iNStr
    str => ls1%Strings(iStr)
    first_vtx => str%Vertices(1)
    ! ASSUMES that LS1_Setup routine created consecutive dipole entries 
    iDP1 = first_vtx%iFDPIndex
    iNDP = str%iNPts-1
    allocate (cDPF(0:iNDP+1,1,1),stat=iStat)
    if (IO_Assert( (iStat==0), "LS1_ComputeCoefficients: Allocation failed",io )) return

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

    call FDP_GetInfluence_ILS(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cDPF(1:iNDP,:,:),io)
    if ( io%lError ) return
    do iVtx=1,iNDP
      iCol = iCol+1
      rARow(iCol) = real( cDPF(iVtx,1,1) )
    end do

    deallocate(cDPF)
  end do
end subroutine LS1_ComputeCoefficients

subroutine LS1_StoreResult(ls1,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
  !! subroutine LS1_StoreResult
  !!
  !! Stores the results of a solution for a single equation associated with
  !! the LS1 module.
  !!
  !! Calling Sequence:
  !!    LS1_StoreResult(ls1,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    real :: rValue
  !!             The new result value from the solution vector
  !!   (in)    integer :: iElementType
  !!             Element type (always ELEM_LS1)
  !!   (in)    integer :: iElementString
  !!             Element string number
  !!   (in)    integer :: iElementVertex
  !!             Element vertex number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  real (kind=ModAEM_Real), intent(in) :: rValue
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_StoreResult: LS1_Create has not been called",io )) return
    if (IO_Assert( (iElementString>=1 .and. iElementString<=ls1%iNStr), &
                    "LS1_StoreResult: Bad element string ID",io)) return
  endif

  str => ls1%Strings(iElementString)

  if ( io%lDebug ) then
    if (IO_Assert( (iElementVertex>=1 .and. iElementVertex<=str%iNPts-1), &
                    "LS1_StoreResult: Bad element vertex ID",io)) return
  endif

  this_vtx => str%Vertices(iElementVertex)
  next_vtx => str%Vertices(iElementVertex+1)
  this_vtx%rStrength = this_vtx%rStrength + rValue 
  next_vtx%rDPStrength = this_vtx%rDPStrength + this_vtx%rStrength * this_vtx%rLength

  return
end subroutine LS1_StoreResult

function rLS1_ComputeRHS(ls1,iElementType,iElementString,iElementVertex, &
                         iElementFlag,rSpecValue,rCheck,io) result(rRHS)
  !! function rLS1_ComputeRHS
  !!
  !! Computes the right-hand side value for the solution
  !!
  !! Calling Sequence:
  !!   rRHS = rLS1_ComputeRHS(ls1,rValue,iElementType,iElementString,iElementVertex, &
  !!                          iElementFlag,rSpecValue,rCheck)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!             AQU_COLLECTION object to be used
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
  type ( LS1_COLLECTION ),pointer :: ls1
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real), intent(in) :: rSpecValue
  real (kind=ModAEM_Real), intent(in) :: rCheck
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rRHS

  ! For LS1, this is easy...
  rRHS = rSpecValue - rCheck

  return
end function rLS1_ComputeRHS

function rLS1_ComputeCheck(ls1,iElementString,iElementVertex,cCPZ,rPot,io) result(rCheck)
  !! function rLS1_ComputeCheck
  !!
  !! Returns the check value for the specified domain and vertex at the point cZ
  !!
  !! Calling Sequence:
  !!    rCheck = rLS1_ComputeCheck(in0,iElementString,iElementVertex,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!             LS1_COLLECTION object to be used
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
  type ( IN0_COLLECTION ),pointer :: ls1
  integer (kind=ModAEM_Integer),intent(in) :: iElementString
  integer (kind=ModAEM_Integer),intent(in) :: iElementVertex
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cCPZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rCheck
  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_ComputeCheck: LS1_Create has not been called",io )) return
  endif

  rCheck = rPot

  return
end function rLS1_ComputeCheck

subroutine LS1_Update(ls1,fwl,fdp,io)
  !! subroutine LS1_Update
  !!
  !! Updates the underlying function objects for the specified layer.
  !!
  !! Calling Sequence:
  !!    LS1_Update(ls1)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer
  !!             FDP_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx


  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Update: LS1_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "LS1_Update: Illegal FDP_COLLECTION object",io )) return
  endif

  do iStr = 1,ls1%iNStr
    str => ls1%Strings(iStr)
    do iVtx = 1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      cRho1 = cmplx(this_vtx%rDPStrength,rZERO,ModAEM_Real)
      cRho3 = cmplx(next_vtx%rDPStrength,rZERO,ModAEM_Real)
      cRho2 = cZERO
      call FDP_Update(fdp,this_vtx%iFDPIndex,(/cRho1,cRho2,cRho3/),io)
      if ( io%lError ) return
    end do
    ! Put a well at the end of the string
    call FWL_Update(fwl,str%iFWLIndex,real(cRho3),io)
    if ( io%lError ) return
  end do

  ! LS1 never forces a regeneration
  ls1%iRegenerate = 0

  return
end subroutine LS1_Update

function lLS1_CheckPoint(ls1,cZ,rTol,io) result(lSing)
  !! logical function lLS1_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lLS1_CheckPoint(ls1,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!             The LS1_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (LS1_COLLECTION),pointer :: ls1
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: vtx

  ! NOTE: Singularities are found at the end points of the BDY elements 
  lSing = .false.
  do iStr = 1,ls1%iNStr
    str => ls1%Strings(iStr)
    do iVtx = 1,str%iNPts
      vtx => str%Vertices(iVtx)
      if ( abs(real(vtx%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(vtx%cZ-cZ))<rTol ) then
        lSing = .true.
        exit
      endif
    end do
  end do
  
  return
end function lLS1_CheckPoint

function lLS1_CheckProximity(ls1,cZ,rProximity,iDir,iElementID,iElementVtx,io) result(lFound)
  !! function lLS1_CheckProximity
  !!
  !! Checks to see if the point specified is near a line-sink, then
  !! moves the point a small distance if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    lFound = LS1_CheckProximity(ls1,cZ,rProximity,iDir,iElementID) result(lFound)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (in)    integer :: iDir
  !!             Flag for checking line-sinks:
  !!               iDir =  0   -  Check without regard to the direction of travel
  !!               iDir = <0   -  Check for back tracing (ignore gaining line-sinks)
  !!               iDir = >0   -  Check for formard tracing (ignore losing line-sinks)
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the line-sink found (if lFound is .true.)
  !!   (out)   integer :: iElementVtx
  !!             The vertex associated with the string segment which was found
  !!   (out)   logical :: lFound
  !!             .true. if a line-sink is found
  !!             .false. if a line-sink is not found
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID,iElementVtx
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lFound
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist
  complex (kind=ModAEM_Real) :: cZL,cZC,cBigZ
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_CheckProximity: LS1_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iStr=1,ls1%iNStr
    str => ls1%Strings(iStr)
    ! Sweep through the string, checking on the basis of iDir
    do iVtx=1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      ! iDir < 0 -- skip gaining line sinks
      if ( iDir < 0 .and. this_vtx%rStrength > rZERO ) cycle
      ! iDir > 0 -- skip losing line sinks
      if ( iDir > 0 .and. this_vtx%rStrength < rZERO ) cycle
      ! Compute the mapped Z-value and then the distance from the linesink
      cZC = rHALF * (this_vtx%cZ + &
                     next_vtx%cZ )
      cZL = rHALF * (this_vtx%cZ - &
                     next_vtx%cZ )
      cBigZ = (cZ-cZC) / cZL
      rDist = abs( aimag(cBigZ) * abs(cZL) )
      if ( rDist < rProximity .and. abs(real(cBigZ)) <= rONE ) then
        if ( rDist < rMinDist ) then
          rMinDist = rDist
          iElementID = str%iID
          iElementVtx = iVtx
          lFound = .true.
        endif
      endif
    end do
  end do

  return
end function lLS1_CheckProximity

subroutine LS1_Read(ls1,io)
  !! subroutine LS1_Read
  !!
  !! Populates an LS1_COLLECTION using data from kIOInputLU
  !!
  !! Calling Sequence:
  !!    call LS1_Read(ls1)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer :: ls1
  !!             LS1_COLLECTION to be populated
  !!
  !! The format of the LS1 section of the input file appears as follows:
  !! LS1
  !! STR NVertices
  !! x y head
  !! ... Up to NVertices
  !! STR NVertices
  !! x y head
  !! ... Up to NVertices
  !! ... Up to NStrings
  !! END
  !!
  !! NOTE: It is assumed that the LS1 line was found by the caller

  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCAL DIRECTIVES ]
  type (DIRECTIVE),dimension(4),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD,dirSTR /)
  ! [ LOCALS ]
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  real (kind=ModAEM_Real) :: rHead
  complex (kind=ModAEM_Real) :: cZ
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMaxStr, iMaxVtx
  integer (kind=ModAEM_Integer) :: iID
  integer (kind=ModAEM_Integer) :: iStr
  logical (kind=ModAEM_Integer) :: lFlag
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,last_vtx

  call IO_MessageText("  Reading LS1 module input",io)

  if (IO_Assert( (associated(ls1)), "LS1_Read: LS1_Create has not been called",io )) return

  ! Use IO_InputRecord to process the model input file.
  nullify ( str,this_vtx,last_vtx )
  do 
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "LS1_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModLS1 "ifXXXRead" routines. 
        if (IO_Assert( .false., "LS1_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(str)), "LS1_Read: No current string",io )) return
        if (IO_Assert( (str%iNPts<size(str%Vertices)), "LS1_Read: Space exhausted",io )) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ,rHead
        if (IO_Assert( (iStat==0), "LS1_Read: I/O Error",io )) return
        str%iNPts = str%iNPts+1
        this_vtx => str%Vertices(str%iNPts)
        this_vtx%cZ = cZ
        this_vtx%rHead  = rHead 
        this_vtx%rHeadCorrection  = rZERO
        this_vtx%rDPStrength = rZERO
        this_vtx%rStrength = rZERO
        if ( str%iNPts > 1 ) then
          last_vtx%rLength = abs(this_vtx%cZ - last_vtx%cZ)
        end if
        this_vtx%iFDPIndex = -1
        last_vtx => this_vtx
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        exit
      case (kOpDBG)
        ! Change the io%lDebug flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetDebug(lFlag,io)
        if ( io%lError ) return
      case (kOpPCD)
        ! Change the IO_Proceed flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetProceed(lFlag,io)
        if ( io%lError ) return
      case (kOpSTR)
        !****************************************************************************
        ! Here for the STR command -- create a new string of line-sinks
        ! the maximum number of vertices is in the input record
        !****************************************************************************
        if (IO_Assert( (ls1%iNStr<size(ls1%Strings)), "LS1_Read: Space exhausted",io )) return
        ! Retrive the number of vertices desired...
        read (unit=sRecord,fmt=*,iostat=iStat) iMaxVtx,iID
        if (IO_Assert((iStat==0), "LS1_Read: I/O Error",io)) return
        if (IO_Assert((iMaxVtx>0), "LS1_Read: Illegal number of vertices",io)) return
        ! OKAY! Allocate the vertices...
        ls1%iNStr = ls1%iNStr+1
        str => ls1%Strings(ls1%iNStr)
        allocate (str%Vertices(iMaxVtx),stat=iStat)
        if (IO_Assert((iStat==0), "LS1_Read: Allocation failed",io)) return
        ! Made it!
        str%iFWLIndex = -1         ! No FWL function yet! 
        str%iID = iID
        str%iNPts = 0      ! Initialize the vertex counter
    end select
  end do  

  call IO_MessageText("  Leaving LS1 module",io)

end subroutine LS1_Read

subroutine LS1_Inquiry(ls1,iLU,io)
  !! subroutine LS1_Inquiry
  !!
  !! Writes an inquiry report for all line-sinks to iLU
  !!
  !! Calling Sequence:
  !!    call LS1_Inquiry(ls1,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rLength
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Inquiry: LS1_Create has not been called",io )) return
  endif

  do iStr=1,ls1%iNStr
    str => ls1%Strings(iStr)
    do iVtx=1,str%iNPts
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      if ( iVtx < str%iNPts ) then
        rLength = abs(next_vtx%cZ - this_vtx%cZ)
      else
        rLength = rZERO
      endif
      write ( unit=iLU, &
              fmt="(""LS1"",3("","",i9),4("","",e14.6))" &
            ) str%iID,iVtx,this_vtx%cZ,rLength,this_vtx%rStrength
    end do
  end do

  return
end subroutine LS1_Inquiry

subroutine LS1_Report(ls1,io)
  !! subroutine LS1_Report
  !!
  !! Writes a debugging report for all line-sinks to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call LS1_Report(ls1)
  !!
  !! Arguments:
  !!   (in)    type ( LS1_COLLECTION ),pointer
  !!             LS1_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS1_COLLECTION ),pointer :: ls1
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( LS1_STRING ),pointer :: str
  type ( LS1_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then 
    if (IO_Assert( (associated(ls1)), &
                    "LS1_Report: LS1_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module LS1'//'Head-specified linesink element information')" )

  if ( .not. associated(ls1%Strings) ) then
    write ( unit=kIOOutputLU, &
            fmt="("" No Strings Allocated""/)" &
          )
  else
    write ( unit=kIOOutputLU, &
            fmt="("" Number of strings: "",i5,""  used: "",i5)" &
          ) size(ls1%Strings),ls1%iNStr

    ! Write the string information for all strings...
    do iStr=1,ls1%iNStr
      str => ls1%Strings(iStr)
      write ( unit=kIOOutputLU, &
              fmt="(" // &
                  "/ "" "",100(""-"") / " // &
                  " "" String: "",i10,"" # Vertices "",i10,"" Well Index"",i10 " // &
                  ")" &
            ) str%iID,str%iNPts,str%iFWLIndex

      ! Write the vertices for this string...
      write ( unit=kIOOutputLU, &
              fmt="("" Vertices:"")" &
            )
      write ( unit=kIOOutputLU, &
              fmt="(""        X"",t15,""       Y"",t30,""   Strength"",t45,""   DP Index"")" &
            )
      do iVtx=1,str%iNPts
        vtx => str%Vertices(iVtx)
        write ( unit=kIOOutputLU, &
                fmt="(1x,d12.5,t15,d12.5,t30,d12.5,t45,i10)" &
              ) vtx%cZ,vtx%rStrength,vtx%iFDPIndex
      end do
    end do

    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iLS1_GetInfo(ls1,SIZE_FWL,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iLS1_GetInfo(ls1,SIZE_FPD,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iLS1_GetInfo(ls1,SIZE_FDP,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iLS1_GetInfo(ls1,SIZE_EQUATIONS,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iLS1_GetInfo(ls1,SIZE_UNKNOWNS,io)
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine LS1_Report

end module m_ls1

