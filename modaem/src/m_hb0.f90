module m_hb0

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

  !! module m_hb0
  !!
  !! Element module for 2-D horizontal no-flow boundaries
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for no-flow strings
  !! (not closed domains at this time). The strings are constructed using
  !! line-doublets (implemented as line-dipoles of imaginary strength).

use u_constants
use u_io
use u_matrix
use f_dipole

implicit none

public

  type :: HB0_VERTEX
    !! type HB0_VERTEX
    !!
    !! Type that holds information for one vertex along a no-flow string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rVertexStrength
    !!     The doublet strength at the vertex
    !!   real :: rCenterStrength
    !!     The doublet strength at the center of the next segment
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rVertexStrength
    real (kind=ModAEM_Real) :: rCenterStrength
    integer (kind=ModAEM_Integer) :: iFDPIndex
  end type HB0_VERTEX

  type :: HB0_STRING
    !! type HB0_STRING
    !!
    !! Type that holds information for one no-flow string
    !!
    !! Members:
    !!   type ( HB0_VERTEX ),dimension(:),pointer :: Vertices
    !!     A vector of HB0_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iID
    !!     The ID number for the string
    !! 
    type ( HB0_VERTEX ),dimension(:),pointer :: Vertices
    integer (kind=ModAEM_Integer) :: iNPts
    integer (kind=ModAEM_Integer) :: iID
  end type HB0_STRING

  type :: HB0_COLLECTION
    !! type HB0_COLLECTION
    !!
    !! Type that holds information for all HB0 elements in a layer
    !!
    !! Members:
    !!   type ( HB0_STRING ),dimension(:),pointer :: Strings
    !!     A vector of HB0_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !! 
    type ( HB0_STRING ),dimension(:),pointer :: Strings    
    integer (kind=ModAEM_Integer) :: iNStr                     
    integer (kind=ModAEM_Integer) :: iRegenerate
  end type HB0_COLLECTION

  ! Matrix generator element flags 
  integer (kind=ModAEM_Integer),private,parameter :: kHB0_Vertex=1
  integer (kind=ModAEM_Integer),private,parameter :: kHB0_Center=2

  real (kind=ModAEM_Real) :: MOVEFACTOR = 1.0001

contains

function HB0_Create(io) result (hb0)
  !! function HB0_Create
  !!
  !! Creates a new HB0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    hb0 => HB0_Create()
  !!
  !! Arguments:
  !!
  ! [ io%lDebug ]
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  type ( HB0_COLLECTION ),pointer :: hb0
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (hb0,stat=iStat)
  if (IO_Assert( (iStat == 0), "HB0_Create: allocation failed",io )) return
  nullify(hb0%strings)
  hb0%iNStr = 0
  hb0%iRegenerate = 1

  return
end function HB0_Create

!**pd New HB0_Alloc subroutine
subroutine HB0_Alloc(hb0, iNStr,io)
  !! Subroutine HB0_Alloc
  !! 
  !! Allocates Strings for the HB0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call HB0_Alloc(hb0, iNStr)
  !!
  !! Arguments:
  !!    (in)    type ( HB0_COLLECTION ),pointer :: ls0
  !!              The HB0_COLLECTION object to be used
  !!    (in)    integer :: iNStr
  !!              The number of strings to make space for
  !!
  !!
  !! Return Value:
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  integer (kind=ModAEM_Integer),intent(in) :: iNStr
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (hb0%Strings(iNStr),stat=iStat)
  if (IO_Assert( (iStat == 0), "HB0_Alloc: allocation failed,",io )) return

end subroutine HB0_Alloc

!**pd  New HB0_Destroy subroutine
subroutine HB0_Destroy(hb0,io)
  !! subroutine HB0_Destroy
  !!
  !! Frees memory allocated for HB0 horizontal no-flow boundaries
  !! and strings of vertices and the HB0 Collection object
  !!
  !! Calling Sequence:
  !!     call HB0_Destroy(ls0)
  !!
  !! Arguments:
  !!  type ( HB0_COLLECTION ),pointer :: hb0
  !!              Pointer to the HB0_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  integer (kind=ModAEM_Integer):: i
  type (hb0_STRING), pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                   "HB0_Destroy: HB0_Create has not been called",io )) return
  endif

  ! First deallocate each string of vertices
  do i = 1, hb0%iNStr
    str => hb0%Strings(i)
    deallocate (str%Vertices,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                   "hb0_Destroy: deallocation of Vertices failed",io )) return
  end do
  ! Then deallocate the strings table
  if ( associated(hb0%Strings) ) then
    deallocate (hb0%Strings,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                   "HB0_Destroy: deallocation of Strings failed",io )) return
  endif
  ! Now deallocate the collection
  deallocate (hb0,stat=iStat)
  if (IO_Assert( (iStat == 0), "HB0_Destroy: deallocation failed",io )) return
  
  return
end subroutine HB0_Destroy

subroutine HB0_New(hb0,Vertices,iNPts,io) 
  !! function HB0_New
  !!
  !! Adds a new HB0_STRING object to the HB0_COLLECTION 'hb0'
  !!
  !! Calling Sequence:
  !!    call HB0_New(hb0,Vertices,iNPt)
  !!
  !! Arguments:
  !!    (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!              The HB0_COLLECTION object to be used
  !!    (in)    type ( HB0_VERTEX ) :: Vertices(:)
  !!              Vector that defines the points along the barrier
  !!    (in)    integer :: iNPt
  !!              The number of vertices in the string
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( HB0_VERTEX ),dimension(:) :: Vertices
  integer (kind=ModAEM_Integer),intent(in) :: iNPts
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  type ( HB0_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_New: HB0_Create has not been called",io )) return
  endif

  if (IO_Assert( (hb0%iNStr < size(hb0%Strings)), &
                  "HB0_New: Space exhausted",io )) return
  if (IO_Assert( (iNPts <= size(Vertices)), &
                  "HB0_New: Size of provided vertices is inconsistent",io )) return

  hb0%iNStr = hb0%iNStr + 1
  str => hb0%Strings(hb0%iNStr)
  allocate ( str%Vertices(iNPts),stat=iStat )
  if (IO_Assert( (iStat==0), "HB0_New: Allocation failed",io )) return
  str%Vertices = Vertices(1:iNPts)
  str%iNPts = iNPts

  return
end subroutine HB0_New

function iHB0_GetInfo(hb0,iOption,io) result(iValue)
  !! function HB0_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iHB0_GetInfo(hb0,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue
  integer (kind=ModAEM_Integer) :: iStr
  type ( HB0_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_GetInfo: HB0_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FDP )
      do iStr=1,hb0%iNStr
        str => hb0%Strings(hb0%iNStr)
        iValue = iValue + str%iNPts-1
      end do
    case ( SIZE_EQUATIONS )
      do iStr=1,hb0%iNStr
        str => hb0%Strings(hb0%iNStr)
        iValue = iValue + 2*str%iNPts-1
      end do
    case ( SIZE_UNKNOWNS )
      do iStr=1,hb0%iNStr
        str => hb0%Strings(hb0%iNStr)
        iValue = iValue + 2*str%iNPts-1
      end do
    case ( INFO_REGENERATE )
      iValue = hb0%iRegenerate
    case default
      iValue = 0
  end select

  return
end function iHB0_GetInfo

subroutine HB0_SetupFunctions(hb0,fdp,io)
  !! subroutine HB0_SetupFunctions
  !!
  !! This routine sets up the functions in f_well and f_dipole for the no-flows
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call HB0_Setup(hb0)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZ3
  complex (kind=ModAEM_Real),dimension(6) :: cCPResult1,cCPResult2
  complex (kind=ModAEM_Real),dimension(3) :: cCPVtx
  character (len=255) :: sBuf
  integer (kind=ModAEM_Integer) :: irv
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_Setup: HB0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "HB0_Setup: Illegal FDP_COLLECTION object",io )) return
  endif

  ! Build dipoles for all segments
  do iStr = 1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx = 1,str%iNPts-1
      vtx => str%Vertices(iVtx)
      cZ1 = vtx%cZ
      cZ2 = str%Vertices(iVtx+1)%cZ
      call FDP_New(fdp,cZ1,cZ2,(/cZERO,cZERO,cZERO/),vtx%iFDPIndex,io)
      if ( io%lError ) return
    end do
  end do

  return
end subroutine HB0_SetupFunctions

subroutine HB0_SetupMatrix(hb0,mat,io)
  !! subroutine HB0_SetupMatrix
  !!
  !! This routine sets up the matrix entries for the module
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call HB0_Setup(hb0)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer :: mat
  !!             MAT_MATRIX object to be used
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZ3
  complex (kind=ModAEM_Real),dimension(6) :: cCPResult1,cCPResult2
  complex (kind=ModAEM_Real),dimension(3) :: cCPVtx
  character (len=255) :: sBuf
  integer (kind=ModAEM_Integer) :: irv
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_Setup: HB0_Create has not been called",io )) return
    if (IO_Assert( (associated(mat)), &
                    "HB0_Setup: Illegal MAT_MATRIX object",io )) return
  endif

  ! Build matrix entries for all segments
  do iStr = 1,hb0%iNStr
    str => hb0%Strings(iStr)
    ! Set up the unknown variables 
    ! Vertex entry -- all vertices except the first and last
    ! No dipole strength at either end of the string
    do iVtx = 1,str%iNPts
      vtx => str%Vertices(iVtx)
      ! Make a doublet vertex variable
      call MAT_CreateVariable(mat,ELEM_HB0,iStr,iVtx,kHB0_Vertex,io)
      if ( io%lError ) return
       ! Make a doublet center variable for all but the last entry...
      if ( iVtx < str%iNPts ) then
        call MAT_CreateVariable(mat,ELEM_HB0,iStr,iVtx,kHB0_Center,io)
        if ( io%lError ) return
      endif
    end do
   
    ! Set up control point sets and equations 
    do iVtx = 1,str%iNPts
      vtx => str%Vertices(iVtx)
      ! Now, create the equation entries...
      ! Two equations only for vertices 1 - #Vertices
      if ( iVtx == 1 ) then
        cZ1 = vtx%cZ
        cZ2 = str%Vertices(iVtx+1)%cZ
        ! Compute control intervals
        call MAT_ComputeControlPoints(cZ1,cZ2,2,cCPResult1,rZero,io)
        if ( io%lError ) return
 
        ! Vertex entry
        call MAT_CreateEquation(mat,(/ cCPResult1(1),cCPResult1(2) /),EQN_FLOW,ELEM_HB0, &
                                iStr,iVtx,kHB0_Vertex,rZERO,cZ2-cZ1,io)
        if ( io%lError ) return
 
        ! Doublet center control interval
        call MAT_CreateEquation(mat,(/ cCPResult1(2),cCPResult1(3) /),EQN_FLOW,ELEM_HB0, &
                                iStr,iVtx,kHB0_Center,rZERO,cZ2-cZ1,io)
        if ( io%lError ) return
 
      else if ( iVtx == str%iNPts ) then
        ! ONLY create a vertex control interval for the last vertex
        cZ1 = str%Vertices(iVtx-1)%cZ
        cZ2 = vtx%cZ

        ! Compute control intervals
        call MAT_ComputeControlPoints(cZ1,cZ2,2,cCPResult1,rZero,io)
        if ( io%lError ) return
 
        ! Vertex entry
        call MAT_CreateEquation(mat,(/ cCPResult1(3),cCPResult1(4) /),EQN_FLOW,ELEM_HB0, &
                               iStr,iVtx,kHB0_Vertex,rZERO,cZ2-cZ1,io)
        if ( io%lError ) return
 
      else
        ! Compute control points for the segment.  There are two unknowns
        ! per segment, and we need two control points per segment (no overspecification
        ! in this version)
        cZ1 = str%Vertices(iVtx-1)%cZ
        cZ2 = vtx%cZ
        cZ3 = str%Vertices(iVtx+1)%cZ
        call MAT_ComputeControlPoints(cZ1,cZ2,2,cCPResult1,rZero,io)
        if ( io%lError ) return
        call MAT_ComputeControlPoints(cZ2,cZ3,2,cCPResult2,rZero,io)
        if ( io%lError ) return
 
        ! Vertex strength control interval is based on neighboring doublets...
        call MAT_CreateEquation(mat,(/ cCPResult1(3),cCPResult1(4),cCPResult2(1),cCPResult2(2) /), &
                               EQN_FLOW,ELEM_HB0,iStr,iVtx,kHB0_Vertex,rZERO,cZ2-cZ1,io)
        if ( io%lError ) return
 
        ! Center equations 
        call MAT_CreateEquation(mat,(/ cCPResult2(2),cCPResult2(3) /),EQN_FLOW,ELEM_HB0, &
                               iStr,iVtx,kHB0_Center,rZERO,cZ2-cZ1,io)
        if ( io%lError ) return
 
      endif 
    end do
  end do
end subroutine HB0_SetupMatrix

subroutine HB0_Prepare(hb0,io) 
  !! subroutine HB0_Prepare
  !! 
  !! Prepares the module for a new iteration
  !!
  !! Do-nothing for m_hb0
  !!
  !! Calling Sequence:
  !!    call HB0_Setup(hb0,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer
  !!             HB0_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( IO_STATUS ),pointer :: io

  return
end subroutine HB0_Prepare

subroutine HB0_ComputeCoefficients(hb0,fdp,cPathZ,iaEqType,iElementType,IElementString, &
                                   iElementVertex,iElementFlag,caOrientation,rARow,io)
  !! subroutine HB0_ComputeCoefficients
  !!
  !! Computes a row of matrix coefficients (with no corrections) for the HB0 
  !! elements in layer iL.
  !!
  !! Calling Sequence:
  !!    call HB0_ComputeCoefficients(hb0,cPathZ,iEqType,cOrientation,rRow)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
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
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( FDP_COLLECTION ),pointer :: fdp
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real), intent(in) :: caOrientation
  integer (kind=ModAEM_Integer), intent(in) :: iaEqType
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real),dimension(:),intent(out) :: rARow
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat,iCol,iStr,iVtx,iDP1,iNDP,iWhich,irv
  complex (kind=ModAEM_Real),dimension(:,:,:),allocatable :: cDPF
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_ComputeCoefficients: HB0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "HB0_ComputeCoefficients: Illegal FDP_COLLECTION object",io )) return
  endif

  iCol = 0
  rARow = rZERO
  do iStr=1,hb0%iNStr
    str => hb0%Strings(iStr)
    ! Assume: the HB0_Setup routine creates consecutive dipole entries 
    iDP1 = str%Vertices(1)%iFDPIndex
    iNDP = str%iNPts-1
    allocate (cDPF(0:iNDP+1,3,1),stat=iStat)
    if (IO_Assert( (iStat==0), "HB0_ComputeCoefficients: Allocation failed",io )) return

    ! Get the appropriate influence functions for the boundary condition type
    select case ( iaEqType )
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

    ! Now, compute the matrix coefficients
    ! Note that xxxGetInfluence returns zero for elements 0 and iNDP+1 of the cDPFx vectors
    cDPF = cZERO
    call FDP_GetInfluence_IDP(fdp,iWhich,iDP1,iNDP,cPathZ,caOrientation,cDPF(1:iNDP,:,:),io)
    if ( io%lError ) return
    do iVtx=1,iNDP+1
      ! Vertex coefficient
      iCol = iCol+1
      rARow(iCol) = -aimag( cDPF(iVtx-1,3,1) + cDPF(iVtx,1,1) )
      ! Center coefficient for all but the last
      if ( iVtx <= iNDP ) then
        iCol = iCol+1
        rARow(iCol) = -aimag( cDPF(iVtx,2,1) )
      endif
    end do

    ! No memory leaks, please!
    deallocate(cDPF)
  end do
end subroutine HB0_ComputeCoefficients

function rHB0_ComputeRHS(hb0,iElementType,iElementString,iElementVertex, &
                         iElementFlag,rSpecValue,rCheck,io) result(rRHS)
  !! function rHB0_ComputeRHS
  !!
  !! Computes the right-hand side value for the solution
  !!
  !! Calling Sequence:
  !!   rRHS = rHB0_ComputeRHS(hb0,rValue,iElementType,iElementString,iElementVertex, &
  !!                          iElementFlag,rSpecValue,rCheck)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
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
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real), intent(in) :: rSpecValue
  real (kind=ModAEM_Real), intent(in) :: rCheck
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rRHS

  ! For LS1, this is easy...
  rRHS = rSpecValue - rCheck

  return
end function rHB0_ComputeRHS

subroutine HB0_StoreResult(hb0,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
  !! subroutine HB0_StoreResult
  !!
  !! Stores the results of a solution for a single equation associated with
  !! the HB0 module.
  !!
  !! Calling Sequence:
  !!    HB0_StoreResult(hb0,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    real :: rValue
  !!             The new result value from the solution vector
  !!   (in)    integer :: iElementType 
  !!             Element type (always ELEM_HB0)
  !!   (in)    integer :: iElementString
  !!             Element string number
  !!   (in)    integer :: iElementVertex
  !!             Element vertex number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!             For HB0, the constants kHB0_Vertex and kHBO_Center are used.
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  real (kind=ModAEM_Real), intent(in) :: rValue
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_StoreResult: HB0_Create has not been called",io )) return
    if (IO_Assert( (iElementString>=1 .and. iElementString<=hb0%iNStr), &
                    "HB0_StoreResult: Bad element string ID",io))  return
  endif

  str => hb0%Strings(iElementString)

  if ( io%lDebug ) then
    if (IO_Assert( (iElementVertex>=1 .and. iElementVertex<=str%iNPts), &
                    "HB0_StoreResult: Bad element vertex ID",io)) return
  endif

  ! All is well.  Store the result...
  vtx => str%Vertices(iElementVertex)
  select case ( iElementFlag ) 
    case ( kHB0_Vertex )
      vtx%rVertexStrength = vtx%rVertexStrength + rValue
    case ( kHB0_Center )
      vtx%rCenterStrength = vtx%rCenterStrength + rValue
  end select

  ! HB0 doesn't (yet) force a regeneration
  hb0%iRegenerate = 0

  return
end subroutine HB0_StoreResult

subroutine HB0_Update(hb0,fdp,io)
  !! subroutine HB0_StoreResult
  !!
  !! Updates the underlying function objects for the specified layer.
  !!
  !! Calling Sequence:
  !!    HB0_Update(hb0)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx,irv
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_Update: HB0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "HB0_Update: Illegal FDP_COLLECTION object",io )) return
  endif

  do iStr = 1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx = 1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      cRho1 = cmplx(rZERO,this_vtx%rVertexStrength,ModAEM_Real)
      cRho2 = cmplx(rZERO,this_vtx%rCenterStrength,ModAEM_Real)
      if ( iVtx == str%iNPts ) then
        cRho3 = cZERO
      else
        cRho3 = cmplx(rZERO,next_vtx%rVertexStrength,ModAEM_Real)
      endif
      call FDP_Update(fdp,this_vtx%iFDPIndex,(/cRho1,cRho2,cRho3/),io)
      if ( io%lError ) return
    end do
  end do
end subroutine HB0_Update

function lHB0_CheckPoint(hb0,cZ,rTol,io) result(lSing)
  !! logical function lHB0_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lHB0_CheckPoint(hb0,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             The HB0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ io%lDebug ]
  type (HB0_COLLECTION),pointer :: hb0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  ! NOTE: Singularities are found at the end points of the BDY elements 
  lSing = .false.
  do iStr = 1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx = 1,str%iNPts
      vtx => str%Vertices(iVtx)
      if ( abs(real(vtx%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(vtx%cZ-cZ))<rTol ) then
        lSing = .true.
        exit
      endif
    end do
  end do
  
  return
end function lHB0_CheckPoint

function lHB0_CheckProximity(hb0,cZ,rProximity,iDir,iElementID,iElementVtx,io) result(lFound)
  !! function lHB0_CheckProximity
  !!
  !! Checks to see if the point specified is near a no-flow, then
  !! moves the point a small distance if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    lFound = HB0_CheckProximity(hb0,cZ,rProximity,iDir,iElementID)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (in)    integer :: iDir
  !!             Flag for checking no-flows:
  !!               iDir =  0   -  Check without regard to the direction of travel
  !!               iDir = <0   -  Check for back tracing (ignore gaining no-flows)
  !!               iDir = >0   -  Check for formard tracing (ignore losing no-flows)
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the no-flow found (if lFound is .true.)
  !!   (out)   integer :: iElementVtx
  !!             The vertex associated with the string segment which was found
  !!   (out)   logical :: lFound
  !!             .true. if a no-flow is found
  !!             .false. if a no-flow is not found
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementVtx
  type ( IO_Status),pointer :: io 

  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lFound
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist
  complex (kind=ModAEM_Real) :: cZL,cZC
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_CheckProximity: HB0_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iStr=1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx=1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      ! Compute the mapped Z-value and then the distance from the linesink
      cZC = rHALF * (next_vtx%cZ + this_vtx%cZ )
      cZL = rHALF * (next_vtx%cZ - this_vtx%cZ )
      rDist = abs( aimag( (cZ-cZC) / cZL ) * abs(cZL) )
      if ( rDist < rProximity ) then
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
end function lHB0_CheckProximity

function lHB0_CheckCrossing(hb0,cZO,cZN,rProximity,iElementID,iElementVtx,cZInt,io) result(lFound)
  !! function lHB0_CheckCrossing
  !!
  !! Checks to see if the line segment cZ0-cZN intersects the barrier. If an 
  !! intersection exists, sets lFound=.true. and also returns the string ID
  !! in iElementID, the vertex in iElementVtx and the location of the crossing 
  !! in cZInt.  If no intersection is found, sets lFound=.false.
  !!
  !! Calling Sequence:
  !!    lFound = lHB0_CheckCrossing(hb0,cZ0,cZN,rProximity,iElementID,iElementVtx,cZInt)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    complex :: cZ0,cZN
  !!             The line segment to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the no-flow found (if lFound is .true.)
  !!   (out)   integer :: iElementVtx
  !!             The vertex associated with the string segment which was found
  !!   (out)   complex :: cZInt
  !!             The location of the intersection (if found)
  !!   (out)   logical :: lFound
  !!             .true. if a barrier is found
  !!             .false. if a barrier is not found
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  complex (kind=ModAEM_Real),intent(in) :: cZO
  complex (kind=ModAEM_Real),intent(in) :: cZN
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementVtx
  complex (kind=ModAEM_Real),intent(out) :: cZInt
  type ( IO_Status),pointer :: io 

  ! [ RETUNR VALUE ]
  logical (kind=ModAEM_Integer) :: lFound
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist,rXInt
  complex (kind=ModAEM_Real) :: cZL,cZC,cBZO,cBZN
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_CheckCrossing: HB0_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iStr=1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx=1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      cZC = rHALF * (next_vtx%cZ + this_vtx%cZ )
      cZL = rHALF * (next_vtx%cZ - this_vtx%cZ )
      cBZO = (cZO-cZC) / cZL 
      cBZN = (cZN-cZC) / cZL 
      ! If product of the imaginary parts < 0, then there is a potential crossing!
      if ( aimag(cBZO) * aimag(cBZN) < rZERO ) then
        rXInt = real(cBZO) - aimag(cBZO) * real(cBZN-cBZO) / aimag(cBZN-cBZO)
        if ( abs(rXInt) <= 1.0_ModAEM_Real ) then
          iElementID = str%iID
          iElementVtx = iVtx
          cZInt = (next_vtx%cZ + this_vtx%cZ ) * rXInt / 2.0
          lFound = .true.
        endif
      endif
    end do
  end do

  return
end function lHB0_CheckCrossing

subroutine HB0_Read(hb0,io)
  !! subroutine HB0_Read
  !!
  !! Populates an HB0_COLLECTION with input from kIOInputLU
  !!
  !! Calling Sequence:
  !!    call HB0_Read(hb0)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ) :: hb0
  !!             HB0_COLLECTION to be populated
  !!
  !! The format of the HB0 section of the input file appears as follows:
  !!
  !! HB0
  !! STR NVertices ID
  !! x y 
  !! ... Up to NVertices
  !! STR NVertices ID
  !! x y s
  !! ... Up to NVertices
  !! ... Up to NStrings
  !!
  !! NOTE: It is assumed that the HB0 line was found by the caller
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( IO_Status),pointer :: io 

  ! Locals -- for Directive parsing
  type (DIRECTIVE),dimension(4),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD,dirSTR /)
  ! Locals -- Input values
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  complex (kind=ModAEM_Real) :: cZ
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMaxStr, iMaxVtx
  integer (kind=ModAEM_Integer) :: iID
  integer (kind=ModAEM_Integer) :: iStr
  logical (kind=ModAEM_Integer) :: lFlag
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  call IO_MessageText("  Reading HB0 module input",io)
  if ( io%lError ) return
 
  if (IO_Assert( (associated(hb0)), "HB0_Read: HB0_Create has not been called",io )) return

  ! Use ifIO_InputRecord to process the model input file.
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "HB0_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModHB0 "ifXXXRead" routines. 
        if (IO_Assert( .false., "HB0_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(str)), "HB0_Read: No current string",io )) return
        if (IO_Assert( (str%iNPts < size(str%Vertices)), "HB0_Read: Space exhausted",io )) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ
        if (IO_Assert( (iStat==0), "HB0_Read: I/O Error",io )) return
        str%iNPts = str%iNPts+1
        vtx => str%Vertices(str%iNPts)
        vtx%cZ = cZ
        vtx%iFDPIndex = -1
        vtx%rVertexStrength = rZERO
        vtx%rCenterStrength = rZERO
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
        ! Here for the STR command -- create a new barrier string
        ! the maximum number of vertices is in the input record
        !****************************************************************************
        if (IO_Assert( (associated(hb0%Strings)), "HB0_Read: No strings allocated",io )) return
        if (IO_Assert( (hb0%iNStr<size(hb0%Strings)), "HB0_Read: Space exhausted",io )) return
        ! Retrive the number of vertices desired...
        read ( unit=sRecord, &
               fmt=*, &
               iostat=iStat &
             ) iMaxVtx,iID
        if (IO_Assert( (iStat==0), "HB0_Read: I/O Error",io )) return
        if (IO_Assert( (iMaxVtx>0), "HB0_Read: Illegal number of vertices",io )) return
        hb0%iNStr = hb0%iNStr+1
        allocate (hb0%Strings(hb0%iNStr)%Vertices(iMaxVtx),stat=iStat)
        if (IO_Assert( (iStat==0), "HB0_Read: Allocation failed",io )) return
        ! Made it!
        hb0%Strings(hb0%iNStr)%iID = iID
        write ( unit=sMessage, &
                fmt="(""HB0_Read: "",i6,"" vertices allocated"")" &
              ) iMaxVtx
        call IO_MessageText(sMessage,io)
        if ( io%lError ) return
        hb0%Strings(hb0%iNStr)%iNPts = 0      ! Initialize the vertex counter
        str => hb0%Strings(hb0%iNStr)
    end select
  end do  

  call IO_MessageText("  Leaving HB0 module",io)

  return
end subroutine HB0_Read

subroutine HB0_Inquiry(hb0,iLU,io)
  !! subroutine HB0_Inquiry
  !!
  !! Writes an inquiry report for all barriers to iLU
  !!
  !! Calling Sequence:
  !!    call HB0_Inquiry(hb0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rLength
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_Inquiry: HB0_Create has not been called",io )) return
  endif

  do iStr=1,hb0%iNStr
    str => hb0%Strings(iStr)
    do iVtx=1,str%iNPts
      vtx => str%Vertices(iVtx)
      if ( iVtx < str%iNPts ) then
        rLength = abs(str%Vertices(iVtx+1)%cZ - vtx%cZ)
      else
        rLength = rZERO
      endif
      write ( unit=iLU, &
              fmt="(""HB0"",3("","",i9),5("","",e14.6))" &
            ) str%iID,iVtx,vtx%cZ,rLength,vtx%rVertexStrength,vtx%rCenterStrength
    end do
  end do

  return
end subroutine HB0_Inquiry

subroutine HB0_Report(hb0,io)
  !! subroutine HB0_Report
  !!
  !! Writes a debugging report for all no-flows to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call HB0_Report(hb0)
  !!
  !! Arguments:
  !!   (in)    type ( HB0_COLLECTION ),pointer :: hb0
  !!             HB0_COLLECTION object to be used
  !!
  ! [ io%lDebug ]
  type ( HB0_COLLECTION ),pointer :: hb0
  type ( IO_Status),pointer :: io 

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( HB0_STRING ),pointer :: str
  type ( HB0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(hb0)), &
                    "HB0_Report: HB0_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module HB0'//'No-flow boundary information')" )

  if ( .not. associated(hb0%Strings) ) then
    write ( unit=kIOOutputLU, &
            fmt="(""No Strings Allocated""/)" &
          )
  else
    write ( unit=kIOOutputLU, &
            fmt="(""Number of strings: "",i5,""  used: "",i5)" &
          ) ubound(hb0%Strings,1),hb0%iNStr

    do iStr=1,hb0%iNStr
    str => hb0%Strings(iStr)
      write ( unit=kIOOutputLU, &
              fmt="(/100(""-"")/"" String: "",i10,"" # Vertices "",i10)" &
            ) str%iID,str%iNPts
      if ( associated(str%Vertices) ) then
        write ( unit=kIOOutputLU, &
                fmt="(""Vertices:"")" &
              )
        write ( unit=kIOOutputLU, &
                fmt="(""       X"",t15,""       Y"",t30,""   FDP Index"",t40," // &
                """   Vertex Str"",t60,""   Center Str"")" &
              )
        do iVtx=1,str%iNPts
          vtx => str%Vertices(iVtx)
          write ( unit=kIOOutputLU, &
                  fmt="(1x,d12.5,t15,d12.5,t30,i10,t40,d12.5,t60,d12.5)" &
                ) vtx%cZ,vtx%iFDPIndex,vtx%rVertexStrength,vtx%rCenterStrength
        end do
      endif
    end do

    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iHB0_GetInfo(hb0,SIZE_FWL,io)
    if ( io%lError ) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iHB0_GetInfo(hb0,SIZE_FPD,io)
    if ( io%lError ) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iHB0_GetInfo(hb0,SIZE_FDP,io)
    if ( io%lError ) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iHB0_GetInfo(hb0,SIZE_EQUATIONS,io)
    if ( io%lError ) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iHB0_GetInfo(hb0,SIZE_UNKNOWNS,io)
    if ( io%lError ) return
   endif

  write ( unit=kIOOutputLU, fmt="(a132)" ) repeat("-",132)

  return
end subroutine HB0_Report

end module m_hb0

