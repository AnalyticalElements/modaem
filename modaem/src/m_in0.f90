module m_in0

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

  !! module m_in0
  !!
  !! Element module for 2-D inhomogeneities
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for inhomogeneity domains
  !! (no common boundaries at this time). The Domains are condomucted using
  !! line-doublets (implemented as line-dipoles of imaginary strength).

use u_constants
use u_io
use u_matrix
use f_dipole

implicit none

public

  private :: lIN0_PointInsideDomain, &
             lIN0_LineIntersectsDomain, &
             lIN0_DomainInsideDomain

  type :: IN0_VERTEX
    !! type IN0_VERTEX
    !!
    !! Type that holds information for one vertex along the edge of an inhomogeneity
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
  end type IN0_VERTEX

  type :: IN0_DOMAIN
    !! type IN0_DOMAIN
    !!
    !! Type that holds information for one inhomogeneity
    !!
    !! Members:
    !!   type ( IN0_VERTEX ),dimension(:),pointer :: Vertices
    !!     A vector of IN0_VERTEX objects
    !!   integer :: iInsideDomain
    !!     Number of the IN0_DOMAIN that this domain is inside
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iID
    !!     The ID number for the string
    !!   real :: rBase
    !!     Base elevation
    !!   real :: rThickness
    !!     Thickness of the in0ifer (set to a large value for unconfined flow)
    !!   real :: rHydCond
    !!     Hydraulic conductivity
    !!   real :: rPorosity
    !!     Porosity
    !!   real :: rConfPot
    !!     Potential at the top of the in0ifer
    !! 
    type ( IN0_VERTEX ),dimension(:),pointer :: Vertices     
    integer (kind=ModAEM_Integer) :: iInsideDomain           
    integer (kind=ModAEM_Integer) :: iNPts                   
    integer (kind=ModAEM_Integer) :: iID
    real (kind=ModAEM_Real) :: rBase
    real (kind=ModAEM_Real) :: rThickness
    real (kind=ModAEM_Real) :: rHydCond
    real (kind=ModAEM_Real) :: rPorosity
    real (kind=ModAEM_Real) :: rConfPot
  end type IN0_DOMAIN

  type :: IN0_COLLECTION
    !! type IN0_COLLECTION
    !!
    !! Type that holds information for all IN0 elements in a layer
    !!
    !! Members:
    !!   type ( IN0_DOMAIN ),dimension(:),pointer :: Domains
    !!     A vector of IN0_DOMAIN objects
    !!   integer :: iNDom
    !!     The number of Domains actually in use
    !! 
    type ( IN0_DOMAIN ),dimension(:),pointer :: Domains
    integer (kind=ModAEM_Integer) :: iNDom                     
    integer (kind=ModAEM_Integer) :: iRegenerate
  end type IN0_COLLECTION

  ! Matrix generator element flags 
  integer (kind=ModAEM_Integer),private,parameter :: kIN0_Vertex=1
  integer (kind=ModAEM_Integer),private,parameter :: kIN0_Center=2
  ! How far do I move the control points off the line segments?
  real (kind=ModAEM_Real),private,parameter :: HB0_NORMAL_OFFSET = rZERO

  real (kind=ModAEM_Real) :: MOVEFACTOR = 1.0001

contains

function IN0_Create(iNDom,rBase,rThickness,rHydCond,rPorosity,io) result (in0)
  !! function IN0_Create
  !!
  !! Creates a new IN0_COLLECTION object, and builds the infinite in0ifer domain
  !!
  !! Calling Sequence:
  !!    IN0 => IN0_Create(iNDom)
  !!
  !! Arguments:
  !!    (in)    integer :: iNDom
  !!              The number of Domains to make space for; this routine will
  !!              add one domain for the "infinite" outside domain
  !!    (in)    real :: rBase
  !!              The base elevation for the infinite domain.  NOTE: Module 
  !!              IN0 does not support base elevation inhomogeneities; this is 
  !!              a placeholder for future expansion.
  !!    (in)    real :: rThickness
  !!              The thickness of the infinite in0ifer
  !!    (in)    real :: rHydCond
  !!              The hydraulic conductivity of the infinite in0ifer
  !!    (in)    real :: rPorosity
  !!              The porosity of the infinite in0ifer
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  integer (kind=ModAEM_Integer),intent(in) :: iNDom
  real (kind=ModAEM_Real),intent(in) :: rBase
  real (kind=ModAEM_Real),intent(in) :: rThickness
  real (kind=ModAEM_Real),intent(in) :: rHydCond
  real (kind=ModAEM_Real),intent(in) :: rPorosity
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  type ( IN0_COLLECTION ),pointer :: in0

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  type ( IN0_DOMAIN ),pointer :: dom

  allocate (in0,stat=iStat)
  if (IO_Assert( (iStat == 0), "IN0_Create: allocation failed",io )) return
  allocate (in0%Domains(iNDom),stat=iStat)
  if (IO_Assert( (iStat == 0), "IN0_Create: allocation failed",io )) return

  in0%iNDom = 1
  in0%iRegenerate = 1

  ! Build the infinite in0ifer
  dom => in0%Domains(1)
  dom%rBase = rBase
  dom%rThickness = rThickness
  dom%rHydCond = rHydCond
  dom%rPorosity = rPorosity
  dom%iInsideDomain = 0
  dom%iNPts = 0
  ! This value is computed by IN0_Setup
  dom%rConfPot = rZERO

  return
end function IN0_Create

subroutine IN0_New(in0,Vertices,iNPts,rBase,rThickness,rHydCond,rPorosity,io) 
  !! function IN0_New
  !!
  !! Adds a new IN0_DOMAIN object to the IN0_COLLECTION 'IN0'
  !!
  !! Calling Sequence:
  !!    call IN0_New(in0,Vertices,iNPt)
  !!
  !! Arguments:
  !!    (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!              The IN0_COLLECTION object to be used
  !!    (in)    type ( IN0_VERTEX ) :: Vertices(:)
  !!              Vector that defines the points along the barrier
  !!    (in)    integer :: iNPt
  !!              The number of vertices in the string
  !!    (in)    real :: rBase
  !!              Base elevation
  !!    (in)    real :: rThickness
  !!              Aquifer thickness
  !!    (in)    real :: rHydCond
  !!              Hydraulic conductivity
  !!    (in)    real :: rPorosity
  !!              Effective porosity
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( IN0_VERTEX ),dimension(:) :: Vertices
  integer (kind=ModAEM_Integer),intent(in) :: iNPts
  real (kind=ModAEM_Real),intent(in) :: rBase
  real (kind=ModAEM_Real),intent(in) :: rThickness
  real (kind=ModAEM_Real),intent(in) :: rHydCond
  real (kind=ModAEM_Real),intent(in) :: rPorosity
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iDom
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_New: IN0_Create has not been called",io )) return
  endif

  if (IO_Assert( (in0%iNDom < size(in0%Domains)), &
                  "IN0_New: Space exhausted",io )) return
  if (IO_Assert( (iNPts <= size(Vertices)), &
                  "IN0_New: Size of provided vertices is inconsistent",io )) return

  in0%iNDom = in0%iNDom + 1
  dom => in0%Domains(in0%iNDom)
  allocate ( dom%Vertices(iNPts),stat=iStat )
  if (IO_Assert( (iStat==0), "IN0_New: Allocation failed",io )) return
  dom%Vertices = Vertices(1:iNPts)
  dom%iNPts = iNPts
  dom%rBase = rBase
  dom%rThickness = rThickness
  dom%rHydCond = rHydCond
  dom%rPorosity = rPorosity
  ! This value is computed in IN0_SetupFunctions
  dom%rConfPot = rZERO

  ! This version of IN0 does not support nested inhomogeneity domains
  ! Check to ensure that no nesting or overlapping of domains has occurred
  ! nesting data structures
  do iDom=2,in0%iNDom-1
    if (IO_Assert( (.not. lIN0_DomainOverlapsDomain(in0,iDom,in0%iNDom,io)), &
                    "IN0_New: New domain overlaps another domain",io )) return
    if (IO_Assert( (.not. lIN0_DomainInsideDomain(in0,iDom,in0%iNDom,io)), &
                   "IN0_New: New domain is inside another domain",io )) return
  end do

  return
end subroutine IN0_New

function iIN0_GetInfo(in0,iOption,io) result(iValue)
  !! function IN0_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iIN0_GetInfo(aqu,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             AQU_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_STATUS ),pointer :: io 
 ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_GetInfo: IN0_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FDP )
      do iDom=1,in0%iNDom
        dom => in0%Domains(iDom)
        iValue = iValue + dom%iNPts
      end do
    case ( SIZE_EQUATIONS )
      do iDom=1,in0%iNDom
        dom => in0%Domains(iDom)
        iValue = iValue + 2*dom%iNPts
      end do
    case ( SIZE_UNKNOWNS )
      do iDom=1,in0%iNDom
        dom => in0%Domains(iDom)
        iValue = iValue + 2*dom%iNPts
      end do
    case ( INFO_REGENERATE )
      iValue = in0%iRegenerate
    case default
      iValue = 0
  end select

  return
end function iIN0_GetInfo

subroutine IN0_SetupFunctions(in0,fdp,io)
  !! subroutine IN0_SetupFunctions
  !!
  !! This routine sets up the functions in f_well and f_dipole for the line-sinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call IN0_SetupFunctions(in0,fdp,io)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!   (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_STATUS ),pointer :: io 
 ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx,i,ii
  real (kind=ModAEM_Real) :: rSum
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZ3
  complex (kind=ModAEM_Real),dimension(6) :: cCPResult1,cCPResult2
  complex (kind=ModAEM_Real),dimension(3) :: cCPVtx
  character (len=255) :: sBuf
  integer (kind=ModAEM_Integer) :: irv
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx
  type ( IN0_VERTEX ) :: temp_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_Setup: IN0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "IN0_Setup: Illegal FDP_COLLECTION object",io )) return
  endif

  ! Compute the necessary constants
  do iDom = 1,in0%iNDom
    dom => in0%Domains(iDom)
    dom%rConfPot = rHALF * dom%rHydCond * dom%rThickness * dom%rThickness
  end do

  ! Build dipoles for all segments
  do iDom = 1,in0%iNDom
    dom => in0%Domains(iDom)

    ! Compute the confined potential for the domain
    dom%rConfPot = rHALF * dom%rHydCond * dom%rThickness * dom%rThickness

    do iVtx = 1,dom%iNPts
      vtx => dom%Vertices(iVtx)
      cZ1 = vtx%cZ
      if ( iVtx < dom%iNPts ) then
        cZ2 = dom%Vertices(iVtx+1)%cZ
      else
        cZ2 = dom%Vertices(1)%cZ
      endif
      call FDP_New(fdp,cZ1,cZ2,(/cZERO,cZERO,cZERO/),vtx%iFDPIndex,io)
      cZ1 = cZ2
    end do
  end do

  return
end subroutine IN0_SetupFunctions

subroutine IN0_SetupMatrix(in0,mat,io)
  !! subroutine IN0_SetupMatrix
  !!
  !! This routine sets up the matrix entries for the module
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call IN0_SetupiMatrix(in0)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!   (in)    type ( MAT_MATRIX ),pointer :: mat
  !!             MAT_MATRIX object to be used
  !!   (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io 
 ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZ3
  complex (kind=ModAEM_Real),dimension(6) :: cCPResult1,cCPResult2
  complex (kind=ModAEM_Real),dimension(3) :: cCPVtx
  character (len=255) :: sBuf
  integer (kind=ModAEM_Integer) :: irv
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_Setup: IN0_Create has not been called",io )) return
    if (IO_Assert( (associated(mat)), &
                    "IN0_Setup: Illegal MAT_MATRIX object",io )) return
  endif

  ! Build matrix entries for all segments
  do iDom = 1,in0%iNDom
    dom => in0%Domains(iDom)
    ! Set up the unknown variables 
    ! Vertex entry -- all vertices except the first and last
    ! No dipole strength at either end of the string
    do iVtx = 1,dom%iNPts
      vtx => dom%Vertices(iVtx)
      ! Make a doublet vertex variable
      call MAT_CreateVariable(mat,ELEM_IN0,iDom,iVtx,kIN0_Vertex,io)
      ! Make a doublet center variable for all but the last entry...
      call MAT_CreateVariable(mat,ELEM_IN0,iDom,iVtx,kIN0_Center,io)
    end do
   
    ! Set up control point sets and equations 
    do iVtx = 1,dom%iNPts
      vtx => dom%Vertices(iVtx)
      ! Now, create the equation entries...
      ! Two equations for each segment of the perimeter path
      cZ1 = vtx%cZ
      if ( iVtx < dom%iNPts ) then
        cZ2 = dom%Vertices(iVtx+1)%cZ
      else
        cZ2 = dom%Vertices(1)%cZ
      endif
      ! Compute control point locations
      call MAT_ComputeControlPoints(cZ1,cZ2,2,cCPResult1,HB0_NORMAL_OFFSET,io)

      ! Vertex entry
      call MAT_CreateEquation(mat,(/ cCPResult1(2) /),EQN_INHO,ELEM_IN0, &
                              iDom,iVtx,kIN0_Vertex,rZERO,cZ2-cZ1,io)

      ! Center entry
      call MAT_CreateEquation(mat,(/ cCPResult1(3) /),EQN_INHO,ELEM_IN0, &
                              iDom,iVtx,kIN0_Center,rZERO,cZ2-cZ1,io)

    end do
  end do
end subroutine IN0_SetupMatrix

subroutine IN0_Prepare(in0,io) 
  !! subroutine IN0_Prepare
  !! 
  !! Prepares the module for a new iteration
  !!
  !! Do-nothing for m_in0
  !!
  !! Calling Sequence:
  !!    call IN0_Setup(in0,aqu,mat)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( IO_STATUS ),pointer :: io

  return
end subroutine IN0_Prepare

subroutine IN0_ComputeCoefficients(in0,fdp,cPathZ,iEqType,iElementType,iElementString, &
                                   iElementVertex,iElementFlag,cOrientation,rARow,io)
  !! subroutine IN0_ComputeCoefficients
  !!
  !! Computes a row of matrix coefficients (with no corrections) for the IN0 
  !! elements in layer iL.
  !!
  !! Calling Sequence:
  !!    call IN0_ComputeCoefficients(in0,cPathZ,iEqType,cOrientation,rRow)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
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
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( FDP_COLLECTION ),pointer :: fdp
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  complex (kind=ModAEM_Real), intent(in) :: cOrientation
  integer (kind=ModAEM_Integer), intent(in) :: iEqType
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real),dimension(:),intent(out) :: rARow
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat,iSkipCol,iDom,iVtx,iDP1,iNDP,iWhich,irv,iThisDP,i,iBaseCol
  complex (kind=ModAEM_Real),dimension(:,:,:),allocatable :: cDPF
  complex (kind=ModAEM_Real),dimension(1,3,1) :: cDPJ
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_ComputeCoefficients: IN0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "IN0_ComputeCoefficients: Illegal FDP_COLLECTION object",io )) return
  endif

  iSkipCol = 0
  do iDom=2,in0%iNDom
    dom => in0%Domains(iDom)
    ! Assume: the IN0_Setup routine creates consecutive dipole entries
    iDP1 = dom%Vertices(1)%iFDPIndex
    iNDP = dom%iNPts
    allocate (cDPF(0:iNDP+1,3,1),stat=iStat)
    if (IO_Assert( (iStat==0), "IN0_ComputeCoefficients: Allocation failed",io )) return

    ! Get the appropriate influence functions for the boundary condition type
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

    ! Now, compute the matrix coefficients
    ! Note that xxxGetInfluence returns zero for elements 0 and iNDP+1 of the cDPFx vectors
    cDPF = cZERO
    call FDP_GetInfluence_IDP(fdp,iWhich,iDP1,iNDP,cPathZ,cOrientation,cDPF(1:iNDP,:,:),io)

    do iVtx=1,iNDP
      iBaseCol = 2*iVtx-1 + iSkipCol

      ! Compute the contributions for all line-doublets
      ! Vertex 1 contribution
      rARow(iBaseCol)     = rARow(iBaseCol)   - aimag(cDPF(iVtx,1,1))
      ! Center contribution
      rARow(iBaseCol+1)   = rARow(iBaseCol+1) - aimag(cDPF(iVtx,2,1))
      ! Vertex 2 contribution
      if ( iVtx < iNDP ) then
        rARow(iBaseCol+2) = rARow(iBaseCol+2) - aimag(cDPF(iVtx,3,1))
      else
        rARow(1+iSkipCol) = rARow(1+iSkipCol) - aimag(cDPF(iVtx,3,1))
      endif

      ! Do I compute the additional term for the inhomogeneity?
      if ( (iEqType == EQN_INHO) .and. &
           (iElementType == ELEM_IN0) .and. &
           (iElementString == iDom) .and. &
           (iElementVertex == iVtx) ) then
        iThisDP = dom%Vertices(iVtx)%iFDPIndex
        call FDP_GetInfluence_IDP(fdp,INFLUENCE_J,iThisDP,1,cPathZ(1:1),cOrientation,cDPJ,io)
        ! Vertex 1 contribution
        rARow(iBaseCol)     = rARow(iBaseCol)   + rIN0_ARecip(in0,iDom,io)*aimag(cDPJ(1,1,1))
        ! Center contribution
        rARow(iBaseCol+1)   = rARow(iBaseCol+1) + rIN0_ARecip(in0,iDom,io)*aimag(cDPJ(1,2,1))
        ! Vertex 2 contribution
        if ( iVtx < iNDP ) then
          rARow(iBaseCol+2) = rARow(iBaseCol+2) + rIN0_ARecip(in0,iDom,io)*aimag(cDPJ(1,3,1))
        else
          rARow(1+iSkipCol) = rARow(1+iSkipCol) + rIN0_ARecip(in0,iDom,io)*aimag(cDPJ(1,3,1))
        endif
      endif
    end do
    iSkipCol = iSkipCol + 2*dom%iNPts

    ! No memory leaks, please!
    deallocate(cDPF)
  end do

  return
end subroutine IN0_ComputeCoefficients

function rIN0_ComputeRHS(in0,iElementType,iElementString,iElementVertex, &
                         iElementFlag,rSpecValue,rCheck) result(rRHS)
  !! function rIN0_ComputeRHS
  !!
  !! Computes the right-hand side value for the solution
  !!
  !! Calling Sequence:
  !!   rRHS = rIN0_ComputeRHS(in0,rValue,iElementType,iElementString,iElementVertex, &
  !!                          iElementFlag,rSpecValue,rCheck)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
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
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real), intent(in) :: rSpecValue
  real (kind=ModAEM_Real), intent(in) :: rCheck
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rRHS

  ! For IN0, this is easy...
  rRHS = rSpecValue - rCheck

  return
end function rIN0_ComputeRHS

function rIN0_ComputeCheck(in0,fdp,iElementString,iElementVertex,cCPZ,rPot,io) result(rCheck)
  !! function rIN0_ComputeCheck
  !!
  !! Returns the check value for the specified domain and vertex at the point cZ
  !!
  !! Calling Sequence:
  !!    rCheck = rIN0_ComputeCheck(in0,iElementString,iElementVertex,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!   (in)    complex,dimension(:) :: cCPZ
  !!             Control point(s) to be used in coefficient calculations
  !!   (in)    integer :: iElementString
  !!             The domain number to be examined
  !!   (in)    integer :: iElementVertex
  !!             The vertex number to be examined
  !!   (in)    complex :: cCPZ(:)
  !!             The location for evaluating the jump
  !!   (in)    real :: rPot
  !!             The potential at cCPZ(1)
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  !! Return Value:
  !!   (out)   real :: rCheck
  !!             The check value
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( FDP_COLLECTION ),pointer :: fdp
  integer (kind=ModAEM_Integer),intent(in) :: iElementString
  integer (kind=ModAEM_Integer),intent(in) :: iElementVertex
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cCPZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_STATUS ),pointer :: io 
 ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rCheck
  ! [ LOCALS ]
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx
  integer (kind=ModAEM_Integer) :: idp

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_ComputeCheck: IN0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "IN0_ComputeCheck: Illegal FDP_COLLECTION object",io )) return
    if (IO_Assert( (iElementString>0 .and. iElementString<=in0%iNDom), &
                    "IN0_ComputeCheck: Illegal domain number",io)) return
  endif
  dom => in0%Domains(iElementString)

  if ( io%lDebug ) then
    if (IO_Assert( (iElementVertex>=1 .and. iElementVertex<=dom%iNPts), &
                    "IN0_Jump: Bad element vertex ID",io)) return
  endif
  vtx => dom%Vertices(iElementVertex)
  rCheck = rPot - rIN0_ARecip(in0,iElementString,io) * rFDP_PotentialJump(fdp,io,vtx%iFDPIndex,cCPZ(1))
  return
end function rIN0_ComputeCheck

subroutine IN0_StoreResult(in0,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
  !! subroutine IN0_StoreResult
  !!
  !! Stores the results of a solution for a single equation associated with
  !! the IN0 module.
  !!
  !! Calling Sequence:
  !!    IN0_StoreResult(in0,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    real :: rValue
  !!             The new result value from the solution vector
  !!   (in)    integer :: iElementType
  !!             Element type (always ELEM_IN0)
  !!   (in)    integer :: iElementString
  !!             Element string number
  !!   (in)    integer :: iElementVertex
  !!             Element vertex number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!             For IN0, the constants kIN0_Vertex and kHBO_Center are used.
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  real (kind=ModAEM_Real), intent(in) :: rValue
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_StoreResult: IN0_Create has not been called",io )) return
    if (IO_Assert( (iElementString>=1 .and. iElementString<=in0%iNDom), &
                    "IN0_StoreResult: Bad element string ID",io)) return
  endif

  dom => in0%Domains(iElementString)

  if ( io%lDebug ) then
    if (IO_Assert( (iElementVertex>=1 .and. iElementVertex<=dom%iNPts), &
                    "IN0_StoreResult: Bad element vertex ID",io)) return
  endif

  ! All is well.  Store the result...
  vtx => dom%Vertices(iElementVertex)
  select case ( iElementFlag ) 
    case ( kIN0_Vertex )
      vtx%rVertexStrength = vtx%rVertexStrength + rValue
    case ( kIN0_Center )
      vtx%rCenterStrength = vtx%rCenterStrength + rValue
  end select

  return
end subroutine IN0_StoreResult

subroutine IN0_Update(in0,fdp,io)
  !! subroutine IN0_StoreResult
  !!
  !! Updates the underlying function objects for the specified layer.
  !!
  !! Calling Sequence:
  !!    IN0_Update(in0)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx,irv
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_Update: IN0_Create has not been called",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "IN0_Update: Illegal FDP_COLLECTION object",io )) return
  endif

  do iDom = 1,in0%iNDom
    dom => in0%Domains(iDom)
    do iVtx = 1,dom%iNPts
      this_vtx => dom%Vertices(iVtx)
      if ( iVtx < dom%iNPts ) then
        next_vtx => dom%Vertices(iVtx+1)
      else
        next_vtx => dom%Vertices(1)
      endif
      cRho1 = cmplx(rZERO,this_vtx%rVertexStrength,ModAEM_Real)
      cRho2 = cmplx(rZERO,this_vtx%rCenterStrength,ModAEM_Real)
      cRho3 = cmplx(rZERO,next_vtx%rVertexStrength,ModAEM_Real)
      call FDP_Update(fdp,this_vtx%iFDPIndex,(/cRho1,cRho2,cRho3/),io)

    end do
  end do

  ! m_in0 doesn't (yet) force a regeneration
  in0%iRegenerate = 0

  return
end subroutine IN0_Update

function lIN0_CheckPoint(in0,cZ,rTol) result(lSing)
  !! logical function lIN0_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lIN0_CheckPoint(in0,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             The IN0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (IN0_COLLECTION),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  ! NOTE: Singularities are found at the end points of the BDY elements 
  lSing = .false.
  do iDom = 1,in0%iNDom
    dom => in0%Domains(iDom)
    do iVtx = 1,dom%iNPts
      vtx => dom%Vertices(iVtx)
      if ( abs(real(vtx%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(vtx%cZ-cZ))<rTol ) then
        lSing = .true.
        exit
      endif
    end do
  end do
  
  return
end function lIN0_CheckPoint

subroutine IN0_CheckProximity(in0,cZ,rProximity,iDir,iElementID,iElementVtx,lFound,io)
  !! subroutine IN0_CheckProximity
  !!
  !! Checks to see if the point specified is near a line-sink, then
  !! moves the point a small distance if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    call IN0_CheckProximity(in0,cZ,rProximity,iDir,iElementID,lFound)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
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
  !!    in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementVtx
  logical (kind=ModAEM_Integer),intent(out) :: lFound
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist
  complex (kind=ModAEM_Real) :: cZL,cZC
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_CheckProximity: IN0_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iDom=1,in0%iNDom
    dom => in0%Domains(iDom)
    do iVtx=1,dom%iNPts-1
      this_vtx => dom%Vertices(iVtx)
      next_vtx => dom%Vertices(iVtx+1)
      ! Compute the mapped Z-value and then the distance from the linesink
      cZC = rHALF * (next_vtx%cZ + this_vtx%cZ )
      cZL = rHALF * (next_vtx%cZ - this_vtx%cZ )
      rDist = abs( aimag( (cZ-cZC) / cZL ) * abs(cZL) )
      if ( rDist < rProximity ) then
        if ( rDist < rMinDist ) then
          rMinDist = rDist
          iElementID = dom%iID
          iElementVtx = iVtx
          lFound = .true.
        endif
      endif
    end do
  end do

  return
end subroutine IN0_CheckProximity

subroutine IN0_CheckCrossing(in0,cZO,cZN,rProximity,iElementID,iElementVtx,cZInt,lFound,io)
  !! subroutine IN0_CheckCrossing
  !!
  !! Checks to see if the line segment cZ0-cZN intersects the barrier. If an 
  !! intersection exists, sets lFound=.true. and also returns the string ID
  !! in iElementID, the vertex in iElementVtx and the location of the crossing 
  !! in cZInt.  If no intersection is found, sets lFound=.false.
  !!
  !! Calling Sequence:
  !!    call IN0_CheckCrossing(in0,cZ0,cZN,rProximity,iElementID,iElementVtx,cZInt,lFound)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    complex :: cZ0,cZN
  !!             The line segment to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the line-sink found (if lFound is .true.)
  !!   (out)   integer :: iElementVtx
  !!             The vertex associated with the string segment which was found
  !!   (out)   complex :: cZInt
  !!             The location of the intersection (if found)
  !!   (out)   logical :: lFound
  !!             .true. if a barrier is found
  !!             .false. if a barrier is not found
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZO
  complex (kind=ModAEM_Real),intent(in) :: cZN
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementVtx
  complex (kind=ModAEM_Real),intent(out) :: cZInt
  logical (kind=ModAEM_Integer),intent(out) :: lFound
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist,rXInt
  complex (kind=ModAEM_Real) :: cZL,cZC,cBZO,cBZN
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_CheckCrossing: IN0_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iDom=1,in0%iNDom
    dom => in0%Domains(iDom)
    do iVtx=1,dom%iNPts-1
      this_vtx => dom%Vertices(iVtx)
      next_vtx => dom%Vertices(iVtx+1)
      cZC = rHALF * (next_vtx%cZ + this_vtx%cZ )
      cZL = rHALF * (next_vtx%cZ - this_vtx%cZ )
      cBZO = (cZO-cZC) / cZL 
      cBZN = (cZN-cZC) / cZL 
      ! If product of the imaginary parts < 0, then there is a potential crossing!
      if ( aimag(cBZO) * aimag(cBZN) < rZERO ) then
        rXInt = real(cBZO) - aimag(cBZO) * real(cBZN-cBZO) / aimag(cBZN-cBZO)
        if ( abs(rXInt) <= 1.0_ModAEM_Real ) then
          iElementID = dom%iID
          iElementVtx = iVtx
          cZInt = (next_vtx%cZ + this_vtx%cZ ) * rXInt / 2.0
          lFound = .true.
        endif
      endif
    end do
  end do

  return
end subroutine IN0_CheckCrossing

subroutine IN0_Inquiry(in0,iLU,io)
  !! subroutine IN0_Inquiry
  !!
  !! Writes an inquiry report for all barriers to iLU
  !!
  !! Calling Sequence:
  !!    call IN0_Inquiry(in0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  real (kind=ModAEM_Real) :: rLength
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_Inquiry: IN0_Create has not been called",io )) return
  endif

  do iDom=1,in0%iNDom
    dom => in0%Domains(iDom)
    do iVtx=1,dom%iNPts
      vtx => dom%Vertices(iVtx)
      if ( iVtx < dom%iNPts ) then
        rLength = abs(dom%Vertices(iVtx+1)%cZ - vtx%cZ)
      else
        rLength = rZERO
      endif
      write ( unit=iLU, &
              fmt="('IN0',3(',',i9),5(',',e14.6))" &
            ) dom%iID,iVtx,vtx%cZ,rLength,vtx%rVertexStrength,vtx%rCenterStrength
    end do
  end do

  return
end subroutine IN0_Inquiry

subroutine IN0_Read(in0,io)
  !! subroutine IN0_Read
  !!
  !! Reads the aquifer information for layer iL from the input LU
  !!
  !! Calling Sequence:
  !!    call IN0_Read(in0)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: aqu
  !!             Layer number to be read
  !!
  !! The format of the aqu section of the input file appears as follows:
  !!
  !! IN0
  !! DOM nvertices base thickness hyd-cond porosity
  !!     (x,y) 
  !!     (x,y)
  !!     ... up to nvertices
  !! DOM ...
  !! ... up to ninho
  !! END
  !! END
  !!
  !! NOTE: It is assumed that the IN0 line was found by the caller
  !!
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( IO_STATUS ),pointer :: io 
  ! Locals -- for Directive parsing
  type (DIRECTIVE),dimension(4),parameter :: dirDirectives = &
                                             (/ dirEND,dirDBG,dirPCD,dirDOM /)
  ! Locals -- Input values
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iError
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMax
  real (kind=ModAEM_Real) :: rBase
  real (kind=ModAEM_Real) :: rThickness
  real (kind=ModAEM_Real) :: rHydCond
  real (kind=ModAEM_Real) :: rPorosity
  real (kind=ModAEM_Real) :: rRefHead
  complex (kind=ModAEM_Real) :: cZ
  logical (kind=ModAEM_Integer) :: lFlag
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  call IO_MessageText("Reading IN0 module input",io)

  if (IO_Assert( (associated(in0)), "IN0_Read: IN0_Create has not been called",io )) return

  do 
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "IN0_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all Mod "ifXXXRead" routines. 
        if (IO_Assert( .false., "IN0_Read: Unexpected EOF",io )) return
      case (kOpData)
        ! A data line was found. If we have a specified perimeter, add the point
        ! to the perimeter.
        if (IO_Assert( (associated(dom%Vertices)), "IN0_Read: No DOM directive",io )) return
        if (IO_Assert( (dom%iNPts < size(dom%Vertices)), &
                        "IN0_Read: Space exhausted",io )) return
        read ( unit=sRecord, &
               fmt=*, &
               iostat=IStat &
             ) cZ
        if (IO_Assert( (iStat==0), "IN0_Read: I/O Error",io )) return
        dom%iNPts = dom%iNPts+1
        vtx => dom%Vertices(dom%iNPts)
        vtx%cZ = cZ
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        return
      case (kOpDBG)
        ! Change the io%lDebug flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetDebug(lFlag,io)
      case (kOpPCD)
        ! Change the IO_Proceed flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetProceed(lFlag,io)
      case (kOpDOM)
        ! Start a new domain
        if (IO_Assert( (associated(in0%Domains)), &
                        "IN0_Read: No domains have been allocated",io )) return
        if (IO_Assert( (in0%iNDom< size(in0%Domains)), &
                        "IN0_Read: Space exhausted",io )) return
        read ( unit=sRecord, &
               fmt=*, &
               iostat=iStat &
             ) iMax,rBase,rThickness,rHydCond,rPorosity
        if (IO_Assert( (iStat==0), "IN0_Read: I/O Error",io )) return
        if (IO_Assert( (iMax>2), "IN0_Read: Illegal dimension",io )) return
        in0%iNDom= in0%iNDom+1
        dom => in0%Domains(in0%iNDom)
        dom%rBase = rBase
        dom%rThickness = rThickness
        dom%rHydCond = rHydCond
        dom%rPorosity = rPorosity
        ! FOR NOW, NO NESTING!
        dom%iInsideDomain = 1
        allocate (dom%Vertices(iMax),stat=iStat)
        if (IO_Assert( (iStat==0), "IN0_Read: Allocation failed",io )) return
        dom%Vertices = IN0_VERTEX(cZERO,rZERO,rZERO,0)
      case default
    end select
  end do

  call IO_MessageText("Leaving IN0 module",io)

  return
end subroutine IN0_Read  

subroutine IN0_Report(in0,io)
  !! subroutine IN0_Report
  !!
  !! Writes a debugging report for all line-sinks to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call IN0_Report(in0)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  type ( IO_STATUS ),pointer :: io 
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom,iVtx
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "IN0_Report: IN0_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module IN0'//'Inhomogeneity equation information')" )

  if ( .not. associated(in0%Domains) ) then
    write ( unit=kIOOutputLU, &
            fmt="('No Domains Allocated'/)" &
          )
  else
    write ( unit=kIOOutputLU, &
            fmt="('Number of Domains: ',i5,'  used: ',i5)" &
          ) ubound(in0%Domains,1),in0%iNDom

    do iDom=1,in0%iNDom
      write ( unit=kIOOutputLU, &
              fmt="(/'Domain ',i10/)" &
            ) iDom
      dom => in0%Domains(iDom)
      write ( unit=kIOOutputLU, &
              fmt="('Base      ',g13.5/'Thickness ',g13.5/'HydCond   ',g13.5/'Porosity  ',g13.5/)" &
            ) dom%rBase,dom%rThickness,dom%rHydCond,dom%rPorosity

      write ( unit=kIOOutputLU, &
              fmt="('Vertices:')" &
            )
      if ( .not. associated(dom%Vertices) ) then
        write ( unit=kIOOutputLU, &
                fmt="(' [None - infinite aquifer]'/)" &
              )
      else
        write ( unit=kIOOutputLU, &
                fmt="('       X',t15,'       Y',t30,'   FDP Index',t40," // &
                "'   Vertex Str',t60,'   Center Str')" &
              )
        do iVtx=1,dom%iNPts
          vtx => dom%Vertices(iVtx)
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
          ) iIN0_GetInfo(in0,SIZE_FWL,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iIN0_GetInfo(in0,SIZE_FPD,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iIN0_GetInfo(in0,SIZE_FDP,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iIN0_GetInfo(in0,SIZE_EQUATIONS,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iIN0_GetInfo(in0,SIZE_UNKNOWNS,io)
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine IN0_Report

function rIN0_ARecip(in0,iDomain,io) result(rARecip)
  !! function IN0_ARecip(in0,iDomain)
  !! 
  !! Returns the factor A*, defined as the reciprocal of (k+ - k-) / k-
  !! where k+ is the inside conductivity and k- is the outside conductivity
  !!
  !! Calling Sequence:
  !!    rARecip = lIN0_ARecip(in0,iDomain)
  !!
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      integer :: iDomain
  !!                Domain index in in0
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!    (out)     real :: rARecip
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iDomain
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rARecip
  ! [ LOCALS ]
  type ( IN0_DOMAIN ),pointer :: inside,outside

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_DomainInsideDomain: IN0_Create has not been called",io )) return
    if (IO_Assert( (iDomain >= 1 .and. iDomain <= in0%iNDom), &
                    "lIN0_DomainInsideDomain: Bad domain index",io )) return
  endif

  inside => in0%Domains(iDomain)
  outside => in0%Domains(inside%iInsideDomain)
  
  rARecip = inside%rHydCond / ( outside%rHydCond - inside%rHydCond )
  
  return
end function rIN0_ARecip

function lIN0_PointInsideDomain(in0,iDomain,cZ,io) result(lInside)
  !! function lIN0_PointInsideDomain
  !!
  !! Tests to see if the specified point is inside the Domain
  !!
  !! Calling Sequence:
  !!    lresult = lIN0_PointInsideDomain(in0,iDomain,cZ)
  !!
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      integer :: iDomain
  !!                Domain index in in0
  !!    (in)      complex :: cZ
  !!                The point to be checked
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           logical :: lInside
  !!             .true. if the point is inside the Domain,
  !!             .false. if the point is outside the Domain
  !!
  !! NOTE:  Uses an algorithm from the comp.graphics FAQ -- THANKS!
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iDomain
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  logical :: lInside
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat,i,ii
  real (kind=ModAEM_Real) :: rSum
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_DomainInsideDomain: IN0_Create has not been called",io )) return
    if (IO_Assert( (iDomain >= 1 .and. iDomain <= in0%iNDom), &
                    "lIN0_DomainInsideDomain: Bad domain index",io )) return
  endif

  ! Use the "arctangent method" for now; it's slow but reliable
  dom => in0%Domains(iDomain)
  if ( size(dom%Vertices) == 0 ) then
    lInside = .true.
  else
    if ( any((cZ==dom%Vertices(1:dom%iNPts)%cZ),1) ) then
      lInside = .true.
    else
      rSum= rZERO
      ii = dom%iNPTs
      do i=1,dom%iNPts
        rSum = rSum + aimag( log( (cZ-dom%Vertices(i)%cZ) / (cZ-dom%Vertices(ii)%cZ) )  )
        ii = i
      end do
      lInside = ( rSum > rONE ) 
    endif
  endif

  return
end function lIN0_PointInsideDomain

function lIN0_LineIntersectsDomain(in0,iDomain,cZ1,cZ2,io) result(lIntersects)
  !! function lIN0_LineIntersectsDomain
  !!
  !! Tests to see if the line segment cZ1-cZ2 intersects the Domain
  !!
  !! Calling Sequence:
  !!    lresult = lIN0_LineIntersectsDomain(in0,iDomain,cZ1,cZ2)
  !!
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      integer :: iDomain
  !!                Domain index in in0
  !!    (in)      complex :: cZ1, cZ2
  !!                End-points of the line to be tested
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           logical :: lIntersects
  !!             .true. if the line segment intersects with any Domain edge
  !!             .false. otherwise
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iDomain
  complex (kind=ModAEM_Real),intent(in) :: cZ1,cZ2
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  logical :: lIntersects
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,j
  complex (kind=ModAEM_Real) :: cMapZ1,cMapZ2
  real (kind=ModAEM_Real) :: rXInt
  type ( IN0_DOMAIN ),pointer :: dom
  type ( IN0_VERTEX ),pointer :: vtx_i,vtx_j
  
  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_DomainInsideDomain: IN0_Create has not been called",io )) return
    if (IO_Assert( (iDomain >= 1 .and. iDomain <= in0%iNDom), &
                    "lIN0_DomainInsideDomain: Bad domain index",io )) return
  endif

  dom => in0%Domains(iDomain)
  lIntersects = .false.
  j = size(dom%Vertices,1)
  do i=1,size(dom%Vertices,1)
    vtx_i => dom%Vertices(i)
    vtx_j => dom%Vertices(j)
    cMapZ1 = ( vtx_i%cZ - rHALF * (cZ2+cZ1) ) / ( rHALF * (cZ2-cZ1) )
    cMapZ2 = ( vtx_j%cZ - rHALF * (cZ2+cZ1) ) / ( rHALF * (cZ2-cZ1) )
    if ( (aimag(cMapZ1) >= rZERO .and. aimag(cMapZ2) < rZERO ) .or. &
         (aimag(cMapZ2) >= rZERO .and. aimag(cMapZ1) < rZERO) ) then
      rXInt = real(cMapZ1) - aimag(cMapZ1) * ((real(cMapZ2)-real(cMapZ1)) / &
                                              (aimag(cMapZ2)-aimag(cMapZ1)))
      if ( abs(rXInt) <= rONE ) then
        lIntersects = .true.
        return
      endif
    endif
    j = i
  end do
  
  return
end function lIN0_LineIntersectsDomain

function lIN0_DomainInsideDomain(in0,iDomain,iDomain2,io) result(lInside)
  !! function lIN0_DomainInsideDomain
  !!
  !! Tests to see if the Domain 'iDomain2' is inside the Domain 'iDomain'
  !!
  !! Calling Sequence:
  !!    lresult = lIN0_DomainInsideDomain(in0,iDomain,iDomain2)
  !!
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      integer :: iDomain
  !!                Domain index in in0
  !!    (in)      integer :: iDomain2
  !!                Domain index in in0
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           logical :: lInside
  !!             .true. if the Domain is completely inside the Domain,
  !!             .false. otherwise
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iDomain
  integer (kind=ModAEM_Integer),intent(in) :: iDomain2
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  logical :: lInside
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,j
  type ( IN0_DOMAIN ),pointer :: dom,dom2
  type ( IN0_VERTEX ),pointer :: vtx_i,vtx_j

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_DomainInsideDomain: IN0_Create has not been called",io )) return
    if (IO_Assert( ((iDomain >= 1 .and. iDomain <= in0%iNDom) .and. &
                     (iDomain2 >= 1 .and. iDomain2 <= in0%iNDom)), &
                    "lIN0_DomainInsideDomain: Bad domain index",io )) return
  endif

  dom => in0%Domains(iDomain)
  dom2 => in0%Domains(iDomain2)
  if ( size(dom%Vertices,1) == 0 ) then
    lInside = .true.
  else
    if ( size(dom2%Vertices) == 0 ) then
      lInside = .false.
    else
      lInside = .true.
      j = size(dom2%Vertices,1)
      do i=1,size(dom2%Vertices,1)
        vtx_i => dom2%Vertices(i)
        vtx_j => dom2%Vertices(j)
        if ( .not. lIN0_PointInsideDomain(in0,iDomain,vtx_i%cZ,io) .or. &
              lIN0_LineIntersectsDomain(in0,iDomain,vtx_i%cZ,vtx_j%cZ,io) ) then
          lInside = .false.
          return
        endif
        j = i
      end do
    endif
  endif

  return
end function lIN0_DomainInsideDomain

function lIN0_DomainOverlapsDomain(in0,iDomain,iDomain2,io) result(lOverlap)
  !! function lIN0_DomainOverlapsDomain
  !!
  !! Tests to see if the Domain 'iDomain2' overlaps the Domain 'iDomain'
  !!
  !! Calling Sequence:
  !!    lresult = lIN0_DomainOverlapsDomain(in0,iDomain,iDomain2)
  !!
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      integer :: iDomain
  !!                Domain index in in0
  !!    (in)      integer :: iDomain2
  !!                Domain index in in0
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           logical :: lInside
  !!             .true. if Domain2 overlaps Domain,
  !!             .false. otherwise
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  integer (kind=ModAEM_Integer),intent(in) :: iDomain
  integer (kind=ModAEM_Integer),intent(in) :: iDomain2
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  logical :: lOverlap
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,j
  type ( IN0_DOMAIN ),pointer :: dom,dom2
  type ( IN0_VERTEX ),pointer :: vtx_i,vtx_j

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_DomainInsideDomain: IN0_Create has not been called",io )) return
    if (IO_Assert( ((iDomain >= 1 .and. iDomain <= in0%iNDom) .and. &
                     (iDomain2 >= 1 .and. iDomain2 <= in0%iNDom)), &
                    "lIN0_DomainInsideDomain: Bad domain index",io )) return
  endif

  dom => in0%Domains(iDomain)
  dom2 => in0%Domains(iDomain2)
  if ( size(dom%Vertices,1) == 0 ) then
    lOverlap = .true.
  else
    if ( size(dom2%Vertices) == 0 ) then
      lOverlap = .false.
    else
      lOverlap = .false.
      j = size(dom2%Vertices,1)
      do i=1,size(dom2%Vertices,1)
        vtx_i => dom2%Vertices(i)
        vtx_j => dom2%Vertices(j)
        if ( lIN0_LineIntersectsDomain(in0,iDomain,vtx_i%cZ,vtx_j%cZ,io) ) then
          lOverlap = .true.
          return
        endif
        j = i
      end do
    endif
  endif

  return
end function lIN0_DomainOverlapsDomain

function IN0_FindDomain(in0,cZ,io) result(dom)
  !! function IN0_FindDomain
  !! 
  !! Returns a pointer to the domain containing the point 'cZ'
  !!
  !! 
  !! Arguments:
  !!    (in)      type ( IN0_COLLECTION ),pointer :: in0
  !!                The IN0_COLLECTION object of interest
  !!    (in)      complex :: cZ
  !!                The point in question
  !!
  !! Return Value:
  !!           type ( IN0_DOMAIN ),pointer :: dom
  !!             Pointer to the domain containing cZ
  !!
  !! Note:
  !!   This version does not support nested inhomogeneities
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  type ( IN0_DOMAIN ),pointer :: dom
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iDom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_FindDomain: IN0_Create has not been called",io )) return
  endif

  ! Default to the infinite domain
  dom => in0%Domains(1)
  do iDom = 2,in0%iNDom
    if ( lIN0_PointInsideDomain(in0,iDom,cZ,io) ) then
      dom => in0%Domains(iDom)
    endif
  end do
 
  return
end function IN0_FindDomain

!! UTILITY ROUTINES
!! These routines provide computational aids for AEM models

function rIN0_HeadToPotential(in0,rHead,cZ,io) result(rPot)
  !! real function rIN0_HeadToPotential
  !!
  !! Converts head to a discharge potential based on the in0ifer properties 
  !! at cZ. Returns the (real) potential.
  !!
  !! Calling Sequence:
  !!    rPot = rIN0_HeadToPotential(iL,rHead,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    real :: rHead
  !!             The head value to be converted
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!   (in)    type ( IO_status ), pointer :: io
  !!             pointer to IO_Status structure

  !!
  !! Return Value:
  !!           real :: rPot
  !!             The discharge potential corresponding to the head 'rHead' 
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  real (kind=ModAEM_Real),intent(in) :: rHead
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rPot
  ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rHd
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_FindDomain: IN0_Create has not been called",io )) return
  endif

  ! For the general case, rK, rBase, rPerm and rThick are functions of position
  dom => IN0_FindDomain(in0,cZ,io)

  rHd = rHead - dom%rBase
  if ( rHd > dom%rThickness ) then
    rPot = dom%rHydCond * dom%rThickness * rHd - &
           rHALF * dom%rHydCond * dom%rThickness * dom%rThickness
  else
    rPot = rHALF * dom%rHydCond * rHd * rHd
  endif

  return
end function rIN0_HeadToPotential

function rIN0_PotentialToHead(in0,rPot,cZ,io) result(rHead)
  !! real function rIN0_HeadToPotential
  !!
  !! Converts head to a discharge potential based on the in0ifer properties 
  !! at cZ. Returns the (real) potential.
  !!
  !! Calling Sequence:
  !!    rPot = rIN0_HeadToPotential(in0,rHead,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    real :: rHead
  !!             The head value to be converted
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!    (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           real :: rHead
  !!             The head corresponding to the discharge potential 'rPot' 
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  real (kind=ModAEM_Real),intent(in) :: rPot
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rHead
  ! [ LOCALS ]
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_FindDomain: IN0_Create has not been called",io )) return
  endif

  ! For the general case, rK, rBase, rPerm and rThick are functions of position
  dom => IN0_FindDomain(in0,cZ,io)

  if ( rPot < rZERO ) then
    rHead = -HUGE(ModAEM_Real)
  else
    if ( rPot > dom%rConfPot ) then                                                  
      rHead = (rPot + dom%rConfPot) / (dom%rHydCond * dom%rThickness) + dom%rBase
    else
      rHead = sqrt(rTWO * rPot / dom%rHydCond) + dom%rBase
    endif
  endif

  return
end function rIN0_PotentialToHead

function cIN0_DischargeToVelocity(in0,cDischarge,cZ,rPot,io) result(cVelocity)
  !! function cIN0_DischargeToVelocity
  !!
  !! Converts Discharge to velocity, using the potential and the location cZ
  !! to compute saturated thickness
  !!
  !! Calling Sequence:
  !!    cV = cIN0_DischargeToVelocity(in0,rDischarge,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    real :: rDischarge
  !!             The discharge value to be converted
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!   (in)    real :: rPot
  !!             The (real) potential at cZ
  !!   (in)    type ( IO_status ), pointer :: io
  !!              pointer to IO_Status structure
  !!
  !! Return Value:
  !!           complex :: cV
  !!             The velocity
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  real (kind=ModAEM_Real),intent(in) :: rPot
  complex (kind=ModAEM_Real),intent(in) :: cDischarge
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cVelocity
  ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rSatdThick
  type ( IN0_DOMAIN ),pointer :: dom

  if ( io%lDebug ) then
    if (IO_Assert( (associated(in0)), &
                    "lIN0_FindDomain: IN0_Create has not been called",io )) return
  endif

  ! For the general case, rK, rBase, rPerm and rThick are functions of position
  dom => IN0_FindDomain(in0,cZ,io)

  ! For the general case, rK, rBase, rPerm and rThick are functions of position, for now 
  ! the position parameters are ignored
  if ( rPot < rZERO ) then
    cVelocity = cZERO
    return
  else
    ! Find the saturated thickness
    if ( rPot > dom%rConfPot ) then                                                  
      rSatdThick = dom%rThickness 
    else
      rSatdThick = rIN0_PotentialToHead(in0,rPot,cZ,io) - dom%rBase
    endif

    ! Compute the velocity
    cVelocity = cDischarge / ( rSatdThick * dom%rPorosity )
  endif

  return
end function cIN0_DischargeToVelocity

function rIN0_SatdThickness(in0,cZ,rPot,io) result (rH)
  !! function cIN0_SatdThickness
  !!
  !! Computes the saturated thickness at point cZ in layer iL where the potential is rPot
  !!
  !! Calling Sequence:
  !!    cV = cIN0_SatdThickness(in0,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!   (in)    real :: rPot
  !!             The (real) potential at cZ
  !!
  !! Return Value:
  !!           complex :: cV
  !!             The velocity
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rH
  ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rHead
  type ( IN0_DOMAIN ),pointer :: dom

!$$ assert : allocated(Layers) : "No layers are allocated"
!$$ assert : (iL>=1 .and. iL<=size(Layers)) : "Wrong layer was specified"

  dom => IN0_FindDomain(in0,cZ,io)
  rHead = rIN0_PotentialToHead(in0,rPot,cZ,io)
  if ( (rHead - dom%rBase) < dom%rThickness ) then
    rH = rHead - dom%rBase
  else
    rH = dom%rThickness
  endif
  
  return
end function rIN0_SatdThickness

function rIN0_Transmissivity(in0,cZ,rPot,io) result (rT)
  !! function rIN0_Transmissivity
  !!
  !! Computes the transmissivity at the point cZ where the potential is rPot
  !!
  !! Calling Sequence:
  !!    cV = cIN0_Transmissivity(in0,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!   (in)    real :: rPot
  !!             The (real) potential at cZ
  !!
  !! Return Value:
  !!           complex :: cV
  !!             The velocity
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rT
  ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rHead
  type ( IN0_DOMAIN ),pointer :: dom

!$$ assert : allocated(Layers) : "No layers are allocated"
!$$ assert : (iL>=1 .and. iL<=size(Layers)) : "Wrong layer was specified"

  dom => IN0_FindDomain(in0,cZ,io)
  rHead = rIN0_PotentialToHead(in0,rPot,cZ,io)
  if ( (rHead - dom%rBase) < dom%rThickness ) then
    rT = dom%rHydCond * (rHead - dom%rBase)
  else
    rT = dom%rHydCond * dom%rThickness
  endif

  return
end function rIN0_Transmissivity

function rIN0_HydCond(in0,cZ,io) result (rHydCond)
  !! function rIN0_HydCond
  !!
  !! Computes the hydraulic conductivity at the point cZ
  !!
  !! Calling Sequence:
  !!    cV = cIN0_HydCond(in0,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( IN0_COLLECTION ),pointer :: in0
  !!             IN0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!
  !! Return Value:
  !!           real :: rHydCond
  !!             The hydraulic conductivity
  !!
  ! [ ARGUMENTS ]
  type ( IN0_COLLECTION ),pointer :: in0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_STATUS ),pointer :: io 
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rHydCond
  ! [ LOCALS ]
  type ( IN0_DOMAIN ),pointer :: dom

!$$ assert : allocated(Layers) : "No layers are allocated"
!$$ assert : (iL>=1 .and. iL<=size(Layers)) : "Wrong layer was specified"

  dom => IN0_FindDomain(in0,cZ,io)
  rHydCond = dom%rHydCond

  return
end function rIN0_HydCond

end module m_in0


