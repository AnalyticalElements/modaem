module m_aqu

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

  !! module m_aqu
  !!
  !! Element module for 2-D uniform (optionally bounded) aquifers
  !!
  !! Module use:
  !!   constants  --  Universal ModAEM constant declarations
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!   f_well     --  Function module for collections of wells
  !!   f_matrix   --  Matrix interfacing routines
  !!
  !! This module provides the necessary functionality for the overlying
  !! 2-D aquifer that works together with the other element modules
  !!
  !! The AQU module is a very special case in ModAEM.  It has three different
  !! "personalities", all in one module:
  !!     (a)  An ELEMENT MODULE that sets up equations to be solved
  !!     (b)  A FUNCTION MODULE that provides the low-level functions (like "f_aqu")
  !!     (c)  A UTILITY MODULE that allows for convenient AEM calculations e.g.
  !!          conversion of heads and potentials, looking up aquifer properties, etc.
  !! In addition, the AQU module must allow the other m_xxx element modules to access
  !! its "global" conversion functions.

use u_constants
use u_io
use u_matrix

use f_well
use f_dipole

use m_in0

implicit none

public

  type,public :: AQU_BDYELEMENT
    !! type AQU_BDYELEMENT
    !! 
    !! Type that holds a vertex for the aquifer perimeter
    !!
    !! Members:
    !!   complex :: cZ1
    !!     The complex coordinate of the first vertex (x,y)
    !!   complex :: cZ2
    !!     The complex coordinate of the second vertex (x,y)
    !!   real :: rSpecHead
    !!     The specified head for the segment
    !!   real :: rSpecFlux
    !!     The specified flux for the segment
    !!   logical :: lFlag
    !!     Flag : if .false., base computations on the specified flux
    !!            if .true.,  base computations on the specified head
    !!   integer :: iFDPIndex
    !!     Pointer in FDP to the dipole 
    !!   real :: rLength
    !!     Length of the line segment
    !!   real :: rStrength
    !!     Sink density of the line segment 
    !!   real :: rDPStrength
    !!     Dipole strength of the line segment
    !!   integer :: iFWLIndex
    !!     Index into FWL for the segment
    !!
    complex (kind=ModAEM_Real) :: cZ1
    complex (kind=ModAEM_Real) :: cZ2
    real (kind=ModAEM_Real) :: rSpecHead
    real (kind=ModAEM_Real) :: rSpecFlux
    real (kind=ModAEM_Real) :: rLength
    real (kind=ModAEM_Real) :: rStrength
    integer (kind=ModAEM_Integer) :: iFDPIndex
    integer (kind=ModAEM_Integer) :: iFWLIndex
    logical (kind=ModAEM_Integer) :: lFlag
  end type AQU_BDYELEMENT

  type,public :: AQU_COLLECTION
    !! type AQU_COLLECTION
    !!
    !! Type that holds information for a layer
    !!
    !! Members:
    !!   complex :: cRefPoint 
    !!     Location of the reference point (if specified)
    !!   real :: rRefHead
    !!     Specified head at the reference point (if specified)
    !!   complex :: cRefUniformFlow
    !!     Infinite aquifer uniform flow discharge vector
    !!   logical :: lReference
    !!     .true. if a reference point has been specified
    !!   real :: rSolConst
    !!     The constant of integration
    !!   type ( AQU_BDYELEMENT ) :: Boundary(:)
    !!     Vector of vertices which make up the perimeter of the aquifer
    !!   integer :: iNBdy
    !!     Number of perimeter points currently in use
    !! 
    complex (kind=ModAEM_Real) :: cRefPoint
    real (kind=ModAEM_Real) :: rRefHead
    complex (kind=ModAEM_Real) :: cRefUniformFlow
    logical (kind=ModAEM_Integer) :: lReference
    real (kind=ModAEM_Real) :: rSolConst
    complex (kind=ModAEM_Real),dimension(:),pointer :: cPerimeter
    integer (kind=ModAEM_Integer) :: iNPerim
    type ( AQU_BDYELEMENT ),dimension(:),pointer :: BdyElements
    integer (kind=ModAEM_Integer) :: iNBdy
    type ( IN0_COLLECTION ),pointer :: in0
    integer (kind=ModAEM_Integer) :: iRegenerate
  end type AQU_COLLECTION

  ! Locally-required constants
  ! Matrix generator flag for reference point equation
  integer (kind=ModAEM_Integer),private,parameter :: kAQUReference = 1
  integer (kind=ModAEM_Integer),private,parameter :: kAQUPerimeter = 2
  ! Fake a small well radius to ensure that the proper fluxes are computed
  ! in FWL_Flow
  real (kind=ModAEM_Real),private,parameter :: BDY_WELL_RADIUS = 1.0e-6_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: BDY_MOVE_FACTOR = 1.0e-4_ModAEM_Real

contains

!! ELEMENT MODULE ROUTINES
!! These routines allow AQU to behave as a ModAEM Element Module

function AQU_Create(iNInho,rBase,rThick,rHydCond,rPorosity,io) result (aqu)
  !! function AQU_Create
  !!
  !! Creates a new AQU_COLLECTION object for an infinite aquifer
  !!
  !! Calling Sequence:
  !!    aqu => AQU_Create(rBase,rThick,rHydCond,rPorosity)
  !!
  !! Arguments:
  !!    (in)    integer :: iNInho
  !!              Maximum number of inhomogeneities allowed
  !!    (in)    real :: rBase
  !!              Base elevation
  !!    (in)    real :: rThick
  !!              Thickness (large value for unconfined flow)
  !!    (in)    real :: rHydCond
  !!              Hydraulic conductivity
  !!    (in)    real :: rPorosity
  !!              Effective porosity
  !!
  ! [ ARGUMENTS ]
  integer (kind=ModAEM_Integer),intent(in) :: iNInho
  real (kind=ModAEM_Real),intent(in) :: rBase
  real (kind=ModAEM_Real),intent(in) :: rThick
  real (kind=ModAEM_Real),intent(in) :: rHydCond
  real (kind=ModAEM_Real),intent(in) :: rPorosity
  type ( IO_Status ),pointer :: io  

! [ RETURN VALUE ]
  type ( AQU_COLLECTION ),pointer :: aqu
  ! [ LOCALS ]
  character (len=132) :: sMessage
  integer (kind=ModAEM_Integer) :: iStat

  allocate (aqu,stat=iStat)
  if (IO_Assert( (iStat == 0), "AQU_Create: allocation failed",io )) return

  aqu%in0 => IN0_Create(iNInho,rBase,rThick,rHydCond,rPorosity,io)
  if ( io%lError ) return

  aqu%cRefPoint = cZERO
  aqu%rRefHead = rZERO
  aqu%cRefUniformFlow = cZERO
  aqu%lReference = .false.
  nullify(aqu%BdyElements)
  aqu%iNBdy = 0
  nullify(aqu%cPerimeter)
  aqu%iNPerim = 0
  ! The following are determined in AQU_Setup
  aqu%rSolConst = rZERO
  aqu%iRegenerate = 1

  return
end function AQU_Create

subroutine AQU_SetReference(aqu,cRefPoint,rRefHead,cRefUniformFlow,io)
  !! subroutine AQU_SetReference
  !!
  !! Sets up the reference point and uniform flow
  !!
  !! Calling Sequence:
  !!    call AQU_SetReference(aqu,cRefPoint,rRefHead,cRefUniformFlow)
  !!
  !! Arguments:
  !!    (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!              The AQU_COLLECTION object to be populated
  !!    (in)    complex :: cRefPoint
  !!              Location of the reference point
  !!    (in)    real :: rRefHead
  !!              The head at the reference point
  !!    (in)    complex :: cRefUniformFlow
  !!              Far-field uniform flow
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cRefPoint
  real (kind=ModAEM_Real),intent(in) :: rRefHead
  complex (kind=ModAEM_Real),intent(in) :: cRefUniformFlow
  type ( IO_Status ),pointer :: io
 
  if ( io%lDebug ) then
    if (IO_Assert( associated(aqu), "AQU_SetReference: No AQU_COLLECTION object",io )) return
  endif

  aqu%cRefPoint = cRefpoint
  aqu%rRefHead = rRefHead
  aqu%cRefUniformFlow = cRefUniformFlow
  aqu%lReference = .true.

  return
end subroutine AQU_SetReference

function iAQU_GetInfo(aqu,iOption,io) result(iValue)
  !! function AQU_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iAQU_GetInfo(aqu,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue
  integer (kind=ModAEM_Integer) :: iStr

  if ( io%lDebug ) then
    if (IO_Assert( (associated(aqu)), &
                    "AQU_GetInfo: AQU_Create has not been called",io )) return
  endif

  iValue = iIN0_GetInfo(aqu%in0,iOption,io)
  if ( io%lError ) return
  select case (iOption)
    case ( SIZE_FWL )
      iValue = iValue + aqu%iNBdy
    case ( SIZE_FDP )
      iValue = iValue + aqu%iNBdy
    case ( SIZE_EQUATIONS )
      iValue = iValue + aqu%iNBdy + 1
    case ( SIZE_UNKNOWNS )
      iValue = iValue + aqu%iNBdy + 1
    case default
      iValue = 0
  end select

  return
end function iAQU_GetInfo

subroutine AQU_SetupFunctions(aqu,fwl,fdp,io)
  !! subroutine AQU_SetupFunction
  !!
  !! This routine sets up the functions in f_well and f_dipole for the perimeter
  !! and sets up the equation for the reference point.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call AQU_SetupFunction(iL)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             FWL_COLLECTION object where functions may be added
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fwl
  !!             FDP_COLLECTION object where functions may be added
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cCPZ1,cCPZ2
  integer (kind=ModAEM_Integer) :: iBdy
  integer (kind=ModAEM_Integer) :: iNWL,iNPD,iNDP,iNEQ,iNUN
  type ( AQU_BDYELEMENT ),pointer :: this_bdy

  if (IO_Assert( (associated(aqu)), "AQU_Setup: AQU_Create has not been called",io )) return

  ! Set the confined potential
  call IN0_SetupFunctions(aqu%in0,fdp,io)
  if ( io%lError ) return

  ! Create the entries in the FWL and FDP modules for the perimeter
  ! Build dipoles for all perimeter segments
  do iBdy = 1,aqu%iNBdy
    this_bdy => aqu%BdyElements(iBdy)
    this_bdy%rLength = abs(this_bdy%cZ2 - this_bdy%cZ1)
    call FDP_New(fdp,this_bdy%cZ1,this_bdy%cZ2,(/cZERO,cZERO,cZERO/), &
                 this_bdy%iFDPIndex,io)
    if ( io%lError ) return
    call FWL_New(fwl,this_bdy%cZ2,rZERO,BDY_WELL_RADIUS,this_bdy%iFWLIndex,io)
    if ( io%lError ) return
  end do

  return
end subroutine AQU_SetupFunctions

subroutine AQU_SetupMatrix(aqu,mat,io)
  !! subroutine AQU_SetupMatrix
  !!
  !! This routine sets up the matrix entries for the module
  !! and sets up the equation for the reference point.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call AQU_SetupMatrix(aqu,mat,io)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer :: mat
  !!             MAT_MATRIX object for the solver
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cCPZ1,cCPZ2
  integer (kind=ModAEM_Integer) :: iBdy
  integer (kind=ModAEM_Integer) :: iNWL,iNPD,iNDP,iNEQ,iNUN
  type ( AQU_BDYELEMENT ),pointer :: this_bdy

  if (IO_Assert( (associated(aqu)), "AQU_Setup: AQU_Create has not been called",io )) return

  ! Set the confined potential
  call IN0_SetupMatrix(aqu%in0,mat,io)
  if ( io%lError ) return

  ! Build the matrix entries for the perimeter (if any)
  do iBdy=1,aqu%iNBdy
    this_bdy => aqu%BdyElements(iBdy)
    call MAT_CreateVariable(mat,ELEM_AQU,kAQUPerimeter,iBdy,0,io)
    if ( io%lError ) return
 
    if ( this_bdy%lFlag ) then
      call MAT_CreateEquation(mat,(/ rHALF*(this_bdy%cZ1+this_bdy%cZ2) /), &
                              EQN_HEAD,ELEM_AQU,kAQUPerimeter,iBdy,0, &
                              rIN0_HeadToPotential(aqu%in0,this_bdy%rSpecHead, &
                                              rHALF*(this_bdy%cZ1+this_bdy%cZ2),io),cZero,io)
      if ( io%lError ) return
    else
      cCPZ1 = this_bdy%cZ1 + (this_bdy%cZ2-this_bdy%cZ1)*BDY_MOVE_FACTOR
      cCPZ2 = this_bdy%cZ2 - (this_bdy%cZ2-this_bdy%cZ1)*BDY_MOVE_FACTOR
      call MAT_CreateEquation(mat,(/ cCPZ1,cCPZ2 /), &
                              EQN_FLOW,ELEM_AQU,kAQUPerimeter,iBdy,0,this_bdy%rSpecFlux,cZero,io)
      if ( io%lError ) return
    endif
  end do

  ! Build a matrix entry for the reference point (if any)
  call MAT_CreateVariable(mat,ELEM_AQU,kAQUReference,0,0,io)
  if ( io%lError ) return
  if ( aqu%lReference ) then
    call MAT_CreateEquation(mat,(/aqu%cRefPoint/),EQN_HEAD,ELEM_AQU,kAQUReference, &
                            0,0,rIN0_HeadToPotential(aqu%in0,aqu%rRefHead,aqu%cRefPoint,io),cZero,io)
    if ( io%lError ) return
  else
    call MAT_CreateEquation(mat,(/ cZERO /),EQN_CONTINUITY,ELEM_AQU,kAQUReference, &
                            0,0,rZERO,cZero,io)
    if ( io%lError ) return
  endif

  return
end subroutine AQU_SetupMatrix

subroutine AQU_Prepare(aqu,io) 
  !! subroutine AQU_Prepare
  !! 
  !! Prepares the module for a new iteration
  !!
  !! Do-nothing for m_aqu
  !!
  !! Calling Sequence:
  !!    call AQU_Prepare(aqu,io)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( MAT_MATRIX ),pointer
  !!             MAT_MATRIX object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( IO_STATUS ),pointer :: io

  return
end subroutine AQU_Prepare

subroutine AQU_ComputeCoefficients(aqu,fwl,fdp,cPathZ,iEqType,iElementType,iElementString, &
                                   iElementVertex,iElementFlag,cOrientation,rARow,io)
  !! subroutine HB0_ComputeCoefficients
  !!
  !! Computes a row of matrix coefficients (with no corrections) for the HB0 
  !! elements in layer iL.
  !!
  !! Calling Sequence:
  !!    call AQU_ComputeCoefficients(aqu,fwl,fdp,mat,cPathZ,iEqType,cOrientation,rRow)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             FWL_COLLECTION object where functions may be added
  !!   (in)    type ( FDP_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object where functions may be added
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
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  integer (kind=ModAEM_Integer),intent(in) :: iEqType
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  real (kind=ModAEM_Real),dimension(:),intent(out) :: rARow
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iIN0NWL,iIN0NPD,iIN0NDP,iIN0NEQ,iIN0NUN
  integer (kind=ModAEM_Integer) :: iStat,iCol,iVtx,iDP1,iNDP,iWhich,iRV
  complex (kind=ModAEM_Real),dimension(:,:,:),allocatable :: cDPF

  if ( io%lDebug ) then
    if (IO_Assert( (associated(aqu)), &
                   "AQU_ComputeCoefficients: No AQU_COLLECTION object",io ) ) return
  endif
  ! Set it up
  rARow = rZERO

  ! Build entries for the inhomogeneities
  iIN0NUN = iIN0_GetInfo(aqu%in0,SIZE_UNKNOWNS,io)
  if ( io%lError ) return
  call IN0_ComputeCoefficients(aqu%in0,fdp,cPathZ,iEqType,iElementType,iElementString, &
                               iElementVertex,iElementFlag,cOrientation,rARow(1:iIN0NUN),io)
  if ( io%lError ) return
  iCol = iIN0NUN + 1

  ! Build entries for the boundary elements
  ! ASSUMES that AQU_Setup routine created consecutive dipole entries 
  if ( associated(aqu%BdyElements) ) then
    iDP1 = aqu%BdyElements(1)%iFDPIndex
    iNDP = aqu%iNBdy
    allocate (cDPF(0:iNDP+1,1,1),stat=iStat)
    if (IO_Assert( (iStat==0), "AQU_ComputeCoefficients: Allocation failed",io )) return

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
      rARow(iCol) = real( cDPF(iVtx,1,1) )
      iCol = iCol+1
    end do
    deallocate(cDPF)
  endif

  ! Set up the closure condition
  select case (iEqType)
    case (EQN_HEAD)
      rARow(iCol) = rONE
    case (EQN_FLOW)
      rARow(iCol) = rZERO
    case (EQN_INHO)
      rARow(iCol) = rONE
    case (EQN_DISCHARGE)
      rARow(iCol) = rZERO
    case (EQN_RECHARGE)
      rARow(iCol) = rZERO
    case (EQN_CONTINUITY)
      rARow(iCol) = rZERO
  end select
  iCol = iCol+1

  return
end subroutine AQU_ComputeCoefficients

function rAQU_ComputeRHS(aqu,iElementType,iElementString,iElementVertex, &
                         iElementFlag,rSpecValue,rCheck,io) result(rRHS)
  !! function rAQU_ComputeRHS
  !!
  !! Computes the right-hand side value for the solution
  !!
  !! Calling Sequence:
  !!   rRHS = rAQU_ComputeRHS(aqu,rValue,iElementType,iElementString,iElementVertex, &
  !!                          iElementFlag,rSpecValue,rCheck)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
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
  type ( AQU_COLLECTION ),pointer :: aqu
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  real (kind=ModAEM_Real), intent(in) :: rSpecValue
  real (kind=ModAEM_Real), intent(in) :: rCheck
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rRHS

  select case ( iElementType )
    case ( ELEM_AQU )
      ! For AQU, this is easy...
      rRHS = rSpecValue - rCheck
    case ( ELEM_IN0 )
      rRHS = rIN0_ComputeRHS(aqu%in0,iElementType,iElementString,iElementVertex,iElementFlag,rSpecValue,rCheck)
  end select

  return
end function rAQU_ComputeRHS

subroutine AQU_StoreResult(aqu,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
  !! subroutine LS1_StoreResult
  !!
  !! Stores the results of a solution for a single equation associated with
  !! the LS1 module.
  !!
  !! Calling Sequence:
  !!    LS1_StoreResult(aqu,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    real :: rValue
  !!             The new result value from the solution vector
  !!   (in)    integer :: iElementType
  !!             Element type (either ELAM_AQU or ELEM_IN0)
  !!   (in)    integer :: iElementString
  !!             Element string number
  !!   (in)    integer :: iElementVertex
  !!             Element vertex number
  !!   (in)    integer :: iElementFlag
  !!             Element flag (e.g. for vertices which yield more than one equation)
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  real (kind=ModAEM_Real), intent(in) :: rValue
  integer (kind=ModAEM_Integer), intent(in) :: iElementType
  integer (kind=ModAEM_Integer), intent(in) :: iElementString
  integer (kind=ModAEM_Integer), intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer), intent(in) :: iElementFlag
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  type ( AQU_BDYELEMENT ),pointer :: this_bdy

  if ( io%lDebug ) then
    if (IO_Assert( (associated(aqu)), "AQU_Update: AQU_Create has not been called",io )) return
  endif

  select case ( iElementType )
    case ( ELEM_AQU )
      select case ( iElementString )
        case ( kAQUReference )
          aqu%rSolConst = aqu%rSolConst + rValue
        case ( kAQUPerimeter )
          this_bdy => aqu%BdyElements(iElementVertex)
          this_bdy%rStrength = this_bdy%rStrength + rValue 
      end select
    case ( ELEM_IN0 )
      call IN0_StoreResult(aqu%in0,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
   if ( io%lError ) return
  end select

  return
end subroutine AQU_StoreResult

subroutine AQU_Update(aqu,fwl,fdp,io)
  !! subroutine AQU_Update
  !!
  !! Updates the underlying function objects for the specified layer.
  !!
  !! Calling Sequence:
  !!    AQU_Update(aqu,fwl,fdp)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             FWL_COLLECTION object to be used
  !!   (in)    type ( fdp_COLLECTION ),pointer :: fdp
  !!             FDP_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iBdy,iRV
  complex (kind=ModAEM_Real) :: cRho1,cRho2,cRho3
  type ( AQU_BDYELEMENT ),pointer :: this_bdy

  if ( io%lDebug ) then
    if (IO_Assert( (associated(aqu)), "AQU_Update: AQU_Create has not been called",io )) return
  endif

  ! Update the boundary elements
  do iBdy = 1,aqu%iNBdy
    this_bdy => aqu%BdyElements(iBdy)
    cRho1 = cZERO
    cRho2 = rHALF * this_bdy%rLength * this_bdy%rStrength
    cRho3 = this_bdy%rLength * this_bdy%rStrength
    call FDP_Update(fdp,this_bdy%iFDPIndex,(/cRho1,cRho2,cRho3/),io)
    if ( io%lError ) return
    call FWL_Update(fwl,this_bdy%iFWLIndex,real(cRho3),io)
    if ( io%lError ) return
  end do

  ! Update the inhomogeneities
  call IN0_Update(aqu%in0,fdp,io)
  if ( io%lError ) return

  ! Never force a regeneration
  aqu%iRegenerate = 0
 
  return
end subroutine AQU_Update

function lAQU_CheckPoint(aqu,cZ,rTol,io) result(lSing)
  !! logical function lAQU_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lAQU_CheckPoint(aqu,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             The AQU_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (AQU_COLLECTION),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iBdy
  type ( AQU_BDYELEMENT ),pointer ::bdy

  ! NOTE: Singularities are found at the end points of the BDY elements 
  lSing = lIN0_CheckPoint(aqu%in0,cZ,rTol)
  do iBdy = 1,aqu%iNBdy
    bdy => aqu%BdyElements(iBdy)
    if ( abs(bdy%cZ1-cZ) < rTol .or. abs(bdy%cZ2-cZ)<rTol ) then
      lSing = .true.
      exit
    endif
  end do
  
  return
end function lAQU_CheckPoint

subroutine AQU_Read(aqu,io)
  !! subroutine AQU_Read
  !!
  !! Reads the aquifer information for layer iL from the input LU
  !!
  !! Calling Sequence:
  !!    call AQU_Read(aqu)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             Layer number to be read
  !!   (out)   integer :: iRV
  !!             Error code or errOK for success
  !!
  !! The format of the aqu section of the input file appears as follows:
  !!
  !! AQU ninho base thickness hyd-cond porosity
  !! REF (z0) refhead (Qx0,Qy0)    [ optional ]
  !! BDY nvertices                 [ optional ] 
  !!     (x1,y1) (x2,y2) head flux flag
  !!     (x1,y1) (x2,y2) head flux flag
  !!     ... up to nvertices
  !! END
  !! IN0 [optional]
  !! DOM nvertices base thickness hyd-cond porosity
  !!     (x,y) 
  !!     (x,y)
  !!     ... up to nvertices
  !! DOM ...
  !! ... up to ninho
  !! END
  !!
  !! NOTE: It is assumed that the AQU line was found by the caller
  !!
  !! The BDY section (perimeter boundary) requests the specified flux or head 
  !! values, based on the 'flag' value. If flag=0, use the flux boundary; if
  !! flag=1, uses the head boundary. Assumes that for vertex 'i', the flux
  !! or head specified is that of segment 'i'; the boundary MUST be closed.
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( IO_Status ),pointer :: io
   ! Locals -- for Directive parsing
  type (DIRECTIVE),dimension(5),parameter :: dirDirectives = &
                                             (/ dirEND,dirREF,dirPRM,dirBDY,dirIN0 /)
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
  complex (kind=ModAEM_Real) :: cRefPoint
  complex (kind=ModAEM_Real) :: cRefUniformFlow
  complex (kind=ModAEM_Real) :: cZ,cZ1,cZ2
  real (kind=ModAEM_Real) :: rSpecHead
  real (kind=ModAEM_Real) :: rSpecFlux
  logical (kind=ModAEM_Integer) :: lFlag
  integer (kind=ModAEM_Integer) :: iNBdy
  type ( AQU_BDYELEMENT ),pointer :: vtx
  ! State variables for the parser
  integer (kind=ModAEM_Integer) :: iParseMode
  integer (kind=ModAEM_Integer),parameter :: PARSE_NONE = 0
  integer (kind=ModAEM_Integer),parameter :: PARSE_BDY = 1
  integer (kind=ModAEM_Integer),parameter :: PARSE_PERIM = 2
  type ( IN0_DOMAIN ),pointer :: dom

  call IO_MessageText("  Reading AQU module input",io)
  if ( io%lError ) return
 
  if (IO_Assert( (associated(aqu)), "AQU_Read: AQU_Create has not been called",io )) return

  iParseMode = PARSE_NONE
  do 
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "AQU_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all Mod "ifXXXRead" routines. 
        if (IO_Assert( .false., "AQU_Read: Unexpected EOF",io )) return
      case (kOpData)
        ! A data line was found. If we have a specified perimeter, add the point
        ! to the perimeter.
        select case (iParseMode)
          case ( PARSE_PERIM )
            if (IO_Assert( (aqu%iNPerim < size(aqu%cPerimeter)), &
                            "AQU_Read: Space exhausted",io )) return
            read ( unit=sRecord, &
                   fmt=*, &
                   iostat=IStat &
                 ) cZ
            if (IO_Assert( (iStat==0), "AQU_Read: I/O Error",io )) return
            aqu%iNPerim = aqu%iNPerim+1
            aqu%cPerimeter(aqu%iNPerim) = cZ
          case ( PARSE_BDY )
            if (IO_Assert( (aqu%iNBdy < size(aqu%BdyElements)), &
                            "AQU_Read: Space exhausted",io )) return
            read ( unit=sRecord, &
                   fmt=*, &
                   iostat=iStat &
                 ) cZ1,cZ2,rSpecHead,rSpecFlux,lFlag
            if (IO_Assert( (iStat==0), "AQU_Read: I/O Error",io )) return
            aqu%iNBdy = aqu%iNBdy+1
            vtx => aqu%BdyElements(aqu%iNBdy)
            vtx%cZ1 = cZ1
            vtx%cZ2 = cZ2
            vtx%rSpecHead = rSpecHead
            vtx%rSpecFlux = rSpecFlux
            vtx%lFlag = lFlag
        end select
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        exit
      case (kOpREF)
        ! Reference record was found. Read the reference point and uniform flow rate
        ! then check for validity.
        read ( unit=sRecord, &
               fmt=*, &
               iostat=iStat &
             ) cRefPoint,rRefHead,cRefUniformFlow
        if (IO_Assert( (iStat==0), "AQU_Read: I/O Error",io )) return
        dom => IN0_FindDomain(aqu%in0,cRefPoint,io)
        if ( io%lError ) return
        if (IO_Assert( ( rRefHead>dom%rBase), "AQU_Read: Base elevation below aquifer top",io )) return
        aqu%cRefPoint = cRefPoint
        aqu%rRefHead = rRefHead
        aqu%cRefUniformFlow = cRefUniformFlow
        aqu%lReference = .true.
      case (kOpPRM)
        ! Start the perimeter
        if (IO_Assert( (.not. associated(aqu%cPerimeter)), &
                        "AQU_Read: Perimteter has already been specified",io )) return
        read ( unit=sRecord, &
               fmt=*, &
               iostat=iStat &
             ) iMax
        if (IO_Assert( (iStat==0), "AQU_Read: I/O Error",io )) return
        if (IO_Assert( (iMax>2), "AQU_Read: Illegal dimension",io )) return
        allocate (aqu%cPerimeter(iMax),stat=iStat)
        if (IO_Assert( (iStat==0), "AQU_Read: Allocation failed",io )) return
        aqu%cPerimeter = cZERO
        iParseMode = PARSE_PERIM
     case (kOpBDY)
        ! Have we already specified boundary elements?
        if (IO_Assert( (.not. associated(aqu%BdyElements)), &
                        "AQU_Read: BDY statement has already been encountered",io )) return
        read ( unit=sRecord, &
               fmt=*, &
               iostat=iStat &
             ) iMax
        if (IO_Assert( (iStat==0), "AQU_Read: I/O Error",io )) return
        if (IO_Assert( (iMax>0), "AQU_Read: Illegal dimension",io )) return
        allocate (aqu%BdyElements(iMax),stat=iStat)
        if (IO_Assert( (iStat==0), "AQU_Read: Allocation failed",io )) return
        aqu%BdyElements = AQU_BDYELEMENT(cZERO,cZERO,rZERO,rZERO,rZERO,rZERO,0,0,.false.)
        aqu%iNBdy = 0
        iParseMode = PARSE_BDY
      case (kOpIN0)
        call IN0_Read(aqu%in0,io)
        if ( io%lError ) return
      case default
    end select
  end do
  call IO_MessageText("  Leaving AQU module",io)

  return
end subroutine AQU_Read  

subroutine AQU_Report(aqu,io)
  ! Reports aquifer information for layer iL
  type ( AQU_COLLECTION ),pointer :: aqu
  type ( IO_Status),pointer :: io
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN,iBdy
  type ( AQU_BDYELEMENT ),pointer :: bdy

!$$ assert : allocated(Layers) : "No layers are allocated"
!$$ assert : (iL>=1 .and. iL<=size(Layers)) : "Wrong layer was specified"

  ! Write out the inhomogeneities for property definitions
  call IN0_Report(aqu%in0,io)
  if ( io%lError ) return
 
  ! Aquifer data from this module
  write ( unit=kIOOutputLU, fmt="('Report for module AQU'//'Aquifer information')" )

  write (unit=kIOOutputLU,fmt="(/""Reference Aquifer Flow Conditions""//)")
  write (unit=kIOOutputLU,fmt="(""    Reference Point:        "",2(e13.6,2x))") aqu%cRefPoint
  write (unit=kIOOutputLU,fmt="(""    Reference Head:         "",e13.6)") aqu%rRefHead
  write (unit=kIOOutputLU,fmt="(""    Uniform Flow Rate:      "",2(e13.6,2x))") aqu%cRefUniformFlow 

  write (unit=kIOOutputLU,fmt="(/""Computed Values""//)")
  write (unit=kIOOutputLU,fmt="(""    Solution Constant:      "",e13.6)") aqu%rSolConst

  ! Report the boundary elements
  write ( unit=kIOOutputLU, fmt="(/'Boundary elements'/)" )
  write ( unit=kIOOutputLU, &
          fmt="('Element#',t11,'Z1',t39,'Z2',t65,'Strength',t78,'FDPIndex',t88,'IWLIndex',t98,'SpecHead'," // &
              "t111,'SpecFlux',t124,'Flag')" &
        )
  do iBdy = 1,aqu%iNBdy
    bdy => aqu%BdyElements(iBdy)
    write ( unit=kIOOutputLU, &
            fmt="(i8,2x,5(g12.5,1x),2(i8,2x),2(g12.5,1x),l1)" &
          ) iBdy,bdy%cZ1,bdy%cZ2,bdy%rStrength,bdy%iFDPIndex,bdy%iFWLIndex,bdy%rSpecHead,bdy%rSpecFlux,bdy%lFlag
  end do

  ! Report the matrix generation values...
  write ( unit=kIOOutputLU, &
          fmt="(/' Function and Matrix values:'/)" &
        )
  write ( unit=kIOOutputLU, &
          fmt="('Number of FWL functions:',t30,i10)" &
        ) iAQU_GetInfo(aqu,SIZE_FWL,io)
  if ( io%lError ) return
  write ( unit=kIOOutputLU, &
          fmt="('Number of FPD functions:',t30,i10)" &
        ) iAQU_GetInfo(aqu,SIZE_FPD,io)
  if ( io%lError ) return
  write ( unit=kIOOutputLU, &
          fmt="('Number of FDP functions:',t30,i10)" &
        ) iAQU_GetInfo(aqu,SIZE_FDP,io)
  if ( io%lError ) return
  write ( unit=kIOOutputLU, &
          fmt="('Number of Equations:',t30,i10)" &
        ) iAQU_GetInfo(aqu,SIZE_EQUATIONS,io)
  if ( io%lError ) return
  write ( unit=kIOOutputLU, &
          fmt="('Number of Unknowns:',t30,i10)" &
        ) iAQU_GetInfo(aqu,SIZE_UNKNOWNS,io)
  if ( io%lError ) return
 
  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine AQU_Report

!! FUNCTION MODULE ROUTINES
!! These routines allow AQU to behave as a ModAEM Function Module

function cAQU_Potential(aqu,cZ,io) result(cOmega)
  !! complex function cAQU_Potential
  !!
  !! Computes the complex potential due to the reference aquifer at cZ
  !!
  !! Calling Sequence:
  !!    cOmega = cFDP_Potential(aqu,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!
  !! Return Value:
  !!           complex :: cOmega
  !!             The complex potential at the point cZ 
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cOmega
  ! [ LOCALS ]

  if (IO_Assert( (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object",io )) return

  if ( aqu%lReference ) then
    cOmega = aqu%rSolConst - conjg(aqu%cRefUniformFlow) * cZ 
  else
    cOmega = aqu%rSolConst
  endif

  return
end function cAQU_Potential

function cAQU_Discharge(aqu,cZ,io) result(cQ)
  !! complex function cAQU_Discharge
  !!
  !! Computes the complex discharge due to the reference aquifer at cZ
  !!
  !! Calling Sequence:
  !!    cOmega = cFDP_Discharge(aqu,cZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the discharge
  !!
  !! Return Value:
  !!           complex :: cW
  !!             The complex discharge at the point cZ in layer iL
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cQ
  ! [ LOCALS ]

  if (IO_Assert( (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object",io )) return

  if ( aqu%lReference ) then
    cQ = aqu%cRefUniformFlow
  else
    cQ = cZERO
  endif

  return
end function cAQU_Discharge

function rAQU_Flow(aqu,cPathZ,io) result(rFlow)
  !! real function cAQU_Flow
  !!
  !! Computes the integrated flow due to the current set of dipoles. 
  !!
  !! Calling Sequence:
  !!    rFlow = fAQU_Flow(aqu,cPathZ,iDP1,iNDP)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    complex :: cPathZ(:)
  !!             The path across which the flow is desired
  !!
  !! Return Value:
  !!           real :: rFlow
  !!             The integrated flux across the path cPathZ
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cPathZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rFlow
  ! [ LOCALS ]

  if (IO_Assert( (associated(aqu)), "rAQU_Flow: No AQU_COLLECTION object",io )) return

  ! The reference point (if any) and the uniform flow function satisfy
  ! Laplace equation, so we need to only look at the ends of the path
  if ( aqu%lReference ) then
    rFlow = aimag( cAQU_Potential(aqu,cPathZ(size(cPathZ)),io) - &
                   cAQU_Potential(aqu,cPathZ(1),io) )
    if ( io%lError ) return
  else
    rFlow = rZERO
  endif

  return
end function rAQU_Flow

function rAQU_Recharge(aqu,cZ,io) result(rG)
  !! complex function rAQU_Recharge
  !!
  !! Computes the recharge rate due to the aquifer
  !!
  !! Calling Sequence:
  !!    rG = rAQU_Recharge(aqu,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             The AQU_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!
  ! [ ARGUMENTS ]
  type (AQU_COLLECTION),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rG
  ! [ LOCALS ]
  
  rG = rZERO

  return
end function rAQU_Recharge

function rAQU_Extraction(aqu,io) result(rQ)
  !! complex function rAQU_Extraction
  !!
  !! Computes the extraction due to the aquifer
  !!
  !! Calling Sequence:
  !!    rQ = rAQU_Extraction(aqu)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             The AQU_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point at which to determine the potential
  !!
  ! [ ARGUMENTS ]
  type (AQU_COLLECTION),pointer :: aqu
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rQ
  ! [ LOCALS ]
  
  if (IO_Assert( (.not. aqu%lReference), "AQU_Extraction: Invalid for problems with uniform flow",io )) return

  rQ = rZERO

  return
end function rAQU_Extraction

!! UTILITY ROUTINES
!! These routines provide computational aids for AEM models
!! They make use of the IN0 equivalent calls

function rAQU_HeadToPotential(aqu,rHead,cZ,io) result(rPot)
  !! real function rAQU_HeadToPotential
  !!
  !! Converts head to a discharge potential based on the in0ifer properties 
  !! at cZ. Returns the (real) potential.
  !!
  !! Calling Sequence:
  !!    rPot = rAQU_HeadToPotential(aqu,rHead,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    real :: rHead
  !!             The head value to be converted
  !!   (in)    complex :: cZ
  !!             The point in the aquifer where the conversion is to take place
  !!
  !! Return Value:
  !!           real :: rPot
  !!             The discharge potential corresponding to the head 'rHead' 
  !! Note:
  !!   If iDP1 is not provided, all dipoles will be used
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  real (kind=ModAEM_Real),intent(in) :: rHead
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rPot

  rPot = rIN0_HeadToPotential(aqu%in0,rHead,cZ,io)
  if ( io%lError ) return
 
  return
end function rAQU_HeadToPotential

function rAQU_PotentialToHead(aqu,rPot,cZ,io) result(rHead)
  !! real function rAQU_HeadToPotential
  !!
  !! Converts head to a discharge potential based on the in0ifer properties 
  !! at cZ. Returns the (real) potential.
  !!
  !! Calling Sequence:
  !!    rPot = rAQU_HeadToPotential(in0,rHead,cZ)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    real :: rHead
  !!             The head value to be converted
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!
  !! Return Value:
  !!           real :: rHead
  !!             The head corresponding to the discharge potential 'rPot' 
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  real (kind=ModAEM_Real),intent(in) :: rPot
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rHead

  rHead = rIN0_PotentialToHead(aqu%in0,rPot,cZ,io)
  if ( io%lError ) return
 
  return
end function rAQU_PotentialToHead

function cAQU_DischargeToVelocity(aqu,cDischarge,cZ,rPot,io) result(cVelocity)
  !! function cAQU_DischargeToVelocity
  !!
  !! Converts Discharge to velocity, using the potential and the location cZ
  !! to compute saturated thickness
  !!
  !! Calling Sequence:
  !!    cV = cAQU_DischargeToVelocity(in0,rDischarge,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
  !!   (in)    real :: rDischarge
  !!             The discharge value to be converted
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
  type ( AQU_COLLECTION ),pointer :: aqu
  real (kind=ModAEM_Real),intent(in) :: rPot
  complex (kind=ModAEM_Real),intent(in) :: cDischarge
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cVelocity

  cVelocity = cIN0_DischargeToVelocity(aqu%in0,cDischarge,cZ,rPot,io)
  if ( io%lError ) return
 
  return
end function cAQU_DischargeToVelocity

function rAQU_SatdThickness(aqu,cZ,rPot,io) result (rH)
  !! function cAQU_SatdThickness
  !!
  !! Computes the saturated thickness at point cZ in layer iL where the potential is rPot
  !!
  !! Calling Sequence:
  !!    cV = cAQU_SatdThickness(in0,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQU_COLLECTION object to be used
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
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rH
  
  rH = rIN0_SatdThickness(aqu%in0,cZ,rPot,io)
  if ( io%lError ) return
 
  return
end function rAQU_SatdThickness

function rAQU_Transmissivity(aqu,cZ,rPot,io) result (rT)
  !! function rAQU_Transmissivity
  !!
  !! Computes the transmissivity at the point cZ where the potential is rPot
  !!
  !! Calling Sequence:
  !!    cV = cAQU_Transmissivity(aqu,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQu_COLLECTION object to be used
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
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real),intent(in) :: rPot
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rT

  rT = rIN0_Transmissivity(aqu%in0,cZ,rPot,io)
  if ( io%lError ) return
 
  return
end function rAQU_Transmissivity

function rAQU_HydCond(aqu,cZ,io) result (rHydCond)
  !! function rAQU_HydCond
  !!
  !! Computes the hydraulic conductivity at the point cZ
  !!
  !! Calling Sequence:
  !!    cV = cAQU_HydCond(aqu,cZ,rPot)
  !!
  !! Arguments:
  !!   (in)    type ( AQU_COLLECTION ),pointer :: aqu
  !!             AQu_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point in the in0ifer where the conversion is to take place
  !!
  !! Return Value:
  !!           complex :: rHydCond
  !!             The hydraulic conductivity
  !!
  ! [ ARGUMENTS ]
  type ( AQU_COLLECTION ),pointer :: aqu
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
   ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rHydCond

  rHydCond = rIN0_HydCond(aqu%in0,cZ,io)
  if ( io%lError ) return
 
  return
end function rAQU_HydCond

end module m_aqu

