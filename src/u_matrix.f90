module u_matrix

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

! This is a generic module which supports the generation and solution of
! a full MxM matrix.

use u_constants
use u_io

implicit none

public

  type,public :: MAT_EQUATION
    !! type MAT_EQUATION
    !!
    !! PUBLIC type that holds information for one equation in the matrix
    !!
    !! Members:
    !!   complex :: cCPZ(:)
    !!     The center of the well
    !! 
    complex (kind=ModAEM_Real),dimension(:),pointer :: cCPZ
    integer (kind=ModAEM_Integer) :: iEqnType
    integer (kind=ModAEM_Integer) :: iElementID
    integer (kind=ModAEM_Integer) :: iElementString
    integer (kind=ModAEM_Integer) :: iElementVertex
    integer (kind=ModAEM_Integer) :: iElementFlag
    complex (kind=ModAEM_Real) :: cOrientation
    real (kind=ModAEM_Real) :: rSpecValue
    real (kind=ModAEM_Real) :: rCheck
  end type MAT_EQUATION

  type,public :: MAT_VARIABLE
    !! type MAT_VARIABLE
    !!
    !! PUBLIC type that holds information for one unknown variable
    !!
    !! Members:
    !!   integer :: iElementID
    !!     The element ID (e.g. WL0, LS0, LS1, HB0)
    !! 
    integer (kind=ModAEM_Integer) :: iElementID
    integer (kind=ModAEM_Integer) :: iElementString
    integer (kind=ModAEM_Integer) :: iElementVertex
    integer (kind=ModAEM_Integer) :: iElementFlag
    real (kind=ModAEM_Real) :: rValue
  end type MAT_VARIABLE

  type,public :: MAT_MATRIX
    !! type MAT_MATRIX
    !!
    !! PUBLIC type that holds information for a matrix
    !!
    !! Members:
    !!   type ( MAT_EQUATION ), pointer :: Equations(:)
    !!     The equation definitions
    !! 
    type (MAT_EQUATION),dimension(:),pointer :: Equations
    type (MAT_VARIABLE),dimension(:),pointer :: Variables
    real (kind=ModAEM_Real) :: rCond
    real (kind=ModAEM_Real),dimension(:,:),pointer :: rMatrix
    real (kind=ModAEM_Real),dimension(:,:),pointer :: rRHS
    integer (kind=ModAEM_Integer),dimension(:),pointer :: iPivot
    integer (kind=ModAEM_Integer) :: iNEqn
    integer (kind=ModAEM_Integer) :: iNVar
  end type MAT_MATRIX

  ! Iterative refinement settings for gesolve()
  integer (kind=ModAEM_Integer),private,save :: iMATMaxIterations = 10
  real (kind=ModAEM_Real),private,save :: rMATIterTolX = 1.0e-6_ModAEM_Real
  real (kind=ModAEM_Real),private,save :: rMATIterTolR = 1.0e-6_ModAEM_Real
  real (kind=ModAEM_Real),private,save :: rMOVEPOINT = 1.0e-4

contains

function MAT_Create(io) result(mat)
  ! [ ARGUMENTS ]
  ! [ RETURN VALUE ]
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iL,iStat

  allocate ( mat,stat=iStat )
  if (IO_Assert( (iStat == 0),"MAT_Create: Allocation failed",io)) return
  
  mat%iNEqn = 0
  mat%iNVar = 0
    nullify ( mat%Equations, &
            mat%Variables, &
            mat%rMatrix, &
            mat%rRHS, &
            mat%iPivot &
          )

  return
end function MAT_Create

subroutine MAT_Alloc(mat,iNEQ,iNUN,io)
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iNEQ
  integer (kind=ModAEM_Integer),intent(in) :: iNUN
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  if ( io%lDebug ) then
    if (IO_Assert( (associated(mat)), &
                    "MAT_Alloc: MAT_Create has not been called",io )) return
    if (IO_Assert( (iNEQ > 0 .and. iNUN > 0), &
                    "MAT_Alloc: Illegal dimensions",io )) return
    if (IO_Assert( (iNEQ == iNUN), &
                    "MAT_Alloc: Overspecification not supported in this version",io)) return
  endif

  allocate ( mat%rMatrix(iNEQ,iNUN), &
             mat%rRHS(iNEQ,1), &
             mat%Equations(iNEQ), &
             mat%Variables(iNUN), & 
             mat%iPivot(iNEQ),stat=iStat &
           )
  if (IO_Assert( (iStat==0), "MAT_Alloc: Allocation failed",io )) return
  
  mat%rMatrix = rZERO
  mat%rRHS = rZERO
  mat%iNEqn = 0
  mat%iNVar = 0
  
  return
end subroutine MAT_Alloc

subroutine MAT_Destroy(mat,io)
  !! This subroutine de-allocates the space prevously allocated for the matrix
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  if ( associated(mat%Equations) ) deallocate (mat%Equations,stat=iStat)
  if ( associated(mat%Variables) ) deallocate (mat%Variables,stat=iStat)
  if ( associated(mat%rMatrix) ) deallocate (mat%rMatrix,stat=iStat)
  if ( associated(mat%rRHS) ) deallocate (mat%rRHS,stat=iStat)
  if ( associated(mat%iPivot) ) deallocate (mat%iPivot,stat=iStat)
  deallocate (mat,stat=iStat)
  if (IO_Assert( (iStat == 0),"MAT_Destroy: Deallocation failed",io)) return

  return
end subroutine MAT_Destroy
  
subroutine MAT_Clear(mat,io)
  !! This subroutine deallocates all the matrix data, leaving the matrix object intact
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  if ( associated(mat%Equations) ) deallocate (mat%Equations,stat=iStat)
  if ( associated(mat%Variables) ) deallocate (mat%Variables,stat=iStat)
  if ( associated(mat%rMatrix) ) deallocate (mat%rMatrix,stat=iStat)
  if ( associated(mat%rRHS) ) deallocate (mat%rRHS,stat=iStat)
  if ( associated(mat%iPivot) ) deallocate (mat%iPivot,stat=iStat)
  if (IO_Assert( (iStat == 0),"MAT_Destroy: Deallocation failed",io)) return

  return
end subroutine MAT_Clear
  
subroutine MAT_CreateVariable(mat,iElementID,iElementString,iElementVertex,iElementFlag,io)
  ! This function is called by the xxxSetup routines to allocate space for an unknown variable
  ! (e.g. an unknown strength coefficient).  This information is stored in the 'variables' buffer.  
  ! The function returns the index into the 'variables' buffer for the unknown value
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iElementID
  integer (kind=ModAEM_Integer),intent(in) :: iElementString
  integer (kind=ModAEM_Integer),intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer),intent(in) :: iElementFlag
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  type ( MAT_VARIABLE ),pointer :: vbl

  if ( io%lDebug ) then
    if (IO_Assert( (associated(mat)), &
                    "MAT_Alloc: MAT_Create has not been called",io )) return
    if (IO_Assert( (mat%iNVar <= size(mat%Variables)), &
                    "MAT_CreateVariable: Space exhausted",io )) return
  endif

  mat%iNVar = mat%iNVar+1
  vbl => mat%Variables(mat%iNVar)
  vbl%iElementID = iElementID
  vbl%iElementString = iElementString
  vbl%iElementVertex = iElementVertex
  vbl%iElementFlag = iElementFlag
  vbl%rValue = rZERO

  return
end subroutine MAT_CreateVariable

subroutine MAT_GetVariable(mat,iVar,rValue,iElementID,iElementString, &
                           iElementVertex,iElementFlag,io)
  ! This function extracts the information about the specified variable
  ! It returns kOK if no error was detected, else an error code
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iVar
  real (kind=ModAEM_Real),intent(out) :: rValue
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementString
  integer (kind=ModAEM_Integer),intent(out) :: iElementVertex
  integer (kind=ModAEM_Integer),intent(out) :: iElementFlag
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  character (len=255) :: sBuf
  type ( MAT_VARIABLE ),pointer :: var

  if ( io%lDebug ) then
    if (IO_Assert( (associated(mat)), &
                    "MAT_Alloc: MAT_Create has not been called",io )) return
    if (IO_Assert( (iVar >= 1 .and. iVar <= mat%iNVar ), &
                    "MAT_GetVariable: Bad variable number",io )) return
  endif

  var => mat%Variables(iVar)
  rValue = var%rValue
  iElementID = var%iElementID 
  iElementString = var%iElementString 
  iElementVertex = var%iElementVertex 
  iElementFlag = var%iElementFlag

  return
end subroutine MAT_GetVariable

subroutine MAT_CreateEquation(mat,cCPZ,iEqnType,iElementID,iElementString, &
                              iElementVertex,iElementFlag,rSpecValue,cOrientation,io)
  ! This function is called by the xxxSetup routines to create entries in the 
  ! matrix generator. It returns the number of the new equation.
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cCPZ
  integer (kind=ModAEM_Integer),intent(in) :: iEqnType 
  integer (kind=ModAEM_Integer),intent(in) :: iElementID
  integer (kind=ModAEM_Integer),intent(in) :: iElementString
  integer (kind=ModAEM_Integer),intent(in) :: iElementVertex
  integer (kind=ModAEM_Integer),intent(in) :: iElementFlag
  real (kind=ModAEM_Real),intent(in) :: rSpecValue
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: ieq,ipt
  real (kind=ModAEM_Real),parameter :: rCPTolerance = 1.0e-2_ModAEM_Real
  logical (kind=ModAEM_Integer) :: lFail
  type ( MAT_EQUATION ),pointer :: eqn

  if ( io%lDebug ) then
    if (IO_Assert( (associated(mat)), &
                    "MAT_CreateEquation: MAT_Create has not been called",io )) return
    if (IO_Assert( (mat%iNEqn < size(mat%Equations)), &
                    "MAT_CreateEquation: Space exhausted",io )) return
  endif
  
  select case ( iEqnType )
    case ( EQN_FLOW )
      if (IO_Assert( ( size(cCPZ) > 1 ), &
                      "MAT_CreateEquation: Path too short for EQN_FLOW",io )) return
    case ( EQN_POTENTIALDIFF )
      if (IO_Assert( ( size(cCPZ) == 2 ), &
                      "MAT_CreateEquation: Exactly two control points required",io )) return
    case default
      if (IO_Assert( ( size(cCPZ) == 1 ), &
                      "MAT_CreateEquation: Exactly one control point required",io )) return
  endselect

  ! Check for overlapping control points
  do ieq = 1,mat%iNEqn
    eqn => mat%Equations(ieq)
    if ( (eqn%iEqnType == iEqnType) .and. (size(cCPZ,1) == size(eqn%cCPZ,1)) ) then
      if (IO_Assert( ( maxval(abs(cCPZ-eqn%cCPZ)) > rCPTolerance ), &
                       "MAT_CreateEquation: Coincident control points",io )) return
    endif
  end do

  mat%iNEqn = mat%iNEqn+1
  eqn => mat%Equations(mat%iNEqn)
  allocate(eqn%cCPZ(ubound(cCPZ,1)))
  eqn%cCPZ = cCPZ
  eqn%iEqnType = iEqnType
  eqn%iElementID = iElementID
  eqn%iElementString = iElementString
  eqn%iElementVertex = iElementVertex
  eqn%iElementFlag = iElementFlag
  eqn%cOrientation = cOrientation
  eqn%rSpecValue = rSpecValue
  eqn%rCheck = rZERO

  return
end subroutine MAT_CreateEquation

subroutine MAT_GetEquation(mat,iEQ,cCPZ,iNCP,iEqnType,iElementID,iElementString, &
                           iElementVertex,iElementFlag,cOrientation,rSpecValue,rCheck,io)
  ! This function extracts the information about the specified equation
  ! It returns kOK if no error was detected, else an error code
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iEQ
  complex (kind=ModAEM_Real),dimension(:),intent(out) :: cCPZ
  integer (kind=ModAEM_Integer),intent(out) :: iNCP
  integer (kind=ModAEM_Integer),intent(out) :: iEqnType
  integer (kind=ModAEM_Integer),intent(out) :: iElementID
  integer (kind=ModAEM_Integer),intent(out) :: iElementString
  integer (kind=ModAEM_Integer),intent(out) :: iElementVertex
  integer (kind=ModAEM_Integer),intent(out) :: iElementFlag
  complex (kind=ModAEM_Real),intent(out) :: cOrientation
  real (kind=ModAEM_Real),intent(out) :: rSpecValue
  real (kind=ModAEM_Real),intent(out) :: rCheck
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  character (len=255) :: sBuf
  type ( MAT_EQUATION ),pointer :: eqn

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_GetEquation: MAT_Create has not been called",io )) return
    if (IO_Assert( ( iEQ <= mat%iNEqn ), &
                    "MAT_GetEquation: Bad equation number",io )) return
  endif

  eqn => mat%Equations(iEQ)
  iNCP = ubound(eqn%cCPZ,1)
  cCPZ(1:iNCP) = eqn%cCPZ
  iEqnType = eqn%iEqnType 
  iElementID = eqn%iElementID 
  iElementString = eqn%iElementString 
  iElementVertex = eqn%iElementVertex 
  iElementFlag = eqn%iElementFlag
  cOrientation = eqn%cOrientation
  rCheck = eqn%rCheck
  rSpecValue = eqn%rSpecValue

  return
end subroutine MAT_GetEquation

subroutine MAT_UpdateEquation(mat,iEQ,rCheck,io)
  ! Updates the Check field in the specified equation
  ! It returns kOK if no error was detected, else an error code
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iEQ
  real (kind=ModAEM_Real),intent(in) :: rCheck
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  character (len=255) :: sBuf
  type ( MAT_EQUATION ),pointer :: eqn

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_UpdateEquation: MAT_Create has not been called",io )) return
    if (IO_Assert( ( iEQ <= mat%iNEqn ), &
                    "MAT_UpdateEquation: Bad equation number",io )) return
  endif

  eqn => mat%Equations(iEQ)
  eqn%rCheck = rCheck

  return
end subroutine MAT_UpdateEquation

subroutine MAT_SetRow(mat,iEQ,rRow,rRHS,io)
  ! Stores the given row data and right-hand side value into the matrix
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iEQ
  real (kind=ModAEM_Real),dimension(:),intent(in) :: rRow
  real (kind=ModAEM_Real),intent(in) :: rRHS
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_SetRow: MAT_Create has not been called",io )) return
    if (IO_Assert( ( iEQ <= mat%iNEqn ), &
                    "MAT_SetRow: Bad equation number",io )) return
    if (IO_Assert( ( size(rRow) == mat%iNVar ), &
                    "MAT_SetRow: Wrong number of terms in matrix row",io)) return
  endif

  ! All's well! Store the row
  mat%rMatrix(iEQ,:) = rRow(1:size(mat%rMatrix,2))
  mat%rRHS(iEQ,1) = rRHS

  return
end subroutine MAT_SetRow

subroutine MAT_SetRHS(mat,iEQ,rRHS,io)
  ! Stores the given row data and right-hand side value into the matrix
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  integer (kind=ModAEM_Integer),intent(in) :: iEQ
  real (kind=ModAEM_Real),intent(in) :: rRHS
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_SetRHS: MAT_Create has not been called",io )) return
    if (IO_Assert( ( iEQ <= mat%iNEqn ), &
                    "MAT_SetRHS: Bad equation number",io )) return
  endif

  mat%rRHS(iEQ,1) = rRHS

  return
end subroutine MAT_SetRHS

subroutine MAT_Decompose(mat,io)
  ! Solves the matrix
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rInfNorm,rCond
  real (kind=ModAEM_Real),dimension(:),allocatable :: rSums,rWork
  integer (kind=ModAEM_Integer) :: iInfo,i
  integer (kind=ModAEM_Integer),dimension(:),allocatable :: iWork
  integer (kind=ModAEM_Integer) :: iRV

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_Decompose: MAT_Create has not been called",io )) return
    if (IO_Assert( ( mat%iNEqn > 0 .and. mat%iNVar > 0 ), &
                    "MAT_Decompose: No rows/columns exist",io)) return
    if (IO_Assert( ( mat%iNVar == mat%iNEqn ), &
                    "MAT_Decompose: Matrix is not square",io )) return
  endif

  ! Compute the 1-Norm of the matrix
  allocate (rSums(size(mat%rMatrix,1)), &
            rWork(4*size(mat%rMatrix,1)), &
            iWork(size(mat%rMatrix,1)) )
  rSums = (/ ( mat%rMatrix(i,:),i=1,size(mat%rMatrix,1) )/)
  rInfNorm = maxval(rSums)
  ! Decompose the matrix
  call DGETRF(size(mat%rMatrix,1), &
              size(mat%rMatrix,2), &
              mat%rMatrix, &
              size(mat%rMatrix,1), &
              mat%iPivot, &
              iRV)
  if (IO_Assert( ( iRV == 0 ), & 
                  "LAPACK solver DGETRF failed",io )) return

  ! Compute the matrix condition number for testing 
  call DGECON('I', &
              size(mat%rMatrix,1), &
              mat%rMatrix, &
              size(mat%rMatrix,2), &
              rInfNorm, &
              rCond, &
              rWork, &
              iWork, &
              iInfo)
  mat%rCond = rCond
  deallocate (rWork,iWork)

  return
end subroutine MAT_Decompose

subroutine MAT_Solve(mat,io)
  ! Solves the matrix
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  real (kind=ModAEM_Real),dimension(:),allocatable :: rX
  integer (kind=ModAEM_Integer) :: iRV
  type ( MAT_VARIABLE ),pointer :: vbl

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_SetRHS: MAT_Create has not been called",io )) return
    if (IO_Assert( ( mat%iNEqn > 0 .and. mat%iNVar > 0 ), &
                    "MAT_SetRHS: No rows/columns exist",io)) return
    if (IO_Assert( ( mat%iNVar == mat%iNEqn ), &
                    "MAT_SetRHS: Matrix is not square",io )) return
  endif

  allocate (rX(ubound(mat%rRHS,1)))
  call DGETRS('N', &
              size(mat%rMatrix,1), &
              1, &
              mat%rMatrix, &
              size(mat%rMatrix,1), &
              mat%iPivot, &
              mat%rRHS, &
              size(mat%rRHS,1), &
              iRV)
  if (IO_Assert( ( iRV == 0 ), & 
                  "LAPACK solver DGETRS failed",io )) return
  
  ! Store the x-vector into the variables
  do i=1,mat%iNVar
    vbl => mat%Variables(i)
    vbl%rValue = mat%rRHS(i,1)
  end do

  deallocate(rX)
  return
end subroutine MAT_Solve

subroutine MAT_ComputeControlPoints(cZ1,cZ2,iNCP,cCPResult,rNormalOffset,io)
  ! This subroutine computes a set of iNCP control points for the linear
  ! element which extends from cZ1 to cZ2.  Points are returned as follows:
  ! cCPResult(1) = cZ1, cCPResult(2:iNCP+1) are computed, cCPResult(iNCP+2) = cCZ2
  complex (kind=ModAEM_Real),intent(in) :: cZ1,cZ2
  integer (kind=ModAEM_Integer),intent(in) :: iNCP
  complex (kind=ModAEM_Real),dimension(:),intent(out) :: cCPResult
  real (kind=ModAEM_Real),intent(in) :: rNormalOffset 
  type ( IO_STATUS ),pointer :: io
  ! Locals
  complex (kind=ModAEM_Real) :: cZC,cZL,cOffset
  real (kind=ModAEM_Real) :: rDTheta
  integer (kind=ModAEM_Integer) :: i


  if ( io%lDebug ) then
    if (IO_Assert( ( abs(cZ2-cZ1) > rMOVEPOINT ), &
                    "MAT_ComputeControlPoints: Endpoints are coincident",io )) return
    if (IO_Assert( ( size(cCPResult) >= iNCP+2 ), &
                    "MAT_ComputeControlPoints: Insufficient space for results",io )) return
  endif

  ! The control points are computed by turning an angle about the center
  ! of the segment.
  cZC = rHALF * (cZ1+cZ2)
  cZL = rHALF * (cZ2-cZ1)
  rDTheta = rPI / (iNCP+1)
  cOffset = conjg(cZ2-cZ1)/abs(cZ2-cZ1) * rNormalOffset
    
  ! Here we go
  cCPResult(1) = cZ1 + rMOVEPOINT*cZL + cOffset
  do i=1,iNCP
    cCPResult(i+1) = cZC - cos(i*rDTheta) * cZL + cOffset 
  end do
  cCPResult(iNCP+2) = cZ2 - rMOVEPOINT*cZL + cOffset
  
  return
end subroutine MAT_ComputeControlPoints

subroutine MAT_Report(mat,label,io)
  ! This routine writes a report of all matrix information
  ! [ ARGUMENTS ]
  type ( MAT_MATRIX ),pointer :: mat
  character (len=*),intent(in) :: label
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i,j,iStat
  type ( MAT_VARIABLE ),pointer :: vbl
  type ( MAT_EQUATION ),pointer :: eqn

  if ( io%lDebug ) then
    if (IO_Assert( ( associated(mat) ), &
                    "MAT_SetRHS: MAT_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module MAT'//'Matrix solver information')" )

  if ( associated(mat%rMatrix) ) then
    write (unit=kIOOutputLU, &
           fmt="('Number of equations:   available: ',i5,'  used: ',i5)") &
           ubound(mat%rMatrix,1),mat%iNEqn
    write (unit=kIOOutputLU, &
           fmt="('Number of unknowns:    available: ',i5,'  used: ',i5)") &
           ubound(mat%rMatrix,2),mat%iNVar
    write (unit=kIOOutputLU,fmt="('Condition number: ',e12.5)" ) mat%rCond

    write (unit=kIOOutputLU,fmt="(/'Matrix equation generator data')")
    do i=1,mat%iNEqn
      eqn => mat%Equations(i)
      write (unit=kIOOutputLU, &
             fmt="(/'Equation',t11,'Type',t21,'Element',t31,'String',t41,'Vertex',t51,'# CPs')")
      write (unit=kIOOutputLU,fmt="(6(i7,3x))") i,eqn%iEqnType,eqn%iElementID, &
                                                eqn%iElementString,eqn%iElementVertex, &
                                                ubound(eqn%cCPZ,1)
      write (unit=kIOOutputLU,fmt="(t11,'Control Points')")
      do j=1,ubound(eqn%cCPZ,1)
        write(unit=kIOOutputLU,fmt="(t11,e14.6,t27,e15.6)") eqn%cCPZ(j)
      end do
    end do
    write (unit=kIOOutputLU,fmt="(/'Matrix unknown variable data'/)")
    write (unit=kIOOutputLU, &
           fmt="('Element',t11,'String',t21,'Vertex',t31,'Flag',t41,'Value')")
    do i=1,mat%iNVar
      vbl => mat%Variables(i)
      write (unit=kIOOutputLU,fmt="(4(i10,1x),e12.6)") vbl%iElementID,vbl%iElementString, &
                                                       vbl%iElementVertex,vbl%iElementFlag, &
                                                       vbl%rValue
    end do

    write (unit=kIOOutputLU,fmt="(/'Matrix values written to grid file mat_report.grd'/)")
    open (unit=kIOGridLU,file="mat_report."//label//".grd",iostat=iStat)
    if (IO_Assert( (iStat==0),"MAT_Report: Could not open grid file mat_report."//label//".grd",io )) return
    write (unit=kIOGridLU,fmt="('DSAA')")
    write (unit=kIOGridLU,fmt=*) 1,size(mat%rMatrix,1)
    write (unit=kIOGridLU,fmt=*) 1,size(mat%rMatrix,2)
    write (unit=kIOGridLU,fmt=*) minval(mat%rMatrix),maxval(mat%rMatrix)
    do i=1,mat%iNEqn
      write (unit=kIOGridLU,fmt=*) (/ (mat%rMatrix(i,j),j=1,size(mat%rMatrix,2)) /)
    end do
    close (unit=kIOGridLU)

  else
    write (unit=kIOOutputLU, &
           fmt="('Matrix is not allocated')")
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine MAT_Report

end module u_matrix
