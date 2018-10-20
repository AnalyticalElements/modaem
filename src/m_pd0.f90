module m_pd0

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

  !! module m_pd0
  !!
  !! Element module for 2-D discharge specified ponds
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_pond     --  Function module for collections of ponds
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! pond elements. Each element requires a location, radius and discharge.

use u_constants
use u_io
use f_pond

implicit none

public

  type :: PD0_POND
    !! type PD0_POND
    !!
    !! Type that holds information for one pond
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the pond
    !!   real :: rRadius
    !!     The radius of the pond
    !!   complex :: rGamma
    !!     Discharge of the pond; positive indicates extraction
    !!     and negative indicates injection.
    !!   integer :: iID
    !!     Identification label for the pond (for interaction with e.g. GUIs)
    !!   integer :: iFPDIndex
    !!     Index for the pond entry in the FPD module
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rGamma
    real (kind=ModAEM_Real) :: rRadius
    integer (kind=ModAEM_Integer) :: iID
    integer (kind=ModAEM_Integer) :: iFPDIndex
  end type PD0_POND

  type :: PD0_COLLECTION
    !! type PD0_COLLECTION
    !!
    !! Type that holds all ponds in a layer
    !!
    !! Members:  type ( IO_Status),pointer :: io


    !!   type ( PD0_POND ),dimension(:),pointer :: Ponds
    !!     Array of PD0_POND objects for the layer; dimensioned for the maximum
    !!     number of ponds according to the input file (see PD0_Read)
    !!   integer :: iCount
    !!     The actual number of ponds in use in the layer
    !! 
    type ( PD0_POND ),dimension(:),pointer :: Ponds
    integer (kind=ModAEM_Integer) :: iCount
  end type PD0_COLLECTION

contains

!**pd Modified to allocate only the collection object
function PD0_Create(io) result (pd0)
  !! function PD0_Create
  !!
  !! Creates a new PD0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    pd0 => PD0_Create()
  !!
  !! Arguments:
  ! [ ARGUMENTS ]
  type ( IO_Status),pointer :: io
  ! [ RETURN VALUE ]
  type ( PD0_COLLECTION ),pointer :: pd0
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (pd0,stat=iStat)
  if (IO_Assert( (iStat == 0), "PD0_Create: allocation failed",io )) return
  nullify (pd0%Ponds)
  pd0%iCount = 0

  return
end function PD0_Create

!**pd New PD0_Alloc subroutine
subroutine PD0_Alloc(pd0, iNPD,io)
  !! Subroutine PD0_Alloc
  !! 
  !! Allocates wells for the PD0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call PD0_Alloc(pd0, iNPD)
  !!
  !! Arguments:
  !!    (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!              The PD0_COLLECTION object to be used
  !!    (in)    integer :: iNPD
  !!              Maximum number of ponds
  !!
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  integer (kind=ModAEM_Integer),intent(in) :: iNPD
  type ( IO_Status),pointer :: io

  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (pd0%Ponds(iNPD),stat=iStat)
  if (IO_Assert( (iStat == 0), "PD0_Alloc: allocation failed",io )) return

end subroutine PD0_Alloc

!**pd  New PD0_Destroy subroutine
subroutine PD0_Destroy(pd0,io)
  !! subroutine PD0_Destroy
  !!
  !! Frees memory allocated for an PD0 Ponds and PD0 Collection object
  !!
  !! Calling Sequence:
  !!     call PD0_Destroy(pd0)
  !!
  !! Arguments:
  !!  type ( PD0_COLLECTION ),pointer :: pd0
  !!              Pointer to the pd0_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( IO_Status),pointer :: io

  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat


  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                   "PD0_Destroy: PD0_Create has not been called",io )) return
  endif

  if (associated(pd0%Ponds)) then
    deallocate (pd0%Ponds,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                     "PD0_Destroy: deallocation of Ponds failed",io )) return
  end if  
  deallocate (pd0,stat=iStat)
  if (IO_Assert( (iStat == 0), "PD0_Destroy: deallocation failed",io )) return
  
  return
end subroutine PD0_Destroy

subroutine PD0_New(pd0,Pond,io) 
  !! function PD0_New
  !!
  !! Adds a new PD0_POND object to the PD0_COLLECTION 'pd0'
  !!
  !! Calling Sequence:
  !!    call PD0_New(pd0,Pond)
  !!
  !! Arguments:
  !!    (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!              The PD0_COLLECTION object to be used
  !!    (in)    type ( PD0_POND ),pointer :: Pond
  !!              Vector that defines the points along the barrier
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( PD0_POND ) :: Pond
  type ( IO_Status),pointer :: io

  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_New: PD0_Create has not been called",io )) return
  endif

  if (IO_Assert( (pd0%iCount < size(pd0%Ponds)), &
                  "PD0_New: Space exhausted",io )) return

  pd0%iCount = pd0%iCount + 1
  pd0%Ponds(pd0%iCount) = Pond

  return
end subroutine PD0_New

function iPD0_GetInfo(pd0,iOption,io) result(iValue)
  !! function PD0_GetInfo
  !!
  !! Returns the following sizing requirements for the PD0module
  !!
  !! Calling Sequence:
  !!    iValue = iPD0_GetInfo(pd0,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_Status),pointer :: io

  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue

  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_GetInfo: PD0_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FPD )
      iValue = pd0%iCount
    case default
      iValue = 0
  end select

  return
end function iPD0_GetInfo

subroutine PD0_SetupFunctions(pd0,fpd,io)
  !! subroutine PD0_Setup
  !!
  !! This routine sets up the functions in f_pond for the pond elements
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_pond by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call PD0_Setup(pd0)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be used
  !!   (in)    type ( FPD_COLLECTION ),pointer :: fpd
  !!             FPD_COLLECTION to be used
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( FPD_COLLECTION ),pointer :: fpd
  type ( IO_Status),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( PD0_POND ),pointer :: pnd

  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_Setup: PD0_Create has not been called",io )) return
    if (IO_Assert( (associated(fpd)), &
                    "PD0_Setup: Illegal FPD_COLLECTION object",io )) return
  endif

  do i = 1,pd0%iCount      
    ! Create a pond function in FPD for each pond
    pnd => pd0%Ponds(i)
    call FPD_New(fpd,pnd%cZ,pnd%rGamma,pnd%rRadius,pnd%iFPDIndex,io)
    if ( io%lError ) return
  end do

  return
end subroutine PD0_SetupFunctions
          
subroutine PD0_ComputeRHS(pd0,fpd,cCPZ,iEqType,cOrientation,rRHS,io)
  !! subroutine PD0_ComputeRHS
  !!
  !! Computes the contribution of the PD0 elements to the right-hand-side
  !! value for the specified equation parameters.
  !!
  !! Calling Sequence:
  !!    PD0_ComputeRHS(pd0,cCPZ,iEqType,cOrientation,rRHS)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be used
  !!   (in)    type ( FPD_COLLECTION ),pointer :: fpd
  !!             FPD_COLLECTION to be used
  !!   (in)    complex,dimension(:) :: cCPZ
  !!             Control point(s) to be used in coefficient calculations
  !!   (in)    integer :: iEqType
  !!             The type of equation to be used
  !!   (in)    complex :: cOrientation
  !!             Orientation unit vector (for discharge-based equations)
  !!   (inout) real :: rRHS
  !!             On return, rRHS is updated by adding the contribution of
  !!             this module to the previous value.
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( FPD_COLLECTION ),pointer :: fpd
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cCPZ
  complex (kind=ModAEM_Real),intent(in) :: cOrientation
  integer (kind=ModAEM_Integer), intent(in) :: iEqType
  real (kind=ModAEM_Real), intent(inout) :: rRHS
  type ( IO_Status),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( PD0_POND ),pointer :: pnd

  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_ComputeRHS: PD0_Create has not been called",io )) return
    if (IO_Assert( (associated(fpd)), &
                    "PD0_ComputeRHS: Illegal FPD_COLLECTION object",io )) return
  endif

!  do i=1,pd0%iCount
!    pnd => pd0%Ponds(i)
!    select case (iEqType)
!      case (EQN_HEAD)
!        rRHS = rRHS - real(cFPD_Potential(fpd,cCPZ(1),pnd%iFPDIndex,1))
!      case (EQN_FLOW)
!        rRHS = rRHS - rFPD_Flow(fpd,cCPZ,pnd%iFPDIndex,1)
!      case (EQN_INHO)
!        rRHS = rRHS - real(cFPD_Potential(fpd,cCPZ(1),pnd%iFPDIndex,1))
!      case (EQN_DISCHARGE)
!        rRHS = rRHS + aimag(cmplx(rZERO,rONE,ModAEM_Real) * &
!                      cFPD_Discharge(fpd,cCPZ(1),pnd%iFPDIndex,1) * &
!                      cOrientation/abs(cOrientation))
!      case (EQN_RECHARGE)
!        rRHS = rRHS
!      case (EQN_CONTINUITY)
!        rRHS = rRHS - rFPD_Extraction(fpd,pnd%iFPDIndex,1)
!    end select
!  end do 

  return
end subroutine PD0_ComputeRHS

function lPD0_CheckPoint(pd0,cZ,rTol,io) result(lSing)
  !! logical function lPD0_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lPD0_CheckPoint(pd0,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             The PD0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (PD0_COLLECTION),pointer :: pd0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_Status),pointer :: io

  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iWel
  type ( PD0_POND ),pointer :: pnd

  lSing = .false.
  do iWel = 1,pd0%iCount
    pnd => pd0%Ponds(iWel)
    if ( abs(real(pnd%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(pnd%cZ-cZ))<rTol ) then
      lSing = .true.
      exit
    endif
  end do
  
  return
end function lPD0_CheckPoint

subroutine PD0_Read(pd0,io)
  !! subroutine PD0_Read
  !!
  !! Reads the ponds for the specified PD0_COLLECTION from kIOInputLU
  !!
  !! Calling Sequence:
  !!    call PD0_Read(pd0)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be populated
  !!
  !! The format of the PD0 section of the input file appears as follows:
  !! PD0
  !! DIM NPonds
  !!     x y rsigma id
  !!     ... Up to NPonds
  !! END
  !!
  !! NOTE: It is assumed that the PD0 line was found by the caller
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( IO_Status),pointer :: io

  ! [ LOCAL DIRECTIVES ]
  type ( DIRECTIVE ),dimension(3),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD /)
  ! [ LOCALS ]
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  real (kind=ModAEM_Real) :: rGamma,rRad
  complex (kind=ModAEM_Real) :: cZ
  integer (kind=ModAEM_Integer) :: iID
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMaxWel
  logical (kind=ModAEM_Integer) :: lFlag
  type ( PD0_POND ),pointer :: pnd

  call IO_MessageText("  Reading PD0 module input",io)
  if ( io%lError ) return

  if (IO_Assert( (associated(pd0)), "PD0_Read: PD0_Create has not been called",io )) return

  ! Process input   
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return

    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "PD0_Read: I/O Error",io)) return
        exit
      case (kOpFileEOF)
        ! EOF is unexpected for all ModPD0 "ifXXXRead" routines. 
        ! Report the condition, but proceed as if EOD was found.
        if (IO_Assert( .false., "PD0_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(pd0%Ponds)), "PD0_Read: No space allocated",io)) return
        if (IO_Assert( (pd0%iCount<size(pd0%Ponds)), "PD0_Read: Space exhausted",io)) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ,rGamma,rRad,iID
        if (IO_Assert( (iStat==0), "PD0_Read: I/O Error",io)) return
        pd0%iCount = pd0%iCount+1
        pnd => pd0%Ponds(pd0%iCount)
        pnd%cZ = cZ
        pnd%rRadius = rRad
        pnd%rGamma = rGamma
        pnd%iID = iID
        ! No FPD is declared here; see PD0_Setup
        pnd%iFPDIndex = -1
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
    end select
  end do  

  call IO_MessageText("  Leaving PD0 module",io)

  return
end subroutine PD0_Read

subroutine PD0_Inquiry(pd0,iLU,io)
  !! subroutine PD0_Inquiry
  !!
  !! Writes an inquiry report for all ponds to iLU
  !!
  !! Calling Sequence:
  !!    call PD0_Inquiry(pd0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_Status),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( PD0_POND ),pointer :: pnd
  
  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_Inquiry: PD0_Create has not been called",io )) return
  endif

  do i=1,pd0%iCount
    pnd => pd0%Ponds(i)
    write ( unit=iLU, &
            fmt="(""PD0"",2("","",i9),4("","",e14.6))" &
          ) pnd%iID,pnd%cZ,pnd%rGamma,pnd%rRadius
  end do

  return
end subroutine PD0_Inquiry

subroutine PD0_Report(pd0,io)
  !! subroutine PD0_Report
  !!
  !! Writes a debugging report for all ponds to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call PD0_Report(pd0)
  !!
  !! Arguments:
  !!   (in)    type ( PD0_COLLECTION ),pointer :: pd0
  !!             PD0_COLLECTION to be used
  !!
  ! [ ARGUMENTS ]
  type ( PD0_COLLECTION ),pointer :: pd0
  type ( IO_Status),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( PD0_POND ),pointer :: pnd

  if ( io%lDebug ) then
    if (IO_Assert( (associated(pd0)), &
                    "PD0_Inquiry: PD0_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module PD0'//'Pond element information')" )

  if ( .not. associated(pd0%Ponds) ) then
    write ( unit=kIOOutputLU, &
            fmt="("" No Ponds Allocated""/)" &
          )
  else
    ! How many strings?
    write ( unit=kIOOutputLU, &
            fmt="("" Number of ponds: "",i5,""  used: "",i5)" &
          ) ubound(pd0%Ponds,1),pd0%iCount

    ! Write the detail information for all ponds...
    write ( unit=kIOOutputLU, &
            fmt="("" Ponds:"")" &
          )
    write ( unit=kIOOutputLU, &
            fmt= "(""    ID"",t15,""       X"",t30,""       Y"",t45, " // &
                 " ""  Discharge"",t60,""  Radius"",t75,""FPD Index"")" &
          )
    do i=1,pd0%iCount
      pnd => pd0%Ponds(i)
      write( unit=kIOOutputLU, &
             fmt="(i10,t15,d12.5,t30,d12.5,t45,d12.5,t60,d12.5,t75,i10)" &
           ) pnd%iID,pnd%cZ,pnd%rGamma,pnd%rRadius,pnd%iFPDIndex
    end do

    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iPD0_GetInfo(pd0,SIZE_FWL,io)
    if ( io%lError ) return

    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iPD0_GetInfo(pd0,SIZE_FPD,io)
    if ( io%lError ) return

    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iPD0_GetInfo(pd0,SIZE_FDP,io)
    if ( io%lError ) return

    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iPD0_GetInfo(pd0,SIZE_EQUATIONS,io)
    if ( io%lError ) return

    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iPD0_GetInfo(pd0,SIZE_UNKNOWNS,io)
  if ( io%lError ) return

  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine PD0_Report

end module m_pd0

