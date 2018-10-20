module m_wl0

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

  !! module m_wl0
  !!
  !! Element module for 2-D discharge specified wells
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! well elements. Each element requires a location, radius and discharge.

use u_constants
use u_io
use f_well

implicit none

public

  type :: WL0_WELL
    !! type WL0_WELL
    !!
    !! Type that holds information for one well
    !!
    !! Members:
    !!   complex :: cZC
    !!     The center of the well
    !!   real :: rRadius
    !!     The radius of the well
    !!   complex :: rDischarge
    !!     Discharge of the well; positive indicates extraction
    !!     and negative indicates injection.
    !!   integer :: iID
    !!     Identification label for the well (for interaction with e.g. GUIs)
    !!   integer :: iFWLIndex
    !!     Index for the well entry in the FWL module
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rDischarge
    real (kind=ModAEM_Real) :: rRadius
    integer (kind=ModAEM_Integer) :: iID
    integer (kind=ModAEM_Integer) :: iFWLIndex
  end type WL0_WELL

  type :: WL0_COLLECTION
    !! type WL0_COLLECTION
    !!
    !! Type that holds all wells in a layer
    !!
    !! Members:
    !!   type ( WL0_WELL ),dimension(:),pointer :: Wells
    !!     Array of WL0_WELL objects for the layer; dimensioned for the maximum
    !!     number of wells according to the input file (see WL0_Read)
    !!   integer :: iCount
    !!     The actual number of wells in use in the layer
    !! 
    type ( WL0_WELL ),dimension(:),pointer :: Wells
    integer (kind=ModAEM_Integer) :: iCount
  end type WL0_COLLECTION

contains

!**pd Modified to allocate ONLY the collection object
function WL0_Create(io) result (wl0)
  !! function WL0_Create
  !!
  !! Creates a new WL0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    wl0 => WL0_Create()
  !!
  !! Arguments:
  !!
  ! [ ARGUMENTS ]
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  type ( WL0_COLLECTION ),pointer :: wl0
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (wl0,stat=iStat)
  if (IO_Assert( (iStat == 0), "WL0_Create: allocation failed",io )) return
  nullify(wl0%Wells)
  wl0%iCount = 0

  return
end function WL0_Create

!**pd New WL0_Alloc subroutine
subroutine WL0_Alloc(wl0, iNWL, io)
  !! Subroutine WL0_Alloc
  !! 
  !! Allocates wells for the WL0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call WL0_Alloc(lw0, iNWL)
  !!
  !! Arguments:
  !!    (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!              The WL0_COLLECTION object to be used
  !!     (in)     integer :: iNWL
  !!                Maximum number of wells
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  integer (kind=ModAEM_Integer),intent(in) :: iNWL
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (wl0%Wells(iNWL),stat=iStat)
  if (IO_Assert( (iStat == 0), "WL0_Alloc: allocation failed",io )) return

end subroutine WL0_Alloc

!**pd  New WL0_Destroy subroutine
subroutine WL0_Destroy(wl0,io)
  !! subroutine WL0_Destroy
  !!
  !! Frees memory allocated for an WL0 Wells and WL0 Collection object
  !!
  !! Calling Sequence:
  !!     call WL0_Destroy(wl0)
  !!
  !! Arguments:
  !!  type ( WL0_COLLECTION ),pointer :: wl0
  !!              Pointer to the WL0_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                   "WL0_Destroy: WL0_Create has not been called",io )) return
  endif

  if (associated(wl0%Wells)) then
    deallocate (wl0%Wells,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                     "WL0_Destroy: deallocation of Wells failed",io )) return
  end if
  deallocate (wl0,stat=iStat)
  if (IO_Assert( (iStat == 0), "WL0_Destroy: deallocation failed",io )) return
  
  return
end subroutine WL0_Destroy


subroutine WL0_New(wl0,Well,io) 
  !! function WL0_New
  !!
  !! Adds a new WL0_WELL object to the WL0_COLLECTION 'wl0'
  !!
  !! Calling Sequence:
  !!    call WL0_New(wl0,Well)
  !!
  !! Arguments:
  !!    (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!              The WL0_COLLECTION object to be used
  !!    (in)    type ( WL0_WELL ),pointer :: Well
  !!              Vector that defines the points along the barrier
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  type ( WL0_WELL ) :: Well
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_New: WL0_Create has not been called",io )) return
  endif

  if (IO_Assert( (wl0%iCount < size(wl0%Wells)), &
                  "WL0_New: Space exhausted",io )) return

  wl0%iCount = wl0%iCount + 1
  wl0%Wells(wl0%iCount) = Well

  return
end subroutine WL0_New

function iWL0_GetInfo(wl0,iOption,io) result(iValue)
  !! function WL0_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iWL0_GetInfo(wl0,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_GetInfo: WL0_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FWL )
      iValue = wl0%iCount
    case default
      iValue = 0
  end select

  return
end function iWL0_GetInfo

subroutine WL0_SetupFunctions(wl0,fwl,io)
  !! subroutine WL0_Setup
  !!
  !! This routine sets up the functions in f_well for the well elements
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call WL0_Setup(wl0)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer :: fwl
  !!             FWL_COLLECTION to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( WL0_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_Setup: WL0_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "WL0_Setup: Illegal FWL_COLLECTION object",io )) return
  endif

  do i = 1,wl0%iCount      
    ! Create a well function in FWL for each well
    wel => wl0%Wells(i)
    call FWL_New(fwl,wel%cZ,wel%rDischarge,wel%rRadius,wel%iFWLIndex,io)
    if ( io%lError ) return
  end do

  return
end subroutine WL0_SetupFunctions
          
function lWL0_CheckPoint(wl0,cZ,rTol,io) result(lSing)
  !! logical function lWL0_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lWL0_CheckPoint(wl0,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             The WL0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (WL0_COLLECTION),pointer :: wl0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iWel
  type ( WL0_WELL ),pointer :: wel

  lSing = .false.
  do iWel = 1,wl0%iCount
    wel => wl0%Wells(iWel)
    if ( abs(real(wel%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(wel%cZ-cZ))<rTol ) then
      lSing = .true.
      exit
    endif
  end do
  
  return
end function lWL0_CheckPoint

subroutine WL0_FindWell(wl0,iWellID,cZWell,rDischarge,rRadius,iFWLIndex,lFound,io)
  !! subroutine WL0_FindWell
  !!
  !! Finds the well specified by the Well ID and returns its parameters
  !!
  !! Calling Sequence:
  !!    lFound = lWL0_FindWell(wl0,iWellID,cZWell,rDischarge,rRadius,iFWLIndex)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!   (in)    integer :: iWellID
  !!             The well ID number
  !!   (out)   complex :: cZWell
  !!             Location of the well
  !!   (out)   real :: rDischarge
  !!             The discharge of the well
  !!   (out)   real :: rRadius
  !!             The radius of the well
  !!   (out)   integer :: iFWLIndex
  !!             Index of the well in f_well
  !!   (out)   logical :: lFound
  !!             .true. if the well was found
  !!             .false. if the well was not found 
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  integer (kind=ModAEM_Integer),intent(in) :: iWellID
  complex (kind=ModAEM_Real),intent(out) :: cZWell
  real (kind=ModAEM_Real),intent(out) :: rRadius,rDischarge
  integer (kind=ModAEM_Integer),intent(out) :: iFWLIndex
  logical (kind=ModAEM_Integer) :: lFound
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( WL0_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_FindWell: WL0_Create has not been called",io )) return
  endif

  lFound = .false.
  do i=1,wl0%iCount
    wel => wl0%Wells(i)
    if ( wel%iID == iWellID ) then
      lFound = .true.
      cZwell = wel%cZ
      rDischarge = wel%rDischarge
      rRadius = wel%rRadius
      iFWLIndex = wel%iFWLIndex
    endif
  end do

  return
end subroutine WL0_FindWell

subroutine WL0_FindWellPointer(wl0,iWellID,well,lfound,io)
  !! subroutine WL0_FindWell
  !!
  !! Finds the well specified by the Well ID and returns a pointer to it
  !!
  !! Calling Sequence:
  !!    lFound = lWL0_FindWellPointer(wl0,iWellID,well,lfound)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!   (in)    integer :: iWellID
  !!             The well ID number
  !!   (out)   type ( WL0_WELL :: well
  !!             Pointer to the well
  !!   (out)   logical :: lFound
  !!             .true. if the well was found
  !!             .false. if the well was not found 
  !!   (in)    type ( IO_STATUS ),pointer :: io
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  integer (kind=ModAEM_Integer),intent(in) :: iWellID
  type ( WL0_WELL ),pointer :: well
  logical (kind=ModAEM_Integer) :: lFound
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i

  if ( io%lDebug ) then
    if ( IO_Assert( (associated(wl0)), &
                    "WL0_FindWell: WL0_Create has not been called",io ) ) return
  endif

  lFound = .false.
  do i=1,wl0%iCount
    well => wl0%Wells(i)
    if ( well%iID == iWellID ) then
      lFound = .true.
      return
    endif
  end do

  return
end subroutine WL0_FindWellPointer

function lWL0_CheckProximity(wl0,cZ,rProximity,iDir,iElementID,io) result(lFound)
  !! subroutine WL0_CheckProximity
  !!
  !! Checks to see if the point specified is within the radius of a well, then
  !! moves the point outside well perimeter if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    lFound = lWL0_CheckProximity(wl0,cZ,rProximity,iDir,iElementID)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
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
  !!             The element ID for the well found (if lFound is .true.)
  !! 
  !! Return Value:
  !!   (out)   logical :: lFound
  !!             .true. if a well is found
  !!             .false. if a well is not found
  !!   Note: See side effects above for iElementID
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
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
  type ( WL0_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_Inquiry: WL0_Create has not been called",io )) return
  endif

  ! Search through the wells, according to the iDir
  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do i=1,wl0%iCount
    wel => wl0%Wells(i)
    if ( iDir < 0 .and. wel%rDischarge >= rZERO ) then
      ! skip pumping wells
      cycle
    endif
    if ( iDir > 0 .and. wel%rDischarge <= rZERO ) then
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
end function lWL0_CheckProximity

subroutine WL0_Read(wl0,io)
  !! subroutine WL0_Read
  !!
  !! Reads the wells for the specified WL0_COLLECTION from kIOInputLU
  !!
  !! Calling Sequence:
  !!    call WL0_Read(wl0)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be populated
  !!
  !! The format of the WL0 section of the input file appears as follows:
  !! WL0
  !! DIM NWells
  !!     x y q r id
  !!     ... Up to NWells
  !!
  !! NOTE: It is assumed that the WL0 line was found by the caller
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  type ( IO_STATUS ),pointer :: io
  ! [ LOCAL DIRECTIVES ]
  type ( DIRECTIVE ),dimension(3),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD /)
  ! [ LOCALS ]
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  real (kind=ModAEM_Real) :: rDischarge,rRad
  complex (kind=ModAEM_Real) :: cZ
  integer (kind=ModAEM_Integer) :: iID
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMaxWel
  logical (kind=ModAEM_Integer) :: lFlag
  type ( WL0_WELL ),pointer :: wel

  call IO_MessageText("  Reading WL0 module input",io)

  if (IO_Assert( (associated(wl0)), "WL0_Read: WL0_Create has not been called",io )) return

  ! Process input   
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "WL0_Read: I/O Error",io)) return
        exit
      case (kOpFileEOF)
        ! EOF is unexpected for all ModWL0 "ifXXXRead" routines. 
        ! Report the condition, but proceed as if EOD was found.
        if (IO_Assert( .false., "WL0_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(wl0%Wells)), "WL0_Read: No space allocated",io)) return
        if (IO_Assert( (wl0%iCount<size(wl0%Wells)), "WL0_Read: Space exhausted",io)) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ,rDischarge,rRad,iID
        if (IO_Assert( (iStat==0), "WL0_Read: I/O Error",io)) return
        wl0%iCount = wl0%iCount+1
        wel => wl0%Wells(wl0%iCount)
        wel%cZ = cZ
        wel%rDischarge = rDischarge
        wel%rRadius = rRad
        wel%iID = iID
        ! No FWL is declared here; see WL0_Setup
        wel%iFWLIndex = -1
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        exit
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
    end select
  end do  

  call IO_MessageText("  Leaving WL0 module",io)

  return
end subroutine WL0_Read

subroutine WL0_Inquiry(wl0,iLU,io)
  !! subroutine WL0_Inquiry
  !!
  !! Writes an inquiry report for all wells to iLU
  !!
  !! Calling Sequence:
  !!    call WL0_Inquiry(wl0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  type ( WL0_WELL ),pointer :: wel
  
  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_Inquiry: WL0_Create has not been called",io )) return
  endif

  do i=1,wl0%iCount
    wel => wl0%Wells(i)
    write ( unit=iLU, &
            fmt="(""WL0"",2("","",i9),4("","",e14.6))" &
          ) wel%iID,wel%cZ,wel%rDischarge,wel%rRadius
  end do

  return
end subroutine WL0_Inquiry

subroutine WL0_Report(wl0,io)
  !! subroutine WL0_Report
  !!
  !! Writes a debugging report for all wells to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call WL0_Report(wl0)
  !!
  !! Arguments:
  !!   (in)    type ( WL0_COLLECTION ),pointer :: wl0
  !!             WL0_COLLECTION to be used
  !!
  ! [ ARGUMENTS ]
  type ( WL0_COLLECTION ),pointer :: wl0
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( WL0_WELL ),pointer :: wel

  if ( io%lDebug ) then
    if (IO_Assert( (associated(wl0)), &
                    "WL0_Inquiry: WL0_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module WL0'//'Well element information')" )

  if ( .not. associated(wl0%Wells) ) then
    write ( unit=kIOOutputLU, &
            fmt="("" No Wells Allocated""/)" &
          )
  else
    ! How many strings?
    write ( unit=kIOOutputLU, &
            fmt="("" Number of wells: "",i5,""  used: "",i5)" &
          ) ubound(wl0%Wells,1),wl0%iCount

    ! Write the detail information for all wells...
    write ( unit=kIOOutputLU, &
            fmt="("" Wells:"")" &
          )
    write ( unit=kIOOutputLU, &
            fmt= "(""    ID"",t15,""       X"",t30,""       Y"",t45, " // &
                 " ""  Discharge"",t60,""  Radius"",t75,""FWL Index"")" &
          )
    do i=1,wl0%iCount
      wel => wl0%Wells(i)
      write( unit=kIOOutputLU, &
             fmt="(i10,t15,d12.5,t30,d12.5,t45,d12.5,t60,d12.5,t75,i10)" &
           ) wel%iID,wel%cZ,wel%rDischarge,wel%rRadius,wel%iFWLIndex
    end do

    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iWL0_GetInfo(wl0,SIZE_FWL,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iWL0_GetInfo(wl0,SIZE_FPD,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iWL0_GetInfo(wl0,SIZE_FDP,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iWL0_GetInfo(wl0,SIZE_EQUATIONS,io)
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iWL0_GetInfo(wl0,SIZE_UNKNOWNS,io)
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine WL0_Report

end module m_wl0
