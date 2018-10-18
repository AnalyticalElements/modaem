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


!! io.f90 -- Common I/O framework for ModAEM
!! 
!! This is a Fortran MODULE for common Input/Output functions
!!
!! Revision History
!!     v0.1  1/20/97  Vic Kelson
!!                    Initial development version
!!     v1.0  6/23/98  Vic Kelson
!!                    Production (WhAEM) version
!!     v1.1  4/02/01  Vic Kelson
!!                    Code cleanup for LGPL release plans
!!
!! This module contains common functions for Input/Output operations,
!! and provides an overall framework for ModAEM I/O.  All ModAEM I/O
!! is handled by use of this framework. As a result, the code can 
!! be easily ported to new I/O frameworks, e.g. GUIs or networks.

module u_io

use u_constants

implicit none

public

  type,public :: IO_STATUS
    logical :: lFilesOpen
    logical :: lDebug
    logical :: lProceed
    character (len=255),dimension(10) :: sFileStack
    integer (kind=ModAEM_Integer),dimension(10) :: iRecordStack
    integer (kind=ModAEM_Integer) :: iStackPtr
    logical (kind=ModAEM_Integer) :: lError
    character (len=80) :: sErrorText
  end type IO_STATUS

  ! I/O Logical Unit Assignments
  ! Only the output LUs are a 'PUBLIC' parameter. It is expected that the
  ! remaining LUs are programmatically accessed by support routines.
  integer (kind=ModAEM_Integer),private,parameter :: kIOInputLU = 20
  integer (kind=ModAEM_Integer),parameter :: kIOScratchLU = 21
  integer (kind=ModAEM_Integer),private,parameter :: kIOErrorLU = 22
  integer (kind=ModAEM_Integer),parameter  :: kIOOutputLU = 23
  integer (kind=ModAEM_Integer),parameter  :: kIOGridLU = 24
  integer (kind=ModAEM_Integer),parameter  :: kIOInqLU = 25
  integer (kind=ModAEM_Integer),parameter  :: kIOCheckLU = 26
  integer (kind=ModAEM_Integer),parameter  :: kIOTraceLU = 27

contains

function IO_Create() result(io)
  ! Creates an IO_STATUS object
  type ( IO_STATUS ),pointer :: io

  allocate (io)
  io%lFilesOpen = .false.
  io%lDebug = .false.
  io%lProceed = .false.
  io%sFileStack = ' '
  io%iRecordStack = 0
  io%iStackPtr = 0
  io%lError = .false.
  io%sErrorText = ' '


end function IO_Create

subroutine IO_OpenAll(sBaseName,io)
  !! subroutine IO_OpenAll
  !!
  !! Opens the standard files for a ModAEM run. The standard I/O channels
  !! opened are: (1) the input LU [kIOInputLU]; (2) the error output LU
  !! [kIOErrorLU]; (3) the report output LU [kIOOutputLU]. All input is
  !! managed using the IO_InputRecord routine -- it is declared PRIVATE.
  !! All error reporting is performed by the IO_ErrorReport and IO_ErrorText
  !! routines -- it is also PRIVATE. General purpose reporting from add-on
  !! modules is allowed -- write to kIOOutputLU. For convenience, LUs for
  !! inquiry and analysis tools are also declared here.
  !!
  !! Arguments:
  !!    character (len=*) :: sBaseName
  !!        Contains the "base" file name for the standard LUs. The
  !!        extensions added are: .aem (input), .out (output) and
  !!        .err (error).

  ! [ ARGUMENTS ]
  character (len=*),intent(in) :: sBaseName
  type ( IO_STATUS ),pointer :: io
  ! [ RESULT VALUE ]
  integer (kind=ModAEM_Integer) :: iRV
  ! [ LOCALS ]
  character (len=255) :: sFName
  integer (kind=ModAEM_Integer) :: iStatus                   

  ! Initialize the input LU.  This will always be a disk file with
  ! the extension ".aem"
  sFName = trim(sBaseName) // ".aem"
  call IO_MessageText("Processing file " // trim(sFName) // "...",io )
  open (unit=kIOInputLU,  file=sFName, action="READ", status="OLD", iostat=iStatus)
  if ( IO_Assert( (iStatus==0), "IO_OpenAll: Could not open input file " // trim(sFName),io) ) return

  ! Initialize the model output LU. This will always be a disk file with
  ! the extension ".out". This is the only LU which is directly written
  ! by the element modules output routines (via PUBLIC PARAMETER kIOOutputLU
  sFName = trim(sBaseName) // ".out"
  open (unit=kIOOutputLU,  file=sFName, action="WRITE", status="REPLACE", iostat=iStatus)
  if ( IO_Assert( (iStatus==0), "IO_OpenAll: Could not open input file " // trim(sFName),io) ) return

  ! Initialize the model output LU. This will always be a disk file with
  ! the extension ".out". This is the only LU which is directly written
  ! by the element modules" output routines (via PUBLIC PARAMETER kIOOutputLU
  sFName = trim(sBaseName) // ".err"
  open (unit=kIOErrorLU,  file=sFName, action="WRITE", status="REPLACE", iostat=iStatus)
  if ( IO_Assert( (iStatus==0), "IO_OpenAll: Could not open input file " // trim(sFName),io) ) return
  
  !**********************************************************************
  ! Files have been successfully opened!
  ! Add any additional I/O initialization code here
  !**********************************************************************
  io%lFilesOpen = .true.
  io%iStackPtr = 1
  io%sFileStack(io%iStackPtr) = sBaseName // ".aem"
  io%iRecordStack(io%iStackPtr) = 0

  return
end subroutine IO_OpenAll

subroutine IO_CloseAll(io)
  !! subroutine IO_CloseAll
  !!
  !! Closes the standard I/O channels. NOTE: Any other open files MUST
  !! be closed by the analysis routines that opened them.
  !!
  !! Arguments: (none)
  type ( IO_STATUS ),pointer :: io

  if ( io%lFilesOpen ) then
    close(unit=kIOInputLU)
    close(unit=kIOOutputLU)
    close(unit=kIOErrorLU)
    io%lFilesOpen = .false.
  endif

end subroutine IO_CloseAll

subroutine IO_InputRecord(dirDirectives,iOpCode,cResult,iKeyCode,io)
  ! iIO_InputRecord reads a record from the input file, then interprets
  ! it.  The valid directives for the current module are provided as the
  ! first parameter.  As the record is read, it is recorded (with the 
  ! current line number) using IO_ErrorReport() and IOOutputReport().
  !
  ! Returns in iOpCode:
  !  kOpError  If a runtime error is detected
  !  kOpEof    If end-of-file is reached
  !  kOpData   If no directive is found
  !  (OpCode)  If read is successful and a directive is found, the OpCode for
  !            the directive is returned (for use in a SELECT CASE structure)
  !
  !  ALSO: For return values kOpPause,kOpData and directive OpCodes
  !        columns 6-end of record are returned in ca255Record
  !  ALSO: For "pause" lines, the function returns the implementation-
  !        specific keycode from ifIOMessageText()
  !
  ! Comment lines: 
  !         #  If a # is found in column 1, the record is considered to be a comment
  !            and is skipped by this routine.
  !         &  If a & is found in column 1, the record is a message to the
  !            user.  It is reported using IO_MessageText() 
  !   
  type (DIRECTIVE),dimension(:),intent(in) :: dirDirectives
  integer (kind=ModAEM_Integer), intent(out) :: iOpCode
  character (len=*),intent(out) :: cResult
  integer (kind=ModAEM_Integer),intent(out) :: iKeyCode
  type ( IO_STATUS ),pointer :: io
  ! Locals
  character (len=132) :: sRawRecord,sRecord
  character (len=255) :: sMessage
  character (len=3) :: sTestDir
  integer (kind=ModAEM_Integer) :: iStatus
  integer (kind=ModAEM_Integer) :: i,ic
  integer (kind=ModAEM_Integer) :: iop

  ! Here we go.  Proceed only if the file is open!
  ! The input loop reads records until valid program input is found,
  ! processing comments, user interaction directives (from "pause" lines)
  ! and scanning for program directives or errors.
  if ( io%lFilesOpen ) then
    do 
      ! Get a record from the input file, check for runtime error
      read ( unit=kIOInputLU,fmt="(a132)",iostat=iStatus ) sRawRecord
      if ( iStatus /= 0 ) then
        if ( iStatus == -1 ) then
          iOpCode = kOpFileEOF
          cResult = ""
          iKeyCode = 0
          exit
        else
          call IO_ErrorReport("IOInputRecord",errRunTime,iStatus,io) 
          iOpCode = kOpError
          cResult = ""
          iKeyCode = 0
          exit
        endif
      else
        ! Here if the file I/O is successful.
        ! Get only the text starting from the first non-blank into the processing record
        if ( len_trim(sRawRecord) == 0 ) then
          cycle
        else
          do ic=1,len_trim(sRawRecord)
            if ( sRawRecord(ic:ic) /= ' ' ) then
              sRecord = trim(sRawRecord(ic:))
              exit
            endif
          end do
        endif
        io%iRecordStack(io%iStackPtr) = io%iRecordStack(io%iStackPtr)+1         
        write (unit=sMessage,fmt="(1x,i5.5,"">> "",a)") io%iRecordStack(io%iStackPtr),trim(sRawRecord)
        call IO_ErrorText(sMessage,io)       ! Echo to error LU
        if ( iStatus /= 0 ) then
          iOpCode = -1
        else if ( sRecord(1:1) == "#" ) then     ! Comment!
          cycle
        else if ( sRecord(1:1) == "&" ) then     ! Comment! Echo as message...
          call IO_MessageText(trim(sRecord(1:)),io)
          cycle
        else
          ! Phew!  Here"s where we scan for program directives. 
          iop = -1                                  ! Flag -- no directive found
          do ic=1,3
            if ( ichar(sRecord(ic:ic))>96 .and. ichar(sRecord(ic:ic))<122 ) then
              sTestDir(ic:ic) = char(ichar(sRecord(ic:ic))-32)
            else
              sTestDir(ic:ic) = sRecord(ic:ic)
            endif
          end do 
          do i=1,ubound(dirDirectives,1)
            if ( sTestDir == dirDirectives(i)%Text ) then
              iop = dirDirectives(i)%OpCode 
              exit
            endif
          end do
          ! Done looking -- was an opcode found??
          if ( iop > 0 ) then
            iOpCode = iop
            cResult = sRecord(4:)
            iKeyCode = 0 
            exit                   ! Exit outer loop -- return to caller
          else
            ! Treat it as a data record
            iOpCode = kOpData
            cResult = sRecord
            iKeyCode = 0
            exit                                      ! Return to caller
          endif
        endif
      end if         ! Read OK block
    end do           ! Outer I/O loop
  else               ! File is open test
  call IO_MessageText("IO_InputRecord: Files are not open",io)
  endif
end subroutine IO_InputRecord

subroutine IO_ScanForEND()
  ! Scans the input file until an END directive is found -- handles unimplemented modules
  ! as placeholders for test codes, and allows for skipping when fatal errors occur
  ! in element module readers.
  type ( IO_STATUS ),pointer :: io
  !
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iKeyCode
  character (len=132) :: sRecord
  type (DIRECTIVE),dimension(1),parameter :: dirDirectives = (/ dirEND /)

  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        exit    
      case (kOpFileEOF)
        exit
      case (kOpEND)
        exit
      case default
        call IO_ErrorText("          [Skipped]",io)
    end select
  end do

  return
end subroutine IO_ScanForEND

subroutine IO_ErrorReport(cRoutine,iErrCode,iRuntimeErrCode,io)
  ! Reports the error according to the specified error code
  character (len=*),intent(in) :: cRoutine
  integer (kind=ModAEM_Integer),intent(in) :: iErrCode
  integer (kind=ModAEM_Integer),intent(in) :: iRuntimeErrCode
  type ( IO_STATUS ),pointer :: io
  ! Locals
  character (len=255) :: sMessage
  integer (kind=ModAEM_Integer) :: i

  if ( iErrCode == errRunTime ) then
    write(unit=sMessage,fmt="(""      ??  "",a,"": runtime error "",i4,"" detected in line "",i5)") cRoutine, &
         iRuntimeErrCode,io%iRecordStack(io%iStackPtr)
  else
    do i=1, ubound(msgErrorMessages,1)
      if ( msgErrorMessages(i)%Code == iErrCode ) then
        write(unit=sMessage,fmt="(""      ??  "",a,"": error "",i4,"" detected in line "",i5,"" - "",a)") cRoutine, &
             iErrCode,io%iRecordStack(io%iStackPtr),trim(msgErrorMessages(i)%Message)
        exit
      endif
    end do
  endif
  call IO_ErrorText(sMessage,io)
  call IO_MessageText(sMessage(7:),io)
end subroutine IO_ErrorReport

subroutine IO_ErrorText(cMessage,io)
  ! Writes a line of text to the error LU
  character (len=*),intent(in) :: cMessage
  type ( IO_STATUS ),pointer :: io

  if ( io%lFilesOpen ) then
    write(unit=kIOErrorLU,fmt="(1x,a)") trim(cMessage)
  else
    write(unit=*,fmt="(1x,a)") trim(cMessage)
  endif
end subroutine IO_ErrorText

subroutine IO_MessageText(sMessage,io)
  ! Writes a run-time message to STDOUT
  character (len=*),intent(in) :: sMessage
  type ( IO_STATUS ),pointer :: io
  ! Locals

  ! Send the message...
  write(unit=*,fmt="(a)") trim(sMessage)
  if ( io%lFilesOpen ) then
    call IO_ErrorText("      **  " // trim(sMessage),io)
  endif

end subroutine IO_MessageText

function IO_Assert(lTest,sText,io) result(lError)
  ! Asserts that lTest is .true.  If .not. lTest, then print the message 
  ! sText and halt execution.  Used for debug error testing.
  logical :: lTest
  character (len=*),intent(in) :: sText
  type ( IO_STATUS ),pointer :: io
  logical :: lError

  if ( .not. lTest ) then
    io%lError = .true.
    io%sErrorText = sText
    lError = .true.
  else
    lError = .false.
  endif

  return
end function IO_Assert

subroutine IO_SetDebug(flag,io)
  ! Sets the IO_Debug flag value
  logical (kind=ModAEM_Integer),intent(in) :: flag
  type ( IO_STATUS ),pointer :: io

  io%lDebug = flag

  return
end subroutine IO_SetDebug

subroutine IO_SetProceed(flag,io)
  ! Sets the IO_Proceed flag value
  logical (kind=ModAEM_Integer),intent(in) :: flag
  type ( IO_STATUS ),pointer :: io

  io%lProceed = flag

  return
end subroutine IO_SetProceed

end module u_io
