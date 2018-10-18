module m_inq

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

  ! Data extraction/inquiry/check module
  ! V.A.Kelson  10/23/98
  !
  ! This module allows for the extraction of results from a ModAEM solution
  !
  ! Commands:
  !
  !   FIL  <filename>               Output file. No extension will be appended.
  !
  !   HEA  (X,Y) L label            Extracts the head at the point (X,Y) in layer L
  !
  !   DIS  (X,Y) L label            Extracts the discharge at the point (X,Y) in layer L
  !
  !   VEL  (X,Y) L label            Extracts the velocity at the point (X,Y) in layer L
  !
  !   FLO  (X1,Y1) (X2,Y2) L label  Extracts the integrated flow across the specified path
  !
  !   WL0  L                        Writes check information for the specified module
  !   WL1  L
  !   LS0  L
  !   LS1  L
  !   HB0  L  
  !
  ! For all extraction commands which requre a label, the label is returned in the extraction file.
  ! All data are written to space-delimeted files with constant column widths.  A header record is
  ! prepended to the file.   

use u_constants
use u_io
use m_wl0
use m_wl1
use m_ls0
use m_ls1
use m_hb0
use m_aem

implicit none

public
  
  ! Global variables for this module

  ! Flags
  logical (kind=ModAEM_Integer),private,save :: fINQFileOpen = .false.

contains

subroutine INQ_Read(aem,io)
  ! This reads and processes the input commands for the INQ package 

  ! Argument list
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io

  ! Locals -- for Directive parsing
  type (DIRECTIVE),dimension(15),parameter :: dirDirectives = &
                                             (/ dirFIL, dirDBG, dirPCD, &
                                                dirHEA, dirDIS, dirPOT, dirVEL, dirFLO, &
                                                dirWL0, dirWL1, dirLS0, &  !**pd
                                                dirLS1, dirHB0, dirAEM, dirEND /)
  ! Locals -- Input values
  character (len=132) :: sRecord   ! Record from input file
  integer (kind=ModAEM_Integer) :: iOpCode         ! Opcode 
  integer (kind=ModAEM_Integer) :: iError          ! Error code
  integer (kind=ModAEM_Integer) :: iStat           ! RunTime error code
  integer (kind=ModAEM_Integer) :: iKeyCode        ! Key code for pause lines
  ! Placeholders for tracing parameters
  real (kind=ModAEM_Real) :: rX1,rX2,rY1,rY2,rValue
  complex (kind=ModAEM_Real) :: cValue
  integer (kind=ModAEM_Integer) :: iLabel
  logical (kind=ModAEM_Integer) :: lFlag
  character (len=132) :: sFile

  call IO_MessageText("Reading INQ module input",io)
  
  ! Clear the status flags
  fINQFileOpen = .false.
  call IO_MessageText("Entering inquiry module INQ",io)
  if (IO_Assert( (associated(aem)), "INQ_Read: No AEM_DOMAIN object",io )) return

  ! The remainder of this routine uses ifIOInputRecord to process 
  ! the model input file.
  iError = errOK
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation.
        if (IO_Assert( .false., "INQ_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModGRI "ifXXXRead" routines. 
        if (IO_Assert( .false., "INQ_Read: Unexpected EOF",io )) return
      case (kOpEND)
        ! END mark was found. Exit the file parser.
        if ( fINQFileOpen ) then
          close(unit=kIOInqLU)
        endif
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
      case (kOpFIL)
        !****************************************************************************
        ! Here for the FIL command -- open the output file
        !****************************************************************************
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        ! If the output file is open, close it now.
        if ( fINQFileOpen ) then
          close (unit=kIOInqLU)
          fINQFileOpen = .false.
        endif
        ! Open the new output file.  If the open fails, write a message.  
        sfile = trim(sfile)
        open (unit=kIOInqLU,  file=sfile, action="WRITE", status="REPLACE", iostat=iStat)
        if (IO_Assert( (iStat==0), "INQ_Read: Open failed",io )) return
        fINQFileOpen = .true.
      case ( kOpHEA )
        !****************************************************************************
        ! Here for the HEA command -- extract the head
        !****************************************************************************
        ! Retrieve the location and label
        read (unit=sRecord,fmt=*,iostat=iStat) rX1,rY1,iLabel
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        if (IO_Assert( fINQFileOpen, "INQ_Read: No output file",io )) return
        rValue = rAEM_Head(aem,cmplx(rX1,rY1,ModAEM_Real),io)
        write (unit=kIOInqLU,fmt="(""HEA"","","",i9,1x,3("","",e14.6))") iLabel,rX1,rY1,rValue
      case ( kOpDIS )
        !****************************************************************************
        ! Here for the DIS command -- extract the discharge
        !****************************************************************************
        ! Retrieve the location and label
        ! Retrieve the maximum time 
        read (unit=sRecord,fmt=*,iostat=iStat) rX1,rY1,iLabel
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        if (IO_Assert( fINQFileOpen, "INQ_Read: No output file",io )) return
        cValue = cAEM_Discharge(aem,cmplx(rX1,rY1,ModAEM_Real),io)
        write (unit=kIOInqLU,fmt="(""DIS"","","",i9,1x,4("","",e14.6))") iLabel,rX1,rY1,cValue
      case ( kOpPOT )
        !****************************************************************************
        ! Here for the POT command -- extract the potential
        !****************************************************************************
        ! Retrieve the location and label
        read (unit=sRecord,fmt=*,iostat=iStat) rX1,rY1,iLabel
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        if (IO_Assert( fINQFileOpen, "INQ_Read: No output file",io )) return
        cValue = cAEM_Potential(aem,cmplx(rX1,rY1,ModAEM_Real),io)
        write (unit=kIOInqLU,fmt="(""POT"","","",i9,1x,4("","",e14.6))") iLabel,rX1,rY1,cValue
      case ( kOpVEL )
        !****************************************************************************
        ! Here for the VEL command -- extract the velocity
        !****************************************************************************
        ! Retrieve the location and label
        ! Retrieve the maximum time 
        read (unit=sRecord,fmt=*,iostat=iStat) rX1,rY1,iLabel
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        if (IO_Assert( fINQFileOpen, "INQ_Read: No output file",io )) return
        cValue = cAEM_Velocity(aem,cmplx(rX1,rY1,ModAEM_Real),io)
        write (unit=kIOInqLU,fmt="(""VEL"","","",i9,1x,4("","",e14.6))") iLabel,rX1,rY1,cValue
      case ( kOpFLO )
        !****************************************************************************
        ! Here for the FLO command -- extract the flow between points
        !****************************************************************************
        ! Retrieve the location and label
        ! Retrieve the maximum time 
        read (unit=sRecord,fmt=*,iostat=iStat) rX1,rY1,rX2,rY2,iLabel
        if (IO_Assert( (iStat==0), "INQ_Read: I/O Error",io )) return
        if (IO_Assert( fINQFileOpen, "INQ_Read: No output file",io )) return
        rValue = rAEM_Flow(aem, (/ cmplx(rX1,rY1,ModAEM_Real),cmplx(rX2,rY2,ModAEM_Real) /) ,io)
        write (unit=kIOInqLU,fmt="(""FLO"","","",i9,1x,5("","",e14.6))") iLabel,rX1,rY1,rX2,rY2,rValue
      case ( kOpWL0 )
        !****************************************************************************
        ! Here for the WL0 command -- extract information about WL0 elements
        !****************************************************************************
        call WL0_Inquiry(aem%wl0,kIOInqLU,io)
      case ( kOpWL1 ) !**pd
        !****************************************************************************
        ! Here for the WL1 command -- extract information about WL1 elements
        !****************************************************************************
        call WL1_Inquiry(aem%wl1,kIOInqLU,io) !**pd
      case ( kOpLS0 )
        !****************************************************************************
        ! Here for the LS0 command -- extract information about LS0 elements
        !****************************************************************************
        call LS0_Inquiry(aem%ls0,kIOInqLU,io)
      case ( kOpLS1 )
        !****************************************************************************
        ! Here for the LS1 command -- extract information about LS1 elements
        !****************************************************************************
        call LS1_Inquiry(aem%ls1,kIOInqLU,io)
      case ( kOpHB0 )
        !****************************************************************************
        ! Here for the HB0 command -- extract information about HB0 elements
        !****************************************************************************
        call HB0_Inquiry(aem%hb0,kIOInqLU,io)
      case ( kOpAEM )
        !****************************************************************************
        ! Here for the AEM command -- extract information about AEM elements
        !****************************************************************************
        call AEM_Inquiry(aem,kIOInqLU,io)
      case default
    end select
  end do

  call IO_MessageText("Leaving INQ module",io)

  return
end subroutine INQ_Read

end module m_inq

