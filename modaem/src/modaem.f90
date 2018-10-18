! ModAEM.f90 -- Main program

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

! 
! This is a Fortran MODULE for general-purpose AEM functions. The regional aquifer
! properties, conversion from head-to-potential, potential/velocity computations, and 
! reporting functions are included here.
!

program modaem

  use u_constants
  use u_io
  use m_aqu
  use m_wl0
  use m_wl1
  use m_ls0
  use m_ls1
  use m_hb0
  use m_aem
  use m_inq
  use a_grid
  use a_tr0

  implicit none


  ! [ RETURN VALUE ]
  type ( AEM_DOMAIN ),pointer :: aem
  ! [ LOCALS ]
  type (DIRECTIVE),dimension(12),parameter :: dirDirectives = &
                   (/ dirEOD,dirAEM,dirRPT,dirSOL,dirGRI, &
                      dirHEA,dirPOT,dirDIS,dirFLO,dirGRA, &
                      dirINQ,dirTR0 /)
  ! Directives
  character (len=255) :: sCommandLine
  character (len=132) :: sBaseName
  character (len=132) :: sFileName
  character (len=132) :: sMessage
  character (len=132) :: sRecord
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iretval
!**pd  integer (kind=ModAEM_Integer) :: iNWL0
!**pd  integer (kind=ModAEM_Integer) :: iNPD0
!**pd  integer (kind=ModAEM_Integer) :: iNLS0
!**pd  integer (kind=ModAEM_Integer) :: iNLS1
!**pd  integer (kind=ModAEM_Integer) :: iNHB0
  integer (kind=ModAEM_Integer) :: iNIter
  integer (kind=ModAEM_Integer) :: i
  ! Placeholders for function module test calls
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZC,cZE1,cZE2
  real (kind=ModAEM_Real) :: rR,rTol,rDPhiDX,rDPhiDY
  type ( IO_Status ),pointer :: io
  
  ! First Create the IO_Status Object
  io => IO_Create()

  ! First, open the names file (with directives), modaem.nam
  open ( unit=kIOScratchLU, file='modaem.nam', status='old', iostat=iStat )
  if (IO_Assert( (iStat==0),"ModAEM: Could not open names file modaem.nam for input", io)) then
    call IO_ErrorReport("ModAEM",errNoData,iStat,io)
    call exit(1)
  endif
  !** FUTURE: Put processing directives (e.g. SQL connections and what-not) here!
  ! Retrieve the AEM file name from modaem.nam
  read (unit=kIOScratchLU,fmt=*,iostat=iStat) sBaseName
  if (IO_Assert( (iStat==0),"ModAEM: Could not read AEM base file name", io)) then
    call IO_ErrorReport("ModAEM",errNoData,iStat,io)
    call exit(1)
  endif
     
  call IO_OpenAll(sBaseName,io)
  nullify( aem )

  ! Here we go!
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation.
        if (IO_Assert( .false., "AEM_Read: I/O Error",io)) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
      case (kOpFileEOF)
        ! EOF is unexpected for all ModAEM "ifXXXRead" routines. 
        if (IO_Assert( .false., "AEM_Read: Unexpected EOF",io)) then
          call IO_ErrorReport("ModAEM",errUnexpectedEOF,0,io)
          call exit(1)
        endif
 
      case (kOpData)
        ! A data line was found. The main ModAEM module has no default
        ! data lines. Report the condition.
        if (IO_Assert( .false., "AEM_Read: Unexpected data record",io)) then
          call IO_ErrorReport("ModAEM",errUnexpectedData,iStat,io)
          call exit(1)
        endif
 
      case (kOpEOD)
        ! EOD mark was found. Exit the file parser.
        call IO_MessageText("EOD Encountered -- Job complete",io)
        exit
      case (kOpAEM)
        ! Populate the AEM_DOMAIN object from the input file according to these sizes...
        if (IO_Assert( (.not. associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has already been created",io)) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        aem => AEM_Create(io)  !**pd
        call AEM_Read(aem,io)
      case (kOpRPT)
        ! Generates a report
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created", io)) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        call AEM_Report(aem,io)
        if ( associated(aem%aqu) ) call AQU_Report(aem%aqu,io)
        if ( associated(aem%wl0) ) call WL0_Report(aem%wl0,io)
        if ( associated(aem%wl1) ) call WL1_Report(aem%wl1,io)
        if ( associated(aem%pd0) ) call PD0_Report(aem%pd0,io)
        if ( associated(aem%ls0) ) call LS0_Report(aem%ls0,io)
        if ( associated(aem%ls1) ) call LS1_Report(aem%ls1,io)
        if ( associated(aem%hb0) ) call HB0_Report(aem%hb0,io)
        if ( associated(aem%ls2) ) call LS2_Report(aem%ls2,io)
        if ( associated(aem%fwl) ) call FWL_Report(aem%fwl,io)
        if ( associated(aem%fpd) ) call FPD_Report(aem%fpd,io)
        if ( associated(aem%fdp) ) call FDP_Report(aem%fdp,io)
        if ( associated(aem%mat) ) call MAT_Report(aem%mat,'final',io)
      case (kOpSOL)
        ! Enter the SOLver module
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created" ,io)) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
          
        read (unit=sRecord,fmt=*,iostat=iStat) iNIter
        call AEM_Solve(aem,iNIter,io)
      case (kOpGRI)
        ! Enter the GRId module
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 

        call GRI_Read(aem,io)
      case (kOpINQ)
        ! Enter the inquiry module
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
          
        call INQ_Read(aem,io)
      case (kOpTR0)
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        ! Enter the TR0 2-D tracing module
        if (IO_Assert( (associated(aem)), "AEM_Read: No DIM statement",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        call TR0_Read(aem,io)
      case (kOpHEA)
        ! Report the head to the error LU
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
         
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) then
          call IO_ErrorReport("ModAEM",errMissingArgument,iStat,io)
          call exit(1)
        endif
  
        write ( unit=sMessage, &
                fmt="(""      >> Head       at "",d13.6,1x,d13.6,""  =  "",d13.6,1x,d13.6)" &
              ) cZ1,rAEM_Head(aem,cZ1,io)
        call IO_ErrorText(sMessage,io)
      case (kOpPOT)
        ! Report the complex potential to the error LU
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) then
          call IO_ErrorReport("ModAEM",errMissingArgument,iStat,io)
          call exit(1)
        endif
 
        write ( unit=sMessage, &
                fmt="(""      >> Potential  at "",d13.6,1x,d13.6,i5,""  =  "",d13.6,1x,d13.6)" &
              ) cZ1,cAEM_Potential(aem,cZ1,io)
        call IO_ErrorText(sMessage,io)
      case (kOpGRA)
        ! Report the estimated gradient in Phi 
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1,rTol
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) then
          call IO_ErrorReport("ModAEM",errMissingArgument,iStat,io)
          call exit(1)
        endif
 
        rDPhiDX = ( real(cAEM_Potential(aem,cZ1-cmplx(rTol,rZERO,ModAEM_Real),io)) - &
                    real(cAEM_Potential(aem,cZ1+cmplx(rTol,rZERO,ModAEM_Real),io)) ) / (rTWO*rTol)
        rDPhiDY = ( real(cAEM_Potential(aem,cZ1-cmplx(rZERO,rTol,ModAEM_Real),io)) - &
                    real(cAEM_Potential(aem,cZ1+cmplx(rZERO,rTol,MOdAEM_Real),io)) ) / (rTWO*rTol)
        write ( unit=sMessage, &
                fmt="(""      >> Gradient   at "",d13.6,1x,d13.6,""  =  "",d13.6,1x,d13.6)" &
              ) cZ1,rDPhiDX,rDPhiDY
        call IO_ErrorText(sMessage,io)
      case (kOpDIS)
        ! Report the complex discharge to the error LU
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) then
          call IO_ErrorReport("ModAEM",errMissingArgument,iStat,io)
          call exit(1)
        endif
 
        write ( unit=sMessage, &
                fmt="(""      >> Discharge  at "",d13.6,1x,d13.6,i5,""  =  "",d13.6,1x,d13.6)" &
              ) cZ1,cAEM_Discharge(aem,cZ1,io)
        call IO_ErrorText(sMessage,io)
      case (kOpFLO)
        ! Report the total flow to the error LU
        if (IO_Assert( (associated(aem)), &
                        "AEM_Read: the AEM_DOMAIN has not been created",io )) then
          call IO_ErrorReport("ModAEM",errInvalidDirective,0,io)
          call exit(1)
        endif
 
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1,cZ2
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) then
          call IO_ErrorReport("ModAEM",errMissingArgument,iStat,io)
          call exit(1)
        endif
 
        write (unit=sMessage, &
               fmt="(""      >> Flow       between "",d13.6,1x,d13.6,"" and "",d13.6,1x,d13.6,""  =  "",d13.6)" &
              ) cZ1,cZ2,rAEM_Flow(aem, (/cZ1,cZ2/), io)
        call IO_ErrorText(sMessage,io)
      case default
    end select
  end do

  !**pd Destroy the aem object and all member objects
  if (associated(aem)) call AEM_Destroy(aem,io)

  call exit(0)
end program ModAEM

