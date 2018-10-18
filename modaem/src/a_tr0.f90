module a_tr0
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

  ! Analysis module for streamline tracing in 2-D.
  ! V.A.Kelson  10/20/98
  !
  ! This module allows for the tracing of streamlines in two-dimensions, and supports
  ! a variety of options for output and selection of release locations
  !
use u_constants
use u_io
use m_wl0
use m_ls0
use m_ls1
use m_hb0
use m_wl1
use m_aem

implicit none

public
  
  ! Global variables for this module
  real (kind=ModAEM_Real),private,save :: rTR0WinX1,rTR0WinY1,rTR0WinX2,rTR0WinY2
  real (kind=ModAEM_Real),private,save :: rTR0MaxTime
  real (kind=ModAEM_Real),private,save :: rTR0Step,rTR0Prox,rTR0Frac,rTR0WellProx
  real (kind=ModAEM_Real),private,save :: rTR0DispL,rTR0Retard,rTR0Decay

  ! Local constants
  integer (kind=ModAEM_Integer),private,parameter :: kTR0Timeout = 0
  integer (kind=ModAEM_Integer),private,parameter :: kTR0LeftWindow = 1
  integer (kind=ModAEM_Integer),private,parameter :: kTR0HitElement = 2
  integer (kind=ModAEM_Integer),private,parameter :: kTR0AquiferDry = 3
  integer (kind=ModAEM_Integer),private,parameter :: kTR0Stagnation = 4
  ! Check factor for termination at HB0 elements
  ! Default step factor -- fraction of tracing window size to be used for stepsize
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_STEP_FACTOR = 0.006_ModAEM_Real
  ! Default proximity factor -- number of step sizes away from elements for reduced step size
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_PROX = 10.0_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_WELLPROX = 10.0_ModAEM_Real
  ! Termination factors -- multiply by step for the termination tolerance at BCs
  real (kind=ModAEM_Real),private,parameter :: TR0_WL0_TERMINATION_FACTOR = 0.2_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: TR0_LS0_TERMINATION_FACTOR = 0.2_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: TR0_LS1_TERMINATION_FACTOR = 0.2_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: TR0_HB0_TERMINATION_FACTOR = 0.1_ModAEM_Real
  real (kind=ModAEM_Real),private,parameter :: TR0_WL1_TERMINATION_FACTOR = 0.2_ModAEM_Real
  ! Default fraction of basic step size to use near elements
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_FRAC = 0.1_ModAEM_Real
  ! Default maximum tracing time
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_MAX_TIME = HUGE(ModAEM_Real)
  ! Default transport parameters (not yet implemented)
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_DISP_L = rZERO
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_RETARDATION = rONE
  real (kind=ModAEM_Real),private,parameter :: TR0_DEFAULT_DECAY = rZERO
  ! Stagnation tolerance (number of step sizes the particle must move in 20 steps)
  real (kind=ModAEM_Real),private,parameter :: TR0_STAGNATION_TOLERANCE = 1.2_ModAEM_Real
  ! Distance factor for releasing particles from a well
  real (kind=ModAEM_Real),private,parameter :: TR0_WL0_DISTANCE_FACTOR = 2.0_ModAEM_Real

contains

subroutine TR0_Read(aem,io)
  ! This reads and processes the input commands for the TR0 package 
  ! Commands:
  !
  !   FIL  <filename>    -        Output file. The extension '.tr0' will be added to the filename.
  !                               For each pathline, the output file will contain:
  !                               A header line of "START X0,Y0,Z0,T0,C0,L,Dir", where:
  !                                 (X0,Y0,Z0) - Starting coordinates
  !                                 (T0) - Starting time (usually 0.0)
  !                                 (C0) - Starting concentration
  !                                 (L)  - Layer for tracing
  !                                 (Dir) - Direction (-1 backward) (1 forward)
  !                               For each point in each pathline, records of X,Y,Z,T,C follow, where
  !                                 (X,Y,Z is the coordinate of a point on the streamline)
  !                                 (T is the time of travel)
  !                                 (C is the concentration)
  !                               At the end of the pathline, "END flag,elem,vertex" is written where
  !                                 (flag) - Termination flag (0=timeout, 1=out of window, 2=hit boundary,
  !                                                            3=aquifer dry,4=stagnation point)
  !                                 (elem) - If at an element, the element type
  !                                 (vertex) - If at an element, the vertex ID  
  !
  !                               NOTE: Since TR0 is a 2-D tracing package, the Z-coordinate remains
  !                               constant throughout all points
  !
  !   TRA  disp-L retard decay    Transport parameters along streamlines. disp-L is the longitudinal
  !                               dispersivity, retard is the retardation factor, and decay is the 
  !                               first-order decay coefficient. Concentrations are computed
  !                               along pathlines using an analytical solution.
  !
  !   WIN  (X1,Y1) (X2,Y2)        Sets the tracing window. Particles that leave the window are terminated.
  !                               Default tuning parameters (see below) are derived from the window size.
  !
  !   TUN  step prox frac small   Tuning parameters for the tracing algorithm
  !                                 (step) - The base step size
  !                                 (prox) - The proximity (in terms of the current step size)
  !                                          to boundary conditions for reducing the step size
  !                                 (frac) - Factor for step size reductions
  !                                 (small) - Smallest allowable step size
  !
  !   TIM  max-time               Specifies the maximum time allowed for particle tracing
  !   
  !   POI  (X0,Y0) Z0 T0 C0 Dir                   Releases a single point at the specified location/
  !                                               layer/direction
  !  
  !   LIN  N (X1,Y1) (X1,Y1) Z0 T0 C0 dir         Releases N particles along the line from (X0,Y0,Z0) to 
  !                                               (X1,Y1,Z0).
  !   
  !   GRI  (X1,Y1) (X2,Y2) N Z0 T0 C0 dir         Releases a grid of particles in the subwindow
  !                                               (X1,Y1) - (X2,Y2) with resolution N
  !  
  !   WL0  id N Z0 T0 C0          Releases N particles in reverse from the well bore of WL0 element #id
  !                               in layer L. The direction of travel is determined by the sign of the 
  !                               well pumping rate (forward for injection wells, backward for pumping
  !                               wells).  


  ! Argument list
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
  ! Locals -- for Directive parsing
  type (DIRECTIVE),dimension(12),parameter :: dirDirectives = &
                                 (/ dirEND, dirDBG, dirPCD, &
                                    dirWIN, dirTRA, dirFIL, dirTUN, dirTIM, &
                                    dirPOI, dirLIN, dirGRI, dirWL0 /)
  ! Locals -- Input values
  character (len=132) :: sRecord   ! Record from input file
  integer (kind=ModAEM_Integer) :: iOpCode         ! Opcode 
  integer (kind=ModAEM_Integer) :: iError          ! Error code
  integer (kind=ModAEM_Integer) :: iStat           ! RunTime error code
  integer (kind=ModAEM_Integer) :: iKeyCode        ! Key code for pause lines
  ! Placeholders for tracing parameters
  real (kind=ModAEM_Real) :: rX0,rY0,rX1,rX2,rY1,rY2,rZ0
  real (kind=ModAEM_Real) :: rStep,rProx,rFrac,rMaxTime,rWellProx
  real (kind=ModAEM_Real) :: rT0,rC0,rDX,rDY
  complex (kind=ModAEM_Real) :: cZ0,cZ1,cZ2,cZWell,cWWell
  real (kind=ModAEM_Real) :: rOrient,rAngle,rDischarge,rRadius
  integer (kind=ModAEM_Integer) :: i,j,iDir,iFlag,iElementType,iElementID,iElementVtx
  integer (kind=ModAEM_Integer) :: iWellID,iNParticles,iFWLIndex,iRes,iResX,iResY
  real (kind=ModAEM_Real) :: rMinX,rMaxX,rMinY,rMaxY,rDelta,rSizeX,rSizeY,rMinWidth
  character (len=132) :: sFile,sMessage
  logical (kind=ModAEM_Integer) :: lFileOpen,lWindowSet,lWellFound,lFlag
  
  ! Clear the status flags
  lFileOpen = .false.
  lWindowSet = .false.
  rTR0MaxTime = TR0_DEFAULT_MAX_TIME
  call IO_MessageText("Entering 2-D trace module TR0",io)

  if (IO_Assert( (associated(aem)), "TR0_Read: No AEM_DOMAIN object",io )) return

  ! The remainder of this routine uses IO_InputRecord to process 
  ! the model input file.
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "TR0_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all xxx_Read routines. 
        if (IO_Assert( .false., "TR0_Read: Unexpected EOF",io )) return
        exit
      case (kOpEND)
        ! END mark was found. Exit the file parser.
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
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        ! If the output file is open, close it now.
        if ( lFileOpen ) then
          close (unit=kIOTraceLU)
          lFileOpen = .false.
        endif
        ! Open the new output file.  If the open fails, write a message.  
        sFile = trim(sFile) // ".tr0"
        open (unit=kIOTraceLU,  file=sFile, action="WRITE", status="REPLACE", iostat=iStat)
        if (IO_Assert( (iStat==0), "TR0_Read: Open failed",io )) return
        lFileOpen = .true.
      case ( kOpWIN )
        !****************************************************************************
        ! Here for the WIN command -- Set the window
        !****************************************************************************
        ! Retrive the window coordinates 
        lWindowSet = .false.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1,cZ2
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        if (IO_Assert( (real(cZ1)<real(cZ2) .and. aimag(cZ1)<aimag(cZ2)), "TR0_Read: Illegal window",io )) return
        rTR0WinX1 = real(cZ1)
        rTR0WinY1 = aimag(cZ1)
        rTR0WinX2 = real(cZ2)
        rTR0WinY2 = aimag(cZ2)
        ! Compute the default tuning parameters
        rTR0Step = TR0_DEFAULT_STEP_FACTOR * abs(cmplx(rX2-rX1,rY2-rY1,ModAEM_Real))
        rTR0Prox = TR0_DEFAULT_PROX
        rTR0Frac = TR0_DEFAULT_FRAC
        rTR0WellProx = TR0_DEFAULT_WELLPROX
        lWindowSet = .true.
        call IO_ErrorText("      >> Default tuning parameters are selected for the new window",io)
      case ( kOpTUN )
        !****************************************************************************
        ! Here for the TUN command -- set the tracing algorithm tuning
        !****************************************************************************
        ! Retrieve the tuning parameters
        read (unit=sRecord,fmt=*,iostat=iStat) rStep,rProx,rFrac,rWellProx
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        ! Note that negative values are ignored.
        if ( rStep > rZERO) rTR0Step = rStep
        if ( rProx > rZERO) rTR0Prox = rProx
        if ( rFrac > rZERO) rTR0Frac = rFrac
        if ( rWellProx > rZERO) rTR0WellProx = rWellProx
        call IO_ErrorText("      >> Specified tuning parameters are in use",io)
      case ( kOpTRA ) 
        !****************************************************************************
        ! Here for the TRA command -- set the transport parameters
        !****************************************************************************
        if (IO_Assert( .false.,"Transport modeling is not yet available",io)) return
      case ( kOpTIM )
        !****************************************************************************
        ! Here for the TIM command -- set the maximum tracing time
        !****************************************************************************
        ! Retrieve the maximum time 
        read (unit=sRecord,fmt=*,iostat=iStat) rMaxTime
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        if (IO_Assert( (rMaxTime>rZERO), "TR0_Read: Illegal maximum time",io )) return
        rTR0MaxTime = rMaxTime
        call IO_ErrorText("      >> Streamline termination time is updated",io)
      case ( kOpPOI )
        !****************************************************************************
        ! Here for the POI command -- release a single point
        !****************************************************************************
        ! Retrieve the maximum time 
        read (unit=sRecord,fmt=*,iostat=iStat) cZ0,rZ0,rT0,rC0,iDir
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        if (IO_Assert( lFileOpen, "TR0_Read: No output file",io )) return
        if (IO_Assert( lWindowSet, "TR0_Read: No window specified",io )) return
        if (IO_Assert( (iDir/=0), "TR0_Read: Illegal direction",io )) return
        ! Fix iDir, then trace away!
        iDir = iDir/abs(iDir)  
        ! For the purpose of the command interpreter, the returned flag information is ignored
        write (unit=sMessage,fmt="("" Starting point: "",2(e13.6,1x),"" Direction: "",i2)") cZ0,iDir
        call IO_MessageText(sMessage,io)
        rX0 = real(cZ0)
        rY0 = aimag(cZ0)
        call TR0_Trace(aem,rX0,rY0,rZ0,rT0,rC0,iDir,iFlag,iElementType,iElementID,iElementVtx,io)
      case ( kOpLIN )
        !****************************************************************************
        ! Here for the LIN command -- release points along a line
        !****************************************************************************
        if (IO_Assert( .false.,"Releases along a line are not yet available",io )) return
      case ( kOpGRI )
        !****************************************************************************
        ! Here for the GRI command -- release points in a grid
        !****************************************************************************
        read (unit=sRecord,fmt=*,iostat=iStat) cZ1,cZ2,iRes,rZ0,rT0,rC0,iDir
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        if (IO_Assert( lFileOpen, "TR0_Read: No output file",io )) return
        if (IO_Assert( lWindowSet, "TR0_Read: No window specified",io )) return
        if (IO_Assert( (iRes>2), "TR0_Read: Bad grid resolution",io )) return
        ! Adjust the grid
        rX1 = real(cZ1)
        rY1 = aimag(cZ1)
        rX2 = real(cZ2)
        rY2 = aimag(cZ2)
        rSizeX = rX2-rX1
        rSizeY = rY2-rY1
        if ( rSizeX >= rSizeY ) then
          iResX = iRes
          rDelta = rSizeX/real(iRes-1,ModAEM_Real)
          iResY = int(rSizeY/rDelta)+1
        else
          iResY = iRes
          rDelta = rSizeY/real(iRes-1,ModAEM_Real)
          iResX = int(rSizeX/rDelta)+1
        endif
        ! Make the real window span the center of the specified window
        rMinX = rHALF * ( rX1+rX2 - (iResX-1)*rDelta )
        rMaxX = rMinX + (iResX-1)*rDelta
        rMinY = rHALF * ( rY1+rY2 - (iResY-1)*rDelta )
        rMaxY = rMinY + (iResY-1)*rDelta
        ! Fix iDir, then trace away for all the points in a grid
        iDir = iDir/abs(iDir)  
        write ( unit=kIOTraceLU, fmt="('GRID ',2i10)" ) iResX,iResY
        do j=1,iResY
          do i=1,iResX
            write ( unit=kIOTraceLU, fmt="('CELL ',2i10)" ) j,i
            rX0 = (i-1)*rDelta + rMinX
            rY0 = rMaxY - (j-1)*rDelta           
            write ( unit=sMessage, &
                    fmt="("" Row: "",i5,"" Col: "",i5,"" Starting point: "",2(e13.6,1x),"" Direction: "",i2)" &
                  ) j,i,rX0,rY0,iDir
            call IO_MessageText(sMessage,io)
            call TR0_Trace(aem,rX0,rY0,rZ0,rT0,rC0,iDir,iFlag,iElementType,iElementID,iElementVtx,io)
          end do
        end do
      case ( kOpWL0 )
        !****************************************************************************
        ! Here for the WL0 command -- release points from a well
        !****************************************************************************
        ! This routine requires that the location of the well element be found (using ifWL0FindWell)
        ! and then using the discharge vector at the well (less the discharge due to the well) to
        ! orient the release points.  This routine does not check to see if the well has pumped the 
        ! aquifer dry; TR0_Trace will report that in the END record.
        ! Retrieve the maximum time 
        !
        ! Modified 6/25/2001 by VAK. Added rMinWidth value for Idaho delineations. If the
        ! computed "ultimate" capture zone width is less than rMinWidth, then ModAEM traces
        ! only one particle back from the wellscreen, in the direction of flow.
        read (unit=sRecord,fmt=*,iostat=iStat) iWellID,iNParticles,rZ0,rT0,rC0,rMinWidth
        if (IO_Assert( (iStat==0), "TR0_Read: I/O Error",io )) return
        if (IO_Assert( lFileOpen, "TR0_Read: No output file",io )) return
        if (IO_Assert( lWindowSet, "TR0_Read: No window specified",io )) return
        call WL0_FindWell(aem%wl0,iWellID,cZWell,rDischarge,rRadius,iFWLIndex,lWellFound,io)
        if (IO_Assert( lWellFound, "TR0_Read: Specified well not found",io )) return
        ! Forward trace from injection wells; backward trace from pumping wells
        if ( rDischarge < 0 ) then
          iDir = 1  
          call IO_ErrorText("      >> Forward tracing is selected from injection well",io)
        else
          iDir = -1
          call IO_ErrorText("      >> Reverse tracing is selected from pumping well",io)
        endif
        ! Now, compute the orientation of the discharge vector, less the well
        cWWell = cAEM_DischargeAtWell(aem,cZWell,iFWLIndex,io)
        if ( cWWell == cZERO ) then 
          rOrient = rZERO
        else
          rOrient = atan2(aimag(cWWell),real(cWWell))
        endif
        ! Distribute particles about the well bore, distributed evenly about the orientation vector
        if ( rDischarge/abs(cWWell) <= rMinWidth ) then
          write ( unit=sMessage, &
                  fmt="('Found small well with ID ',i5,' -- tracing one particle')" &
                ) iWellID
          call IO_MessageText(sMessage,io)
          iNParticles = 1
        endif
        write (unit=sMessage,fmt="(""Releasing "",i5,"" particles from well ID "",i5)") iNParticles,iWellID
        write ( unit=kIOTraceLU, fmt="('WL0   ',i10)" ) iWellID
        call IO_MessageText(sMessage,io)
        do i=1,iNParticles
          rAngle = rOrient + (float(i-1)+rHALF) * 2*rPI / float(iNParticles)
          rX1 = real(cZWell) + (rRadius * TR0_WL0_DISTANCE_FACTOR) * cos(rAngle)
          rY1 = aimag(cZWell) + (rRadius * TR0_WL0_DISTANCE_FACTOR) * sin(rAngle)
          write (unit=sMessage,fmt="("" Starting point: "",2(e13.6,1x),"" Direction: "",i2)") rX1,rY1,iDir
          call IO_MessageText(sMessage,io)
          call TR0_Trace(aem,rX1,rY1,rZ0,rT0,rC0,iDir,iFlag,iElementType,iElementID,iElementVtx,io)
        end do
      case default
    end select
  end do
  ! Return the most recent error code found.
  return
end subroutine TR0_Read

subroutine TR0_SetWindow(rWinX1,rWinY1,rWinX2,rWinY2,io)
  ! Sets the tracing window
  real (kind=ModAEM_Real),intent(in) :: rWinX1
  real (kind=ModAEM_Real),intent(in) :: rWinY1
  real (kind=ModAEM_Real),intent(in) :: rWinX2
  real (kind=ModAEM_Real),intent(in) :: rWinY2
  type ( IO_Status ),pointer :: io

  if (IO_Assert( ( (rWinX2 > rWinX1) .and. (rWinY2 > rWinY1) ), &
                  "TR0_SetWindow: Illegal window specified",io ) ) return

  rTR0WinX1 = rWinX1
  rTR0WinY1 = rWinY1
  rTR0WinX2 = rWinX2
  rTR0WinY2 = rWinY2

  return
end subroutine TR0_SetWindow

subroutine TR0_SetMaxTime(rMaxTime,io)
  ! Sets the maximum time for tracing
  real (kind=ModAEM_Real),intent(in) :: rMaxTime
  type ( IO_Status ),pointer :: io

  if (IO_Assert( (rMaxTime > rZERO),"TR0_MaxTime: Illegal time specified",io )) return
  rTR0MaxTime = rMaxTime

  return
end subroutine TR0_SetMaxTime

subroutine TR0_SetTuning(rStep,rProx,rFrac,rWellProx,io)
  ! Sets tuning parameters for tracing
  real (kind=ModAEM_Real),intent(in),optional :: rStep
  real (kind=ModAEM_Real),intent(in),optional :: rProx
  real (kind=ModAEM_Real),intent(in),optional :: rFrac
  real (kind=ModAEM_Real),intent(in),optional :: rWellProx
  type ( IO_Status ),pointer :: io

  if ( present(rStep) ) then
    rTR0Step = rStep
  else
    rTR0Step = TR0_DEFAULT_STEP_FACTOR * abs(cmplx(rTR0WinX2-rTR0WinX1,rTR0WinY2-rTR0WinY1,ModAEM_Real))
  endif

  if ( present(rProx) ) then
    rTR0Prox = rProx
  else
    rTR0Prox = TR0_DEFAULT_PROX
  endif

  if ( present(rFrac) ) then
    rTR0Frac = rFrac
  else
    rTR0Frac = TR0_DEFAULT_FRAC
  endif

  if ( present(rWellProx) ) then
    rTR0WellProx = rWellProx
  else
    rTR0WellProx = TR0_DEFAULT_WELLPROX
  endif

  return
end subroutine TR0_SetTuning

subroutine TR0_SetTransport(rDispL,rRetard,rDecay,io)
  ! Sets parameters for transport along pathlines
  real (kind=ModAEM_Real),intent(in) :: rDispL
  real (kind=ModAEM_Real),intent(in) :: rRetard
  real (kind=ModAEM_Real),intent(in) :: rDecay
  type ( IO_Status ),pointer :: io

  rTR0DispL = rDispL
  rTR0Retard = rRetard
  rTR0Decay = rDecay

  return
end subroutine TR0_SetTransport

subroutine TR0_Trace(aem,rX0,rY0,rZ0,rT0,rC0,iDir,iFlag,iElementType,iElementID,iElementVtx,io)
  ! Computes a single pathline starting at cZO, using the predictor-corrector method
  type ( AEM_DOMAIN ),pointer :: aem
  real (kind=ModAEM_Real),intent(inout) :: rX0,rY0,rZ0,rT0,rC0
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  integer (kind=ModAEM_Integer),intent(out) :: iFlag,iElementType,iElementID,iElementVtx
  type ( IO_Status ),pointer :: io
   ! Locals
  integer (kind=ModAEM_Integer) :: iPrev,iNSteps,i1,i2,iE,iV
  real (kind=ModAEM_Real) :: rBaseTol,rStagTol,rDZ20,rStep,rDeltaTime,rTime,rConc
  complex (kind=ModAEM_Real) :: cZO,cQO,cZN,cQN,cQ,cZInt
  complex (kind=ModAEM_Real),dimension(20) :: cPrevZ
  logical (kind=ModAEM_Integer) :: lTimeout,lFlag

  ! Tolerance used to determine whether an element has been reached
  rBaseTol = 1.1_ModAEM_Real * rTR0Step
  cZN = cmplx(rX0,rY0,ModAEM_Real)
  cQN = cAEM_Velocity(aem,CZN,io)
  rTime = rT0
  rConc = rC0
  iPrev = 1
  iNSteps = 0
  rStagTol = TR0_STAGNATION_TOLERANCE * rTR0Step
  lTimeout = .false.

  ! Write the start record 
  write (unit=kIOTraceLU,fmt="(""START "",5(e13.6,1x),2i10)") rX0,rY0,rZ0,rT0,rC0,iDir

  do 
    ! Write out the results of the previous step
    write (unit=kIOTraceLU,fmt="(""      "",5(e13.6,1x))") cZN,rZ0,rTime,rConc

    ! Prepare to take a step along the streamline
    cZO = cZN
    cQO = cQN
    cPrevZ(iPrev) = cZO
    iPrev = iPrev+1
    if ( iPrev > 20 ) iPrev = mod(iPrev,20)

    ! Abort if we are out of the window
    if ( real(cZN)  < rTR0WinX1 .or. real(cZN)  > rTR0WinX2 .or. &
         aimag(cZN) < rTR0WinY1 .or. aimag(cZN) > rTR0WinY2 ) then
      iFlag = kTR0LeftWindow
      iElementType = 0
      iElementID = 0
      iElementVtx = 0
      exit
    endif

    ! Abort if we have not moved in the past 20 steps
    if ( iNSteps > 20 ) then
      i1 = iPrev
      i2 = iPrev-1
      if ( i2 < 1 ) i2 = 20
      rDZ20 = abs(cPrevZ(i1) - cPrevZ(i2))
      if ( rDZ20 < rStagTol ) then
        iFlag = kTR0Stagnation
        iElementType = 0
        iElementID = 0
        iElementVtx = 0
        exit  
      endif
    endif
    
    ! Abort if the maximum travel time is achieved
    if ( lTimeout) then
      iFlag = kTR0Timeout
      iElementType = 0
      iElementID = 0
      iElementVtx = 0
      exit
    endif

    ! Abort if the head is below the aquifer bottom
    if ( real(cAEM_Potential(aem,cZO,io)) <= rZERO ) then
      iFlag = kTR0AquiferDry
      iElementType = 0
      iElementID = 0
      iElementVtx = 0
      exit
    endif
    
    ! Abort if at an element -- check all element types 
    ! Check the WL0 elements
    if ( lWL0_CheckProximity(aem%wl0,cZO,rTR0Step*TR0_WL0_TERMINATION_FACTOR,iDir,iE,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_WL0
      iElementID = iE
      iElementVtx = 0
      exit
    endif
    ! Check the LS0 elements
    if ( lLS0_CheckProximity(aem%ls0,cZO,rTR0Step*TR0_LS0_TERMINATION_FACTOR,iDir,iE,iV,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_LS0
      iElementID = iE
      iElementVtx = iV
      exit
    endif
    ! Check the LS1 elements
    if ( lLS1_CheckProximity(aem%ls1,cZO,rTR0Step*TR0_LS1_TERMINATION_FACTOR,iDir,iE,iV,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_LS1
      iElementID = iE
      iElementVtx = iV
      exit
    endif
    ! Check the HB0 elements.  NOTE: In THEORY, the no-flow boundary cannot terminate a streamline, since
    ! it removes no water.  In practice, it is possible for a streamline to "get lost" at a no-flow boundary 
    ! due to inaccuracies in the solution.  Check at a smaller tolerance near no-flow boundaries.
    if ( lHB0_CheckProximity(aem%hb0,cZO,rTR0Step*TR0_HB0_TERMINATION_FACTOR,iDir,iE,iV,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_HB0
      iElementID = iE
      iElementVtx = iV
      exit
    endif
    ! Check the WL1 elements
    if ( lWL1_CheckProximity(aem%wl1,cZO,rTR0Step*TR0_WL1_TERMINATION_FACTOR,iDir,iE,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_WL1
      iElementID = iE
      iElementVtx = iV
      exit
    endif
    
    ! PHEW. Now, proceed with the tracing...
    ! Adjust step size if necessary...  Note: all the xxxCheckProximity calls are made without regard to iDir
    if ( lWL0_CheckProximity(aem%wl0,cZO,rTR0WellProx*rTR0Step,0,iE,io) .or. &
         lLS0_CheckProximity(aem%ls0,cZO,rTR0Prox*rTR0Step,0,iE,iV,io) .or. &
         lLS1_CheckProximity(aem%ls1,cZO,rTR0Prox*rTR0Step,0,iE,iV,io) .or. &
         lHB0_CheckProximity(aem%hb0,cZO,rTR0Prox*rTR0Step,0,iE,iV,io) .or. &
         lWL1_CheckProximity(aem%wl1,cZO,rTR0Prox*rTR0Step,0,iE,io) ) then
       rStep = rTR0Frac * iDir * rTR0Step
    else
       rStep = iDir * rTR0Step
    endif

    ! Determine the complex coordinate of the next point on the streamline
    cZN = cZO + rStep*cQO/abs(cQO)
    ! Check to see if we have entered a well.
    if ( lWL0_CheckProximity(aem%wl0,cZN,rStep,iDir,iE,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_WL0
      iElementID = iE
      iElementVtx = 0
      exit
    endif

    ! Compute the complex discharge Qx-iQy at this point.
    cQN = cAEM_Velocity(aem,cZN,io)
    ! Estimate the average value of the complex discharge on the interval.
    cQ = .5*(cQO+cQN)
    ! Obtain the second approximation of the point on the streamline.
    cZN = cZO + rStep*cQ/abs(cQ)

    ! Check again to see if we have entered a well, to correct the location
    if ( lWL0_CheckProximity(aem%wl0,cZN,rStep,iDir,iE,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_WL0
      iElementID = iE
      iElementVtx = 0
      exit
    endif

    ! Check to see if we have crossed a no-flow boundary
    if ( lHB0_CheckCrossing(aem%hb0,cZO,cZN,rTR0Prox*rTR0Step,iE,iV,cZInt,io) ) then
      iFlag = kTR0HitElement
      iElementType = ELEM_HB0
      iElementID = iE
      iElementVtx = 0
      exit
    endif

    ! Compute the time-of-travel and report the results for this point
    ! NOTE: The transport solution is not yet implemented.
    rDeltaTime = abs(cZN-cZO)/abs(cQ)
    if ( rTime+rDeltaTime > rTR0MaxTime ) then
      if ( abs(cZN-cZO) > rZERO ) then
        cZN = cZO + (cZN-cZO)*(rTR0MaxTime-rTime)/rDeltaTime
        rTime = rTR0MaxTime 
        lTimeout = .true.
      end if
    else
      rTime = rTime + rDeltaTime
    end if

    ! Repeat the procedure for the next point.
    iNSteps = iNSteps + 1

  end do

  ! Write the end record 
  write (unit=kIOTraceLU,fmt="(""END   "",5(e13.6,1x),4i10)") cZN,rZ0,rTime,rConc, &
                                                              iFlag,iElementType,iElementID,iElementVtx

  return
end subroutine TR0_Trace

end module a_tr0

