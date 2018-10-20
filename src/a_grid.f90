module a_grid

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

use u_constants
use u_io
use m_aem
use m_aqu

implicit none

public

  integer (kind=ModAEM_Integer),parameter :: GRID_HEAD=1
  integer (kind=ModAEM_Integer),parameter :: GRID_POTENTIAL=2
  integer (kind=ModAEM_Integer),parameter :: GRID_STREAMFUNCTION=3
  integer (kind=ModAEM_Integer),parameter :: GRID_QX=4
  integer (kind=ModAEM_Integer),parameter :: GRID_QY=5
  integer (kind=ModAEM_Integer),parameter :: GRID_VX=6
  integer (kind=ModAEM_Integer),parameter :: GRID_VY=7

  integer (kind=ModAEM_Integer),parameter :: FILE_SURFER=1
  integer (kind=ModAEM_Integer),parameter :: FILE_MATLAB=2
  
  real (kind=ModAEM_Real),parameter :: rSMALL_GRID = 1.0e-1

contains

subroutine GRI_Read(aem,io)
  ! This routine handles grid generation and I/O commands

  ! Argument list
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! Locals -- for Directive parsing
  type ( DIRECTIVE ),dimension(13),parameter :: dirDirectives = &
                   (/ dirEND, dirDBG, dirPCD, dirOPT, &
                      dirWIN, dirDIM, dirHEA, dirPOT, dirPSI, dirDIS, dirVEL, dirQ_X, dirQ_Y /)
  ! Locals -- Input values
  character (len=132) :: sRecord   ! Record from input file
  integer (kind=ModAEM_Integer) :: iOpCode         ! Opcode 
  integer (kind=ModAEM_Integer) :: iError          ! Error code
  integer (kind=ModAEM_Integer) :: iStat           ! RunTime error code
  integer (kind=ModAEM_Integer) :: iKeyCode        ! Key code for pause lines
  ! Placeholders for grid parameters
  complex (kind=ModAEM_Real) :: cLL, cUR 
  integer (kind=ModAEM_Integer) :: iRes,iOption,ic
  complex (kind=ModAEM_Real) :: cZ
  logical (kind=ModAEM_Integer) :: lFlag
  character (len=132) :: sFile,sOption,sFixOption
  
  ! The remainder of this routine uses ifIO_InputRecord to process 
  ! the model input file.
  iOption = FILE_SURFER
  call IO_MessageText("Entering GRI module",io)
  if (IO_Assert( (associated(aem)), "GRI_Read: No AEM_DOMAIN object",io )) return
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "GRI_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModGRI "ifXXXRead" routines. 
        ! Report the condition, but proceed as if EOD was found.
        if (IO_Assert( .false., "GRI_Read: Unexpected EOF",io )) return
      case (kOpEND)
        ! END mark was found. Exit the file parser.
        call IO_MessageText("Leaving GRI module",io)
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
      case (kOpDIM)
        !****************************************************************************
        ! Here for the DIM command -- dimension a grid
        !****************************************************************************
        ! Retrive the grid dimension 
        read (unit=sRecord,fmt=*,iostat=iStat) iRes
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        if (IO_Assert( (iRes>2), "GRI_Read: Illegal resolution",io )) return
      case ( kOpWIN )
        !****************************************************************************
        ! Here for the WIN command -- Set the window
        !****************************************************************************
        ! Retrive the window coordinates 
        read (unit=sRecord,fmt=*,iostat=iStat) cLL,cUR
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        if (IO_Assert( (real(cLL)<real(cUR)), "GRI_Read: Illegal window",io )) return
        if (IO_Assert( (aimag(cLL)<aimag(cUR)), "GRI_Read: Illegal window",io )) return
      case ( kOpHEA )
        !****************************************************************************
        ! Here for the HEA command -- generate and write a grid of heads
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        call GRI_MakeGrid(aem,GRID_HEAD,cLL,cUR,iRes,sFile,iOption,io)
      case ( kOpPOT )
        !****************************************************************************
        ! Here for the POT command -- generate and write a grid of potentials
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        call GRI_MakeGrid(aem,GRID_POTENTIAL,cLL,cUR,iRes,sFile,iOption,io)
      case ( kOpPSI )
        !****************************************************************************
        ! Here for the PSI command -- generate and write a grid of streamfunctions
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        call GRI_MakeGrid(aem,GRID_STREAMFUNCTION,cLL,cUR,iRes,sFile,iOption,io)
      case ( kOpQ_X )
        !****************************************************************************
        ! Here for the QX command -- generate and write a grid of discharges
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        call GRI_MakeGrid(aem,GRID_QX,cLL,cUR,iRes,sFile,iOption,io)
      case ( kOpQ_Y )
        !****************************************************************************
        ! Here for the QY command -- generate and write a grid of discharges
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sFile
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        call GRI_MakeGrid(aem,GRID_QY,cLL,cUR,iRes,sFile,iOption,io)
      case ( kOpOPT )
        !****************************************************************************
        ! Here for the QY command -- generate and write a grid of discharges
        !****************************************************************************
        ! Retrive the filename
        read (unit=sRecord,fmt=*,iostat=iStat) sOption
        if (IO_Assert( (iStat==0), "GRI_Read: I/O Error",io)) return
        do ic=1,len(sOption)
          if ( ichar(sOption(ic:ic))>96 .and. ichar(sOption(ic:ic))<122 ) then
            sFixOption(ic:ic) = char(ichar(sOption(ic:ic))-32)
          else
            sFixOption(ic:ic) = sOption(ic:ic)
          endif
        end do
        if ( trim(sFixOption) == "SURFER" ) then
           iOption = FILE_SURFER
           call IO_MessageText("  Selecting SURFER output",io)
        elseif ( trim(sFixOption) == "MATLAB" ) then
           iOption = FILE_MATLAB
           call IO_MessageText("  Selecting Matlab output",io)
        else
          if ( IO_Assert( (.false.),"GRI_Read: Illegal file format option",io ) ) return
        endif
    end select
  end do
  ! Return the most recent error code found.
  return
end subroutine GRI_Read

subroutine GRI_MakeGrid(aem,iFunc,cLL,cUR,iRes,sFile,iOption,io)
  !! subroutine MakeGrid
  !!
  !! Makes a grid of function 'func' for the AEM_DOMAIN 'aem', writing the
  !! results to sFile
  !!
  !! Arguments
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              AEM_DOMAIN to be used
  !!    (in)    integer :: iFunc
  !!              Function to be gridded (see constants defined above)
  !!                GRID_HEAD            Head (makes '.head.grd' file)
  !!                GRID_POTENTIAL       Potential (makes '.pot.grd' file)
  !!                GRID_STREAMFUNCTION  Streamfunction (makes '.psi.grd' file)
  !!                GRID_GRID_QX              X-discharge (makes '.GRID_QX.grd' file)
  !!                GRID_GRID_QY              Y-discharge (makes '.GRID_QY.grd' file)
  !!                GRID_GRID_VX              X-Velocity (makes 'GRID_VX.grd' file)
  !!                GRID_GRID_VY              Y-Velocity (makes 'GRID_VY.grd' file)
  !!    (in)    complex :: cLL
  !!              Lower-left corner of the grid
  !!    (in)    complex :: cUR
  !!              Upper-rigth corner of the grid
  !!    (in)    integer :: iRes
  !!              Number of points along the longer axis
  !!    (in)    character :: sFile
  !!              File to be written, without extensions (they will be added
  !!              based on the function gridded; see above)
  !!    (in)    integer :: iOption
  !!              Option for file format
  !!
  !! Note: Uses the kIOGridLU LU for output
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  integer (kind=ModAEM_Integer),intent(in) :: iFunc
  complex (kind=ModAEM_Real),intent(in) :: cLL
  complex (kind=ModAEM_Real),intent(in) :: cUR
  integer (kind=ModAEM_Integer),intent(in) :: iRes
  character (len=*),intent(in) :: sFile
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_Status ),pointer :: io
  ! [ LOCALS ]
  character (len=255) :: sFileName
  character (len=4) :: sExt
  integer (kind=ModAEM_Integer) :: iStat,iResX,iResY,i,j
  real (kind=ModAEM_Real) :: rSizeX,rSizeY,rDelta,rX,rY
  real (kind=ModAEM_Real) :: rMinX,rMaxX,rMinY,rMaxY
  real (kind=ModAEM_Real),dimension(:,:),allocatable :: rGrid

  if (IO_Assert( (associated(aem)), 'MakeGrid: No AEM_DOMAIN object',io )) return
  if (IO_Assert( (abs(cUR-cLL)>rSMALL_GRID), 'MakeGrid: Bad window',io )) return
  if (IO_Assert( (iRes>2), 'MakeGrid: Bad resolution',io )) return

  rSizeX = real(cUR-cLL)
  rSizeY = aimag(cUR-cLL)
  if ( rSizeX >= rSizeY ) then
    iResX = iRes
    rDelta = rSizeX/real(iRes-1,ModAEM_Real)
    iResY = int(rSizeY/rDelta)+1
  else
    iResY = iRes
    rDelta = rSizeY/real(iRes-1,ModAEM_Real)
    iResX = int(rSizeX/rDelta)+1
  endif
  allocate ( rGrid(iResY,iResX), stat=iStat )
  if (IO_Assert( (iStat==0), 'MakeGrid: Allocation failed',io )) return

  ! Make the real window span the center of the specified window
  rMinX = rHALF * ( real(cLL+cUR) - (iResX-1)*rDelta )
  rMaxX = rMinX + (iResX-1)*rDelta
  rMinY = rHALF * ( aimag(cLL+cUR) - (iResY-1)*rDelta )
  rMaxY = rMinY + (iResY-1)*rDelta

  ! Generate the grid
  select case ( iFunc )
    case ( GRID_HEAD )
      call IO_MessageText("  Generating a grid of heads",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = rAEM_Head(aem,cmplx(rX,rY,ModAEM_Real),io)
        end do
      end do
    case ( GRID_POTENTIAL )
      call IO_MessageText("  Generating a grid of potentials",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = real(cAEM_Potential(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
    case ( GRID_STREAMFUNCTION )
      call IO_MessageText("  Generating a grid of the streamfunction",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = aimag(cAEM_Potential(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
    case ( GRID_QX )
      call IO_MessageText("  Generating a grid of Qx",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = real(cAEM_Discharge(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
    case ( GRID_QY )
      call IO_MessageText("  Generating a grid of Qy",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = aimag(cAEM_Discharge(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
    case ( GRID_VX )
      call IO_MessageText("  Generating a grid of Vx",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = real(cAEM_Velocity(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
    case ( GRID_VY )
      call IO_MessageText("  Generating a grid of Vy",io)
      do j=1,iResY
        rY = rMinY + (j-1)*rDelta
        do i=1,iResX
          rX = rMinX + (i-1)*rDelta
          rGrid(j,i) = aimag(cAEM_Velocity(aem,cmplx(rX,rY,ModAEM_Real),io))
        end do
      end do
  end select

  ! Write the grid file
  if ( iOption == FILE_SURFER ) then
    sExt = ".grd"
  elseif ( iOption == FILE_MATLAB ) then
    sExt = ".m"
  endif
  select case ( iFunc )
    case ( GRID_HEAD )
      sFileName =  trim(sFile) // "_head" // sExt
    case ( GRID_POTENTIAL )
      sFileName =  trim(sFile) // "_potential" // sExt
    case ( GRID_STREAMFUNCTION )
      sFileName =  trim(sFile) // "_psi" // sExt
    case ( GRID_QX )
      sFileName =  trim(sFile) // "_qx" // sExt
    case ( GRID_QY )
      sFileName =  trim(sFile) // "_qy" // sExt
    case ( GRID_VX )
      sFileName =  trim(sFile) // "_vx" // sExt
    case ( GRID_VY )
      sFileName =  trim(sFile) // "_vy" // sExt
  end select

  open ( unit=kIOGridLU, file=trim(sFileName), iostat=iStat)
  if (IO_Assert( (iStat==0), 'MakeGrid: Could not open file',io )) return
  call IO_MessageText("  Writing grid to " // trim(sFileName) ,io)

  if ( iOption == FILE_SURFER ) then
    write ( unit=kIOGridLU, &
            fmt="('DSAA')" &
          )
    write ( unit=kIOGridLU, &
            fmt="(2(1x,i5))" &
          ) iResX,iResY
    write ( unit=kIOGridLU, &
            fmt="(2(1x,g12.5))" &
          ) rMinX,rMaxX
    write ( unit=kIOGridLU, &
            fmt="(2(1x,g12.5))" &
          ) rMinY,rMaxY
    write ( unit=kIOGridLU, &
            fmt="(2(1x,g12.5))" &
          ) minval(rGrid),maxval(rGrid)
    do j=1,iResY
      write ( unit=kIOGridLU, &
              fmt=* &
            ) rGrid(j,:)
    end do
    ! That's it
    close ( unit=kIOGridLU )
  
  elseif ( iOption == FILE_MATLAB ) then
    write ( unit=kIOGridLU, &
            fmt="('x=[')" &
          ) 
    do i=1,iResX
      write ( unit=kIOGridLU, &
              fmt=* &
            ) rMinX + (i-1)*rDelta
    end do
    write ( unit=kIOGridLU, &
            fmt="(']')" &
          ) 
    write ( unit=kIOGridLU, &
            fmt="('y=[')" &
          ) 
    do j=1,iResY
      write ( unit=kIOGridLU, &
              fmt=* &
            ) rMinY + (j-1)*rDelta
    end do
    write ( unit=kIOGridLU, &
            fmt="(']')" &
          ) 


    write ( unit=kIOGridLU, &
            fmt="('f=[')" &
          ) 
    do j=1,iResY
      write ( unit=kIOGridLU, &
              fmt=* &
            ) rGrid(:,j)
      write ( unit=kIOGridLU, &
              fmt="(';')" &
             ) 
    end do
    write ( unit=kIOGridLU, &
            fmt="(']')" &
          ) 
    ! That's it
    close ( unit=kIOGridLU )

  endif

  return
end subroutine GRI_MakeGrid

end module a_grid

