module u_grid

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

  ! module u_grdd (GRD)
  ! Data structures and convenient methods for regularly grddded data 

use u_constants
use u_io

implicit none

public

  type :: GRD_GRID
    !! type GRD_GRID
    !!
    !! PUBLIC type that holds the information for regularly-grddded data
    !!
    !! Members:
    !!   real :: rMinX
    !!     Minimum X value
    !!   real :: rMaxX
    !!     Maximum X value
    !!   real :: rMinY
    !!     Minimum Y value
    !!   real :: rMaxY
    !!     Maximum Y value
    !!   real :: rMinZ
    !!     Minimum Z value
    !!   real :: rMaxZ
    !!     Maximum Z value
    !!   integer :: iResX
    !!     Number of points along the X axis
    !!   integer :: iResY
    !!     Number of points along the Y axis
    !!   real :: rDelta
    !!     The spacing between grdd points
    !!   real :: rValues(:,:)
    !!     The values of points in the grd
    !!
    real (kind=ModAEM_Real) :: rMinX
    real (kind=ModAEM_Real) :: rMaxX
    real (kind=ModAEM_Real) :: rMinY
    real (kind=ModAEM_Real) :: rMaxY
    real (kind=ModAEM_Real) :: rMinZ
    real (kind=ModAEM_Real) :: rMaxZ
    integer (kind=ModAEM_Integer) :: iResX
    integer (kind=ModAEM_Integer) :: iResY
    real (kind=ModAEM_Real) :: rDelta
    real (kind=ModAEM_Real),dimension(:,:),pointer :: rValues
  end type GRD_GRID

  ! Define the smalles legal grid window here
  real (kind=ModAEM_Real),private,parameter :: rSMALL_GRID = 1.0e-3_ModAEM_Real

contains

function GRD_Create(cLL,cUR,iRes,io) result(grd)
  !! function GRD_Create
  !!
  !! Creates a new GRD_GRID object
  !!
  !! Calling Sequence:
  !!    grd => GRD_Create(cLL,cUR,iRes)
  !!
  !! Arguments:
  !!    (in)    complex :: cLL
  !!              Lower-left corner of the grdd
  !!    (in)    complex :: cUR
  !!              Upper-rigth corner of the grdd
  !!    (in)    integer :: iRes
  !!              Number of points along the longer axis
  !!
  ! [ ARGUMENTS ]
  complex (kind=ModAEM_Real),intent(in) :: cLL
  complex (kind=ModAEM_Real),intent(in) :: cUR
  integer (kind=ModAEM_Integer),intent(in) :: iRes
  type ( IO_Status ),pointer :: io

  ! [ RETURN VALUE ]
  type ( GRD_GRID ),pointer :: grd

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  real (kind=ModAEM_Real) :: rSizeX,rSizeY

  if (IO_Assert( (abs(cUR-cLL)>rSMALL_GRID), 'GRD_Create: Bad window',io )) return
  if (IO_Assert( (iRes>2), 'GRD_Create: Bad resolution',io )) return

  allocate ( grd,stat=iStat )
  if (IO_Assert( (iStat==0), 'GRD_Create: Allocation failed',io )) return

  ! Compute the actual data sizes for the grdd and allocate space
  rSizeX = real(cUR-cLL)
  rSizeY = aimag(cUR-cLL)
  if ( rSizeX >= rSizeY ) then
    grd%iResX = iRes
    grd%rDelta = rSizeX/real(iRes-1,ModAEM_Real)
    grd%iResY = int(rSizeY/grd%rDelta)+1
  else
    grd%iResY = iRes
    grd%rDelta = rSizeY/real(iRes-1,ModAEM_Real)
    grd%iResX = int(rSizeX/grd%rDelta)+1
  endif
  allocate ( grd%rValues(grd%iResY,grd%iResX), stat=iStat )
  if (IO_Assert( (iStat==0), 'MakeGrid: Allocation failed',io )) return

  ! Make the real window span the center of the specified window
  grd%rMinX = rHALF * ( real(cLL+cUR) - (grd%iResX-1)*grd%rDelta )
  grd%rMaxX = grd%rMinX + (grd%iResX-1)*grd%rDelta
  grd%rMinY = rHALF * ( aimag(cLL+cUR) - (grd%iResY-1)*grd%rDelta )
  grd%rMaxY = grd%rMinY + (grd%iResY-1)*grd%rDelta
  grd%rValues = rZERO

  return
end function GRD_Create

function GRD_CreateFromGRD(sFileName,io) result(grd)
  !! function GRD_CreateFromGRD
  !!
  !! Creates a new GRD_GRID object from a SURFER(tm)-compatible ASCII .GRD file
  !!
  !! Calling Sequence:
  !!    grd => GRD_CreateFromGRD(sFileName)
  !!
  !! Arguments:
  !!    (in)    character(len=*) :: sFileName
  !!              Name of the SURFER(tm) GRD file
  !!
  ! [ ARGUMENTS ]
  character (len=*),intent(in) :: sFileName
  type ( IO_Status ),pointer :: io

  ! [ RETURN VALUE ]
  type ( GRD_GRID ),pointer :: grd
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  character (len=4) :: sDSAA
  integer (kind=ModAEM_Integer) :: iNX,iNY
  real (kind=ModAEM_Real) :: rDX,rDY,rMinX,rMaxX,rMinY,rMaxY,rMinZ,rMaxZ

  if ( io%lDebug ) then
    if (IO_Assert( (.not. associated(grd)), "GRD_CreateFromGRD: GRD_GRID object is already in use",io )) return
  endif
 
  ! Open up the GRD file
  open ( unit=kIOGridLU, file=trim(sFileName), status="OLD", iostat=iStat )
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: Could not open " // sFileName // " for input",io )) return
 
  ! Read the grid header
  read ( unit=kIOGridLU, fmt="(a4)", iostat=iStat ) sDSAA
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: I/O Error",io )) return
  if (IO_Assert( (sDSAA=="DSAA"), "GRD_CreateFromGRD: Illegal header record",io )) return
  ! Dimensions and limits
  read ( unit=kIOGridLU, fmt=*, iostat=iStat ) iNX,iNY
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: I/O Error",io )) return
  read ( unit=kIOGridLU, fmt=*, iostat=iStat ) rMinX,rMaxX
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: I/O Error",io )) return
  read ( unit=kIOGridLU, fmt=*, iostat=iStat ) rMinY,rMaxY
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: I/O Error",io )) return
  read ( unit=kIOGridLU, fmt=*, iostat=iStat ) rMinZ,rMaxZ
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: I/O Error",io )) return
  ! Made it this far... build the object and populate it
  allocate ( grd,stat=iStat )
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: Allocation failed",io )) return
  grd%iResX = iNX
  grd%iResY = iNY
  grd%rMinX = rMinX
  grd%rMaxX = rMaxX
  grd%rMinY = rMinY
  grd%rMaxY = rMaxY
  grd%rMinZ = rMinZ
  grd%rMaxZ = rMaxZ
  allocate ( grd%rValues(iNY,iNX), stat=iStat )
  if (IO_Assert( (iStat==0), "GRD_CreateFromGRD: Allocation failed",io )) return
  read ( unit=kIOGridLU, fmt=*, iostat=iStat ) grd%rValues
  
  return
end function GRD_CreateFromGrd

subroutine GRD_WriteToGRD(grd,sFileName,io)
  !! function GRD_CreateFromGRD
  !!
  !! Creates a new GRD_GRID object from a SURFER(tm)-compatible ASCII .GRD file
  !!
  !! Calling Sequence:
  !!    grd => GRD_CreateFromGRD(grd,sFileName)
  !!
  !! Arguments:
  !!    (in)    type ( GRD_GRID ),pointer :: grd
  !!              GRD_GRID object to be written
  !!    (in)    character(len=*) :: sFileName
  !!              Name of the SURFER(tm) GRD file
  !!
  ! [ ARGUMENTS ]
  type ( GRD_GRID ),pointer :: grd
  character (len=*),intent(in) :: sFileName
  type ( IO_Status ),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: j,iStat

  open ( unit=kIOGridLU, file=trim(sFileName), status="NEW", iostat=iStat)
  if (IO_Assert( (iStat==0), "MakeGrid: Could not open " // sFileName // " for output",io )) return

  call IO_MessageText("  Writing grid to " // trim(sFileName), io)
  write ( unit=kIOGridLU, &
          fmt="('DSAA')" &
        )
  write ( unit=kIOGridLU, &
          fmt="(2(1x,i5))" &
        ) grd%iResX,grd%iResY
  write ( unit=kIOGridLU, &
          fmt="(2(1x,g12.5))" &
        ) grd%rMinX,grd%rMaxX
  write ( unit=kIOGridLU, &
          fmt="(2(1x,g12.5))" &
        ) grd%rMinY,grd%rMaxY
  write ( unit=kIOGridLU, &
          fmt="(2(1x,g12.5))" &
        ) minval(grd%rValues),maxval(grd%rValues)
  do j=1,grd%iResY
    write ( unit=kIOGridLU, &
            fmt=* &
          ) grd%rValues(j,:)
  end do
  ! That's it
  close ( unit=kIOGridLU )

  return
end subroutine GRD_WriteToGRD

subroutine GRD_WriteToMatlab(grd,sFileName,io)
  !! function GRD_WriteToMatlab
  !!
  !! Creates a matlab '.m' file for contouring
  !!
  !! Calling Sequence:
  !!    grd => GRD_CreateFromMatlab(grd,sFileName)
  !!
  !! Arguments:
  !!    (in)    type ( GRD_GRID ),pointer :: grd
  !!              GRD_GRID object to be written
  !!    (in)    character(len=*) :: sFileName
  !!              Name of the Matlab(tm) '.m' file
  !!
  ! [ ARGUMENTS ]
  type ( GRD_GRID ),pointer :: grd
  character (len=*),intent(in) :: sFileName
  type ( IO_Status ),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: j,iStat

  open ( unit=kIOGridLU, file=trim(sFileName), status="NEW", iostat=iStat)
  if (IO_Assert( (iStat==0), "MakeGrid: Could not open " // sFileName // " for output",io )) return

  call IO_MessageText("  Writing grid to " // trim(sFileName), io)
  write ( unit=kIOGridLU, &
          fmt="('[x,y]=meshgrid(',g12.5,':',g12.5,':',g12.5,',',g12.5,':',g12.5,':',g12.5,')')" &
        ) grd%rMinX,(grd%rMaxX-grd%rMinX)/float(grd%iResX-1),grd%rMaxX, &
          grd%rMinY,(grd%rMaxY-grd%rMinY)/float(grd%iResY-1),grd%rMaxY
  write ( unit=kIOGridLU, &
          fmt="('f=[')" &
        )
  do j=1,grd%iResY
    write ( unit=kIOGridLU, &
            fmt=* &
          ) grd%rValues(j,:)
  end do
  write ( unit=kIOOutputLU, &
          fmt="(']')" &
        )
  ! That's it
  close ( unit=kIOGridLU )

  return
end subroutine GRD_WriteToMatlab

subroutine GRD_Destroy(grd,io)
  !! subroutine GRD_Destroy
  !! 
  !! Destroys a GRD_GRID object
  !!
  !! Usage:
  !!        call GRD_Destroy(grd)
  !!
  !! Arguments:
  !!   (in)     type(GRD_GRID),pointer :: grd
  !!              The grid to be destroyed
  !!
  ! [ ARGUMENTS ]
  type ( GRD_GRID ),pointer :: grd
  type ( IO_Status ),pointer :: io

  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  if ( io%lDebug ) then
    if (IO_Assert( (associated(grd)), "GRD_Destroy: No GRD_GRID object",io )) return
    if (IO_Assert( (associated(grd%rValues)), "GRD_Destroy: GRD_GRID object has no data",io )) return
  endif

  deallocate( grd%rValues, stat=iStat )
  if (IO_Assert( (iStat==0), "GRD_Destroy: Deallocation error",io )) return
  deallocate( grd, stat=iStat )
  if (IO_Assert( (iStat==0), "GRD_Destroy: Deallocation error",io )) return

  return
end subroutine GRD_Destroy

function cGRD_GetZ(grd,iRow,iCol,io) result(cZ)
  !! function cGRD_GetZ
  !!
  !! Gets the complex coordinate of the center of the cell (iRow,iCol)
  !!
  !! Calling Sequence:
  !!    cZ = cGRD_GetZ(grd,iRow,iCol)
  !!
  !! Arguments:
  !!     (in)     type ( GRD_GRID ),pointer :: grd
  !!                GRD_GRID object to be queried
  !!     (in)     integer :: iRow
  !!                Row number
  !!     (in)     integer :: iCol
  !!                Column number
  !!
  !! Return Value:
  !!              complex :: cZ
  !!                Complex coordinate of the center of cell (iRow,iCol)
  !!                returns HUGE if outside the grid.
  !!
  ! [ ARGUMENTS ]
  type ( GRD_GRID ),pointer :: grd
  integer (kind=ModAEM_Integer),intent(in) :: iRow
  integer (kind=ModAEM_Integer),intent(in) :: iCol
  type ( IO_Status ),pointer :: io

  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cZ

  if ( io%lDebug ) then
    if (IO_Assert( (associated(grd)), "GRD_GetZ: No GRD_GRID object",io )) return
  endif

  if ( iRow < 1 .or. &
       iRow > grd%iResY .or. &
       iCol < 1 .or. &
       iCol > grd%iResX ) then 
    cZ = cmplx(HUGE(ModAEM_Real),HUGE(ModAEM_Real),ModAEM_Real)
  else
    cZ = cmplx( grd%rMinY+(grd%iResY-iRow-1)*grd%rDelta, &
                grd%rMinX+(iCol-1)*grd%rDelta, &
                ModAEM_Real )
  endif

  return
end function cGRD_GetZ

subroutine GRD_Lookup(grd,cZ,iRow,iCol,io)
  !! subroutine GRD_Lookup
  !!
  !! Looks up the row,col coordinate of the point cZ
  !!
  !! Calling Sequence:
  !!    call GRD_Lookup(grd,cZ,iRow,iCol)
  !!
  !! Arguments:
  !!     (in)     type ( GRD_GRID ),pointer :: grd
  !!                GRD_GRID object to be queried
  !!     (in)     complex :: cZ
  !!                Complex coordinate 
  !!     (out)    integer :: iRow
  !!                Row number or -1 if outside the grid
  !!     (out)    integer :: iCol
  !!                Column number or -1 if outside the grid
  !!
  !!
  ! [ ARGUMENTS ]
  type ( GRD_GRID ),pointer :: grd
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(out) :: iRow
  integer (kind=ModAEM_Integer),intent(out) :: iCol
  type ( IO_Status ),pointer :: io
  
  if ( io%lDebug ) then
    if (IO_Assert( (associated(grd)), "GRD_Lookup: No GRD_GRID object",io )) return
  endif

  if ( real(cZ) < grd%rMinX .or. &
       real(cZ) > grd%rMaxX .or. &
       aimag(cZ) < grd%rMinY .or. &
       aimag(cZ) > grd%rMaxY ) then
    iRow = -1
    iCol = -1
  else
    ! Compute iRow and iCol (adjust for off-by one exactly on boundaries)
    iRow = int(grd%rMaxY-aimag(cZ))/grd%rDelta + 1
    if ( iRow > grd%iResY ) iRow = grd%iResY
    iCol = int(real(cZ)-grd%rMinX)/grd%rDelta + 1
    if ( iCol > grd%iResX ) iRow = grd%iResX
  endif
  
  return
end subroutine GRD_Lookup

end module u_grid

