program contour

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


! Generates a contour plot for WhAEM, using a SURFER-format grid file (.GRD).
!
! Usage:
!         contour -f <file> -min <minimum> -max <maximum> -inc <increment> -o <file>
! or
!         contour -f <file> -d <number of contours> -o <file>
!
! contour will add the extension .GRD to the input filename, then

use u_constants
use contouring

implicit none
    
  ! Constants
  integer (kind=ModAEM_Integer),parameter :: kInLU = 20
  integer (kind=ModAEM_Integer),parameter :: kOutLU = 21
  ! Locals
  integer (kind=ModAEM_Integer) :: iStat
  integer (2) :: iArg
  character (len=80) :: s80filename,sTemp
  real (kind=ModAEM_Real) :: rZMin,rZMax,rZIncr

  call GETARG(1, sTemp)
  read ( unit=sTemp, fmt=*, iostat=iStat ) s80filename,rZMin,rZIncr,rZMax
  if ( iStat /= 0 ) then
    print *,"Usage: contour <filename> <minimum> <interval> <maximum>"
    stop
  endif

  ! Open the files
  open (unit=kInLU,  file=trim(s80filename)//".grd", action="READ", status="OLD", iostat=iStat)
  if ( iStat /= 0 ) then
    print *,"Could not open GRD file ",trim(s80filename)//".grd"
    stop
  endif
  !
  open (unit=kOutLU,  file=trim(s80filename)//".ctr", action="WRITE", status="REPLACE", iostat=iStat)
  if ( iStat /= 0 ) then
    print *,"Could not open CTR file ",trim(s80filename)//".ctr"
    stop
  endif

  ! Read the grid file
  call CTR_ReadGRD(kInLU)

  ! Now, draw the contours
  call CTR_DrawContours(kOutLU,rZMin,rZIncr,rZMax)

  close (unit=kInLU)
  close (unit=kOutLU)

  stop
end program contour


