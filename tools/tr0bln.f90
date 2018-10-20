program tr0bln

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

! Converts a ModAEM trace file from the TR0 module to a SURFER boundary-line file

use u_constants
use u_io

implicit none
    
  ! Constants
  integer (kind=ModAEM_Integer),parameter :: kNamLU = 20
  integer (kind=ModAEM_Integer),parameter :: kInLU = 21
  integer (kind=ModAEM_Integer),parameter :: kOutLU = 22
  ! Locals
  integer (kind=ModAEM_Integer) :: i,j,iStat
  integer (2) :: iArg
  character (len=80) :: s80filename
  character (len=132) :: s132buf
  real (kind=ModAEM_Real),dimension(10000) :: rX,rY

  iArg = 1
  open ( unit=kNamLU, file='modaem.nam' )
  read ( unit=kNamLU, fmt=* ) s80filename
  close ( unit=kNamLU )
  open (unit=kInLU,  file=trim(s80filename)//".tr0", action="READ", status="OLD", iostat=iStat)
  if (IO_Assert( (iStat==0), "Could not open TR0 file" )) return
  open (unit=kOutLU,  file=trim(s80filename)//".bln", action="WRITE", status="REPLACE", iostat=iStat)
  if (IO_Assert( (iStat==0), "Could not open BLN file" )) return

  ! Read TR0 file and write BLN file
  do
    ! Read a record and process it.
    read (unit=kInLU,fmt="(a132)",iostat=iStat) s132buf
    if ( iStat /= 0 ) stop "Execution complete"
    if ( s132buf(1:3) == "STA" ) then
      i = 0
    elseif ( s132buf(1:3) == "END" ) then
      i = i+1
      read ( unit=s132buf(6:),fmt=* ) rX(i),rY(i)
      write (unit=kOutLU,fmt=*) i,1
      do j=1,i
        write (unit=kOutLU,fmt=*) rX(j),rY(j)
      end do
    elseif ( s132buf(1:3) == "GRI" ) then
      cycle
    elseif ( s132buf(1:3) == "CEL" ) then
      cycle
    elseif ( s132buf(1:3) == "WL0" ) then
      cycle
    else
      i = i+1
      read ( unit=s132buf,fmt=* ) rX(i),rY(i)
    endif
  end do

  close (unit=kInLU)
  close (unit=kOutLU)

  stop
end program tr0bln

