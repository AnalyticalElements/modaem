module m_ls0

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

  !! module m_ls0
  !!
  !! Element module for 2-D discharge specified line-sinks
  !!
  !! Module use:
  !!   u_constants  --  Universal ModAEM constant declarations
  !!   f_well     --  Function module for collections of wells
  !!   f_dipole   --  Function module for collections of line-dipoles
  !!
  !! This module provides the necessary functionality for discharge-specified
  !! line-sink elements. Elements are defined as strings of points, with the
  !! discharge rate (per length) specified for each segment.
  !!
  !! Note: ModAEM uses the "traditional" line-sink function (Strack, 1989)
  !! for matrix generation, but uses strings of dipoles terminated by a well
  !! for computational performance, once a solution is achieved.

use u_constants
use u_io
use u_matrix
use f_well
use f_dipole

implicit none

public

  type :: LS0_VERTEX
    !! type LS0_VERTEX
    !!
    !! Type that holds information for one vertex along a line-sink string
    !!
    !! Members:
    !!   complex :: cZC
    !!     The complex coordinate of the vertex
    !!   real :: rSigma
    !!     The sink density at the vertex
    !!   integer :: iFDPIndex
    !!     Index for the vertex entry in the FDP module. Note: set to -1 for the
    !!     last vertex of the string; an element is considered to extend from vertex
    !!     'i' to vertex 'i+1'.
    !! 
    complex (kind=ModAEM_Real) :: cZ
    real (kind=ModAEM_Real) :: rSigma
    integer (kind=ModAEM_Integer) :: iFDPIndex
  end type LS0_VERTEX

  type :: LS0_STRING
    !! type LS0_STRING
    !!
    !! Type that holds information for one line-sink string
    !!
    !! Members:
    !!   type ( LS0_VERTEX ),dimension(:),pointer :: Vertices
    !!     A vector of LS0_VERTEX objects
    !!   integer :: iNPts
    !!     The number of vertices actually in use
    !!   integer :: iFWLIndex
    !!     Index for the string entry in the FWL module. The well extracts the
    !!     total extraction rate for the string.
    !!   integer :: iID
    !!     The ID number for the string
    !! 
    type ( LS0_VERTEX ),dimension(:),pointer :: Vertices
    integer (kind=ModAEM_Integer) :: iNPts
    integer (kind=ModAEM_Integer) :: iFWLIndex
    integer (kind=ModAEM_Integer) :: iID
  end type LS0_STRING

  type :: LS0_COLLECTION
    !! type LS0_COLLECTION
    !!
    !! Type that holds information for all LS0 elements in a layer
    !!
    !! Members:
    !!   type ( LS0_STRING ),dimension(:),pointer :: Strings
    !!     A vector of LS0_STRING objects
    !!   integer :: iNStr
    !!     The number of strings actually in use
    !! 
    type ( LS0_STRING ),dimension(:),pointer :: Strings
    integer (kind=ModAEM_Integer) :: iNStr
  end type LS0_COLLECTION

contains

!**pd Modified to allocate only the collection object
function LS0_Create(io) result (ls0)
  !! function LS0_Create
  !!
  !! Creates a new LS0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    ls0 => LS0_Create()
  !!
  !! Arguments:
  !!
  ! [ ARGUMENTS ]
  ! [ RETURN VALUE ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat

  allocate (ls0,stat=iStat)
  if (IO_Assert( (iStat == 0), "LS0_Create: allocation failed",io )) return
  nullify(ls0%Strings)
  ls0%iNStr = 0

  return
end function LS0_Create


!**pd New LS0_Alloc subroutine
subroutine LS0_Alloc(ls0, iNStr,io)
  !! Subroutine LS0_Alloc
  !! 
  !! Allocates Strings for the LS0_COLLECTION object
  !!
  !! Calling Sequence:
  !!    call LS0_Alloc(ls0, iNStr)
  !!
  !! Arguments:
  !!    (in)    type ( LS0_COLLECTION ),pointer :: ls0
  !!              The LS0_COLLECTION object to be used
  !!    (in)    integer :: iNStr
  !!              The number of strings to make space for
  !!
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  integer (kind=ModAEM_Integer),intent(in) :: iNStr
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  
  allocate (ls0%Strings(iNStr),stat=iStat)
  if (IO_Assert( (iStat == 0), "LS0_Alloc: allocation failed",io )) return

end subroutine LS0_Alloc


!**pd  New LS0_Destroy subroutine
subroutine LS0_Destroy(ls0,io)
  !! subroutine LS0_Destroy
  !!
  !! Frees memory allocated for LS0 Linesinks and strings of vertices 
  !! and the LS0 Collection object
  !!
  !! Calling Sequence:
  !!     call LS0_Destroy(ls0)
  !!
  !! Arguments:
  !!  type ( LS0_COLLECTION ),pointer :: ls0
  !!              Pointer to the LS0_COLLECTION object to be used
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat
  integer (kind=ModAEM_Integer) :: i
  type (LS0_STRING), pointer :: str


  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                   "LS0_Destroy: LS0_Create has not been called",io )) return
  endif
  
  ! First deallocate each string of vertices
  do i = 1, ls0%iNStr
    str => ls0%Strings(i)
    deallocate (str%Vertices,stat=iStat)
    if (IO_Assert( (iStat == 0), "LS0_Destroy: deallocation of Vertices failed",io )) return
  end do
  ! Then deallocate the strings table
  if ( associated(ls0%Strings) ) then
    deallocate (ls0%Strings,stat=iStat)
    if (IO_Assert( (iStat == 0), &
                     "LS0_Destroy: deallocation of Strings failed",io )) return
  endif     
  ! Now deallocate the collection
  deallocate (ls0,stat=iStat)
  if (IO_Assert( (iStat == 0), "LS0_Destroy: deallocation failed",io )) return
  
  return
end subroutine LS0_Destroy


subroutine LS0_New(ls0,Vertices,iNPts,io) 
  !! function LS0_New
  !!
  !! Adds a new LS0_STRING object to the LS0_COLLECTION 'ls0'
  !!
  !! Calling Sequence:
  !!    call LS0_New(ls0,Vertices,iNPt)
  !!
  !! Arguments:
  !!    (in)    type ( LS0_COLLECTION ) :: ls0
  !!              The LS0_COLLECTION object to be used
  !!    (in)    type ( LS0_VERTEX ) :: Vertices(:)
  !!              Vector that defines the points along the barrier
  !!    (in)    integer :: iNPt
  !!              The number of vertices in the string
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( LS0_VERTEX ),dimension(:) :: Vertices
  integer (kind=ModAEM_Integer),intent(in) :: iNPts
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStat
  type ( LS0_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_New: LS0_Create has not been called",io )) return
  endif

  if (IO_Assert( (ls0%iNStr < size(ls0%Strings)), &
                  "LS0_New: Space exhausted",io )) return
  if (IO_Assert( (iNPts <= size(Vertices)), &
                  "LS0_New: Size of provided vertices is inconsistent",io )) return

  ls0%iNStr = ls0%iNStr + 1
  str => ls0%Strings(ls0%iNStr)
  allocate ( str%Vertices(iNPts),stat=iStat )
  if (IO_Assert( (iStat==0),"LS0_New: Allocation failed",io )) return
  str%Vertices = Vertices(1:iNPts)
  str%iNPts = iNPts

  return
end subroutine LS0_New

function iLS0_GetInfo(ls0,iOption,io) result(iValue)
  !! function LS0_GetInfo
  !!
  !! Returns the following sizing requirements for the WL0module
  !!
  !! Calling Sequence:
  !!    iValue = iLS0_GetInfo(wl0,iOption)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer :: ls0
  !!             LS0_COLLECTION to be used
  !!   (out)   integer :: iOption
  !!             The (see u_constants.f90) to be retrieved
  !!
  !! Return Value:
  !!   integer :: iOption
  !!     The requested information for the object. Note: Unrecognized options
  !!     should always return zero; (via 'case default' in 'select' structure)
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  integer (kind=ModAEM_Integer),intent(in) :: iOption
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  integer (kind=ModAEM_Integer) :: iValue
  integer (kind=ModAEM_Integer) :: iStr
  type ( LS0_STRING ),pointer :: str

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_GetInfo: LS0_Create has not been called",io )) return
  endif

  iValue = 0
  select case (iOption)
    case ( SIZE_FWL )
      iValue = ls0%iNStr
    case ( SIZE_FDP )
      do iStr=1,ls0%iNStr
        str => ls0%Strings(ls0%iNStr)
        iValue = iValue + str%iNPts-1
      end do
    case default
      iValue = 0
  end select

  return
end function iLS0_GetInfo

subroutine LS0_SetupFunctions(ls0,fwl,fdp,io)
  !! subroutine LS0_Setup
  !!
  !! This routine sets up the functions in f_well and f_dipole for the line-sinks
  !! Since this module creates given-strength elements, the strengths of 
  !! all functions are computed at set-up time.
  !!
  !! Note: This routine assumes that sufficient space has been allocated
  !! in f_well and in f_dipole by SOL_Alloc.
  !!
  !! Calling Sequence:
  !!    call LS0_Setup(ls0,fwl,fdp)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer
  !!             LS0_COLLECTION object to be used
  !!   (in)    type ( FWL_COLLECTION ),pointer
  !!             FWL_COLLECTION function object 
  !!   (in)    type ( LS0_COLLECTION ),pointer
  !!             FDP_COLLECTION function object 
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( FWL_COLLECTION ),pointer :: fwl
  type ( FDP_COLLECTION ),pointer :: fdp
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rSigma,rDisch
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cRho1,cRho2,cRho3
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: this_vtx,next_vtx,last_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_Setup: LS0_Create has not been called",io )) return
    if (IO_Assert( (associated(fwl)), &
                    "LS0_Setup: Illegal FWL_COLLECTION object",io )) return
    if (IO_Assert( (associated(fdp)), &
                    "LS0_Setup: Illegal FDP_COLLECTION object",io )) return
  endif

  do iStr = 1,ls0%iNStr
    str => ls0%Strings(iStr)
    cRho1 = cZERO
    ! Build dipoles for all segments
    do iVtx = 1,str%iNPts-1  ! Set up nVertices-1 dipoles...
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      cZ1 = this_vtx%cZ
      cZ2 = next_vtx%cZ

      ! Compute the dipole strengths in terms of the given sink density.  NOTE: these
      ! computations may be easily adjusted for linear-strength line-sinks.
      rSigma = rHALF * (this_vtx%rSigma + next_vtx%rSigma)
      rDisch = rSigma*abs(cZ2-cZ1)
      cRho2 = cZERO
      cRho3 = cRho1 + rDisch
      call FDP_New(fdp,cZ1,cZ2,(/cRho1,cRho2,cRho3/),this_vtx%iFDPIndex,io)
      if ( io%lError ) return
      cRho1 = cRho3                               ! Move on to the next one with this Rho value
    end do

    ! Put a well at the end of the string
    last_vtx => str%Vertices(str%iNPts)
    cZ1 = last_vtx%cZ
    call FWL_New(fwl,cZ1,real(cRho3),rZERO,str%iFWLIndex,io)
    if ( io%lError ) return
  end do

  return
end subroutine LS0_SetupFunctions

function lLS0_CheckPoint(ls0,cZ,rTol,io) result(lSing)
  !! logical function lLS0_CheckPoint
  !!
  !! Checks for a singularity near cZ
  !!
  !! Calling Sequence:
  !!    lSing = lLS0_CheckPoint(ls0,cZ,rTol)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer :: ls0
  !!             The LS0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rTol
  !!             Tolerance about cZ to be checked
  !!
  ! [ ARGUMENTS ]
  type (LS0_COLLECTION),pointer :: ls0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  real (kind=ModAEM_Real), intent(in) :: rTol
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lSing
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: vtx

  ! NOTE: Singularities are found at the end points of the BDY elements 
  lSing = .false.
  do iStr = 1,ls0%iNStr
    str => ls0%Strings(iStr)
    do iVtx = 1,str%iNPts
      vtx => str%Vertices(iVtx)
      if ( abs(real(vtx%cZ-cZ,ModAEM_Real)) < rTol .and. abs(aimag(vtx%cZ-cZ))<rTol ) then
        lSing = .true.
        exit
      endif
    end do
  end do
  
  return
end function lLS0_CheckPoint

function lLS0_CheckProximity(ls0,cZ,rProximity,iDir,iElementID,iElementVtx,io) result(lFound)
  !! function lLS0_CheckProximity
  !!
  !! Checks to see if the point specified is near a line-sink, then
  !! moves the point a small distance if it is.  Sets lChange=.true. if a
  !! change is made to cZ.  NOTE: lChange is preset by the caller, so this 
  !! routine only changes it to .true. if necessary.
  !!
  !! An example application is found in cAEM_Potential()
  !!
  !! Calling Sequence:
  !!    lFound = lLS0_CheckProximity(ls0,cZ,rProximity,iDir,iElementID)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer
  !!             LS0_COLLECTION object to be used
  !!   (in)    complex :: cZ
  !!             The point to be checked
  !!   (in)    real :: rProximity
  !!             The tolerance about cZ to be examined
  !!   (in)    integer :: iDir
  !!             Flag for checking line-sinks:
  !!               iDir =  0   -  Check without regard to the direction of travel
  !!               iDir = <0   -  Check for back tracing (ignore gaining line-sinks)
  !!               iDir = >0   -  Check for formard tracing (ignore losing line-sinks)
  !!   (out)   integer :: iElementID
  !!             The element (string) ID for the line-sink found (if lFound is .true.)
  !!   (out)   integer :: iElementVtx
  !!             The vertex associated with the string segment which was found
  !!
  !! Return Value:
  !!   (out)   logical :: lFound
  !!             .true. if a line-sink is found
  !!             .false. if a line-sink is not found
  !!   Note side effects above for iElementID and iElementVtx
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iDir
  real (kind=ModAEM_Real),intent(in) :: rProximity
  integer (kind=ModAEM_Integer),intent(out) :: iElementID,iElementVtx
  type ( IO_STATUS ),pointer :: io
  ! [ RETURN VALUE ]
  logical (kind=ModAEM_Integer) :: lFound
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rDist,rMinDist
  complex (kind=ModAEM_Real) :: cZL,cZC,cBigZ
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_CheckProximity: LS0_Create has not been called",io )) return
  endif

  iElementID = 0
  lFound = .false.
  rMinDist = HUGE(ModAEM_Real)
  do iStr=1,ls0%iNStr
    str => ls0%Strings(iStr)
    ! Sweep through the string, checking on the basis of iDir
    do iVtx=1,str%iNPts-1
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      ! iDir < 0 -- skip gaining line sinks
      if ( iDir < 0 .and. this_vtx%rSigma > rZERO ) cycle
      ! iDir > 0 -- skip losing line sinks
      if ( iDir > 0 .and. this_vtx%rSigma < rZERO ) cycle
      ! Compute the mapped Z-value and then the distance from the linesink
      cZC = rHALF * (this_vtx%cZ + &
                     next_vtx%cZ )
      cZL = rHALF * (this_vtx%cZ - &
                     next_vtx%cZ )
      cBigZ = (cZ-cZC) / cZL
      rDist = abs( aimag(cBigZ) * abs(cZL) )
      if ( rDist < rProximity .and. abs(real(cBigZ)) <= rONE ) then
        if ( rDist < rMinDist ) then
          rMinDist = rDist
          iElementID = str%iID
          iElementVtx = iVtx
          lFound = .true.
        endif
      endif
    end do
  end do

  return
end function lLS0_CheckProximity

subroutine LS0_FindStringPointer(ls0,iLSID,LSString,lfound,io)
  !! subroutine LS0_FindStringPointer
  !!
  !! Finds the linesink string specified by the ID and returns a pointer to it
  !!
  !! Calling Sequence:
  !!    call LS0_FindStringPointer(ls0,iLSID,LSString,lfound)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer :: ls0
  !!             LS0_COLLECTION to be used
  !!   (in)    integer :: iLSID
  !!             The linesink string ID number
  !!   (out)   type ( LS0_STRING ) :: LSString
  !!             Pointer to the linesink string
  !!   (out)   logical :: lFound
  !!             .true. if the well was found
  !!             .false. if the well was not found 
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  integer (kind=ModAEM_Integer),intent(in) :: iLSID
  type ( LS0_STRING ),pointer :: LSString
  logical (kind=ModAEM_Integer) :: lFound
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i

  if ( io%lDebug ) then
    if ( IO_Assert( (associated(ls0)), &
                    "LS0_FindStringPointer: LS0_Create has not been called",io ) ) return
  endif

  lFound = .false.
  do i=1,ls0%iNStr
    LSString => ls0%Strings(i)
    if ( LSString%iID == iLSID ) then
      lFound = .true.
      return
    endif
  end do

  return
end subroutine LS0_FindStringPointer

subroutine LS0_Read(ls0,io)
  !! subroutine LS0_Read
  !!
  !! Reads the line-sinks for the specified LS0_COLLECTIOn from kIO_InputLU
  !!
  !! Calling Sequence:
  !!    call LS0_Read(ls0)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer :: ls0
  !!             LS0_COLLECTION to be populated
  !!
  !! The format of the LS0 section of the input file appears as follows:
  !! LS0
  !! DIM Nstrings
  !! STR NVertices
  !! x y strength
  !! ... Up to NVertices
  !! STRING NVertices
  !! x y strength
  !! ... Up to NVertices
  !! ... Up to NStrings
  !!
  !! NOTE: It is assumed that the LS0 line was found by the caller

  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( IO_STATUS ),pointer :: io
  ! [ LOCAL DIRECTIVES ]
  type (DIRECTIVE),dimension(4),parameter :: dirDirectives = (/ dirEND,dirDBG,dirPCD,dirSTR /)
  ! [ LOCALS ]
  character (len=132) :: sRecord
  character (len=132) :: sMessage
  real (kind=ModAEM_Real) :: rSigma
  complex (kind=ModAEM_Real) :: cZ
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iMaxStr, iMaxVtx
  integer (kind=ModAEM_Integer) :: iID
  integer (kind=ModAEM_Integer) :: iStr
  logical (kind=ModAEM_Integer) :: lFlag
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: vtx

  call IO_MessageText("  Reading LS0 module input",io)
  if ( io%lError ) return

  if (IO_Assert( (associated(ls0)), "LS0_Read: LS0_Create was not called",io )) return

  ! Use ifIO_InputRecord to process the model input file.
  nullify ( str,vtx )
  do 
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation. This
        ! condition is fatal; warn the user, and exit.
        if (IO_Assert( .false., "LS0_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModLS0 "ifXXXRead" routines. 
        if (IO_Assert( .false., "LS0_Read: Unexpected EOF",io )) return
      case (kOpData)
        !****************************************************************************
        ! Here for data records
        !****************************************************************************
        if (IO_Assert( (associated(str)), "LS0_Read: No STR command was found",io )) return
        if (IO_Assert( (str%iNPts<size(str%Vertices)), "LS0_Read: Space exhausted",io )) return
        ! Oh, all right...  Retrieve the data and put it in there.
        read (unit=sRecord,fmt=*,iostat=iStat) cZ,rSigma
        if (IO_Assert( (iStat==0), "LS0_Read: I/O Error",io )) return
        str%iNPts = str%iNPts+1
        vtx => str%Vertices(str%iNPts)
        vtx%cZ = cZ
        vtx%rSigma = rSigma
        vtx%iFDPIndex = -1
      case (kOpEND)
        ! EOD mark was found. Exit the file parser.
        exit
      case (kOpDBG)
        ! Change the io%lDebug flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetDebug(lFlag,io)
        if ( IO%lError ) return
      case (kOpPCD)
        ! Change the IO_Proceed flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetProceed(lFlag,io)
        if ( IO%lError ) return
      case (kOpSTR)
        !****************************************************************************
        ! Here for the STR command -- create a new string of line-sinks
        ! the maximum number of vertices is in the input record
        !****************************************************************************
        if (IO_Assert( (associated(ls0%Strings)), "LS0_Read: No space is allocated",io )) return
        if (IO_Assert( (ls0%iNStr < size(ls0%Strings)), "LS0_Read: Space exhausted",io )) return
        ! Retrive the number of vertices desired...
        read (unit=sRecord,fmt=*,iostat=iStat) iMaxVtx,iID
        if (IO_Assert( (iStat==0), "LS0_Read: I/O Error",io )) return
        if (IO_Assert( (iMaxVtx>0), "LS0_Read: Illegal number of vertices",io )) return
        ! OKAY! Allocate the vertices...
        ls0%iNStr = ls0%iNStr+1
        str => ls0%Strings(ls0%iNStr)
        allocate (str%Vertices(iMaxVtx),stat=iStat)
        if (IO_Assert( (iStat==0), "LS0_Read: Allocation failed",io )) return
        ! Made it!
        str%iFWLIndex = -1         ! No FWL function yet! 
        str%iID = iID
        write ( unit=sMessage, &
                fmt="("" LS0_Read: "",i6,"" vertices allocated"")" &
              ) iMaxVtx
        call IO_MessageText(sMessage,io)
        if ( IO%lError ) return
        str%iNPts = 0      ! Initialize the vertex counter
    end select
  end do  

  call IO_MessageText("  Leaving LS0 module",io)

end subroutine LS0_Read

subroutine LS0_Inquiry(ls0,iLU,io)
  !! subroutine LS0_Inquiry
  !!
  !! Writes an inquiry report for all line-sinks to iLU
  !!
  !! Calling Sequence:
  !!    call LS0_Inquiry(ls0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer
  !!             LS0_COLLECTION object to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  real (kind=ModAEM_Real) :: rLength
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: this_vtx,next_vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_Inquiry: LS0_Create has not been called",io )) return
  endif

  do iStr=1,ls0%iNStr
    str => ls0%Strings(iStr)
    do iVtx=1,str%iNPts
      this_vtx => str%Vertices(iVtx)
      next_vtx => str%Vertices(iVtx+1)
      if ( iVtx < str%iNPts ) then
        rLength = abs(next_vtx%cZ - this_vtx%cZ)
      else
        rLength = rZERO
      endif
      write ( unit=iLU, &
              fmt="(""LS0"",3("","",i9),4("","",e14.6))" &
            ) str%iID, iVtx,this_vtx%cZ,rLength,this_vtx%rSigma
    end do
  end do

  return
end subroutine LS0_Inquiry

subroutine LS0_Report(ls0,io)
  !! subroutine LS0_Report
  !!
  !! Writes a debugging report for all line-sinks to kIOOutputLU
  !!
  !! Calling Sequence:
  !!    call LS0_Report(ls0)
  !!
  !! Arguments:
  !!   (in)    type ( LS0_COLLECTION ),pointer
  !!             LS0_COLLECTION object to be used
  !!
  ! [ ARGUMENTS ]
  type ( LS0_COLLECTION ),pointer :: ls0
  type ( IO_STATUS ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iStr,iVtx
  integer (kind=ModAEM_Integer) :: nWL,nPD,nDP,nEQ,nUN
  type ( LS0_STRING ),pointer :: str
  type ( LS0_VERTEX ),pointer :: vtx

  if ( io%lDebug ) then
    if (IO_Assert( (associated(ls0)), &
                    "LS0_Inquiry: LS0_Create has not been called",io )) return
  endif

  write ( unit=kIOOutputLU, fmt="('Report for module LS0'//'Discharge-specified linesink element information')" )

  if ( .not. associated(ls0%Strings) ) then
    write ( unit=kIOOutputLU, &
            fmt="("" No Strings Allocated""/)" &
          )
  else
    write ( unit=kIOOutputLU, &
            fmt="("" Number of strings: "",i5,""  used: "",i5)" &
          ) size(ls0%Strings),ls0%iNStr

    ! Write the string information for all strings...
    do iStr=1,ls0%iNStr
      str => ls0%Strings(iStr)
      write ( unit=kIOOutputLU, &
              fmt="(" // &
                  "/ "" "",100(""-"") / " // &
                  " "" String: "",i10,"" # Vertices "",i10,"" Well Index"",i10 " // &
                  ")" &
            ) str%iID,str%iNPts,str%iFWLIndex

      ! Write the vertices for this string...
      write ( unit=kIOOutputLU, &
              fmt="("" Vertices:"")" &
            )
      write ( unit=kIOOutputLU, &
              fmt="(""        X"",t15,""       Y"",t30,""   Strength"",t45,""   DP Index"")" &
            )
      do iVtx=1,str%iNPts
        vtx => str%Vertices(iVtx)
        write ( unit=kIOOutputLU, &
                fmt="(1x,d12.5,t15,d12.5,t30,d12.5,t45,i10)" &
              ) vtx%cZ,vtx%rSigma,vtx%iFDPIndex
      end do
    end do

    ! Report the matrix generation values...
    write ( unit=kIOOutputLU, &
            fmt="(/' Function and Matrix values:'/)" &
          )
    write ( unit=kIOOutputLU, &
            fmt="('Number of FWL functions:',t30,i10)" &
          ) iLS0_GetInfo(ls0,SIZE_FWL,io)
    if (io%lError) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of FPD functions:',t30,i10)" &
          ) iLS0_GetInfo(ls0,SIZE_FPD,io)
    if (io%lError) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of FDP functions:',t30,i10)" &
          ) iLS0_GetInfo(ls0,SIZE_FDP,io)
    if (io%lError) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of Equations:',t30,i10)" &
          ) iLS0_GetInfo(ls0,SIZE_EQUATIONS,io)
    if (io%lError) return
    write ( unit=kIOOutputLU, &
            fmt="('Number of Unknowns:',t30,i10)" &
          ) iLS0_GetInfo(ls0,SIZE_UNKNOWNS,io)
    if (io%lError) return
  endif

  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

  return
end subroutine LS0_Report

end module m_ls0

