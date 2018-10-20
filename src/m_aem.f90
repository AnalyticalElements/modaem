module m_aem

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

! AEMver control module for ModAEM. Allows selection of AEMver types, handles
! matrix generation and AEMution.

use u_constants
use u_io
use u_matrix
use f_well
use f_dipole
use f_pond
use m_aqu
use m_wl0
use m_wl1
use m_pd0
use m_ls0
use m_ls1
use m_ls2
use m_hb0

implicit none

public

  type :: AEM_DOMAIN
    !! type AEM_DOMAIN
    !!
    !! Type that holds a complete AEM domain
    !!
    !! Members:
    !!   type ( FWL_COLLECTION ),pointer :: fwl
    !!     Collection of well functions
    !!   type ( FPD_COLLECTION ),pointer :: fpd
    !!     Collection of pond functions
    !!   type ( FDP_COLLECTION ),pointer :: fdp
    !!     Collection of dipole functions
    !!   type ( WL0_COLLECTION ),pointer :: wl0
    !!     WL0 elements -- discharge wells
    !!   type ( WL1_COLLECTION ),pointer :: wl1
    !!     WL1 elements -- head-specified wells
    !!   type ( PD0_COLLECTION ),pointer :: pd0
    !!     PD0 elements -- discharge wells
    !!   type ( LS0_COLLECTION ),pointer :: ls0
    !!     LS0 elements -- discharge linesinks
    !!   type ( LS1_COLLECTION ),pointer :: ls1
    !!     LS1 elements -- head linesinks
    !!   type ( LS2_COLLECTION ),pointer :: ls2
    !!     LS2 elements -- head linesinks with resistance
    !!   type ( HB0_COLLECTION ),pointer :: hb0
    !!     HB0 elements -- no-flow barriers
    !!   type ( MAT_MATRIX ),pointer :: mat
    !!     Solution matrix object
    !!   integer :: iAQUStart
    !!     Starting equation number for the AQU module
    !!   integer :: iAQUNUnk
    !!     Number of unknowns for the AQU module
    !!   integer :: iLS1Start
    !!     Starting equation number for the LS1 module
    !!   integer :: iLS1NUnk
    !!     Number of unknowns for the LS1 module
    !!   integer :: iHB0Start
    !!     Starting equation number for the HB0 module
    !!   integer :: iHB0NUnk
    !!     Number of unknowns for the HB0 module
    !!   integer :: iWL1Start
    !!     Starting equation number for the WL1 module
    !!   integer :: iWL1NUnk
    !!     Number of unknowns for the WL1 module
    !! 
    integer (kind=ModAEM_Integer) :: iAQUStart
    integer (kind=ModAEM_Integer) :: iAQUNUnk
    integer (kind=ModAEM_Integer) :: iLS1Start
    integer (kind=ModAEM_Integer) :: iLS1NUnk
    integer (kind=ModAEM_Integer) :: iLS2Start
    integer (kind=ModAEM_Integer) :: iLS2NUnk
    integer (kind=ModAEM_Integer) :: iHB0Start
    integer (kind=ModAEM_Integer) :: iHB0NUnk
    integer (kind=ModAEM_Integer) :: iWL1Start
    integer (kind=ModAEM_Integer) :: iWL1NUnk
    ! Pointers to other collection objects
    type ( AQU_COLLECTION ),pointer :: aqu
    type ( FWL_COLLECTION ),pointer :: fwl
    type ( FPD_COLLECTION ),pointer :: fpd
    type ( FDP_COLLECTION ),pointer :: fdp
    type ( WL0_COLLECTION ),pointer :: wl0
    type ( PD0_COLLECTION ),pointer :: pd0
    type ( LS0_COLLECTION ),pointer :: ls0
    type ( LS1_COLLECTION ),pointer :: ls1
    type ( LS2_COLLECTION ),pointer :: ls2
    type ( HB0_COLLECTION ),pointer :: hb0
    type ( WL1_COLLECTION ),pointer :: wl1
    type ( MAT_MATRIX ),pointer :: mat
  end type AEM_DOMAIN

  real (kind=ModAEM_Real),private,parameter :: MOVEPOINT=1.0e-6_ModAEM_Real

contains

!**pd Modified
!**pd function AEM_Create(iNWL0,iNPD0,iNLS0,iNLS1,iNHB0) result(aem)
function AEM_Create(io) result(aem)
  !! function AEM_Create
  !!
  !! Creates a new AEM_DOMAIN object
  !!
  !! Calling Sequence:
  !!    AEM => AEM_Create()
  !!
  !! Arguments:
  !!
  !! Return Value:
  !!   On success, AEM points to a new AEM_DOMAIN object
  !!   On failure (allocation error), fatal error
  !!
  ! [ ARGUMENTS ]
  ! [ RETURN VALUE ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat

  allocate (aem,stat=iStat)
  if (IO_Assert( (iStat == 0), "AEM_Create: allocation failed",io )) return
  
  ! Initialize
  aem%iAQUStart = 0
  aem%iAQUNUnk = 0
  aem%iLS1Start = 0
  aem%iLS1NUnk = 0
  aem%iLS2NUnk = 0
  aem%iLS2NUnk = 0
  aem%iHB0Start = 0
  aem%iHB0NUnk = 0
  aem%iWL1Start = 0
  aem%iWL1Start = 0
  ! Create the module objects
  nullify(aem%fwl)
  nullify(aem%fpd)
  nullify(aem%fdp)
  nullify(aem%aqu)
  !**pd Init pointers to NULL for now
  aem%wl0 => WL0_Create(io)
  if ( io%lError ) return
  aem%pd0 => PD0_Create(io)
  if ( io%lError ) return
  aem%ls0 => LS0_Create(io)
  if ( io%lError ) return
  aem%ls1 => LS1_Create(io)
  if ( io%lError ) return
  aem%ls2 => LS2_Create(io)
  if ( io%lError ) return
  aem%hb0 => HB0_Create(io)
  if ( io%lError ) return
  aem%wl1 => WL1_Create(io)
  if ( io%lError ) return

  aem%mat => MAT_Create(io)
  if ( io%lError ) return

  return
end function AEM_Create


subroutine AEM_Destroy(aem,io)
  !! subroutine AEM_Destroy
  !!
  !! Frees memory allocated for an AEM_DOMAIN object
  !!
  !! Calling Sequence:
  !!     call AEM_Destroy(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  !!
  !! Return Value:
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  ! [ LOCALS ]
  integer (kind=ModAEM_Integer):: iStat

  call WL0_Destroy(aem%wl0,io)
  if ( io%lError ) return
  call PD0_Destroy(aem%pd0,io)  
  if ( io%lError ) return
  call LS0_Destroy(aem%ls0,io)  
  if ( io%lError ) return
  call LS1_Destroy(aem%ls1,io)  
  if ( io%lError ) return
  call HB0_Destroy(aem%hb0,io)  
  if ( io%lError ) return
  call WL1_Destroy(aem%wl1,io)  
  if ( io%lError ) return
  call MAT_Destroy(aem%mat,io)
  if ( io%lError ) return
  deallocate (aem,stat=iStat)
  if (IO_Assert( (iStat == 0), "AEM_Destroy: deallocation failed",io )) return
  
  return
end subroutine AEM_Destroy

function cAEM_Potential(aem,cZ,io) result(cOmega)
  !! function cAEM_Potential
  !!
  !! Returns the complex potential at cZ
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Potential(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cOmega
  ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cZArg

  cZArg = cZ
  do
    if ( lAQU_CheckPoint(aem%aqu,cZArg,MOVEPOINT,io) .or. &
         lWL0_CheckPoint(aem%wl0,cZArg,MOVEPOINT,io) .or. &
         lPD0_CheckPoint(aem%pd0,cZArg,MOVEPOINT,io) .or. &
         lLS0_CheckPoint(aem%ls0,cZArg,MOVEPOINT,io) .or. &
         lLS1_CheckPoint(aem%ls1,cZArg,MOVEPOINT,io) .or. &
         lHB0_CheckPoint(aem%hb0,cZArg,MOVEPOINT,io) .or. &
         lWL1_CheckPoint(aem%wl1,cZArg,MOVEPOINT,io) &
       ) then
      cZArg = cZArg + MOVEPOINT
    else
      exit
    endif
  end do

  cOmega = cAQU_Potential(aem%aqu,cZArg,io) + &
           cFWL_Potential(aem%fwl,cZArg,io) + &
           cFPD_Potential(aem%fpd,cZArg,io) + &
           cFDP_Potential(aem%fdp,cZArg,io) 

  return
end function cAEM_Potential

function cAEM_Discharge(aem,cZ,io) result(cQ)
  !! function cAEM_Discharge
  !!
  !! Returns the discharge vector at cZ
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Discharge(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cQ
  ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cZArg

  cZArg = cZ
  do
    if ( lAQU_CheckPoint(aem%aqu,cZArg,MOVEPOINT,io) .or. &
         lWL0_CheckPoint(aem%wl0,cZArg,MOVEPOINT,io) .or. &
         lPD0_CheckPoint(aem%pd0,cZArg,MOVEPOINT,io) .or. &
         lLS0_CheckPoint(aem%ls0,cZArg,MOVEPOINT,io) .or. &
         lLS1_CheckPoint(aem%ls1,cZArg,MOVEPOINT,io) .or. &
         lHB0_CheckPoint(aem%hb0,cZArg,MOVEPOINT,io) .or. &
         lWL1_CheckPoint(aem%wl1,cZArg,MOVEPOINT,io) &
       ) then
      cZArg = cZArg + MOVEPOINT
    else
      exit
    endif
  end do

  cQ = cAQU_Discharge(aem%aqu,cZArg,io) + &
       cFWL_Discharge(aem%fwl,cZArg,io) + &
       cFPD_Discharge(aem%fpd,cZArg,io) + &
       cFDP_Discharge(aem%fdp,cZArg,io) 

  return
end function cAEM_Discharge

function rAEM_Recharge(aem,cZ,io) result(rGamma)
  !! function rAEM_Recharge
  !!
  !! Returns the recharge rate at cZ
  !!
  !! Calling Sequence:
  !!    rGamma = cAEM_Discharge(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rGamma
  ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cZArg

  cZArg = cZ
  ! NOTE: No singularities in recharge (in ModAEM)

  rGamma = rAQU_Recharge(aem%aqu,cZArg,io) + &
       rFWL_Recharge(aem%fwl,cZArg,io) + &
       rFPD_Recharge(aem%fpd,cZArg,io) + &
       rFDP_Recharge(aem%fdp,cZArg,io) 

  return
end function rAEM_Recharge

function rAEM_Extraction(aem,io) result(rQ)
  !! function rAEM_Extraction
  !!
  !! Returns the net extraction rate of the model (should be zero!)
  !!
  !! Calling Sequence:
  !!    rQ = cAEM_Extraction(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rQ
  ! [ LOCALS ]

  rQ = rAQU_Extraction(aem%aqu,io) + &
       rFWL_Extraction(aem%fwl,io) + &
       rFPD_Extraction(aem%fpd,io) + &
       rFDP_Extraction(aem%fdp,io) 

  return
end function rAEM_Extraction

function cAEM_DischargeAtWell(aem,cZ,iFWLIndex,io) result(cQ)
  !! function cAEM_DischargeAtWell
  !!
  !! Returns the discharge vector at cZ, excluding the well with index iFWLIndex
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Discharge(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  integer (kind=ModAEM_Integer),intent(in) :: iFWLIndex
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cQ
  ! [ LOCALS ]
  complex (kind=ModAEM_Real) :: cZArg

  cZArg = cZ
  do
    if ( lAQU_CheckPoint(aem%aqu,cZArg,MOVEPOINT,io) .or. &
         lWL0_CheckPoint(aem%wl0,cZArg,MOVEPOINT,io) .or. &
         lPD0_CheckPoint(aem%pd0,cZArg,MOVEPOINT,io) .or. &
         lLS0_CheckPoint(aem%ls0,cZArg,MOVEPOINT,io) .or. &
         lLS1_CheckPoint(aem%ls1,cZArg,MOVEPOINT,io) .or. &
         lHB0_CheckPoint(aem%hb0,cZArg,MOVEPOINT,io) .or. &
         lWL1_CheckPoint(aem%wl1,cZArg,MOVEPOINT,io) &
       ) then
      cZArg = cZArg + MOVEPOINT
    else
      exit
    endif
  end do

  cQ = cAEM_Discharge(aem,cZArg,io) - &
       cFWL_Discharge(aem%fwl,cZArg,io,iFWLIndex,1) 

  return
end function cAEM_DischargeAtWell

function rAEM_Head(aem,cZ,io) result(rHead)
  !! function rAEM_Head
  !!
  !! Returns the head at cZ
  !!
  !! Calling Sequence:
  !!    rHead = rAEM_Head(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rHead
  ! [ LOCALS ]

  rHead = rAQU_PotentialToHead(aem%aqu,real(cAEM_Potential(aem,cZ,io)),cZ,io)

  return
end function rAEM_Head

function rAEM_Flow(aem,cZ,io) result(rFlow)
  !! function rAEM_Flow
  !!
  !! Returns the integrated flow across the path cZ
  !!
  !! Calling Sequence:
  !!    rFlow = rAEM_Flow(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),dimension(:),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  real (kind=ModAEM_Real) :: rFlow
  ! [ LOCALS ]
  complex (kind=ModAEM_Real),dimension(:),allocatable :: cZArg
  integer (kind=ModAEM_Integer) :: iStat,i

  allocate ( cZArg(size(cZ)),stat=iStat )
  if (IO_Assert( (iStat==0),"rAEM_Flow: Allocation failed",io )) return
  cZArg = cZ
  do i=1,size(cZArg)
    do
      if ( lAQU_CheckPoint(aem%aqu,cZArg(i),MOVEPOINT,io) .or. &
           lWL0_CheckPoint(aem%wl0,cZArg(i),MOVEPOINT,io) .or. &
           lPD0_CheckPoint(aem%pd0,cZArg(i),MOVEPOINT,io) .or. &
           lLS0_CheckPoint(aem%ls0,cZArg(i),MOVEPOINT,io) .or. &
           lLS1_CheckPoint(aem%ls1,cZArg(i),MOVEPOINT,io) .or. &
           lHB0_CheckPoint(aem%hb0,cZArg(i),MOVEPOINT,io) .or. &
           lWL1_CheckPoint(aem%wl1,cZArg(i),MOVEPOINT,io) &
         ) then
        cZArg = cZArg + MOVEPOINT
      else
        exit
      endif
    end do
  end do

  rFlow = rAQU_Flow(aem%aqu,cZArg,io) + &
          rFWL_Flow(aem%fwl,cZArg,io) + &
          rFPD_Flow(aem%fpd,cZArg,io) + &
          rFDP_Flow(aem%fdp,cZArg,io) 
  deallocate( cZArg )

  return
end function rAEM_Flow

function cAEM_Velocity(aem,cZ,io) result(cV)
  !! function cAEM_DischargeAtWell
  !!
  !! Returns the head at cZ
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Discharge(aem,cZ)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    complex :: cZ
  !!              Complex coordinate in question
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  complex (kind=ModAEM_Real),intent(in) :: cZ
  type ( IO_Status ),pointer :: io
  ! [ RETURN VALUE ]
  complex (kind=ModAEM_Real) :: cV
  ! [ LOCALS ]

  cV = cAQU_DischargeToVelocity(aem%aqu,cAEM_Discharge(aem,cZ,io), &
                                cZ,real(cAEM_Potential(aem,cZ,io)),io)

  return
end function cAEM_Velocity

subroutine AEM_AllocateFunctions(aem,io)
  !! subroutine AEM_Alloc
  !!
  !! Allocates the space in the MAT, FWL, FPD, and FDP modules
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Discharge(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
 ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iNFWL,iNFPD,iNFDP
  integer (kind=ModAEM_Integer) :: i

  !*******************************************************************************************
  ! Compute the required space for all function modules
  !*******************************************************************************************
  iNFWL = 0
  iNFPD = 0
  iNFDP = 0

  !*******************************************************************************************
  ! Given-strength (no equations necessary) go here
  !*******************************************************************************************
  ! WL0 Module (given wells)
  iNFWL = iNFWL + iWL0_GetInfo(aem%wl0,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iWL0_GetInfo(aem%wl0,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iWL0_GetInfo(aem%wl0,SIZE_FDP,io)
  if ( io%lError ) return
  ! PD0 Module (given ponds)
  iNFWL = iNFWL + iPD0_GetInfo(aem%pd0,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iPD0_GetInfo(aem%pd0,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iPD0_GetInfo(aem%pd0,SIZE_FDP,io)
  if ( io%lError ) return
  ! LS0 Module (given linesinks)
  iNFWL = iNFWL + iLS0_GetInfo(aem%ls0,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iLS0_GetInfo(aem%ls0,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iLS0_GetInfo(aem%ls0,SIZE_FDP,io)
  if ( io%lError ) return
  !*******************************************************************************************
  ! Modify for new element types: 
  ! Add additional given-strength modules above...
  !*******************************************************************************************

  !*******************************************************************************************
  ! Unknown-strength (main matrix equations necessary) go here
  !*******************************************************************************************
  ! AQU Module (reference point, etc.)
  iNFWL = iNFWL + iAQU_GetInfo(aem%aqu,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iAQU_GetInfo(aem%aqu,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iAQU_GetInfo(aem%aqu,SIZE_FDP,io)
  if ( io%lError ) return
  ! LS1 Module (head-specified linesinks)
  iNFWL = iNFWL + iLS1_GetInfo(aem%ls1,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iLS1_GetInfo(aem%ls1,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iLS1_GetInfo(aem%ls1,SIZE_FDP,io)
  if ( io%lError ) return
  ! HB0 Module (no-flow boundaries)
  iNFWL = iNFWL + iHB0_GetInfo(aem%hb0,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iHB0_GetInfo(aem%hb0,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iHB0_GetInfo(aem%hb0,SIZE_FDP,io)
  if ( io%lError ) return
  ! WL1 Module (head-specified wells)
  iNFWL = iNFWL + iWL1_GetInfo(aem%wl1,SIZE_FWL,io)
  if ( io%lError ) return
  iNFPD = iNFPD + iWL1_GetInfo(aem%wl1,SIZE_FPD,io)
  if ( io%lError ) return
  iNFDP = iNFDP + iWL1_GetInfo(aem%wl1,SIZE_FDP,io)
  if ( io%lError ) return

  !*******************************************************************************************
  ! Modify for new element types: 
  ! Add additional unknown-strength modules above...
  !*******************************************************************************************

  !*******************************************************************************************
  ! Now, allocate the space for all the functions and the matrix generator 
  !*******************************************************************************************
  ! Well functions...
  aem%fwl => FWL_Create(iNFWL,io)
  !*******************************************************************************************
  ! Pond functions...
  aem%fpd => FPD_Create(iNFPD,io)
  !*******************************************************************************************
  ! dipole functions...
  aem%fdp => FDP_Create(iNFDP,io)
  !*******************************************************************************************
  ! Modify for new function types: 
  ! Add additional function modules above...
  !*******************************************************************************************

  !*******************************************************************************************
  ! If we made it here, allocation is complete. Now use the element module setup
  ! routines to set up the functions and matrix generator
  !*******************************************************************************************
  call AQU_SetupFunctions ( aem%aqu, aem%fwl, aem%fdp,io)
  if ( io%lError ) return
  call WL0_SetupFunctions ( aem%wl0, aem%fwl, io)
  if ( io%lError ) return
  call PD0_SetupFunctions ( aem%pd0, aem%fpd, io)
  if ( io%lError ) return
  call LS0_SetupFunctions ( aem%ls0, aem%fwl, aem%fdp, io)
  if ( io%lError ) return
  call LS1_SetupFunctions ( aem%ls1, aem%fwl, aem%fdp, io)
  if ( io%lError ) return
  call HB0_SetupFunctions ( aem%hb0, aem%fdp, io)
  if ( io%lError ) return
  call WL1_SetupFunctions ( aem%wl1, aem%fwl, io)
  if ( io%lError ) return

  return
end subroutine AEM_AllocateFunctions

subroutine AEM_AllocateMatrix(aem,io)
  !! subroutine AEM_Alloc
  !!
  !! Allocates the space in the MAT module
  !!
  !! Calling Sequence:
  !!    cOmega = cAEM_Discharge(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
 ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iNEQ,iNUN
  integer (kind=ModAEM_Integer) :: i

  !*******************************************************************************************
  ! Compute the required space for all function modules
  !*******************************************************************************************
  iNEQ = 0
  iNUN = 0

  !*******************************************************************************************
  ! Unknown-strength (main matrix equations necessary) go here
  !*******************************************************************************************
  ! AQU Module (reference point, etc.)
  iNEQ = iNEQ + iAQU_GetInfo(aem%aqu,SIZE_EQUATIONS,io)
  if ( io%lError ) return
  aem%iAQUStart = 1
  aem%iAQUNUnk = iAQU_GetInfo(aem%aqu,SIZE_UNKNOWNS,io)
  if ( io%lError ) return
  iNUN = iNUN + aem%iAQUNUnk
  ! LS1 Module (head-specified linesinks)
  iNEQ = iNEQ + iLS1_GetInfo(aem%ls1,SIZE_EQUATIONS,io)
  if ( io%lError ) return
  aem%iLS1Start = aem%iAQUNUnk + 1
  aem%iLS1NUnk = iLS1_GetInfo(aem%ls1,SIZE_UNKNOWNS,io)
  if ( io%lError ) return
  iNUN = iNUN + aem%iLS1NUnk 
  ! HB0 Module (no-flow boundaries)
  iNEQ = iNEQ + iHB0_GetInfo(aem%hb0,SIZE_EQUATIONS,io)
  if ( io%lError ) return
  aem%iHB0Start = aem%iAQUNUnk + aem%iLS1NUnk + 1
  aem%iHB0NUnk = iHB0_GetInfo(aem%hb0,SIZE_UNKNOWNS,io)
  if ( io%lError ) return
  iNUN = iNUN + aem%iHB0NUnk
  ! WL1 Module (head-specified wells)
  iNEQ = iNEQ + iWL1_GetInfo(aem%wl1,SIZE_EQUATIONS,io)
  if ( io%lError ) return
  aem%iWL1Start = aem%iAQUNUnk + aem%iLS1NUnk + aem%iHB0NUnk + 1
  aem%iWL1NUnk = iWL1_GetInfo(aem%wl1,SIZE_UNKNOWNS,io)
  iNUN = iNUN + aem%iWL1NUnk
  if ( io%lError ) return

  !*******************************************************************************************
  ! Modify for new element types: 
  ! Add additional unknown-strength modules above...
  !*******************************************************************************************

  !*******************************************************************************************
  ! Now, allocate the matrix and control-point arrays. 
  call MAT_Alloc(aem%mat,iNEQ,iNUN,io)
  if ( io%lError ) return

  !*******************************************************************************************
  ! If we made it here, allocation is complete. Now use the element module setup
  ! routines to set up the functions and matrix generator
  !*******************************************************************************************
  call AQU_SetupMatrix ( aem%aqu, aem%mat,io)
  if ( io%lError ) return
  call LS1_SetupMatrix ( aem%ls1, aem%aqu,aem%mat, io)
  if ( io%lError ) return
  call HB0_SetupMatrix ( aem%hb0, aem%mat, io)
  if ( io%lError ) return
  call WL1_SetupMatrix ( aem%wl1, aem%aqu, aem%mat, io)
  if ( io%lError ) return

  return
end subroutine AEM_AllocateMatrix

subroutine AEM_Solve(aem,iNIter,io)
  !! subroutine AEM_Solve
  !!
  !! Performs the solution process
  !!
  !! Calling Sequence:
  !!    call AEM_Solve(aem,iNIter)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!    (in)    integer :: iNIter
  !!              Number of iterations to be performed
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  integer (kind=ModAEM_Integer),intent(inout) :: iNIter
  type ( IO_Status ),pointer :: io
  ! Locals
  integer (kind=ModAEM_Integer) :: iIter
  character (len=132) :: sMessage

  ! Heeeere we go!
  call IO_MessageText("Allocating space for functions and matrix",io)
  call AEM_AllocateFunctions(aem,io)

  do iIter = 1,iNIter
    if ( iWL1_GetInfo(aem%wl1,INFO_REGENERATE,io) /= 0 .or. &
         iLS1_GetInfo(aem%ls1,INFO_REGENERATE,io) /= 0 .or. &
         iLS2_GetInfo(aem%ls2,INFO_REGENERATE,io) /= 0 .or. &
         iHB0_GetInfo(aem%hb0,INFO_REGENERATE,io) /= 0 &
       ) then
      if ( io%lError ) return
      call IO_MessageText("Allocating matrix",io)
      call MAT_Clear(aem%mat,io)
      call AEM_AllocateMatrix(aem,io)
      if ( io%lError ) return
    endif

    call IO_MessageText("Generating matrix",io)
    if ( io%lError ) return
    call AEM_GenerateMatrix(aem,io)
    if ( io%lError ) return

!    if ( io%lDebug ) then
      call MAT_Report(aem%mat,'pre',io)
      if ( io%lError ) return
!    endif

    call IO_MessageText("Decomposing matrix",io)
    call MAT_Decompose(aem%mat,io)
    if ( io%lError ) return

    ! Initialize the CHECK values
    call AEM_ComputeCheck(aem,io)
    if ( io%lError ) return

    ! Iterative Solution scheme
    call IO_MessageText("Generating solution...",io)
    if ( io%lError ) return
    write (unit=sMessage,fmt="('  Iteration: ',i5)") iIter
    call IO_MessageText(sMessage,io)
    call AEM_GenerateRHS(aem,io)
    if ( io%lError ) return
    call AEM_SolveMatrix(aem,io)
    if ( io%lError ) return
    call AEM_Update(aem,io)
    if ( io%lError ) return
  end do
  call IO_MessageText("Solution complete",io)

  return
end subroutine AEM_Solve

subroutine AEM_GenerateMatrix(aem,io)
  !! subroutine AEM_GenerateMatrix
  !!
  !! Generates the matrix coefficients
  !!
  !! Calling Sequence:
  !!    call AEM_GenerateMatrix(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iRow
  real (kind=ModAEM_Real),dimension(:),allocatable :: rARow
  real (kind=ModAEM_Real) :: rRHS,rCheck,rSpecValue
  integer (kind=ModAEM_Integer) :: iStat
  complex (kind=ModAEM_Real),dimension(10) :: cCPZ
  integer (kind=ModAEM_Integer) :: iNCP,iEqType,ic
  integer (kind=ModAEM_Integer) :: iElementType,iElementString,iElementVertex,iElementFlag
  complex (kind=ModAEM_Real) :: cOrientation
  
  ! First allocate the storage for the matrix row being generated
  allocate (rARow(aem%mat%iNVar),stat=iStat)
  if (IO_Assert( (iStat==0), "AEM_GenerateMatrix: Allocation failed",io )) return

  ! This is the main matrix generator loop. It generates the matrix row-by-row, and
  ! is parallelizable.
  do iRow = 1, aem%mat%iNEqn
    ! Get the matrix generator information for the row
    call MAT_GetEquation(aem%mat,iRow,cCPZ,iNCP,iEqType,iElementType,iElementString, &
                         iElementVertex,iElementFlag,cOrientation,rSpecValue,rCheck,io) 
    if ( io%lError ) return

    !******************************************************************************************************
    ! Correct the matrix coefficients and compute the right-hand side constant for the equation; each
    ! element module with unknown strength coefficients should have a call here.
    !******************************************************************************************************
    ! The array constructor below is used to ensure that the proper array size is received
    if ( aem%iAQUNUnk > 0 ) then
      call AQU_ComputeCoefficients(aem%aqu,aem%fwl,aem%fdp,(/(cCPZ(ic),ic=1,iNCP)/), &
                      iEqType,iElementType,iElementString,iElementVertex,iElementFlag, &
                      cOrientation,rARow(aem%iAQUStart:aem%iAQUStart+aem%iAQUNUnk-1),io)
      if ( io%lError ) return
    endif

    ! Generate the coefficients for the LS1 element module and store them in the matrix row
    if ( aem%iLS1NUnk > 0 ) then
      call LS1_ComputeCoefficients(aem%ls1,aem%fwl,aem%fdp,(/(cCPZ(ic),ic=1,iNCP)/), &
                      iEqType,iElementType,iElementString,iElementVertex,iElementFlag, &
                      cOrientation,rARow(aem%iLS1Start:aem%iLS1Start+aem%iLS1NUnk-1),io)
      if ( io%lError ) return
    endif

    ! Generate the coefficients for the HB0 element module and store them in the matrix row
    if ( aem%iHB0NUnk > 0 ) then
      call HB0_ComputeCoefficients(aem%hb0,aem%fdp,(/(cCPZ(ic),ic=1,iNCP)/), &
                      iEqType,iElementType,iElementString,iElementVertex,iElementFlag, &
                      cOrientation,rARow(aem%iHB0Start:aem%iHB0Start+aem%iHB0NUnk-1),io)
      if ( io%lError ) return
    endif

    ! Generate the coefficients for the LS1 element module and store them in the matrix row
    if ( aem%iWL1NUnk > 0 ) then
      call WL1_ComputeCoefficients(aem%wl1,aem%fwl,(/(cCPZ(ic),ic=1,iNCP)/), &
                      iEqType,iElementType,iElementString,iElementVertex,iElementFlag, &
                      cOrientation,rARow(aem%iWL1Start:aem%iWL1Start+aem%iWL1NUnk-1),io)
      if ( io%lError ) return
    endif

    !******************************************************************************************************
    ! Add more element modules here
    !******************************************************************************************************

    ! Store the row into the matrix
    call MAT_SetRow(aem%mat,iRow,rARow(1:aem%mat%iNVar),rRHS,io)
    if ( io%lError ) return
  end do

  deallocate (rARow)

  return
end subroutine AEM_GenerateMatrix

subroutine AEM_GenerateRHS(aem,io)
  !! subroutine AEM_GenerateRHS
  !!
  !! Generates the right-hand side vector
  !!
  !! Calling Sequence:
  !!    call AEM_GenerateMatrix(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iRow
  real (kind=ModAEM_Real) :: rRHS,rCheck,rSpecValue
  integer (kind=ModAEM_Integer) :: iStat
  complex (kind=ModAEM_Real),dimension(10) :: cCPZ
  integer (kind=ModAEM_Integer) :: iNCP,iEqType,iElementType,ic
  integer (kind=ModAEM_Integer) :: iElementString,iElementVertex,iElementFlag
  complex (kind=ModAEM_Real) :: cOrientation
  
  ! This is the main matrix generator loop. It generates the matrix row-by-row, and
  ! is parallelizable.
  do iRow = 1, aem%mat%iNEqn
    ! Get the matrix generator information for the row
    call MAT_GetEquation(aem%mat,iRow,cCPZ,iNCP,iEqType,iElementType,iElementString, &
                         iElementVertex,iElementFlag,cOrientation,rSpecValue,rCheck,io)

    ! Call the proper module to compute the right-hand side
    select case ( iElementType )
      case ( ELEM_AQU )
        rRHS = rAQU_ComputeRHS(aem%aqu,iElementType,iElementString,iElementVertex, &
                            iElementFlag,rSpecValue,rCheck,io)
        if ( io%lError ) return
      case ( ELEM_IN0 )
        rRHS = rAQU_ComputeRHS(aem%aqu,iElementType,iElementString,iElementVertex, &
                            iElementFlag,rSpecValue,rCheck,io)
        if ( io%lError ) return
      case ( ELEM_LS1 )
        rRHS = rLS1_ComputeRHS(aem%ls1,iElementType,iElementString,iElementVertex, &
                            iElementFlag,rSpecValue,rCheck,io)
        if ( io%lError ) return
      case ( ELEM_HB0 )
        rRHS = rHB0_ComputeRHS(aem%hb0,iElementType,iElementString,iElementVertex, &
                            iElementFlag,rSpecValue,rCheck,io)
        if ( io%lError ) return
      case ( ELEM_WL1 )
        rRHS = rWL1_ComputeRHS(aem%wl1,iElementType,iElementString,iElementVertex, &
                            iElementFlag,rSpecValue,rCheck,io)
        if ( io%lError ) return
    end select

    call MAT_SetRHS(aem%mat,iRow,rRHS,io)
    if ( io%lError ) return
  end do

  return
end subroutine AEM_GenerateRHS

subroutine AEM_SolveMatrix(aem,io)
  !! subroutine AEM_SolveMatrix
  !!
  !! Generates the right-hand side vector
  !!
  !! Calling Sequence:
  !!    call AEM_SolveMatrix(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: i
  real (kind=ModAEM_Real) :: rValue
  integer (kind=ModAEM_Integer) :: iElementType,iElementString,iElementVertex,iElementFlag

  call MAT_Solve(aem%mat,io)

  ! Now, finish the job by setting the final results for all elements...
  do i=1,aem%mat%iNEqn
    ! Get the result from the matrix AEMver
    call MAT_GetVariable(aem%mat,i,rValue,iElementType, &
                        iElementString,iElementVertex,iElementFlag,io)

    ! Store it in the proper element data structures
    select case ( iElementType )
      case ( ELEM_AQU )
        call AQU_StoreResult(aem%aqu,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
        if ( io%lError ) return
      case ( ELEM_IN0 )
        call AQU_StoreResult(aem%aqu,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
        if ( io%lError ) return
      case ( ELEM_LS1 )
        call LS1_StoreResult(aem%ls1,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
        if ( io%lError ) return
      case ( ELEM_HB0 )
        call HB0_StoreResult(aem%hb0,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
        if ( io%lError ) return
      case ( ELEM_WL1 )
        call WL1_StoreResult(aem%wl1,rValue,iElementType,iElementString,iElementVertex,iElementFlag,io)
        if ( io%lError ) return
    end select
  end do

  return
end subroutine AEM_SolveMatrix

subroutine AEM_Update(aem,io)
  !! subroutine AEM_Update
  !!
  !! Update all element modules.  The update procedures should:
  !!   1)  Update the strength parameters for all functional modules
  !!   2)  Update any element module-dependent features (e.g. turning
  !!       off "percolating" line-sinks.
  !!   3)  Report any "check" information about the AEMution to this
  !!       point
  !!
  !! Calling Sequence:
  !!    call AEM_Update(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]

  ! Aquifer module
  call AQU_Update(aem%aqu,aem%fwl,aem%fdp,io)
  if ( io%lError ) return
  ! LS1 module
  call LS1_Update(aem%ls1,aem%fwl,aem%fdp,io)
  if ( io%lError ) return
  ! HB0 module
  call HB0_Update(aem%hb0,aem%fdp,io)
  if ( io%lError ) return
  ! WL1 module
  call WL1_Update(aem%wl1,aem%fwl,io)
  if ( io%lError ) return
 ! Add new elements here!

  call AEM_ComputeCheck(aem,io)

  return
end subroutine AEM_Update

subroutine AEM_ComputeCheck(aem,io)
  !! subroutine AEM_ComputeCheck
  !!
  !! Updates the check information for all equations
  !!
  !! Calling Sequence:
  !!    call AEM_Update(aem)
  !!
  !! Arguments:
  !!    (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!              The AEM_DOMAIN object to be used
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  real (kind=ModAEM_Real) :: rValue
  integer (kind=ModAEM_Integer) :: iRow
  real (kind=ModAEM_Real),dimension(:),allocatable :: rARow
  real (kind=ModAEM_Real) :: rRHS,rCheck,rSpecValue
  integer (kind=ModAEM_Integer) :: iStat,ipt
  complex (kind=ModAEM_Real),dimension(10) :: cCPZ
  integer (kind=ModAEM_Integer) :: iNCP,iEqType,ic
  integer (kind=ModAEM_Integer) :: iElementType,iElementString,iElementVertex,iElementFlag
  complex (kind=ModAEM_Real) :: cOrientation

  ! Update the check fields for all equations
  do iRow = 1, aem%mat%iNEqn
    ! Get the matrix generator information for the row
    call MAT_GetEquation(aem%mat,iRow,cCPZ,iNCP,iEqType,iElementType,iElementString, &
                         iElementVertex,iElementFlag,cOrientation,rSpecValue,rCheck,io) 

    ! Look up the modeled result for the equation
    select case (iEqType) 
      case (EQN_HEAD)
        rCheck = real(cAEM_Potential(aem,cCPZ(1),io))
        if ( io%lError ) return
      case (EQN_FLOW)
        rCheck = rAEM_Flow(aem,(/ (cCPZ(ipt),ipt=1,iNCP) /),io)
        if ( io%lError ) return
      case (EQN_INHO)
        rCheck = rIN0_ComputeCheck(aem%aqu%in0,aem%fdp,iElementString,iElementVertex,cCPZ, &
                                   real(cAEM_Potential(aem,cCPZ(1),io)),io)
        if ( io%lError ) return
      case (EQN_DISCHARGE)
        rCheck = real(cAEM_Discharge(aem,cCPZ(1),io))
        if ( io%lError ) return
      case (EQN_RECHARGE)
        rCheck = rAEM_Recharge(aem,cCPZ(1),io)
        if ( io%lError ) return
      case (EQN_CONTINUITY)
        rCheck = rAEM_Extraction(aem,io)
        if ( io%lError ) return
      case (EQN_POTENTIALDIFF)
        rCheck = real(cAEM_Potential(aem,cCPZ(1),io)) - real(cAEM_Potential(aem,cCPZ(2),io))
        if ( io%lError ) return
    end select
    
    call MAT_UpdateEquation(aem%mat,iRow,rCheck,io) 
  end do

  return
end subroutine AEM_ComputeCheck

subroutine AEM_Inquiry(aem,iLU,io)
  !! subroutine AEM_Inquiry
  !!
  !! Writes an inquiry report for all wells to iLU
  !!
  !! Calling Sequence:
  !!    call AEM_Inquiry(wl0,iLU)
  !!
  !! Arguments:
  !!   (in)    type ( AEM_DOMAIN ),pointer :: aem
  !!             AEM_DOMAIN to be used
  !!   (in)    integer :: iLU
  !!             The output LU to receive output
  !!
  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  integer (kind=ModAEM_Integer),intent(in) :: iLU
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  integer (kind=ModAEM_Integer) :: iRow,ipt
  complex (kind=ModAEM_Real),dimension(10) :: cCPZ
  integer (kind=ModAEM_Integer) :: iNCP,iEqType,iElementType,iID
  integer (kind=ModAEM_Integer) :: iElementString,iElementVertex,iElementFlag
  complex (kind=ModAEM_Real) :: cOrientation
  real (kind=ModAEM_Real) :: rSpecVal,rCheck

  if ( io%lDebug ) then
    if (IO_Assert( (associated(aem)), &
                    "AEM_Inquiry: AEM_Create has not been called",io ) ) return
  endif

  do iRow = 1,aem%mat%iNEqn
    ! Get the matrix generator information for the row
    call MAT_GetEquation(aem%mat,iRow,cCPZ,iNCP,iEqType,iElementType,iElementString, &
                        iElementVertex,iElementFlag,cOrientation,rSpecVal,rCheck,io) 
    select case (iEqType)
      case (EQN_HEAD)
        select case (iElementType)
          case (ELEM_AQU)
            write ( unit=iLU, &
                    fmt="(""AEM"",5("","",i9),2("","",e14.6))" &
                  ) iEqType,iElementType,iElementVertex,iElementFlag, &
                            rAEM_Head(aem,cCPZ(1),io),rSpecVal
          case (ELEM_LS1)
            write ( unit=iLU, &
                    fmt="(""AEM"",5("","",i9),2("","",e14.6))" &
                  ) iEqType,iElementType,iElementVertex,iElementFlag, &
                            rAEM_Head(aem,cCPZ(1),io),rSpecVal
          case (ELEM_WL1)
            write ( unit=iLU, &
                    fmt="(""AEM"",5("","",i9),2("","",e14.6))" &
                  ) iEqType,iElementType,iElementVertex,iElementFlag, &
                            rAEM_Head(aem,cCPZ(1),io),rSpecVal
        end select
      case (EQN_FLOW)
        select case (iElementType)
          case (ELEM_HB0)
            write ( unit=iLU, &
                    fmt="(""AEM"",5("","",i9),2("","",e14.6))" &
                  ) iEqType,iElementType,iElementVertex,iElementFlag, &
                            rAEM_Flow(aem,(/ (cCPZ(ipt),ipt=1,iNCP) /),io),rSpecVal
        end select
      case (EQN_INHO)
        ! not yet implemented
      case (EQN_DISCHARGE)
        ! not yet implemented
    end select
  end do
  return
end subroutine AEM_Inquiry

subroutine AEM_Read(aem,io)
  !! Reads an input file (including processing directives) and returns
  !! a populated AEM_DOMAIN object.

  ! [ ARGUMENTS ]
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
   ! [ LOCALS ]
  type (DIRECTIVE),dimension(10),parameter :: dirDirectives = &
                   (/ dirEND,dirDBG,dirPCD,dirAQU,dirWL0,dirPD0, &
                      dirWL1,dirLS0,dirLS1,dirHB0/)
  ! Directives
  character (len=132) :: sMessage
  character (len=132) :: sRecord
  integer (kind=ModAEM_Integer) :: iOpCode
  integer (kind=ModAEM_Integer) :: iStat
  integer (kind=ModAEM_Integer) :: iKeyCode
  integer (kind=ModAEM_Integer) :: iretval
  integer (kind=ModAEM_Integer) :: iNAQU
  integer (kind=ModAEM_Integer) :: iNWL0, iNPD0, iNLS0, iNLS1, iNHB0, iNWL1, iNLS2
  ! Placeholders for function module test calls
  complex (kind=ModAEM_Real) :: cZ1,cZ2,cZC,cZE1,cZE2
  real (kind=ModAEM_Real) :: rR,rTol,rDPhiDX,rDPhiDY
  real (kind=ModAEM_Real) :: rBase,rThick,rHydCond,rPorosity
  logical (kind=ModAEM_Integer) :: lFlag

  call IO_MessageText("Reading AEM module input",io)

  if (IO_Assert( (associated(aem)), &
                  "AEM_Read: the AEM_DOMAIN has not been created",io ) ) return

  ! Here we go!
  do
    call IO_InputRecord(dirDirectives,iOpCode,sRecord,iKeyCode,io)
    if ( io%lError ) return
    select case (iOpCode)
      case (kOpError)
        ! A RunTime error was found during a file read operation.
        if (IO_Assert( .false., "AEM_Read: I/O Error",io )) return
      case (kOpFileEOF)
        ! EOF is unexpected for all ModAEM "ifXXXRead" routines. 
        if (IO_Assert( .false., "AEM_Read: Unexpected EOF",io )) return
      case (kOpData)
        ! A data line was found. The main ModAEM module has no default
        ! data lines. Report the condition.
        if (IO_Assert( .false., "AEM_Read: Unexpected data record",io )) return
      case (kOpEND)
        ! END mark was found. Exit the file parser.
        exit
      case (kOpDBG)
        ! Change the io%lDebug flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetDebug(lFlag,io)
        if ( io%lError ) return
      case (kOpPCD)
        ! Change the IO_Proceed flag
        read (unit=sRecord,fmt=*,iostat=iStat) lFlag
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        call IO_SetProceed(lFlag,io)
        if ( io%lError ) return
      case (kOpAQU)
        ! Read the infinite aquifer properties and then enter the AQU element module
        read ( unit=sRecord,fmt=*,iostat=iStat ) iNAqu,rBase,rThick,rHydCond,rPorosity
        if (IO_Assert( (iStat==0), "AEM_Read: I/O Error",io )) return
        aem%aqu => AQU_Create(iNAqu,rBase,rThick,rHydCond,rPorosity,io)
        call AQU_Read(aem%aqu,io)
        if ( io%lError ) return
      case (kOpWL0)
        ! Enter the WL0 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNWL0      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call WL0_Alloc(aem%wl0, iNWL0,io)                    !**pd
        if ( io%lError ) return
        call WL0_Read(aem%wl0,io)
        if ( io%lError ) return
      case (kOpWL1)
        ! Enter the WL1 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNWL1      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call WL1_Alloc(aem%wl1, iNWL1,io)                    !**pd
        if ( io%lError ) return
        call WL1_Read(aem%wl1,io)
        if ( io%lError ) return
      case (kOpPD0)
        ! Enter the PD0 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNPD0      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call PD0_Alloc(aem%pd0, iNPD0,io)                    !**pd
        if ( io%lError ) return
        call PD0_Read(aem%pd0,io)
        if ( io%lError ) return
      case (kOpLS0)
        ! Enter the LS0 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNLS0      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call LS0_Alloc(aem%ls0, iNLS0,io)                    !**pd
        if ( io%lError ) return
        call LS0_Read(aem%ls0,io)
        if ( io%lError ) return
      case (kOpLS1)
        ! Enter the LS1 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNLS1      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call LS1_Alloc(aem%ls1, iNLS1,io)                    !**pd
        if ( io%lError ) return
        call LS1_Read(aem%ls1,io)
        if ( io%lError ) return
      case (kOpHB0)
        ! Enter the HB0 element module
        read (unit=sRecord,fmt=*,iostat=iStat) iNHB0      !**pd
        if (IO_Assert( (iStat==0),"AEM_Read: I/O error",io)) return !**pd
        call HB0_Alloc(aem%hb0, iNHB0,io)                    !**pd
        if ( io%lError ) return
        call HB0_Read(aem%hb0,io)
        if ( io%lError ) return
    end select
  end do

  call IO_MessageText("Leaving AEM module",io)

  return
end subroutine AEM_Read  

subroutine AEM_Report(aem,io)
  ! The AEMver report shows the error at the boundary conditions specified in the matrix
  ! Arguments
  type ( AEM_DOMAIN ),pointer :: aem
  type ( IO_Status ),pointer :: io
  ! Locals
  integer (kind=ModAEM_Integer) :: iRow
  complex (kind=ModAEM_Real),dimension(10) :: cCPZ
  integer (kind=ModAEM_Integer) :: iNCP,iEqType,iElementType
  integer (kind=ModAEM_Integer) :: iElementString,iElementVertex,iElementFlag
  complex (kind=ModAEM_Real) :: cOrientation
  real (kind=ModAEM_Real) :: rCheck,rSpecVal

  write ( unit=kIOOutputLU, fmt="(a132)" ) repeat("-",132)
  write ( unit=kIOOutputLU, fmt="('Report for module AEM'//'Matrix equation information')" )

  write ( unit=kIOOutputLU, &
          fmt="('Row',t16,'EqType',t31,'EType',t46,'EString',t61,'EVertex',t76,'EFlag',t91,'SpecVal',t106,'Check')" )
  do iRow = 1,aem%mat%iNEqn
    ! Get the matrix generator information for the row
    call MAT_GetEquation(aem%mat,iRow,cCPZ,iNCP,iEqType,iElementType,iElementString, &
                        iElementVertex,iElementFlag,cOrientation,rSpecVal,rCheck,io) 
    write ( unit=kIOOutputLU, &
            fmt="(6(i10,5x),2(g13.5,2x))" &
          ) iRow,iEqType,iElementType,iElementString,iElementVertex,iElementFlag,rSpecVal,rCheck
  end do
  write ( unit=kIOOutputLU, fmt="(/a132)" ) repeat("-",132)

end subroutine AEM_Report

end module m_aem

