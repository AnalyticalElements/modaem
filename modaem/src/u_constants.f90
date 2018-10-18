module u_constants

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

  implicit none

  integer,public,parameter :: ModAEM_Integer = selected_int_kind(8)
  integer,public,parameter :: ModAEM_Logical = selected_int_kind(8)
  integer,public,parameter :: ModAEM_Real    = selected_real_kind(10,8)
  integer,public,parameter :: ModAEM_Double  = selected_real_kind(10,8)

  ! Progam-defined types for input file interpreter
  type,public :: DIRECTIVE
    integer (kind=ModAEM_Integer) :: OpCode
    character (len=3) :: Text
  end type DIRECTIVE

  ! OpCode named constants
  ! Error and message codes
  ! Data record found
  integer (kind=ModAEM_Integer),public,parameter :: kOpData = 0
  ! Error detected
  integer (kind=ModAEM_Integer),public,parameter :: kOpError = -1
  ! End-of-file detected
  integer (kind=ModAEM_Integer),public,parameter :: kOpFileEOF = -2
  ! Program control
  ! Begin definition of an AEM model
  integer (kind=ModAEM_Integer),public,parameter :: kOpAEM = 1001
  ! End of data for module
  integer (kind=ModAEM_Integer),public,parameter :: kOpEND = 1002
  ! End of model input data
  integer (kind=ModAEM_Integer),public,parameter :: kOpEOD = 1003
  ! Dimensions an element module
  integer (kind=ModAEM_Integer),public,parameter :: kOpDIM = 1004
  ! Dimensions a limear element item
  integer (kind=ModAEM_Integer),public,parameter :: kOpSTR = 1005
  ! Tell a module to do its thing
  integer (kind=ModAEM_Integer),public,parameter :: kOpRUN = 1006
  ! Switch the IO_Debug flag
  integer (kind=ModAEM_Integer),public,parameter :: kOpDBG = 1007
  ! Switch the IO_Proceed flag
  integer (kind=ModAEM_Integer),public,parameter :: kOpPCD = 1008
  ! Generate an overall report
  integer (kind=ModAEM_Integer),public,parameter :: kOpRPT = 1009
  ! Create a well
  integer (kind=ModAEM_Integer),public,parameter :: kOpWEL = 1010
  ! Standard model modules
  ! Reference information follows
  integer (kind=ModAEM_Integer),public,parameter :: kOpREF = 2001
  ! Regional aquifer property data follows
  integer (kind=ModAEM_Integer),public,parameter :: kOpPRM = 2002
  ! Aquifer boundary follows
  integer (kind=ModAEM_Integer),public,parameter :: kOpBDY = 2003
  ! Solution module
  integer (kind=ModAEM_Integer),public,parameter :: kOpSOL = 2004
  ! Grid module
  integer (kind=ModAEM_Integer),public,parameter :: kOpGRI = 2005
  ! Data inquiry module
  integer (kind=ModAEM_Integer),public,parameter :: kOpINQ = 2006
  ! Analytic element modules
  ! Aquifer data module
  integer (kind=ModAEM_Integer),public,parameter :: kOpAQU = 3000
  ! Strength-specified wells
  integer (kind=ModAEM_Integer),public,parameter :: kOpWL0 = 3001
  ! head-specified wells
  integer (kind=ModAEM_Integer),public,parameter :: kOpWL1 = 3002
  ! Strength-specified line-sinks
  integer (kind=ModAEM_Integer),public,parameter :: kOpLS0 = 3003
  ! head-specified line-sinks
  integer (kind=ModAEM_Integer),public,parameter :: kOpLS1 = 3004
  ! No-flow horizontal barrier
  integer (kind=ModAEM_Integer),public,parameter :: kOpHB0 = 3005
  ! Inhomogeneity module
  integer (kind=ModAEM_Integer),public,parameter :: kOpIN0 = 3007
  ! Inhomogeneity domain
  integer (kind=ModAEM_Integer),public,parameter :: kOpDOM = 3008
  ! Discharge-specifed pond module
  integer (kind=ModAEM_Integer),public,parameter :: kOpPD0 = 3009
  ! Linesinks with resistance, routing
  integer (kind=ModAEM_Integer),public,parameter :: kOpLS2 = 3010
  ! Analysis Modules
  ! Pathline trace module
  integer (kind=ModAEM_Integer),public,parameter :: kOpTR0 = 3506
  ! Other commands
  ! Set the window
  integer (kind=ModAEM_Integer),public,parameter :: kOpWIN = 5001
  ! Select/inquire heads
  integer (kind=ModAEM_Integer),public,parameter :: kOpHEA = 5002
  ! Select/inquire potentials
  integer (kind=ModAEM_Integer),public,parameter :: kOpPOT = 5003
  ! Select/inquire streamfunction
  integer (kind=ModAEM_Integer),public,parameter :: kOpPSI = 5004
  ! Select/inquire discharge
  integer (kind=ModAEM_Integer),public,parameter :: kOpDIS = 5005
  ! Select/inquire velocity
  integer (kind=ModAEM_Integer),public,parameter :: kOpVEL = 5006
  ! Select/inquire total flow
  integer (kind=ModAEM_Integer),public,parameter :: kOpFLO = 5007
  ! Select/inquire gradient in Phi
  integer (kind=ModAEM_Integer),public,parameter :: kOpGRA = 5008
  ! Create a Gaussian RBF in LKG
  integer (kind=ModAEM_Integer),public,parameter :: kOpGAU = 5009
  ! Create a Bessel RBF in LKG
  integer (kind=ModAEM_Integer),public,parameter :: kOpBES = 5010
  ! Read a SURFER grid of leakances in LKG
  integer (kind=ModAEM_Integer),public,parameter :: kOpSUR = 5011
  ! Set a layer leakance in LKG
  integer (kind=ModAEM_Integer),public,parameter :: kOpLAY = 5012
  ! Grid of leakage errors
  integer (kind=ModAEM_Integer),public,parameter :: kOpLER = 5013
  ! Specify an analysis output file
  integer (kind=ModAEM_Integer),public,parameter :: kOpFIL = 5014
  ! Transport parameters
  integer (kind=ModAEM_Integer),public,parameter :: kOpTRA = 5015
  ! Tuning parameters
  integer (kind=ModAEM_Integer),public,parameter :: kOpTUN = 5016
  ! Timing parameters
  integer (kind=ModAEM_Integer),public,parameter :: kOpTIM = 5017 
  ! Specify a single point
  integer (kind=ModAEM_Integer),public,parameter :: kOpPOI = 5018
  ! Specify a line
  integer (kind=ModAEM_Integer),public,parameter :: kOpLIN = 5019
  ! Select/inquire Qx
  integer (kind=ModAEM_Integer),public,parameter :: kOpQ_X = 5020
  ! Select/inquire Qy
  integer (kind=ModAEM_Integer),public,parameter :: kOpQ_Y = 5021
  ! Set a module option
  integer (kind=ModAEM_Integer),public,parameter :: kOpOPT = 5022
  
  ! Pre-packaged standard directives
  type (DIRECTIVE),public,parameter :: dirAEM = DIRECTIVE(kOpAEM,"AEM")
  type (DIRECTIVE),public,parameter :: dirEND = DIRECTIVE(kOpEND,"END")
  type (DIRECTIVE),public,parameter :: dirEOD = DIRECTIVE(kOpEOD,"EOD")
  type (DIRECTIVE),public,parameter :: dirDIM = DIRECTIVE(kOpDIM,"DIM")
  type (DIRECTIVE),public,parameter :: dirSTR = DIRECTIVE(kOpSTR,"STR")
  type (DIRECTIVE),public,parameter :: dirRUN = DIRECTIVE(kOpRUN,"RUN")
  type (DIRECTIVE),public,parameter :: dirDBG = DIRECTIVE(kOpRUN,"DBG")
  type (DIRECTIVE),public,parameter :: dirPCD = DIRECTIVE(kOpRUN,"PCD")
  type (DIRECTIVE),public,parameter :: dirRPT = DIRECTIVE(kOpRPT,"RPT")
  type (DIRECTIVE),public,parameter :: dirWEL = DIRECTIVE(kOpWEL,"WEL")
  !
  type (DIRECTIVE),public,parameter :: dirREF = DIRECTIVE(kOpREF,"REF")
  type (DIRECTIVE),public,parameter :: dirPRM = DIRECTIVE(kOpPRM,"PRM")
  type (DIRECTIVE),public,parameter :: dirBDY = DIRECTIVE(kOpBDY,"BDY")
  type (DIRECTIVE),public,parameter :: dirSOL = DIRECTIVE(kOpSOL,"SOL")
  type (DIRECTIVE),public,parameter :: dirGRI = DIRECTIVE(kOpGRI,"GRI")
  type (DIRECTIVE),public,parameter :: dirINQ = DIRECTIVE(kOpINQ,"INQ")
  type (DIRECTIVE),public,parameter :: dirTR0 = DIRECTIVE(kOpTR0,"TR0")
  !
  type (DIRECTIVE),public,parameter :: dirAQU = DIRECTIVE(kOpAQU,"AQU")
  type (DIRECTIVE),public,parameter :: dirWL0 = DIRECTIVE(kOpWL0,"WL0")
  type (DIRECTIVE),public,parameter :: dirWL1 = DIRECTIVE(kOpWL1,"WL1")
  type (DIRECTIVE),public,parameter :: dirLS0 = DIRECTIVE(kOpLS0,"LS0")
  type (DIRECTIVE),public,parameter :: dirLS1 = DIRECTIVE(kOpLS1,"LS1")
  type (DIRECTIVE),public,parameter :: dirHB0 = DIRECTIVE(kOpHB0,"HB0")
  type (DIRECTIVE),public,parameter :: dirIN0 = DIRECTIVE(kOpIN0,"IN0")
  type (DIRECTIVE),public,parameter :: dirDOM = DIRECTIVE(kOpDOM,"DOM")
  type (DIRECTIVE),public,parameter :: dirPD0 = DIRECTIVE(kOpPD0,"PD0")
  type (DIRECTIVE),public,parameter :: dirLS2 = DIRECTIVE(kOpLS2,"LS2")
  !
  type (DIRECTIVE),public,parameter :: dirWIN = DIRECTIVE(kOpWIN,"WIN")
  type (DIRECTIVE),public,parameter :: dirHEA = DIRECTIVE(kOpHEA,"HEA")
  type (DIRECTIVE),public,parameter :: dirPOT = DIRECTIVE(kOpPOT,"POT")
  type (DIRECTIVE),public,parameter :: dirPSI = DIRECTIVE(kOpPSI,"PSI")
  type (DIRECTIVE),public,parameter :: dirDIS = DIRECTIVE(kOpDIS,"DIS")
  type (DIRECTIVE),public,parameter :: dirVEL = DIRECTIVE(kOpVEL,"VEL")
  type (DIRECTIVE),public,parameter :: dirFLO = DIRECTIVE(kOpFLO,"FLO")
  type (DIRECTIVE),public,parameter :: dirGRA = DIRECTIVE(kOpGRA,"GRA")
  type (DIRECTIVE),public,parameter :: dirGAU = DIRECTIVE(kOpGAU,"GAU")
  type (DIRECTIVE),public,parameter :: dirBES = DIRECTIVE(kOpBES,"BES")
  type (DIRECTIVE),public,parameter :: dirSUR = DIRECTIVE(kOpSUR,"SUR")
  type (DIRECTIVE),public,parameter :: dirLAY = DIRECTIVE(kOpLAY,"LAY")
  type (DIRECTIVE),public,parameter :: dirLER = DIRECTIVE(kOpLER,"LER")
  type (DIRECTIVE),public,parameter :: dirFIL = DIRECTIVE(kOpFIL,"FIL")
  type (DIRECTIVE),public,parameter :: dirTRA = DIRECTIVE(kOpTRA,"TRA")
  type (DIRECTIVE),public,parameter :: dirTUN = DIRECTIVE(kOpTUN,"TUN")
  type (DIRECTIVE),public,parameter :: dirTIM = DIRECTIVE(kOpTIM,"TIM")
  type (DIRECTIVE),public,parameter :: dirPOI = DIRECTIVE(kOpPOI,"POI")
  type (DIRECTIVE),public,parameter :: dirLIN = DIRECTIVE(kOpLIN,"LIN")
  type (DIRECTIVE),public,parameter :: dirQ_X = DIRECTIVE(kOpQ_X,"Q_X")
  type (DIRECTIVE),public,parameter :: dirQ_Y = DIRECTIVE(kOpQ_Y,"Q_Y")
  type (DIRECTIVE),public,parameter :: dirOPT = DIRECTIVE(kOpOPT,"OPT")

  ! Other constants...

  ! Error Constants
  integer (kind=ModAEM_Integer),public,parameter :: errOK = 0
  integer (kind=ModAEM_Integer),public,parameter :: errRunTime = -1
  integer (kind=ModAEM_Integer),public,parameter :: errFatal = -2
  integer (kind=ModAEM_Integer),public,parameter :: errNonFatal = -3
  integer (kind=ModAEM_Integer),public,parameter :: errInternal = -4
  integer (kind=ModAEM_Integer),public,parameter :: errAssertion = -5
  integer (kind=ModAEM_Integer),public,parameter :: errInvalidDirective    = -101
  integer (kind=ModAEM_Integer),public,parameter :: errIllegalValue        = -102
  integer (kind=ModAEM_Integer),public,parameter :: errAquiferNotSet       = -103
  integer (kind=ModAEM_Integer),public,parameter :: errNoSolution          = -104
  integer (kind=ModAEM_Integer),public,parameter :: errNotAllocated        = -105
  integer (kind=ModAEM_Integer),public,parameter :: errPreviouslyAllocated = -106
  integer (kind=ModAEM_Integer),public,parameter :: errSpaceExhausted      = -107
  integer (kind=ModAEM_Integer),public,parameter :: errUnexpectedEOF       = -108
  integer (kind=ModAEM_Integer),public,parameter :: errUnexpectedData      = -109
  integer (kind=ModAEM_Integer),public,parameter :: errMissingArgument     = -110
  integer (kind=ModAEM_Integer),public,parameter :: errNoData              = -111
  integer (kind=ModAEM_Integer),public,parameter :: errNotImplemented      = -112
  integer (kind=ModAEM_Integer),public,parameter :: errNoOutputFile        = -113


  ! Error messages
  type,public :: ERRORMSG
    integer (kind=ModAEM_Integer) :: Code
    character (len=80) :: Message
  end type ERRORMSG 
  !
  type (ERRORMSG),public,parameter,dimension(12) :: msgErrorMessages = &
    (/ ERRORMSG(errInvalidDirective,    "Invalid program directive"), &
       ERRORMSG(errIllegalValue,        "Illegal data value"), &
       ERRORMSG(errAquiferNotSet,       "Aquifer parameters are not set"), &
       ERRORMSG(errNoSolution,          "No model solution is present"), &
       ERRORMSG(errNotAllocated,        "Array has not been allocated"), &
       ERRORMSG(errPreviouslyAllocated, "Array has already been allocated"), &
       ERRORMSG(errSpaceExhausted,      "Buffer space exhausted"), &
       ERRORMSG(errUnexpectedEOF,       "Unexpected end-of-file detected"), &
       ERRORMSG(errUnexpectedData,      "Unexpected data record found"), &
       ERRORMSG(errMissingArgument,     "Missing argument in parameter list"), &
       ERRORMSG(errNoData,              "No data available"), &
       ERRORMSG(errAssertion,           "Assertion error") &
    /)

  ! Constants for matrix generator
  ! Equation types
  integer (kind=ModAEM_Integer),public,parameter :: EQN_HEAD=0
  integer (kind=ModAEM_Integer),public,parameter :: EQN_FLOW=1
  integer (kind=ModAEM_Integer),public,parameter :: EQN_INHO=2
  integer (kind=ModAEM_Integer),public,parameter :: EQN_DISCHARGE=3
  integer (kind=ModAEM_Integer),public,parameter :: EQN_RECHARGE=4
  integer (kind=ModAEM_Integer),public,parameter :: EQN_CONTINUITY=5
  integer (kind=ModAEM_Integer),public,parameter :: EQN_POTENTIALDIFF=6
  integer (kind=ModAEM_Integer),public,parameter :: EQN_TOTALFLOW=7
  ! Element types
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_AQU=1
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_WL0=2
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_SD0=3
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_LS0=4
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_LS1=5
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_HB0=6
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_IN0=7
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_RCH=8
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_WL1=9
  integer (kind=ModAEM_Integer),public,parameter :: ELEM_LS2=10
  ! Influence Function types
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_P = 1
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_W = 2
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_F = 3
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_G = 4
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_Q = 5
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_J = 6
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_D = 7
  integer (kind=ModAEM_Integer),public,parameter :: INFLUENCE_Z = 8
  ! Size and other inquiry options for XXX_GetInfo()
  integer (kind=ModAEM_Integer),public,parameter :: SIZE_FWL = 1
  integer (kind=ModAEM_Integer),public,parameter :: SIZE_FDP = 2
  integer (kind=ModAEM_Integer),public,parameter :: SIZE_FPD = 3
  integer (kind=ModAEM_Integer),public,parameter :: SIZE_UNKNOWNS = 4
  integer (kind=ModAEM_Integer),public,parameter :: SIZE_EQUATIONS = 5
  integer (kind=ModAEM_Integer),public,parameter :: INFO_REGENERATE = 6

  ! Computational constants
  real (kind=ModAEM_Real),public,parameter :: rPI = 3.14159265359_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rTINY = 1.0e-20_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rZERO = 0.0_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rHALF = 0.5_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rONE_FOURTH = 0.25_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rTHREE_FOURTHS = 0.75_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rONE_EIGHTH = 0.125_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rONE = 1.0_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rTWO = 2.0_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rFOUR = 4.0_ModAEM_Real
  real (kind=ModAEM_Real),public,parameter :: rEIGHT = 8.0_ModAEM_Real
  complex (kind=ModAEM_Real),public,parameter :: cZERO = (0.0_ModAEM_Real,0.0_ModAEM_Real)
  complex (kind=ModAEM_Real),public,parameter :: cONE = (1.0_ModAEM_Real,0.0_ModAEM_Real)
  complex (kind=ModAEM_Real),public,parameter :: cI = (0.0_ModAEM_Real,1.0_ModAEM_Real)
  ! Tolerance near a vertex along a linear feature for point moves
  real (kind=ModAEM_Real),public,parameter :: rVERTEXTOL = 0.00001_ModAEM_Real

! [No code is present in this module]

end module u_constants

