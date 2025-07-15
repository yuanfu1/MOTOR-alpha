! Description:
!> @file
!!   Stores pre-calculated data for bpr calculation (rttov_bpr_init/calc)
!
!> @brief
!!   Stores pre-calculated data for bpr calculation (rttov_bpr_init/calc)
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE RTTOV_BPR_MOD

USE PARKIND1,    ONLY : JPRB, JPIM

IMPLICIT NONE

PRIVATE
PUBLIC :: NPHANGLE, PHANGIND90, PHASE_INIT, VALS_PER_DEG, &
          ARX, XARR0, CXARR0, MUX, CPHI, XE, &
          TACOS1, TACOS2, TACOS3, TACOS4, TACOS5, &
          OFFACOS1, OFFACOS2, OFFACOS3, OFFACOS4, OFFACOS5

INTEGER(KIND=JPIM) :: NPHANGLE 
INTEGER(KIND=JPIM) :: PHANGIND90 ! INDEX OF PHASE ANGLE 90 DEGREES

LOGICAL :: PHASE_INIT=.FALSE.

INTEGER(KIND=JPIM),PARAMETER :: VALS_PER_DEG=100

REAL(KIND=JPRB), POINTER :: ARX(:)
REAL(KIND=JPRB), POINTER :: XARR0(:)
REAL(KIND=JPRB), POINTER :: CXARR0(:)
REAL(KIND=JPRB), POINTER :: MUX(:)
REAL(KIND=JPRB) :: CPHI(1:360_JPIM)
REAL(KIND=JPRB) :: XE(0:180*VALS_PER_DEG)

REAL(KIND=JPRB) :: TACOS1(0:6400)
REAL(KIND=JPRB) :: TACOS2(0:35250)
REAL(KIND=JPRB) :: TACOS3(0:7420)
REAL(KIND=JPRB) :: TACOS4(0:790)
REAL(KIND=JPRB) :: TACOS5(0:100 )
REAL(KIND=JPRB) :: OFFACOS1, OFFACOS2, OFFACOS3, OFFACOS4, OFFACOS5

END MODULE
