! Description:
!> @file
!!   Initialise look-up tables for bpr calculations.
!
!> @brief
!!   Initialise look-up tables for bpr (back-scattering
!!   parameter) calculations.
!!
!! @details
!!   This subroutine is used when calculating bpr values
!!   from phase functions to pass into RTTOV in the
!!   rttov_opt_param structure for IR cloud/aerosol
!!   scattering calculations.
!!
!!   This should be called once before calling rttov_bpr_calc.
!!   If the phase function angle array changes you must call
!!   rttov_bpr_dealloc and then call this subroutine again.
!!
!! @param[out]  err       status on exit
!! @param[in]   phangle   array of angles on which phase functions are defined
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
SUBROUTINE RTTOV_BPR_INIT( ERR, PHANGLE )

!INTF_OFF
  USE RTTOV_CONST, ONLY : DEG2RAD
  USE RTTOV_BPR_MOD
#include "throw.h"
!INTF_ON
  USE PARKIND1, ONLY : JPRB, JPIM

  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  REAL(KIND=JPRB),    INTENT(IN)  :: PHANGLE(:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=JPIM) :: IANG, I
  REAL(KIND=JPRB) :: V

  TRY
  NPHANGLE = SIZE(PHANGLE)

  ALLOCATE(ARX(NPHANGLE), STAT = ERR)
  THROWM(err.NE.0,"Allocation of ARX failed")
  ALLOCATE(XARR0(NPHANGLE), STAT = ERR)
  THROWM(err.NE.0,"Allocation of XARR0 failed")
  ALLOCATE(CXARR0(NPHANGLE), STAT = ERR)
  THROWM(err.NE.0,"Allocation of CXARR0 failed")
  ALLOCATE(MUX(NPHANGLE), STAT = ERR)
  THROWM(err.NE.0,"Allocation of MUX failed")

  IF( PHANGLE(2_JPIM) .LE. PHANGLE(1_JPIM) ) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"Phase angles should be increasing order")
  ENDIF

  DO IANG = 1_JPIM, NPHANGLE
    IF(PHANGLE(IANG) .GE. 90._JPRB) THEN
       PHANGIND90 = IANG
       EXIT
    ENDIF
  ENDDO

  ! ordering cosinus from -1 to +1;  so increasing sec(mu)
  DO IANG = 1_JPIM, NPHANGLE
    ARX(IANG) = COS(PHANGLE(NPHANGLE+1_JPIM - IANG) * DEG2RAD)
  ENDDO

  DO IANG = 0, 180*VALS_PER_DEG
    V = IANG*1.0_JPRB/VALS_PER_DEG * DEG2RAD
    XE(IANG) = COS(V)
  ENDDO

  XARR0(:)  = PHANGLE(:) * DEG2RAD
  CXARR0(:) = COS(XARR0(:))
  CPHI      = COS( (/ (I, I=1,360) /)  * DEG2RAD)
  MUX       = SQRT(1._JPRB - CXARR0(:)**2_JPIM)
  !TACOS      = ACOS( (/ (I, I=0,1000000) /) /1000000.0_JPRB)

  ! arcosinus
  ! pour conserver la precision de 0.01 degre il faut une precision de :
  ! 10-4 cos=0.640000, 0.000000  ~angle 50.0, 90.0 degree
  ! 10-5 cos=0.992500, 0.640000  ~angle 7.00, 50.0 degree
  ! 10-6 cos=0.999920, 0.992500  ~angle 0.70, 7.00 degree
  ! 10-7 cos=0.999999, 0.999920  ~angle 0.08, 0.70 degree
  ! 10-8 cos=1.000000, 0.999999  ~angle 0.00, 0.08 degree
  TACOS1      = ACOS( (/ (I, I=       0,     6400) /)     /10000.0_JPRB)
  TACOS2      = ACOS( (/ (I, I=   64000,    99250) /)    /100000.0_JPRB)
  TACOS3      = ACOS( (/ (I, I=  992500,   999920) /)   /1000000.0_JPRB)
  TACOS4      = ACOS( (/ (I, I= 9999200,  9999990) /)  /10000000.0_JPRB)
  TACOS5      = ACOS( (/ (I, I=99999900,100000000) /) /100000000.0_JPRB)

  OFFACOS1 = 0.0_JPRB
  OFFACOS2 = 0.640000_JPRB
  OFFACOS3 = 0.9925000_JPRB
  OFFACOS4 = 0.9999200_JPRB
  OFFACOS5 = 0.9999990_JPRB

  PHASE_INIT = .TRUE.

  CATCH
END SUBROUTINE RTTOV_BPR_INIT
