! Description:
!> @file
!!   Compute phase function for ice cloud from Baran parameterisation
!
!> @brief
!!   Compute phase function for ice cloud from Baran parameterisation
!!
!! @param[in]     g                asymmetry parameter
!! @param[in]     phangle          angles on which to calculate phase function (degrees)
!! @param[out]    phase            phase function
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE RTTOV_BARAN_CALC_PHASE( G, PHANGLE, PHASE)

  USE PARKIND1, ONLY : JPRB
!INTF_OFF
  USE PARKIND1, ONLY : JPIM
  USE RTTOV_CONST, ONLY : DEG2RAD
!INTF_ON
  IMPLICIT NONE

  REAL(KIND=JPRB),INTENT(IN)  :: G
  REAL(KIND=JPRB),INTENT(IN)  :: PHANGLE(:)
  REAL(KIND=JPRB),INTENT(OUT) :: PHASE(:)
!INTF_END

  REAL(KIND=JPRB) :: v  ! angle (rd)
  REAL(KIND=JPRB) :: AA
  REAL(KIND=JPRB) :: BB
  REAL(KIND=JPRB) :: P
  REAL(KIND=JPRB) :: X
  REAL(KIND=JPRB) :: ALPHA
  REAL(KIND=JPRB) :: NORM
  REAL(KIND=JPRB) :: BETA
  REAL(KIND=JPRB) :: DD

  INTEGER(KIND=JPIM) :: NPHANGLE
  INTEGER(KIND=JPIM) :: J

  NPHANGLE = SIZE(PHANGLE)

  AA=1.-G*G
  
  IF ( ( G .LT. 0.2_JPRB ) .AND. ( G .GE. 0.0_JPRB ) ) THEN
    DO J=1,NPHANGLE
      V=PHANGLE(J)
      X=V* DEG2RAD
      BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB
      P=(AA/BB)
      PHASE(J)=P
    ENDDO
  ENDIF
  
  IF ( ( G .LT. 0.7_JPRB ) .AND. ( G .GE. 0.2_JPRB ) ) THEN
    IF ( ( G .LT. 0.7_JPRB ) .AND. ( G .GE. 0.6_JPRB ) ) THEN
      ALPHA=(1.0_JPRB/(SQRT(1.095_JPRB*G)))
    ENDIF
    IF ( ( G .LT. 0.6_JPRB ) .AND. ( G .GE. 0.45_JPRB ) ) THEN
      ALPHA=(1.0_JPRB/(SQRT(1.23_JPRB*G)))
    ENDIF
    IF ( ( G .LT. 0.45_JPRB ) .AND. ( G .GE. 0.3_JPRB ) ) THEN
      ALPHA=(1.0_JPRB/(SQRT(1.5_JPRB*G)))   
    ENDIF
    IF ( ( G .LT. 0.3_JPRB ) .AND. ( G .GE. 0.2_JPRB ) ) THEN
      ALPHA=((1.0_JPRB/(SQRT((1.0_JPRB-G)))))*1.25_JPRB
    ENDIF   
    DO J=1_JPIM,NPHANGLE
      V=PHANGLE(J)
      X=V* DEG2RAD
      IF ( V .LT. 54.8_JPRB ) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB
        P=(AA/BB)*COS(X)*ALPHA
        PHASE(J)=P
      ENDIF
      IF ( V .GE. 54.8_JPRB ) THEN
        BB=(1.0_JPRB+G*G-1.8_JPRB*G*COS(X)*SIN(X))**1.5_JPRB
        P=AA/BB
        PHASE(J)=P
      ENDIF
    ENDDO
  ENDIF
  
  IF ( G .GE. 0.7_JPRB ) THEN
    IF ( ( G .LE. 0.8_JPRB ) .AND. ( G .GE. 0.7_JPRB ) ) THEN
      NORM=0.1481E+03-0.2025E+03*G+0.4949E+02*G*G
      ALPHA=NORM/SQRT(G)
    ENDIF
    IF ( ( G .LE. 0.9_JPRB ) .AND. ( G .GT. 0.8_JPRB ) ) THEN
      NORM=0.2771E+03-0.5102E+03*G+0.2329E+03*G*G
      ALPHA=NORM/SQRT(G)
    ENDIF
    IF ( G .GT. 0.90_JPRB ) THEN
      NORM=0.4219E+03-0.8271E+03*G+0.4063E+03*G*G
      ALPHA=NORM/SQRT(G)
    ENDIF
    DO J = 1_JPIM,NPHANGLE
      V=PHANGLE(J)
      X=V* DEG2RAD
      IF ( V .LE. 3.0_JPRB) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB
        P=(AA/BB)*(COS(X))**128.0_JPRB*ALPHA
        PHASE(J)=P
      ENDIF
      IF ( ( V .GT. 3.0_JPRB ) .AND. ( V .LT. 30.0_JPRB ) ) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(1.3_JPRB*X))**1.2_JPRB
        P=(AA/BB)*COS(X)
        PHASE(J)=P
      ENDIF
      IF ( ( V .GE. 30.0_JPRB ) .AND. ( V .LT. 54.8_JPRB ) ) THEN
        DD=((1.0_JPRB-G)/4.6_JPRB)+G
        IF ( ( G .GE. 0.7_JPRB ) .AND. ( G .LT. 0.9_JPRB ) ) THEN
          BETA=0.68_JPRB
        ENDIF
        IF ( G .GE. 0.9_JPRB ) THEN
          BETA=0.71_JPRB
        ENDIF
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(DD*X))**BETA
        P=(AA/BB)*COS(X)
        PHASE(J)=P
      ENDIF
      IF ( ( V .GE. 54.8_JPRB ) .AND. ( V .LE. 95.0_JPRB ) ) THEN
        BB=(1.0_JPRB+G*G-1.5_JPRB*G*COS(X)*SIN(X))**1.5_JPRB
        P=(AA/BB)
        PHASE(J)=P
      ENDIF
      IF ( V .GT. 95.0_JPRB ) THEN
        PHASE(J)=P ! take last valid value
      ENDIF
    ENDDO ! ANGLE
  ENDIF ! G

END SUBROUTINE
