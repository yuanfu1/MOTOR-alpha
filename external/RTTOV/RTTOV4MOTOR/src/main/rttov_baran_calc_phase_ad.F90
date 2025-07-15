! Description:
!> @file
!!   AD of phase function for ice cloud from Baran parameterisation
!
!> @brief
!!   AD of phase function for ice cloud from Baran parameterisation
!!
!! @param[in]     g                asymmetry parameter
!! @param[in,out] g_ad             asymmetry parameter increment
!! @param[in]     phangle          angles on which to calculate phase function (degrees)
!! @param[in,out] phase_ad         phase function increments
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
SUBROUTINE RTTOV_BARAN_CALC_PHASE_AD( G, G_AD, PHANGLE, PHASE_AD)

  USE PARKIND1, ONLY : JPRB
!INTF_OFF
  USE PARKIND1, ONLY : JPIM
  USE RTTOV_CONST, ONLY : DEG2RAD
!INTF_ON
  IMPLICIT NONE

  REAL(KIND=JPRB),INTENT(IN)   :: G
  REAL(KIND=JPRB),INTENT(INOUT):: G_AD
  REAL(KIND=JPRB),INTENT(IN)   :: PHANGLE(:)
  REAL(KIND=JPRB),INTENT(INOUT):: PHASE_AD(:)
!INTF_END

  REAL(KIND=JPRB) :: v  ! angle (rd)
  REAL(KIND=JPRB) :: AA, AA_AD
  REAL(KIND=JPRB) :: BB, BB_AD
  REAL(KIND=JPRB) :: P_AD
  REAL(KIND=JPRB) :: X
  REAL(KIND=JPRB) :: ALPHA, ALPHA_AD
  REAL(KIND=JPRB) :: NORM,  NORM_AD
  REAL(KIND=JPRB) :: BETA
  REAL(KIND=JPRB) :: DD,   DD_AD

  INTEGER(KIND=JPIM) :: NPHANGLE
  INTEGER(KIND=JPIM) :: J

  NPHANGLE = SIZE(PHANGLE)
  AA    = 1.-G*G

  AA_AD = 0.0_JPRB
  BB_AD = 0.0_JPRB
  P_AD  = 0.0_JPRB
  DD_AD = 0.0_JPRB
  ALPHA_AD = 0.0_JPRB
  NORM_AD  = 0.0_JPRB

  IF ( ( G .LT. 0.2_JPRB ) .AND. ( G .GE. 0.0_JPRB ) ) THEN
    DO J=NPHANGLE,1_JPIM, -1_JPIM
      V=PHANGLE(J)
      X=V* DEG2RAD

      BB    = (1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB

      P_AD  = P_AD + PHASE_AD(J)
      PHASE_AD(J) = 0.0_JPRB

      AA_AD = AA_AD + P_AD / BB
      BB_AD = BB_AD + (-AA/(BB**2_JPIM))*P_AD
      P_AD  = 0.0_JPRB

      G_AD =  G_AD + 1.5_JPRB * (1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**0.5_JPRB *&
           & ( 2.0_JPRB*G - 2.0_JPRB*COS(X)) * BB_AD
      BB_AD = 0.0_JPRB
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

    DO J=NPHANGLE,1_JPIM, -1_JPIM
      V=PHANGLE(J)
      X=V* DEG2RAD

      P_AD  = P_AD + PHASE_AD(J)
      PHASE_AD(J) = 0.0_JPRB

      IF ( V .LT. 54.8_JPRB ) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB

        ALPHA_AD = ALPHA_AD + (AA/BB)*COS(X)*P_AD
        AA_AD = AA_AD +  P_AD / BB * COS(X)*ALPHA
        BB_AD = BB_AD - (AA * P_AD / BB**2_JPIM )* COS(X)*ALPHA
        P_AD  = 0.0_JPRB

        G_AD =  G_AD + &
           & 1.5_JPRB * (1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**0.5_JPRB *&
           & ( G - COS(X))* ( 2.0_JPRB*BB_AD)
        BB_AD = 0.0_JPRB

      ENDIF
      IF ( V .GE. 54.8_JPRB ) THEN
        BB=(1.0_JPRB+G*G-1.8_JPRB*G*COS(X)*SIN(X))**1.5_JPRB

        AA_AD = AA_AD +  P_AD / BB
        BB_AD = BB_AD - (AA * P_AD / BB**2_JPIM )
        P_AD  = 0.0_JPRB

        G_AD =  G_AD + &
           & 1.5_JPRB * (1.0_JPRB+G*G-1.8_JPRB*G*COS(X)*SIN(X))**0.5_JPRB *&
           & ( 2.0_JPRB*G -1.8_JPRB*COS(X)*SIN(X) )*BB_AD
        BB_AD = 0.0_JPRB

      ENDIF
    ENDDO

    IF ( ( G .LT. 0.7_JPRB ) .AND. ( G .GE. 0.6_JPRB ) ) THEN
      G_AD = G_AD - (ALPHA_AD*ALPHA) / (2.0_JPRB*G)
      ALPHA_AD = 0.0_JPRB
    ENDIF
    IF ( ( G .LT. 0.6_JPRB ) .AND. ( G .GE. 0.45_JPRB ) ) THEN
      G_AD = G_AD - (ALPHA_AD*ALPHA) / (2.0_JPRB*G)
      ALPHA_AD = 0.0_JPRB
    ENDIF
    IF ( ( G .LT. 0.45_JPRB ) .AND. ( G .GE. 0.3_JPRB ) ) THEN
      G_AD = G_AD - (ALPHA_AD*ALPHA) / (2.0_JPRB*G)
      ALPHA_AD = 0.0_JPRB
    ENDIF
    IF ( ( G .LT. 0.3_JPRB ) .AND. ( G .GE. 0.2_JPRB ) ) THEN
      G_AD = G_AD + (ALPHA_AD*ALPHA) / (2.0_JPRB*(1.0_JPRB-G))
      ALPHA_AD = 0.0_JPRB
    ENDIF
  ENDIF



  IF ( G .GE. 0.7_JPRB ) THEN
    IF ( ( G .LE. 0.8_JPRB ) .AND. ( G .GE. 0.7_JPRB ) ) THEN
      NORM=0.1481E+03-0.2025E+03*G+0.4949E+02*G*G
    ENDIF
    IF ( ( G .LE. 0.9_JPRB ) .AND. ( G .GT. 0.8_JPRB ) ) THEN
      NORM=0.2771E+03-0.5102E+03*G+0.2329E+03*G*G
    ENDIF
    IF ( G .GT. 0.90_JPRB ) THEN
      NORM=0.4219E+03-0.8271E+03*G+0.4063E+03*G*G
    ENDIF
    ALPHA=NORM/SQRT(G)

    DO J=NPHANGLE,1_JPIM, -1_JPIM

      V=PHANGLE(J)
      X=V* DEG2RAD

      P_AD = P_AD + PHASE_AD(J)
      PHASE_AD(J) = 0.0_JPRB

      IF ( V .LE. 3.0_JPRB) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**1.5_JPRB

        ALPHA_AD = ALPHA_AD + (AA/BB)*P_AD*(COS(X))**128.0_JPRB
        AA_AD = AA_AD + (COS(X))**128.0_JPRB * ( P_AD / BB ) * ALPHA
        BB_AD = BB_AD - (COS(X))**128.0_JPRB * (( AA * P_AD ) / BB**2_JPIM ) * ALPHA
        P_AD = 0.0_JPRB

        G_AD =  G_AD + 1.5_JPRB * (1.0_JPRB+G*G-2.0_JPRB*G*COS(X))**0.5_JPRB *&
                & ( G - COS(X)) * 2.0_JPRB*BB_AD
        BB_AD = 0.0_JPRB

      ENDIF
      IF ( ( V .GT. 3.0_JPRB ) .AND. ( V .LT. 30.0_JPRB ) ) THEN
        BB=(1.0_JPRB+G*G-2.0_JPRB*G*COS(1.3_JPRB*X))**1.2_JPRB

        AA_AD = AA_AD + ( P_AD / BB ) * COS(X)
        BB_AD = BB_AD - (( AA * P_AD ) / BB**2_JPIM ) * COS(X)
        P_AD = 0.0_JPRB

        G_AD =  G_AD + 1.2_JPRB * (1.0_JPRB+G*G-2.0_JPRB*G*COS(1.3_JPRB*X))**0.2_JPRB *&
                & ( G - COS(1.3_JPRB*X)) * 2.0_JPRB*BB_AD
        BB_AD = 0.0_JPRB

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

        AA_AD = AA_AD +  COS(X) * P_AD / BB
        BB_AD = BB_AD - COS(X) * (AA * P_AD / BB**2_JPIM )
        P_AD  = 0.0_JPRB

        G_AD =  G_AD + BETA * (1.0_JPRB+G*G-2.0_JPRB*G*COS(DD*X))**(BETA-1.0_JPRB) *&
            & (G - COS(DD*X)) * 2.0_JPRB*BB_AD
        DD_AD = DD_AD + BETA * (1.0_JPRB+G*G-2.0_JPRB*G*COS(DD*X))**(BETA-1.0_JPRB) *&
            & ( 2.0_JPRB*G*SIN(DD*X)*BB_AD*X )
        BB_AD = 0.0_JPRB

        G_AD =  G_AD + (1.0_JPRB - 1.0_JPRB/4.6_JPRB)*DD_AD
        DD_AD = 0.0_JPRB
      ENDIF

      IF ( ( V .GE. 54.8_JPRB ) .AND. ( V .LE. 95.0_JPRB ) ) THEN
        BB=(1.0_JPRB+G*G-1.5_JPRB*G*COS(X)*SIN(X))**1.5_JPRB

        AA_AD = AA_AD +  P_AD / BB
        BB_AD = BB_AD - (AA * P_AD / BB**2_JPIM )
        P_AD  = 0.0_JPRB

        G_AD = G_AD + 1.5_JPRB * (1.0_JPRB+G*G-1.5_JPRB*G*COS(X)*SIN(X))**0.5_JPRB *&
              & ( 2.0_JPRB*G -1.5_JPRB*COS(X)*SIN(X)) * BB_AD
        BB_AD = 0.0_JPRB

      ENDIF

    ENDDO ! ANGLE

    NORM_AD = NORM_AD + (ALPHA_AD*SQRT(G) ) / G
    G_AD = G_AD - (0.5_JPRB*ALPHA*ALPHA_AD ) / G
    ALPHA_AD = 0.0_JPRB

    IF ( ( G .LE. 0.8_JPRB ) .AND. ( G .GE. 0.7_JPRB ) ) THEN
      G_AD = G_AD + (-0.2025E+03 + 2.0_JPRB*0.4949E+02*G) * NORM_AD
    ENDIF
    IF ( ( G .LE. 0.9_JPRB ) .AND. ( G .GT. 0.8_JPRB ) ) THEN
      G_AD = G_AD + (-0.5102E+03 + 2.0_JPRB*0.2329E+03*G) * NORM_AD

    ENDIF
    IF ( G .GT. 0.90_JPRB ) THEN
      G_AD = G_AD + (-0.8271E+03 + 2.0_JPRB*0.4063E+03*G) * NORM_AD

    ENDIF
    NORM_AD = 0.0_JPRB

  ENDIF ! G

  G_AD = G_AD - 2.0_JPRB * AA_AD * G
  AA_AD = 0.0_JPRB

END SUBROUTINE
