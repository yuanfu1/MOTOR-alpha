! Description:
!> @file
!!   Calculate bpr for a given phase function.
!
!> @brief
!!   Calculate bpr (back-scattering parameter) for a given phase function.
!!
!! @details
!!   Calculate the pseudo back-scattering parameter bpr using the
!!   look-up tables calculated by rttov_bpr_init and stored in
!!   rttov_scattering_mod.
!!
!!   The rttov_bpr_init subroutine prepares the look-up tables to speed up
!!   the calculations. These include:
!!   - all arrays related to the phase angles
!!   - necessary angles for the integration over 1-360
!!   - arc-cosine function keeping 0.01 degree precision.
!!
!!   This improves the speed by a factor of 10 over the original code
!!   commented out below: in the original subroutine most of time is spent
!!   in the interpolation routine which is called in the inner loop
!!
!! @param[out]  err       status on exit
!! @param[in]   pha       phase function values at angles in phangle array
!! @param[in]   phangle   array of angles on which phase function is defined
!! @param[out]  bpr       calculated back-scattering parameter
!! @param[in]   nthreads  number of OpenMP threads to use, optional, default 1
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
SUBROUTINE RTTOV_BPR_CALC( &
             & ERR,     &
             & PHA,     &
             & PHANGLE, &
             & BPR,     &
             & NTHREADS)

  USE PARKIND1, ONLY : JPRB, JPIM
!INTF_OFF
#include "throw.h"
  USE PARKIND1, ONLY : JPLM
  USE RTTOV_CONST, ONLY : PI, DEG2RAD
  USE RTTOV_SCATTERING_MOD, ONLY : INTER, INTEGRATE
  USE RTTOV_BPR_MOD
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(OUT)          :: ERR
  REAL(KIND=JPRB),    INTENT(IN)           :: PHA(:)
  REAL(KIND=JPRB),    INTENT(IN)           :: PHANGLE(SIZE(PHA))
  REAL(KIND=JPRB),    INTENT(OUT)          :: BPR
  INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: NTHREADS
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=JPIM) :: IANG,I,J,K
  REAL(KIND=JPRB) :: ANG,ANG1,ARY(NPHANGLE)
  REAL(KIND=JPRB) :: PHAS(NPHANGLE,NPHANGLE)
  REAL(KIND=JPRB) :: MU,MU1,COST,PHASINT
  REAL(KIND=JPRB) :: HR,ECOEF,ERROR
  REAL(KIND=JPRB) :: XARR1(360_JPIM),XARR2(PHANGIND90),XARR3(NPHANGLE-PHANGIND90+1_JPIM)
  REAL(KIND=JPRB) :: YARR1(360_JPIM),YARR2(PHANGIND90),YARR3(NPHANGLE-PHANGIND90+1_JPIM)
  REAL(KIND=JPRB) :: AV

  REAL(KIND=JPRB) :: YE(0:180*VALS_PER_DEG)
  INTEGER(KIND=JPIM) :: IV, NYE
  REAL(KIND=JPRB) :: V, MU2, MU3, MU4
  INTEGER(KIND=JPIM) :: NTHREADS1
  LOGICAL(KIND=JPLM) :: ERRORFLAG

  TRY

  ! The arrays declared above require NPHANGLE and PHANGIND90 to be defined
  ! on entry which requires rttov_bpr_init to have been called already
  IF (.NOT. PHASE_INIT) THEN
    ERR = ERRORSTATUS_FATAL
    THROWM(ERR.NE.0_JPIM, "Error: call rttov_bpr_init first")
  ENDIF

  NTHREADS1 = 1
  IF (PRESENT(NTHREADS)) THEN
    IF (NTHREADS > 1) NTHREADS1 = NTHREADS
  ENDIF

  DO IANG = 1_JPIM, NPHANGLE
    ARY(IANG) = PHA(NPHANGLE+1_JPIM - IANG)
  ENDDO

  NYE =  SIZE(YE)
  DO IANG = 0, 180*VALS_PER_DEG
    CALL INTER(NPHANGLE, NPHANGLE, 2_JPIM, XE(IANG), ARX, ARY, YE(IANG), HR)
  ENDDO

  ERRORFLAG = .FALSE.
  ! CALCULATE AZIMUTHALLY-AVERAGED PHASE FUNCTION

!$OMP PARALLEL DO NUM_THREADS(NTHREADS1) DEFAULT(PRIVATE) SCHEDULE(DYNAMIC) &
!$OMP             SHARED(NPHANGLE, CXARR0, MUX, CPHI, OFFACOS1, OFFACOS2, OFFACOS3, OFFACOS4, OFFACOS5, &
!$OMP                    TACOS1, TACOS2, TACOS3, TACOS4, TACOS5, NYE, XE, YE, PHAS, ERRORFLAG)
  DO I = 1_JPIM, NPHANGLE
    MU  = CXARR0(I)
    MU4 = MUX(I)
    DO J = 1_JPIM, NPHANGLE
      MU1 = CXARR0(J)
      MU2 = MU * MU1
      MU3 = MU4 * MUX(J)
      DO K = 1_JPIM, 360_JPIM

        COST = MU2 + MU3* CPHI(K)

!              ! exact calculation (as of original subroutine):
!              !V = ACOS(COST)

        ! approximation by arc-cosinus tables
        AV = ABS(COST)
        IF( AV .LT. OFFACOS2 ) THEN
          V = TACOS1( NINT((AV-OFFACOS1)*10000.0_JPRB) )
        ELSEIF( AV .LT. OFFACOS3 ) THEN
          V = TACOS2( NINT((AV-OFFACOS2)*100000.0_JPRB) )
        ELSEIF( AV .LT. OFFACOS4 ) THEN
          V = TACOS3( NINT((AV-OFFACOS3)*1000000.0_JPRB) )
        ELSEIF( AV .LT. OFFACOS5 ) THEN
          V = TACOS4( NINT((AV-OFFACOS4)*10000000.0_JPRB) )
        ELSEIF( AV .LE. 1.0_JPRB ) THEN
          V = TACOS5( NINT((AV-OFFACOS5)*100000000.0_JPRB) )
        ELSE
          V = 0.0_JPRB
        ENDIF
!              V = TACOS( NINT(AV*1000000) )
        IF( COST .LE. 0.0_JPRB) V = PI -V
        ! approximation by arc-cosinus tables; end

        IV = INT(V * VALS_PER_DEG / DEG2RAD)
        IF( IV .LE. NYE-2_JPIM) THEN
          PHASINT = YE(IV) + (COST - XE(IV))*(YE(IV+1)-YE(IV)) / &
                                & (XE(IV+1) - XE(IV))
        ELSE
          PHASINT = YE(IV)
        ENDIF
        XARR1(K) = K * DEG2RAD
        YARR1(K) = PHASINT
      ENDDO
      CALL INTEGRATE(360_JPIM,XARR1,YARR1,ECOEF,ERROR,ERR)
      IF (ERR.NE.0_JPIM) ERRORFLAG = .TRUE.
      PHAS(I,J) = ECOEF / (2._JPRB * PI)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  THROW(ERRORFLAG)

  ! perform double integration as in eqn 32 in tm474 to obtain bpr
  DO I = 1, PHANGIND90
    ANG      = XARR0(I)
    XARR2(I) = ANG
    DO J = PHANGIND90, NPHANGLE
      ANG1 = XARR0(J)
      XARR3(J - (PHANGIND90-1)) = ANG1
      YARR3(J - (PHANGIND90-1)) = PHAS(I,J) * SIN(ANG1)
    ENDDO
    CALL INTEGRATE(NPHANGLE-PHANGIND90+1_JPIM,XARR3(1:NPHANGLE-PHANGIND90+1),YARR3(1:NPHANGLE-PHANGIND90+1),ECOEF,ERROR,ERR)
    THROW(ERR.NE.0_JPIM)
    YARR2(I) = ECOEF * SIN(ANG)
  ENDDO

  CALL INTEGRATE(PHANGIND90,XARR2(1:PHANGIND90),YARR2(1:PHANGIND90),ECOEF,ERROR,ERR)
  THROW(ERR.NE.0_JPIM)
  BPR = ECOEF / 2._JPRB

  CATCH
END SUBROUTINE RTTOV_BPR_CALC

! SUBROUTINE rttov_calc_bpr0( &
!              & pha, &
!              & phangle, &
!              & bpr) 
! ! original subroutine provided by J. Vidot based on rttov_mie_params
! 
! use parkind1, only : jprb, jpim
! 
! IMPLICIT NONE
! integer(kind=jpim), parameter :: nphangle = 208_JPIM                  , &
!                                      phangind90 = 118_JPIM
! real(kind=jprb), parameter    :: deg2rad = 0.0174532925199432957_JPRB , &
!                                      pi  = 3.1415926535897932385_JPRB
! Real(kind=jprb),intent(in)  :: pha(nphangle)
! Real(kind=jprb),intent(in)  :: phangle(nphangle)
! Real(kind=jprb),intent(out) :: bpr
! 
! Integer(kind=jpim) :: iang,i,j,k,ifail
! Real(kind=jprb) :: ang,ang1,arx(nphangle),ary(nphangle)
! Real(kind=jprb) :: xarr0(nphangle),intg(nphangle)
! Real(kind=jprb) :: phas(nphangle,nphangle)
! Real(kind=jprb) :: mu,mu1,phi,cosT,phasint
! Real(kind=jprb) :: hr,ecoef,error
! Real(kind=jprb) :: xarr1(360_JPIM),xarr2(phangind90),xarr3(nphangle-phangind90+1_JPIM)
! Real(kind=jprb) :: yarr1(360_JPIM),yarr2(phangind90),yarr3(nphangle-phangind90+1_JPIM)
! 
!         do iang = 1_JPIM, nphangle
!           arx(iang) = cos(phangle(nphangle+1_JPIM - iang) * deg2rad)
!           ary(iang) = pha(nphangle+1_JPIM - iang)
!         enddo
! 
!         ! Calculate azimuthally-averaged phase function
!         xarr0(:) = phangle(:) * deg2rad
!         phas = 0_JPIM
!         do i = 1_JPIM, nphangle
!           ang = xarr0(i)
!           do j = 1_JPIM, nphangle
!             ang1 = xarr0(j)
!               do k = 1_JPIM, 360_JPIM
!                 mu  = cos(ang)
!                 mu1 = cos(ang1)
!                 phi = cos(k * deg2rad)
!                 cosT = mu * mu1 + sqrt(1._jprb - mu**2_JPIM) * sqrt(1._jprb - mu1**2_JPIM) * phi
!                 call INTER(nphangle, nphangle, 2_JPIM, cosT, arx, ary, phasint, HR)
!                 xarr1(k) = k * deg2rad
!                 yarr1(k) = phasint
!               enddo
!             call integrate(360_JPIM,xarr1(1:360),yarr1(1:360),ecoef,error,ifail)
!             phas(i,j) = ecoef / (2._jprb * pi) 
!           enddo
!         enddo
! 
!         ! Perform double integration as in eqn 32 in TM474 to obtain bpr
!         intg(:) = 0
!         do i = 1, phangind90
!           ang = xarr0(i)
!           xarr2(i) = ang
!           do j = phangind90, nphangle
!             ang1 = xarr0(j)
!             xarr3(j - (phangind90-1)) = ang1
!             yarr3(j - (phangind90-1)) = phas(i,j) * sin(ang1)
!           enddo
!           call integrate(nphangle-phangind90+1,xarr3(1:nphangle-phangind90+1),yarr3(1:nphangle-phangind90+1),ecoef,error,ifail)
!           intg(i) = ecoef
!           yarr2(i) = intg(i) * sin(ang)
!         enddo
! 
!         call integrate(phangind90,xarr2(1:phangind90),yarr2(1:phangind90),ecoef,error,ifail)
!         bpr = ecoef / 2._jprb
! 
! END SUBROUTINE
