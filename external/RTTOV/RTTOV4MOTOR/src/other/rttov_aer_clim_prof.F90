! Description:
!> @file
!!   Generates climatological aerosol profiles.
!
!> @brief
!!   Generates climatological aerosol profiles.
!!
!! @details
!!   The inputs are profiles of pressure, temperature
!!   and water vapour.
!!
!!   The outputs are 10 aerosol profiles which can be passed
!!   to RTTOV IR aerosol scattering calculations. They correspond
!!   to the following climatological compositions:
!!
!!       \li  J=1  -->Continental clean
!!       \li  J=2  -->Continental average
!!       \li  J=3  -->Continental polluted
!!       \li  J=4  -->Urban
!!       \li  J=5  -->Desert
!!       \li  J=6  -->Maritime clean
!!       \li  J=7  -->Maritime polluted
!!       \li  J=8  -->Maritime tropical
!!       \li  J=9  -->Arctic
!!       \li  J=10 -->Antarctic
!!
!!   Each climatological composition comprises profiles of each
!!   of the pre-defined RTTOV aerosol particle types.
!!
!! @param[in]   p            input pressure profile (hPa)
!! @param[in]   t            input temperature profile on levels in p (K)
!! @param[in]   q            input water vapour profile on levels in p
!! @param[in]   gas_units    input water vapour units: 0=>ppmv dry; 1=>kg/kg; 2=>ppmv moist
!! @param[in]   mmr_aer      output aerosol units: .FALSE.=>cm-3; .TRUE.=>kg/kg
!! @param[in]   levsurf      level number corresponding to the surface (index into p)
!! @param[in]   latitude     latitude of profile (-90 to +90 degrees)
!! @param[in]   elevation    surface elevation (km)
!! @param[in]   scalefactor  factor by which to scale the aerosol profile values
!! @param[out]  aerprof      generated aerosol profiles, dimensions: (nlayers,10,naer_opac)
!!                              where nlayers=SIZE(p)-1, naer_opac is the number of RTTOV aerosol
!!                              particle types in OPAC scaercoef files, units: cm-3
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
SUBROUTINE rttov_aer_clim_prof( &
        P,                      &
        T,                      &
        Q,                      &
        GAS_UNITS,              &
        MMR_AER,                &
        LEVSURF,                &
        LATITUDE,               &
        ELEVATION,              &
        SCALEFACTOR,            &
        AERPROF)

  USE parkind1, ONLY : jprb, jpim, jplm
  USE rttov_const, ONLY: naer_opac
!INTF_OFF
  USE rttov_const, ONLY: &
    flatt,   &
    omega,   &
    eqrad,   &
    grave,   &
    rgc,     &
    mair,    &
    mh2o,    &
    deg2rad, &
    aer_opac_confac, &
    ngases_unit
!INTF_ON
  IMPLICIT NONE

  REAL(jprb),    INTENT(IN)  :: P(:)                            ! Input pressure profile (hPa)
  REAL(jprb),    INTENT(IN)  :: T(SIZE(p))                      ! Input temperature profile (K)
  REAL(jprb),    INTENT(IN)  :: Q(SIZE(p))                      ! Input water vapour profile (ppmv)
  INTEGER(jpim), INTENT(IN)  :: GAS_UNITS                       ! Input q units: 0=>ppmv dry; 1=>kg/kg; 2=>ppmv moist
  LOGICAL(jplm), INTENT(IN)  :: MMR_AER                         ! Output aerosol units: .FALSE.=>cm-3; .TRUE.=>kg/kg
  INTEGER(jpim), INTENT(IN)  :: LEVSURF                         ! Surface level
  REAL(jprb),    INTENT(IN)  :: LATITUDE                        ! Latitude for profiles
  REAL(jprb),    INTENT(IN)  :: ELEVATION                       ! Surface elevation (km)
  REAL(jprb),    INTENT(IN)  :: SCALEFACTOR                     ! Factor by which to scale profiles
  REAL(jprb),    INTENT(OUT) :: AERPROF(SIZE(p)-1,10,naer_opac) ! Output aerosol profiles (nlayer, nclim, naer_opac)
!INTF_END

! Note: input/output profiles conform to RTTOV standard (descending levels).
! Internal profiles are in the reverse direction (similar to rttov_locpat).

  INTEGER(jpim) :: NLEV
  INTEGER(jpim) :: RLEVSURF             ! Surface level on reversed levels
  INTEGER(jpim) :: I, II, J, IR, ILAY, RILAY
  REAL(jprb)    :: RLH
  REAL(jprb)    :: HFACSL1
  REAL(jprb)    :: HFACSL2
  REAL(jprb)    :: HFACFT1
  REAL(jprb)    :: HFACFT2
  REAL(jprb)    :: DFLAT
  REAL(jprb)    :: FAC
  REAL(jprb)    :: ETA
  REAL(jprb)    :: BETA
  REAL(jprb)    :: GRAVL
  REAL(jprb)    :: GRAVH
  REAL(jprb)    :: LATR
  REAL(jprb)    :: NDENSST
  REAL(jprb)    :: REARTH
  REAL(jprb)    :: TEMPP(SIZE(P))
  REAL(jprb)    :: PRES(SIZE(P))
  REAL(jprb)    :: WATER(SIZE(P))
  REAL(jprb)    :: DAIR(SIZE(P))
  REAL(jprb)    :: DMAIR(SIZE(P))
  REAL(jprb)    :: HGPL(SIZE(P))
  REAL(jprb)    :: INT(SIZE(P))
  REAL(jprb)    :: HL(SIZE(P))
  REAL(jprb)    :: LTICK(SIZE(P))
  REAL(jprb)    :: HFAC(SIZE(P))
  REAL(jprb)    :: PLAY(SIZE(P)-1)
  REAL(jprb)    :: TLAY(SIZE(P)-1)
  REAL(jprb)    :: QLAY(SIZE(P)-1)
  REAL(jprb)    :: QMR(SIZE(P))
  REAL(jprb)    :: FACTOR(SIZE(P)-1)
  REAL(jprb)    :: CHEIG   (10,3)
  REAL(jprb)    :: ZFAC    (10,4)
  REAL(jprb)    :: NDENS   (5,10)
  REAL(jprb)    :: NDENSFT (3)
  INTEGER(jpim) :: AFAC    (3)
  INTEGER(jpim) :: CCOMP   (5,10)
  INTEGER(jpim) :: CNUM    (10)

!-----End of header-------------------------------------------------------------

  DATA CHEIG/02._jprb,02._jprb,02._jprb,02._jprb,06._jprb,02._jprb,02._jprb,02._jprb,02._jprb,10._jprb,  &
             12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,12._jprb,  &
             35._jprb,35._jprb,35._jprb,35._jprb,35._jprb,35._jprb,35._jprb,35._jprb,35._jprb,35._jprb/

  DATA ZFAC /8._jprb ,8._jprb ,8._jprb ,8._jprb ,2._jprb ,1._jprb ,1._jprb ,1._jprb ,99._jprb,8._jprb,   &
             8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb ,8._jprb,   &
             99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,  &
             99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb,99._jprb/

  DATA CCOMP/2,1,0,0,0,                                                    &
             2,1,3,0,0,                                                    &
             2,1,3,0,0,                                                    &
             2,1,3,0,0,                                                    &
             2,6,7,8,0,                                                    &
             2,4,5,0,0,                                                    &
             2,4,5,3,0,                                                    &
             2,4,5,0,0,                                                    &
             2,1,4,3,0,                                                    &
             10,4,9,0,0/

  DATA NDENS/2600._jprb , 0.15_jprb   , 0._jprb     , 0._jprb   , 0._jprb,            &
             7000._jprb , 0.4_jprb    , 8300._jprb  , 0._jprb   , 0._jprb,            &
             15700._jprb, 0.6_jprb    , 34300._jprb , 0._jprb   , 0._jprb,            &
             28000._jprb, 1.5_jprb    , 130000._jprb, 0._jprb   , 0._jprb,            &
             2000._jprb , 269.5_jprb  , 30.5_jprb   , 0.142_jprb, 0._jprb,            &
             1500._jprb , 20._jprb    , 3.2e-3_jprb , 0._jprb   , 0._jprb,            &
             3800._jprb , 20._jprb    , 3.2e-3_jprb , 5180._jprb, 0._jprb,            &
             590._jprb  , 10._jprb    , 1.3e-3_jprb , 0._jprb   , 0._jprb,            &
             1300._jprb , 0.01_jprb   , 1.9_jprb    , 5300._jprb, 0._jprb,            &
             42.9_jprb  , 0.47e-1_jprb, 0.53e-2_jprb, 0._jprb   , 0._jprb/

  DATA NDENSFT/438._jprb,0.1241_jprb,292._jprb/

  DATA NDENSST/3._jprb/

  DATA CNUM /2,3,3,3,4,3,4,3,4,3/

!-------------------------------------------------------------------------------


  NLEV = SIZE(P)
  RLEVSURF = NLEV - LEVSURF + 1

  AERPROF(:,:,:)=0.

  IF (GAS_UNITS < 0 .OR. GAS_UNITS > ngases_unit) THEN
    PRINT *, 'Invalid gas_units'
    RETURN
  ENDIF

!---------Compute tickness of the layer----------------------------------------

!-------1.  Caculate the earth's radius at given latitude assuming the---------
!       Earth is an ellipsoid of revolution                                    |
!------------------------------------------------------------------------------

  DFLAT  = (1.0 - FLATT)**2
  LATR   = DEG2RAD*ABS(LATITUDE)
  REARTH = SQRT(EQRAD**2*DFLAT/(SIN(LATR)**2 + DFLAT*COS(LATR)**2))

!-------2.  The valute of earth's gravity at surface at given latitude is------
!       computed using the international gravity formula.                      |
!------------------------------------------------------------------------------

  FAC   = (OMEGA**2*(EQRAD*1000))/(GRAVE)
  BETA  = 5*FAC/2-FLATT-17*FAC*FLATT/14
  ETA   = FLATT*(5*FAC-FLATT)/8
  GRAVL = GRAVE*(1+BETA*(SIN(DEG2RAD*LATITUDE))**2+ETA*(SIN(2*DEG2RAD*LATITUDE))**2)

!-------3.  The value of the gravity as a function of altitude H can be--------
!       expressed using the inverse-square law of gravitation:                 |
!                                                                              |
!       GRAVL(H)=GRAVL*(REARTH/(H+REARTH))**2=                                 |
!                GRAVL*(1-2*H/REARTH+3*H**2/REARTH**2+terms of higher order)   |
!                                                                              |
!       If we eliminate the second and higher order terms we can write:        |
!                                                                              |
!       GRAVL(H)=GRAVL-2*GRAVL*H/REARTH=GRAVL-GRAVH*H                          |
!------------------------------------------------------------------------------

  GRAVH = 2.0E-3*GRAVL/REARTH

!-------4.  Unpack input profile------------------------------------------------

  DO I=1,NLEV
    PRES(I)  = P(NLEV-I+1)
    TEMPP(I) = T(NLEV-I+1)

    ! Q contains input water vapour profile in units of GAS_UNITS
    ! QMR contains Q in units of kg/kg
    ! WATER contains Q on reversed levels in units of ppmv moist
    IF (GAS_UNITS == 2) THEN
      ! ppmv moist => ppmv dry
      QMR(I) =  Q(I) / (1._jprb - 1.E-06_jprb * Q(I))
      ! ppmv dry => kg/kg
      QMR(I) = QMR(I) / (Mair * 1.E06_jprb / Mh2o + QMR(I))

      ! ppmv moist
      WATER(I) = Q(NLEV-I+1)
    ELSEIF (GAS_UNITS == 1) THEN
      ! kg/kg
      QMR(I) = Q(I)

      ! kg/kg => ppmv dry
      WATER(I) = (1.E06_jprb * Mair / Mh2o) * Q(NLEV-I+1) / (1._jprb - Q(NLEV-I+1))
      ! ppmv dry => ppmv moist
      WATER(I) = WATER(I) / (1._jprb + 1.E-06_jprb * WATER(I))
    ELSE
      ! ppmv dry => kg/kg
      QMR(I) = Q(I) / (Mair * 1.E06_jprb / Mh2o + Q(I))

      ! ppmv dry => ppmv moist
      WATER(I) = Q(NLEV-I+1) / (1._jprb + 1.E-06_jprb * Q(NLEV-I+1))
    ENDIF
  ENDDO

!-------5.  Calculate density for dry air (ideal gas)---------------------------

  DO IR=1,NLEV
     DAIR(IR)=1.E+02*PRES(IR)*MAIR/(1000*RGC*TEMPP(IR))
  END DO


!-------6.  The density of dry air is adjusted to account for the presence ----
!       of water vapour by introducing a scaling factor (i.e. the              |
!       temperature is replaced by the virtual temperature).The density        |
!       for moist air is thus obtained. It is assumed that the partial         |
!       pressure for water vapour can be written as PRES*WATER*1.E-6           |
!------------------------------------------------------------------------------

  DO IR=1,NLEV
    DMAIR(IR)= DAIR(IR)*(MAIR*(1.0E6 - WATER(IR))+MH2O*WATER(IR))/(1.0E6*MAIR)
  ENDDO


!-------7.  Set up the height of the surface pressure level---------------------

  HGPL(RLEVSURF) = ELEVATION


!-------8  Compute height of pressure levels:levels above surface. ------------
!       The height of pressure levels H is obtained by integrating             |
!       the hydrostatic equation dPRES=-GRAVL(H)*DMAIR*dH between              |
!       two adjacent pressure levels                                           |
!                                                                              |
!           -P2         - H2                                                   |
!          |           |                                                       |
!          | dP/DMAIR= | GRAVL(H)*dH                                           |
!          |           |                                                       |
!         -  P1       -   H1                                                   |
!                                                                              |
!       The integration id DMAIR is carried out assuming this                  |
!       quantity varies linearly with pressure. The integral is                |
!       then computed applying the trapezoidal formula.                        |
!------------------------------------------------------------------------------

  DO IR = RLEVSURF+1,NLEV
    INT(IR)     =- 0.5*(PRES(IR)-PRES(IR-1))*                          &
                   (1/DMAIR(IR)+1/DMAIR(IR-1))
    INT(IR)     =INT(IR)*1.0E2            ! Integral in units [m2.s-1]
    RLH     =GRAVL/GRAVH
    HL(IR)=1.0E3*HGPL(IR-1)
    HGPL(IR)=(RLH -  SQRT(RLH**2 - HL(IR)*(2.0*RLH - HL(IR))     - &
               2.0*INT(IR)/GRAVH))*1.0E-3
  ENDDO


!-------8.1  Compute height of pressure levels:levels below surface-------------

  DO IR =RLEVSURF-1,1,-1
    INT(IR)     =0.5*(PRES(IR+1)-PRES(IR-1+1))*                        &
                 (1/DMAIR(IR+1)+1/DMAIR(IR))
    INT(IR)     =INT(IR)*1.0E2            ! Integral in units [m2.s-1]
    RLH     =GRAVL/GRAVH
    HL(IR)      =1.0E3*HGPL(IR+1)
    HGPL(IR)=(RLH -  SQRT(RLH**2 - HL(IR)*(2.0*RLH - HL(IR))     - &
               2.0*INT(IR)/GRAVH))*1.0E-3
  ENDDO



!---------Compute tickness of the layer-----------------------------------------

  DO ILAY=1,NLEV-1
    LTICK(ILAY)=HGPL(ILAY+1)-HGPL(ILAY)
  ENDDO

!-------------------------------------------------------------------------------

!---------Compute aerosol number density in each layer for a given -------------
!         climatological composition:
!
!         J=1  -->Continental clean
!         J=2  -->Continental average
!         J=3  -->Continental polluted
!         J=4  -->Urban
!         J=5  -->Desert
!         J=6  -->Maritime clean
!         J=7  -->Maritime polluted
!         J=8  -->Maritime tropical
!         J=9  -->Arctic
!         J=10 -->Antarctic

  DO J=1,10

!-----------Nearest layer-------------------------------------------------------

    DO I=1,3
      DO ILAY=1,NLEV
        IF(HGPL(ILAY)>CHEIG(J,I))THEN
          AFAC(I)=ILAY
          EXIT
        ENDIF
      ENDDO
    ENDDO

!-----------Compute the average number density in each layer:planetary layer----

    DO ILAY=1,AFAC(1)-2
      RILAY = NLEV - ILAY
      IF(ZFAC(J,1) /= 99._jprb)THEN
        HFAC(ILAY)=ZFAC(J,1)                                         * &
                   (EXP(-HGPL(ILAY)/ZFAC(J,1))                       - &
                   EXP(-HGPL(ILAY+1)/ZFAC(J,1)))
      ELSE IF(ZFAC(J,1) == 99._jprb)THEN
        HFAC(ILAY)=LTICK(ILAY)
      ENDIF
      DO I=1,CNUM(J)
        II   =CCOMP(I,J)
        AERPROF(RILAY,J,II)=SCALEFACTOR*NDENS(I,J)*HFAC(ILAY)
      ENDDO
    ENDDO

    RILAY = NLEV - ILAY

    IF(ZFAC(J,1) /= 99._jprb)THEN
      HFACSL1     =ZFAC(J,1)                                         * &
                  (EXP(-HGPL(ILAY)/ZFAC(J,1))                        - &
                  EXP(-CHEIG(J,1)/ZFAC(J,1)))
    ELSE IF(ZFAC(J,1) == 99._jprb)THEN
      HFACSL1     =CHEIG(J,1)-HGPL(ILAY)
    ENDIF
    DO I=1,CNUM(J)
      II   =CCOMP(I,J)
      AERPROF(RILAY,J,II)=SCALEFACTOR*NDENS(I,J)*HFACSL1
    ENDDO

    IF(ZFAC(J,2) /= 99._jprb)THEN
      HFACSL2=ZFAC(J,2)                                              * &
              (EXP(-CHEIG(J,1)/ZFAC(J,2))                            - &
              EXP(-HGPL(ILAY+1)/ZFAC(J,2)))
    ELSE IF(ZFAC(J,2) == 99._jprb)THEN
      HFACSL2     =(HGPL(ILAY+1)-CHEIG(J,1))
    ENDIF
    DO I=1,CNUM(3)
      II   =CCOMP(I,3)
      AERPROF(RILAY,J,II)=SCALEFACTOR*NDENSFT(I)*HFACSL2
    ENDDO

!-----------Compute the average number density in each layer:free troposphere---

    DO ILAY=AFAC(1),AFAC(2)-2
      RILAY = NLEV - ILAY
      IF(ZFAC(J,2) /= 99._jprb)THEN
        HFAC(ILAY)=ZFAC(J,2)                                         * &
                   (EXP(-HGPL(ILAY)/ZFAC(J,2))                       - &
                   EXP(-HGPL(ILAY+1)/ZFAC(J,2   )))
      ELSE IF(ZFAC(J,2) == 99._jprb)THEN
        HFAC(ILAY)=LTICK(ILAY)
      ENDIF
      DO I=1,CNUM(3)
        II   =CCOMP(I,3)
        AERPROF(RILAY,J,II)=SCALEFACTOR*NDENSFT(I)*HFAC(ILAY)
      ENDDO
    ENDDO

    RILAY = NLEV - ILAY

    IF(ZFAC(J,2) /= 99._jprb)THEN
      HFACFT1     =ZFAC(J,2)                                         * &
                   (EXP(-HGPL(ILAY)/ZFAC(J,2))                       - &
                   EXP(-CHEIG(J,2)/ZFAC(J,2)))
    ELSE IF(ZFAC(J,2) == 99._jprb)THEN
      HFACFT1     =CHEIG(J,2)-HGPL(ILAY)
    ENDIF
    DO I=1,CNUM(3)
      II   =CCOMP(I,3)
      AERPROF(RILAY,J,II)=SCALEFACTOR*NDENSFT(I)*HFACFT1
    ENDDO

    IF(ZFAC(J,3) /= 99._jprb)THEN
      HFACFT2  =ZFAC(J,3)                                            * &
                (EXP(-CHEIG(J,2)/ZFAC(J,3))                          - &
                EXP(-HGPL(ILAY+1)/ZFAC(J,3)))
    ELSE IF(ZFAC(J,3) == 99._jprb)THEN
      HFACFT2     =(HGPL(ILAY+1)-CHEIG(J,2))
    ENDIF
    AERPROF(RILAY,J,10)=SCALEFACTOR*NDENSST*HFACFT2

!-----------Compute the average number density in each layer:stratosphere-------

    DO ILAY=AFAC(2),AFAC(3)-1
      RILAY = NLEV - ILAY
      IF(ZFAC(J,3) /= 99._jprb)THEN
        HFAC(ILAY)=ZFAC(J,3)                                         * &
                   (EXP(-HGPL(ILAY)/ZFAC(J,3))                       - &
                   EXP(-HGPL(ILAY+1)/ZFAC(J,3   )))
      ELSE IF(ZFAC(J,3) == 99._jprb)THEN
        HFAC(ILAY)=LTICK(ILAY)
      ENDIF
      AERPROF(RILAY,J,10)=SCALEFACTOR*NDENSST*HFAC(ILAY) !HFACFT2
    ENDDO

    DO ILAY=1,NLEV-1
      RILAY = NLEV - ILAY
      AERPROF(RILAY,J,:) = AERPROF(RILAY,J,:)/LTICK(ILAY)
    ENDDO

    IF (MMR_AER) THEN
      PLAY = 0.5_jprb * 100._jprb * (P(1:NLEV-1) + P(2:NLEV)) ! Units: Pa
      TLAY = 0.5_jprb * (T(1:NLEV-1) + T(2:NLEV))
      QLAY = 0.5_jprb * (QMR(1:NLEV-1) + QMR(2:NLEV))
      FACTOR = (1._jprb + QLAY * (Mair - Mh2o) / Mh2o) * rgc / Mair
      FACTOR = PLAY / (TLAY * FACTOR)
      DO ILAY=1,NLEV-1
        AERPROF(ILAY,J,:) = AERPROF(ILAY,J,:) * AER_OPAC_CONFAC(:) / FACTOR(ILAY)
      ENDDO
    ENDIF

  ENDDO ! J

END SUBROUTINE rttov_aer_clim_prof

