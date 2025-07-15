! Description:
!> @file
!!   K of MW surface emissivity calculation
!
!> @brief
!!   K of MW surface emissivity calculation
!!
!! @param[in]     opts                options to configure the simulations
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in,out] profiles_k          profile increments
!! @param[in]     geometry            internal geometry structure
!! @param[in]     coef                optical depth coefficients structure
!! @param[in]     chanprof            specifies channels and profiles to simulate
!! @param[in]     transmission_aux    RTTOV internal auxiliary transmission structure
!! @param[in,out] transmission_aux_k  auxiliary transmission increments
!! @param[in]     calcemis            flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity_k        Jacobians wrt surface emissivities
!! @param[in,out] reflectivity_k      Jacobians wrt surface reflectances
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
SUBROUTINE rttov_calcemis_mw_k( &
            & opts,               &
            & profiles,           &
            & profiles_k,         &
            & geometry,           &
            & coef,               &
            & chanprof,           &
            & transmission_aux,   &
            & transmission_aux_k, &
            & calcemis,           &
            & emissivity_k,       &
            & reflectivity_k)

  USE rttov_types, ONLY :  &
       & rttov_options,          &
       & rttov_chanprof,         &
       & rttov_coef,             &
       & rttov_profile,          &
       & rttov_transmission_aux, &
       & rttov_geometry
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : &
       & pi,              &
       & surftype_sea,    &
       & pol_v,           &
       & pol_h,           &
       & pol_s3,          &
       & max_fastem_version
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim
  USE mod_rttov_fastem3_coef, ONLY : fastem3_coef
  USE rttov_tessem_mod, ONLY : rttov_tessem_ad
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options),          INTENT(IN)            :: opts
  TYPE(rttov_chanprof),         INTENT(IN)            :: chanprof(:)
  TYPE(rttov_profile),          INTENT(IN)   , TARGET :: profiles(:)
  TYPE(rttov_geometry),         INTENT(IN)   , TARGET :: geometry(size(profiles))
  TYPE(rttov_coef),             INTENT(IN)            :: coef
  TYPE(rttov_transmission_aux), INTENT(IN)            :: transmission_aux
  LOGICAL(KIND=jplm),           INTENT(IN)            :: calcemis  (size(chanprof))
  TYPE(rttov_profile),          INTENT(INOUT), TARGET :: profiles_k(size(chanprof))
  TYPE(rttov_transmission_aux), INTENT(INOUT)         :: transmission_aux_k
  REAL(KIND=jprb),              INTENT(INOUT)         :: emissivity_k  (size(chanprof))
  REAL(KIND=jprb),              INTENT(INOUT)         :: reflectivity_k(size(chanprof))
!INTF_END

#include "rttov_fastem5_k.interface"

!lddocal constants:
  REAL(KIND=jprb)            , PARAMETER             :: windscale = 999999.0_JPRB
  REAL(KIND=jprb)            , PARAMETER             :: windlimit = 0.0001_JPRB
  REAL(KIND=jprb)            , PARAMETER             :: quadcof  (4, 2  ) =      &
    & Reshape((/0.0_JPRB, 1.0_JPRB, 1.0_JPRB, 2.0_JPRB, 1.0_JPRB,  - 1.0_JPRB, 1.0_JPRB,  - 1.0_JPRB/), (/4, 2/))
  REAL(KIND=jprb)            , PARAMETER             ::      &
    & freqfixed(4)      = Reshape((/7.0_JPRB, 10.0_JPRB, 19.0_JPRB, 37.0_JPRB/), (/4/))

!local variables:
  REAL(KIND=jprb) :: tcelsius
  REAL(KIND=jprb) :: tcelsius_sq
  REAL(KIND=jprb) :: tcelsius_cu
  REAL(KIND=jprb) :: f1, f2
  REAL(KIND=jprb) :: del1          , del2
  REAL(KIND=jprb) :: einf
  REAL(KIND=jprb) :: fen, fen_sq
  REAL(KIND=jprb) :: den1          , den2
  REAL(KIND=jprb) :: perm_free
  REAL(KIND=jprb) :: sigma
  REAL(KIND=jprb) :: perm_real1    , perm_real2
  REAL(KIND=jprb) :: perm_imag1    , perm_imag2    , perm_imag3
  REAL(KIND=jprb) :: perm_Real     , perm_imag
  REAL(KIND=jprb) :: perm_static   , perm_infinite
  REAL(KIND=jprb) :: freq_ghz      , freq_ghz_sq
  REAL(KIND=jprb) :: fresnel_v_Real, fresnel_v_imag
  REAL(KIND=jprb) :: fresnel_h_Real, fresnel_h_imag
  REAL(KIND=jprb) :: fresnel(4)
  REAL(KIND=jprb) :: small_rough_cor, foam_cor(4)
  REAL(KIND=jprb) :: large_rough_cor(4)
  REAL(KIND=jprb) :: small_rough, large_rough
  REAL(KIND=jprb) :: emiss_save(4)
  REAL(KIND=jprb) :: variance        , varm
  REAL(KIND=jprb) :: wind10
  REAL(KIND=jprb) :: wind10_sq       , windsec
  REAL(KIND=jprb) :: wind10_direction, windangle    , windratio ! Note wind azimuth is in radians
  REAL(KIND=jprb) :: emissstokes    (size(chanprof), 4)
  REAL(KIND=jprb) :: emissstokes_k  (size(chanprof), 4)
  REAL(KIND=jprb) :: reflectstokes_k(size(chanprof), 4)
  REAL(KIND=jprb) :: u19, phi, dfreq
  REAL(KIND=jprb) :: tbfixed      (4, 4, 3)                     ! Surface brightness temperature azimuthal
                                                                !   variation terms for 37, 19, 10, 7 GHz
  REAL(KIND=jprb) :: efixed       (4, 4, 3)                     ! Emissivity azimuthal variation terms for
                                                                !   7, 10, 19, 37 GHz
  REAL(KIND=jprb) :: einterpolated(4, 3   )                     ! Emissivity azimuthal variation terms for
                                                                !   interpolated to required frequency
  REAL(KIND=jprb) :: a1e, a2e, a3e                              ! coefficients used in azimuthal emissivity model
  REAL(KIND=jprb) :: zrough_v     , zrough_h
  REAL(KIND=jprb) :: zreflmod_v   , zreflmod_h
  REAL(KIND=jprb) :: delta        , delta2
  REAL(KIND=jprb) :: qdepol       , emissfactor
  REAL(KIND=jprb) :: emissfactor_v, emissfactor_h
  REAL(KIND=jprb) :: zc(12), zx(9)
  REAL(KIND=jprb) :: opdpsfc, freqr
  REAL(KIND=jprb), POINTER :: c(:)
  COMPLEX(KIND=jprb) :: perm1       , perm2
  COMPLEX(KIND=jprb) :: rhth        , rvth
  COMPLEX(KIND=jprb) :: permittivity
  INTEGER(KIND=jpim) :: i, j, chan, istokes, ifreq, m
  INTEGER(KIND=jpim) :: iquadrant                            ! Determines which quadrant (NE, SE, SW, NW) the wind
                                                             !   is blowing to
  INTEGER(KIND=jpim) :: pol_id                               ! polarisation indice
  INTEGER(KIND=jpim) :: i_freq      , j_stokes, ich          ! indices used in azimuthal emissivity model
  INTEGER(KIND=jpim) :: jcof        , jcofm1
  TYPE(rttov_profile),  POINTER :: prof
  TYPE(rttov_profile),  POINTER :: prof_k
  TYPE(rttov_geometry), POINTER :: geom
  REAL   (KIND=jprb) :: tcelsius_k
  REAL   (KIND=jprb) :: tcelsius_sq_k
  REAL   (KIND=jprb) :: tcelsius_cu_k
  REAL   (KIND=jprb) :: f1_k, f2_k
  REAL   (KIND=jprb) :: del1_k           , del2_k
  REAL   (KIND=jprb) :: einf_k
  REAL   (KIND=jprb) :: fen_k            , fen_sq_k
  REAL   (KIND=jprb) :: den1_k           , den2_k
  REAL   (KIND=jprb) :: sigma_k
  REAL   (KIND=jprb) :: perm_real1_k     , perm_real2_k
  REAL   (KIND=jprb) :: perm_imag1_k     , perm_imag2_k    , perm_imag3_k
  REAL   (KIND=jprb) :: perm_Real_k      , perm_imag_k
  REAL   (KIND=jprb) :: perm_static_k    , perm_infinite_k
  REAL   (KIND=jprb) :: fresnel_v_Real_k , fresnel_v_imag_k
  REAL   (KIND=jprb) :: fresnel_h_Real_k , fresnel_h_imag_k
  REAL   (KIND=jprb) :: fresnel_v_k      , fresnel_h_k
  REAL   (KIND=jprb) :: small_rough_cor_k, foam_cor_k
  REAL   (KIND=jprb) :: large_rough_cor_k(2)
  REAL   (KIND=jprb) :: small_rough_k     , large_rough_k
  REAL   (KIND=jprb) :: variance_k        , varm_k
  REAL   (KIND=jprb) :: wind10_k
  REAL   (KIND=jprb) :: wind10_sq_k       , windsec_k
  REAL   (KIND=jprb) :: wind10_direction_k, windangle_k, windratio_k  ! Note wind azimuth is in radians
  REAL   (KIND=jprb) :: azimuthal_emiss_k , azimuthal_emiss, u19_k, phi_k
  REAL   (KIND=jprb) :: tbfixed_k      (4, 4, 3)             ! Surface brightness temperature azimuthal
                                                             !   variation terms for 37, 19, 10, 7 GHz
  REAL   (KIND=jprb) :: efixed_k       (4, 4, 3)             ! Emissivity azimuthal variation terms for
                                                             !   7, 10, 19, 37 GHz
  REAL   (KIND=jprb) :: einterpolated_k(4, 3   )             ! Emissivity azimuthal variation terms for
                                                             !   interpolated to required frequency
  REAL   (KIND=jprb) :: a1e_k          , a2e_k, a3e_k        ! coefficients used in azimuthal emissivity model
  REAL   (KIND=jprb) :: opdpsfc_k      , freqr_k
  REAL   (KIND=jprb) :: zrough_v_k     , zrough_h_k
  REAL   (KIND=jprb) :: zreflmod_v_k   , zreflmod_h_k
  REAL   (KIND=jprb) :: delta_k        , delta2_k
  REAL   (KIND=jprb) :: qdepol_k       , emissfactor_k
  REAL   (KIND=jprb) :: emissfactor_v_k, emissfactor_h_k
  REAL   (KIND=jprb) :: zx_k(9)
  COMPLEX(KIND=jprb) :: perm1_k          , perm2_k
  COMPLEX(KIND=jprb) :: rhth_k           , rvth_k
  COMPLEX(KIND=jprb) :: permittivity_k
  REAL   (KIND=jprb) :: test_variance
! rttov_fastem4
  REAL   (KIND=jprb) :: Zenith_Angle     , Salinity        , Transmittance, Rel_Azimuth, Wind_Speed
  REAL(KIND=jprb), DIMENSION(4) :: jemissivity, jreflectivity
  REAL(KIND=jprb) :: Salinity_k   , Transmittance_k, Rel_Azimuth_k
  REAL(KIND=jprb) :: Temperature_k, Wind_Speed_k
  REAL(KIND=jprb), DIMENSION(4) :: jemissivity_k, jreflectivity_k
  INTEGER(KIND=jpim) :: nchannels   ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! The rttov_check_options subroutine guarantees the selected FASTEM version is valid
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW_K', 0_jpim, ZHOOK_HANDLE)
  nchannels          = size(chanprof)

!Loop over channels
  DO i = 1, nchannels
    IF (.NOT. calcemis(i)) CYCLE
    chan  = chanprof(i)%chan
    prof => profiles(chanprof(i)%prof)
    prof_k => profiles_k(i)
    geom => geometry(chanprof(i)%prof)
!-------------------------------
!0. Point to fastem coefficients
!-------------------------------
    c => fastem3_coef
    pol_id                = coef%fastem_polar(chan) + 1
    reflectstokes_k(i, :) = 0.0_JPRB
    emissstokes_k(i, :)   = 0.0_JPRB
    wind10_k              = 0.0_JPRB
    wind10_direction_k    = 0.0_JPRB
    tcelsius_k            = 0.0_JPRB
!
! Now calc channel emissivity after mixing v and h pol
!
    emissfactor_v_k       = pol_v(1, pol_id) * coef%pol_fac_v(chan) + &
                            pol_v(2, pol_id) * geom%sinview_sq + pol_v(3, pol_id) * geom%cosview_sq
    emissfactor_h_k       = pol_h(1, pol_id) * coef%pol_fac_h(chan) + &
                            pol_h(2, pol_id) * geom%sinview_sq + pol_h(3, pol_id) * geom%cosview_sq
    reflectstokes_k(i, 1) = emissfactor_v_k * reflectivity_k(i)
    reflectstokes_k(i, 2) = emissfactor_h_k * reflectivity_k(i)
    reflectstokes_k(i, 3) = pol_s3(0, pol_id) * reflectivity_k(i)
    reflectstokes_k(i, 4) = pol_s3(1, pol_id) * reflectivity_k(i)
    reflectivity_k(i)     = 0._jprb
    emissstokes_k(i, 1)   = emissfactor_v_k * emissivity_k(i)
    emissstokes_k(i, 2)   = emissfactor_h_k * emissivity_k(i)
    emissstokes_k(i, 3)   = pol_s3(0, pol_id) * emissivity_k(i)
    emissstokes_k(i, 4)   = pol_s3(1, pol_id) * emissivity_k(i)
    emissivity_k(i)       = 0._jprb
!---------------
!1. Sea surfaces
!---------------
    IF (prof%skin%surftype == surftype_sea) THEN
!-------------------------------------------
!1.1 Calculate channel independent variables
!-------------------------------------------
      wind10_sq = prof%s2m%u * prof%s2m%u + prof%s2m%v * prof%s2m%v
      wind10    = Sqrt(wind10_sq)
      freq_ghz       = coef%frequency_ghz(chan)
      freq_ghz_sq    = freq_ghz * freq_ghz
      windsec   = wind10 * geom%seczen
      IF (prof%s2m%u >= 0.0_JPRB .AND. prof%s2m%v >= 0.0_JPRB) iquadrant = 1
      IF (prof%s2m%u >= 0.0_JPRB .AND. prof%s2m%v < 0.0_JPRB ) iquadrant = 2
      IF (prof%s2m%u < 0.0_JPRB .AND. prof%s2m%v >= 0.0_JPRB ) iquadrant = 4
      IF (prof%s2m%u < 0.0_JPRB .AND. prof%s2m%v < 0.0_JPRB  ) iquadrant = 3
      IF (abs(prof%s2m%v) >= windlimit) THEN
        windratio = prof%s2m%u / prof%s2m%v
      ELSE
        windratio = 0.0
        IF (abs(prof%s2m%u) > windlimit) THEN
          windratio = windscale * prof%s2m%u
        ENDIF
      ENDIF
      windangle        = atan(abs(windratio))
      wind10_direction = quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2)
      SELECT CASE (opts%rt_mw%fastem_version)
      CASE (0)
        emissstokes_k(i,:) = emissstokes_k(i,:) - reflectstokes_k(i,:)
        emissstokes_k(i,3:4) = 0._jprb
        CALL rttov_tessem_ad(freq_ghz, prof%zenangle, wind10, prof%skin%t, prof%skin%salinity, &
                             wind10_k, prof_k%skin%t, prof_k%skin%salinity, &
                             emissstokes_k(i,2), emissstokes_k(i,1))
        wind10_sq_k = 0._jprb
        wind10_direction_k = 0._jprb
      CASE (1:3)
!Set values for temperature polynomials (convert from kelvin to celsius)
        tcelsius       = prof%skin%t - 273.15_JPRB
        tcelsius_sq    = tcelsius * tcelsius                                                 !quadratic
        tcelsius_cu    = tcelsius_sq * tcelsius                                              !cubic
!Define two relaxation frequencies, f1 and f2
        f1 = c(1) + c(2) * tcelsius + c(3) * tcelsius_sq
        f2 = c(4) + c(5) * tcelsius + c(6) * tcelsius_sq + c(7) * tcelsius_cu
!Static permittivity estatic = del1+del2+einf
        del1           = c(8) + c(9) * tcelsius + c(10) * tcelsius_sq + c(11) * tcelsius_cu
        del2           = c(12) + c(13) * tcelsius + c(14) * tcelsius_sq + c(15) * tcelsius_cu
        einf           = c(18) + c(19) * tcelsius
!-----------------------------------------------------
!1.2 calculate permittivity using double-debye formula
!-----------------------------------------------------
        fen = 2.0_JPRB * c(20) * freq_ghz * 0.001_JPRB
        fen_sq         = fen * fen
        den1           = 1.0_JPRB + fen_sq * f1 * f1
        den2           = 1.0_JPRB + fen_sq * f2 * f2
        perm_real1     = del1 / den1
        perm_real2     = del2 / den2
        perm_imag1     = del1 * fen * f1 / den1
        perm_imag2     = del2 * fen * f2 / den2
        perm_free      = 8.854E-03_JPRB
        sigma          = 2.906_JPRB + 0.09437_JPRB * tcelsius
        perm_imag3     = sigma / (2.0_JPRB * c(20) * perm_free * freq_ghz)
        perm_Real      = perm_real1 + perm_real2 + einf
!        perm_imag    = perm_imag1 + perm_imag2 + perm_imag3 + perm_imag3
        perm_imag      = perm_imag1 + perm_imag2 + perm_imag3
        permittivity   = Cmplx(perm_Real, perm_imag, jprb)
!-------------------------------------------------------------
!1.3 calculate complex reflection coefficients and corrections
!-------------------------------------------------------------
!1.3.1) Fresnel reflection coefficients
!------
        perm1          = sqrt(permittivity - geom%sinzen_sq)
        perm2          = permittivity * geom%coszen
        rhth           = (geom%coszen - perm1) / (geom%coszen + perm1)
        rvth           = (perm2 - perm1) / (perm2 + perm1)
!    fresnel_v_real = dble(rvth)
        fresnel_v_Real = Real(rvth)
        fresnel_v_imag = Aimag(rvth)
        fresnel(1)     = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
!    fresnel_h_real = dble(rhth)
        fresnel_h_Real = Real(rhth)
        fresnel_h_imag = Aimag(rhth)
        fresnel(2)     = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
        fresnel(3)     = 0.0_JPRB
        fresnel(4)     = 0.0_JPRB
!1.3.2) Small scale correction to reflection coefficients
!------
        IF (freq_ghz >= 15.0) THEN
          small_rough_cor = Exp(c(21) * wind10 * geom%coszen_sq / (freq_ghz_sq))
        ELSE
          small_rough_cor = 1.0
        ENDIF
!1.3.3) Large scale geometric correction
!------
!Point to correct coefficients for this version. There are 36 altogether.
!Those for FASTEM-2/3 are stored in section 24:59 of the array, those for
!FASTEM1 in section 60:95.
        IF (opts%rt_mw%fastem_version == 2 .OR. opts%rt_mw%fastem_version == 3) THEN
          c => fastem3_coef(24:59)
        ELSE
          c => fastem3_coef(60:95)
        ENDIF
        DO j = 1, 12
          zc(j) = c(j * 3 - 2) + c(j * 3 - 1) * freq_ghz + c(j * 3) * freq_ghz_sq
        ENDDO
!Point back to all coefficients again
        c => fastem3_coef
        large_rough_cor(1) =                                                                                                 &
          & (zc(1) + zc(2) * geom%seczen + zc(3) * geom%seczen_sq + zc(4) * wind10 + zc(5) * wind10_sq + zc(6) * windsec) /  &
          & 100._JPRB
        large_rough_cor(2) =                                                                                                 &
          & (zc(7) + zc(8) * geom%seczen + zc(9) * geom%seczen_sq + zc(10) * wind10 + zc(11) * wind10_sq + zc(12) * windsec) &
          &  / 100._JPRB
        large_rough_cor(3) = 0.0_JPRB
        large_rough_cor(4) = 0.0_JPRB
! Introduce emiss_v_save and emiss_h_save arrays to be able
! to simplify further AD code
        emiss_save(:)      = 1.0 - fresnel(:) * small_rough_cor + large_rough_cor(:)
!Apply foam correction
        IF (.NOT. opts%rt_mw%supply_foam_fraction) THEN
          foam_cor(1)        = c(22) * (wind10 ** c(23))
          foam_cor(2)        = c(22) * (wind10 ** c(23))
        ELSE
          foam_cor(1)        = prof%skin%foam_fraction
          foam_cor(2)        = prof%skin%foam_fraction
        ENDIF
!Currently ignore foam effects on 3rd and 4th elements.
        foam_cor(3)        = 0.0_JPRB
        foam_cor(4)        = 0.0_JPRB
        emissstokes(i, :)  = emiss_save(:) - foam_cor(:) * emiss_save(:) + foam_cor(:)
        emissstokes(i, 3)  = 0.0
        emissstokes(i, 4)  = 0.0
        IF (opts%rt_mw%fastem_version == 3) THEN
! Add azimuthal component from Fuzhong Weng (NOAA/NESDIS) based on work by Dr. Gene Poe (NRL)
! Angle between wind direction and satellite azimuthal view angle
! Assume 19m wind = 10m wind for now (fix later).
          phi = pi - wind10_direction + prof%azangle * pi / 180.0_JPRB
          u19 = wind10
          DO ich = 0, 15
            a1e = c(141 + ich * 12) + u19 * (c(142 + ich * 12) + u19 * (c(143 + ich * 12) + u19 * c(144 + ich * 12)))
            a2e = c(145 + ich * 12) + u19 * (c(146 + ich * 12) + u19 * (c(147 + ich * 12) + u19 * c(148 + ich * 12)))
            a3e = c(149 + ich * 12) + u19 * (c(150 + ich * 12) + u19 * (c(151 + ich * 12) + u19 * c(152 + ich * 12)))
            i_freq = int(ich / 4_jpim) + 1                       ! 37, 19, 10, 7 GHz
            j_stokes                     = mod(ich, 4_jpim) + 1
            tbfixed(j_stokes, i_freq, 1) = a1e
            tbfixed(j_stokes, i_freq, 2) = a2e
            tbfixed(j_stokes, i_freq, 3) = a3e
          ENDDO
          DO M = 1, 3
            DO istokes = 1, 4
              efixed(1, istokes, M) = tbfixed(istokes, 4, M)! 7   GHz
              efixed(2, istokes, M) = tbfixed(istokes, 3, M)! 10  GHz
              efixed(3, istokes, M) = tbfixed(istokes, 2, M)! 19  GHz
              efixed(4, istokes, M) = tbfixed(istokes, 1, M)! 37  GHz
            ENDDO
! Interpolate results to required frequency based on 7, 10, 19, 37 GHz
            IF (freq_ghz .LE. freqfixed(1)) THEN
              einterpolated(:, M) = efixed(1, :, M)
            ELSE IF (freq_ghz .GE. freqfixed(4)) THEN
              einterpolated(:, M) = efixed(4, :, M)
            ELSE
              IF (freq_ghz .LT. freqfixed(2)                                 ) ifreq               = 2
              IF (freq_ghz .LT. freqfixed(3) .AND. freq_ghz .GE. freqfixed(2)) ifreq               = 3
              IF (freq_ghz .GE. freqfixed(3)                                 ) ifreq               = 4
              dfreq               = (freq_ghz - freqfixed(ifreq - 1)) / (freqfixed(ifreq) - freqfixed(ifreq - 1))
              einterpolated(:, M) = efixed(ifreq - 1, :, M) + dfreq * (efixed(ifreq, :, M) - efixed(ifreq - 1, :, M))
            ENDIF
          ENDDO
          DO istokes = 1, 4
            azimuthal_emiss = 0.0_JPRB
            DO M = 1, 3
              IF (istokes .LE. 2) THEN
                azimuthal_emiss = azimuthal_emiss +      &
                  & einterpolated(istokes, M) * cos(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
              ELSE
                azimuthal_emiss = azimuthal_emiss +      &
                  & einterpolated(istokes, M) * sin(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
              ENDIF
            ENDDO
            emissstokes(i, istokes) = emissstokes(i, istokes) + azimuthal_emiss
          ENDDO
        ENDIF
        IF ((opts%rt_mw%fastem_version == 2 .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) < 0.9999_JPRB .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) > 0.00001_JPRB) THEN
!Convert windspeed to slope variance using the Cox and Munk model
          variance      = 0.00512_JPRB * wind10 + 0.0030_JPRB
          varm          = variance * c(138)
          variance      = varm * (c(139) * freq_ghz + c(140))
          test_variance = variance
          IF (variance > varm) THEN
            variance = varm
          ENDIF
          IF (variance < 0.0_JPRB) THEN
            variance = 0.0_JPRB
          ENDIF
!Compute surface to space optical depth
          opdpsfc  =  - log(transmission_aux%thermal_path1%tau_surf(0, i)) / geom%seczen
!Define nine predictors for the effective angle calculation
          zx(1)    = 1.0_JPRB
          zx(2)    = variance
          zx(4)    = 1.0_JPRB / geom%coszen
          zx(3)    = zx(2) * zx(4)
          zx(5)    = zx(3) * zx(3)
          zx(6)    = zx(4) * zx(4)
          zx(7)    = zx(2) * zx(2)
          zx(8)    = log(opdpsfc)
          zx(9)    = zx(8) * zx(8)
          zrough_v = 1.0_JPRB
          zrough_h = 1.0_JPRB
          DO jcof = 1, 7
            jcofm1   = jcof - 1
!Switched h to v Deblonde SSMIS june 7, 2001
            zrough_h =      &
              & zrough_h + zx(jcof) * (c(96 + jcofm1 * 3) + zx(8) * c(97 + jcofm1 * 3) + zx(9) * c(98 + jcofm1 * 3))
            zrough_v =      &
              & zrough_v + zx(jcof) * (c(117 + jcofm1 * 3) + zx(8) * c(118 + jcofm1 * 3) + zx(9) * c(119 + jcofm1 * 3))
          ENDDO
          zreflmod_v =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v) / &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          zreflmod_h =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h) / &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
        ENDIF
!.......end of forward part....................................
!
! * Now run K code of fastem
!
! Only apply non-specular correction for Fastem-3 if theta < 60 degrees
        IF ((opts%rt_mw%fastem_version == 2 .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) < 0.9999_JPRB .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) > 0.00001_JPRB) THEN
          zreflmod_v_k                      = reflectstokes_k(i, 1) * (1.0_JPRB - emissstokes(i, 1))
          zreflmod_h_k                      = reflectstokes_k(i, 2) * (1.0_JPRB - emissstokes(i, 2))
          zreflmod_v_k                      = zreflmod_v_k - 0.5_JPRB * reflectstokes_k(i, 3) * emissstokes(i, 3) -      &
            & 0.5_JPRB * reflectstokes_k(i, 4) * emissstokes(i, 4)
          zreflmod_h_k                      = zreflmod_h_k - 0.5_JPRB * reflectstokes_k(i, 3) * emissstokes(i, 3) -      &
            & 0.5_JPRB * reflectstokes_k(i, 4) * emissstokes(i, 4)
          emissstokes_k(i, 4)               =      &
            & emissstokes_k(i, 4) - 0.5_JPRB * (zreflmod_v + zreflmod_h) * reflectstokes_k(i, 4)
          emissstokes_k(i, 3)               =      &
            & emissstokes_k(i, 3) - 0.5_JPRB * (zreflmod_v + zreflmod_h) * reflectstokes_k(i, 3)
          emissstokes_k(i, 2)               = emissstokes_k(i, 2) - reflectstokes_k(i, 2) * zreflmod_h
          emissstokes_k(i, 1)               = emissstokes_k(i, 1) - reflectstokes_k(i, 1) * zreflmod_v
          zrough_v_k                        =                                                          &
            & - zreflmod_v_k * (transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v *            &
            & Log(transmission_aux%thermal_path1%tau_surf(0, i))) /                                    &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          transmission_aux_k%thermal_path1%tau_surf(0, i) = &
            & transmission_aux_k%thermal_path1%tau_surf(0, i) + zreflmod_v_k * ( -                     &
            & zrough_v * transmission_aux%thermal_path1%tau_surf(0, i) ** (zrough_v - 1.0_JPRB) *      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) +                             &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v)) /                &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) ** 2
          zrough_h_k                        =                                                          &
            & - zreflmod_h_k * (transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h *            &
            & Log(transmission_aux%thermal_path1%tau_surf(0, i))) /  &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          transmission_aux_k%thermal_path1%tau_surf(0, i) = &
            & transmission_aux_k%thermal_path1%tau_surf(0, i) + zreflmod_h_k * ( -                     &
            & zrough_h * transmission_aux%thermal_path1%tau_surf(0, i) ** (zrough_h - 1.0_JPRB) *      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) +                             &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h)) /                &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) ** 2
          zx_k(:) = 0._JPRB
          DO jcof = 1, 7
            jcofm1     = jcof - 1
!Switched h to v Deblonde SSMIS june 7, 2001
            zx_k(9)    = zx_k(9) + zrough_v_k * zx(jcof) * c(119 + jcofm1 * 3)
            zx_k(8)    = zx_k(8) + zrough_v_k * zx(jcof) * c(118 + jcofm1 * 3)
            zx_k(jcof) = zrough_v_k * (c(117 + jcofm1 * 3) + zx(8) * c(118 + jcofm1 * 3) + zx(9) * c(119 + jcofm1 * 3))
            zx_k(9)    = zx_k(9) + zrough_h_k * zx(jcof) * c(98 + jcofm1 * 3)
            zx_k(8)    = zx_k(8) + zrough_h_k * zx(jcof) * c(97 + jcofm1 * 3)
            zx_k(jcof) =      &
              & zx_k(jcof) + zrough_h_k * (c(96 + jcofm1 * 3) + zx(8) * c(97 + jcofm1 * 3) + zx(9) * c(98 + jcofm1 * 3))
          ENDDO
          zrough_v_k                        = 0._JPRB
          zrough_h_k                        = 0._JPRB
!Define nine predictors for the effective angle calculation
          zx_k(8) = zx_k(8) + zx_k(9) * 2 * zx(8)
          opdpsfc_k                         = zx_k(8) / opdpsfc
          zx_k(2) = zx_k(2) + zx_k(7) * 2 * zx(2)
          zx_k(4) = zx_k(4) + zx_k(6) * 2 * zx(4)
          zx_k(3) = zx_k(3) + zx_k(5) * 2 * zx(3)
          zx_k(2) = zx_k(2) + zx_k(3) * zx(4)
          zx_k(4) = 0._JPRB
          variance_k                        = zx_k(2)
          zx_k(1) = 0._JPRB
!Compute surface to space optical depth
          transmission_aux_k%thermal_path1%tau_surf(0, i) = &
            & transmission_aux_k%thermal_path1%tau_surf(0, i) - &
            & opdpsfc_k / (transmission_aux%thermal_path1%tau_surf(0, i) * geom%seczen)
          IF (test_variance < varm) THEN
            varm_k = variance_k * (c(139) * freq_ghz + c(140))
          ELSE
            varm_k = variance_k
          ENDIF
          variance_k = varm_k * c(138)
          wind10_k   = wind10_k + variance_k * 0.00512_JPRB
        ELSE
          emissstokes_k(i, :) = emissstokes_k(i, :) - reflectstokes_k(i, :)
        ENDIF
        IF (opts%rt_mw%fastem_version == 3) THEN
          azimuthal_emiss_k = 0.0_JPRB
          phi_k             = 0.0_JPRB
          DO istokes = 1, 4
            azimuthal_emiss_k = emissstokes_k(i, istokes)
            DO M = 1, 3
              IF (istokes .LE. 2) THEN
                einterpolated_k(istokes, M) =      &
                  & azimuthal_emiss_k * cos(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
                phi_k = phi_k -                                                                                    &
                  & azimuthal_emiss_k * einterpolated(istokes, M) * m * sin(m * phi) * (1.0_JPRB - geom%coszen) /  &
                  & (1.0_JPRB - 0.6018_JPRB)
              ELSE
                einterpolated_k(istokes, M) =      &
                  & azimuthal_emiss_k * sin(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
                phi_k = phi_k +                                                                                    &
                  & azimuthal_emiss_k * einterpolated(istokes, M) * m * cos(m * phi) * (1.0_JPRB - geom%coszen) /  &
                  & (1.0_JPRB - 0.6018_JPRB)
              ENDIF
            ENDDO
          ENDDO
          efixed_k(:,:,:) = 0.0_JPRB
          DO M = 1, 3
            IF (freq_ghz .LE. freqfixed(1)) THEN
              efixed_k(1, :, M) = efixed_k(1, :, M) + einterpolated_k(:, M)
            ELSE IF (freq_ghz .GE. freqfixed(4)) THEN
              efixed_k(4, :, M) = efixed_k(4, :, M) + einterpolated_k(:, M)
            ELSE
              IF (freq_ghz .LT. freqfixed(2)                                 ) ifreq = 2
              IF (freq_ghz .LT. freqfixed(3) .AND. freq_ghz .GE. freqfixed(2)) ifreq = 3
              IF (freq_ghz .GE. freqfixed(3)                                 ) ifreq = 4
              dfreq = (freq_ghz - freqfixed(ifreq - 1)) / (freqfixed(ifreq) - freqfixed(ifreq - 1))
              efixed_k(ifreq, :, M)     = efixed_k(ifreq, :, M) + einterpolated_k(:, M) * dfreq
              efixed_k(ifreq - 1, :, M) = efixed_k(ifreq - 1, :, M) + einterpolated_k(:, M) * (1.0 - dfreq)
            ENDIF
            DO istokes = 1, 4
              tbfixed_k(istokes, 4, M) = efixed_k(1, istokes, M)! 7   GHz
              tbfixed_k(istokes, 3, M) = efixed_k(2, istokes, M)! 10  GHz
              tbfixed_k(istokes, 2, M) = efixed_k(3, istokes, M)! 19  GHz
              tbfixed_k(istokes, 1, M) = efixed_k(4, istokes, M)! 37  GHz
            ENDDO
          ENDDO
          u19_k = 0.0_JPRB
          DO ich = 0, 15
            i_freq   = int(ich / 4) + 1! 37, 19, 10, 7 GHz
            j_stokes = mod(ich, 4_jpim) + 1
            a3e_k    = tbfixed_k(j_stokes, i_freq, 3)
            a2e_k    = tbfixed_k(j_stokes, i_freq, 2)
            a1e_k    = tbfixed_k(j_stokes, i_freq, 1)
            u19_k    =      &
              & u19_k + a3e_k * (c(150 + ich * 12) + u19 * (2.0 * c(151 + ich * 12) + 3.0 * u19 * c(152 + ich * 12)))
            u19_k    =      &
              & u19_k + a2e_k * (c(146 + ich * 12) + u19 * (2.0 * c(147 + ich * 12) + 3.0 * u19 * c(148 + ich * 12)))
            u19_k    =      &
              & u19_k + a1e_k * (c(142 + ich * 12) + u19 * (2.0 * c(143 + ich * 12) + 3.0 * u19 * c(144 + ich * 12)))
          ENDDO
          wind10_k           = wind10_k + u19_k
          wind10_direction_k =  - 1.0_JPRB * phi_k
        ENDIF
! Be careful do TL first because the next 2 lines of the direct model
! have variables in input/output of the statement
        foam_cor_k = 0.0_JPRB
        DO Ich = 1, 4
          foam_cor_k            = foam_cor_k + emissstokes_k(i, ich) * (1.0_JPRB - emiss_save(ich))
          emissstokes_k(i, Ich) = emissstokes_k(i, ich) * (1.0_JPRB - foam_cor(ich))
        ENDDO
!Apply foam correction
        IF (.NOT. opts%rt_mw%supply_foam_fraction) THEN
          wind10_k             = wind10_k + foam_cor_k * c(22) * c(23) * (wind10 ** (c(23) - 1.0_JPRB))
        ELSE
          prof_k%skin%foam_fraction = prof_k%skin%foam_fraction + foam_cor_k
        ENDIF
!1.3.3) Large scale geometric correction
!------
        fresnel_v_k          =  - emissstokes_k(i, 1) * small_rough_cor
        small_rough_cor_k    =  - emissstokes_k(i, 1) * fresnel(1)
        large_rough_cor_k(1) = emissstokes_k(i, 1)
        fresnel_h_k          =  - emissstokes_k(i, 2) * small_rough_cor
        small_rough_cor_k    = small_rough_cor_k - emissstokes_k(i, 2) * fresnel(2)
        large_rough_cor_k(2) = emissstokes_k(i, 2)
        windsec_k            = large_rough_cor_k(2) * zc(12) / 100._JPRB
        wind10_sq_k          = large_rough_cor_k(2) * zc(11) / 100._JPRB
        wind10_k             = wind10_k + large_rough_cor_k(2) * zc(10) / 100._JPRB
        windsec_k            = windsec_k + large_rough_cor_k(1) * zc(6) / 100._JPRB
        wind10_sq_k          = wind10_sq_k + large_rough_cor_k(1) * zc(5) / 100._JPRB
        wind10_k             = wind10_k + large_rough_cor_k(1) * zc(4) / 100._JPRB
        wind10_k             = wind10_k + windsec_k * geom%seczen
!1.3.2) Small scale correction to reflection coefficients
!------
        IF (freq_ghz >= 15.0) THEN
          wind10_k = wind10_k + small_rough_cor_k * small_rough_cor * c(21) * geom%coszen_sq / (freq_ghz_sq)
        ENDIF
!1.3.1) Fresnel reflection coefficients
!------
        fresnel_h_real_k = fresnel_h_k * 2 * fresnel_h_real
        fresnel_h_imag_k = fresnel_h_k * 2 * fresnel_h_imag
        rhth_k           = CMPLX(fresnel_h_real_k,  - fresnel_h_imag_k, jprb)
        fresnel_v_real_k = fresnel_v_k * 2 * fresnel_v_real
        fresnel_v_imag_k = fresnel_v_k * 2 * fresnel_v_imag
        rvth_k           = CMPLX(fresnel_v_real_k,  - fresnel_v_imag_k, jprb)
        perm1_k          =  - rvth_k * 2 * perm2 / (perm2 + perm1) ** 2
        perm2_k          = rvth_k * 2 * perm1 / (perm2 + perm1) ** 2
        perm1_k          = perm1_k - rhth_k * 2 * geom%coszen / (geom%coszen + perm1) ** 2
        permittivity_k   = perm2_k * geom%coszen
        permittivity_k   = permittivity_k + perm1_k * 0.5_JPRB / perm1
!-----------------------------------------------------
!1.2 calculate permittivity using double-debye formula
!-----------------------------------------------------
        perm_Real_k      = Real(permittivity_k)
        perm_imag_k      =  - Aimag(permittivity_k)
        perm_imag1_k     = perm_imag_k
        perm_imag2_k     = perm_imag_k
!        perm_imag3_k = perm_imag_k
        perm_imag3_k     = perm_imag_k
        einf_k           = perm_real_k
        perm_real1_k     = perm_real_k
        perm_real2_k     = perm_real_k
        sigma_k          = perm_imag3_k / (2.0_JPRB * c(20) * perm_free * freq_ghz)
        tcelsius_k       = 0.09437_JPRB * sigma_k
        del2_k           = perm_imag2_k * fen * den2 * f2 / (den2 * den2)
        den2_k           =  - perm_imag2_k * fen * del2 * f2 / (den2 * den2)
        f2_k             = perm_imag2_k * fen * den2 * del2 / (den2 * den2)
        del1_k           = perm_imag1_k * fen * den1 * f1 / (den1 * den1)
        den1_k           =  - perm_imag1_k * fen * del1 * f1 / (den1 * den1)
        f1_k             = perm_imag1_k * fen * den1 * del1 / (den1 * den1)
        del2_k           = del2_k + perm_real2_k * den2 / (den2 * den2)
        den2_k           = den2_k - perm_real2_k * del2 / (den2 * den2)
        del1_k           = del1_k + perm_real1_k * den1 / (den1 * den1)
        den1_k           = den1_k - perm_real1_k * del1 / (den1 * den1)
        f2_k             = f2_k + den2_k * 2 * fen_sq * f2
        f1_k             = f1_k + den1_k * 2 * fen_sq * f1
!Static permittivity estatic = del1+del2+einf
        tcelsius_k       = tcelsius_k + c(19) * einf_k
        tcelsius_k       = tcelsius_k + del2_k * c(13)
        tcelsius_sq_k    = del2_k * c(14)
        tcelsius_cu_k    = del2_k * c(15)
        tcelsius_k       = tcelsius_k + del1_k * c(9)
        tcelsius_sq_k    = tcelsius_sq_k + del1_k * c(10)
        tcelsius_cu_k    = tcelsius_cu_k + del1_k * c(11)
!Define two relaxation frequencies, f1 and f2
        tcelsius_k       = tcelsius_k + f2_k * c(5)
        tcelsius_sq_k    = tcelsius_sq_k + f2_k * c(6)
        tcelsius_cu_k    = tcelsius_cu_k + f2_k * c(7)
        tcelsius_k       = tcelsius_k + f1_k * c(2)
        tcelsius_sq_k    = tcelsius_sq_k + f1_k * c(3)
!Set values for temperature polynomials (convert from kelvin to celsius)
        tcelsius_k       = tcelsius_k + tcelsius_cu_k * 3 * tcelsius_sq
        tcelsius_k       = tcelsius_k + tcelsius_sq_k * 2 * tcelsius
      CASE (4:max_fastem_version)
        jreflectivity_k(:)  = reflectstokes_k(i, :)
        jemissivity_k(:)    = emissstokes_k(i, :)
        Zenith_Angle        = 1.0_JPRB / geom%seczen
        Zenith_Angle        = acos(Zenith_Angle) * 180.0_JPRB / pi
        Salinity            = prof%skin%salinity
        Transmittance       = transmission_aux%thermal_path1%tau_surf(0, i)
        Wind_Speed          = wind10
! relative azimuth angle in degree
        Rel_Azimuth         = (wind10_direction * 180.0_JPRB / pi - prof%azangle)
        Rel_Azimuth_k       = 0.0_JPRB
        Transmittance_k     = 0.0_JPRB
        Wind_Speed_k        = 0.0_JPRB
        Temperature_k       = 0.0_JPRB
        Salinity_k          = 0.0_JPRB
        CALL rttov_fastem5_k( &
              & opts%rt_mw%fastem_version,&
              & freq_ghz,         &
              & Zenith_Angle,     &
              & prof%skin%t,      &
              & Salinity,         &
              & Wind_Speed,       &
              & jemissivity_k,    &
              & jreflectivity_k,  &
              & Temperature_k,    &
              & Salinity_k,       &
              & Wind_Speed_k,     &
              & jemissivity,      &
              & jreflectivity,    &
              & Transmittance,    &
              & Rel_Azimuth,      &
              & Transmittance_k,  &
              & Rel_Azimuth_k,    & ! Output
              & Supply_Foam_Fraction = opts%rt_mw%supply_foam_fraction, &
              & Foam_Fraction = prof%skin%foam_fraction, &
              & Foam_Fraction_k = prof_k%skin%foam_fraction)
        prof_k%skin%salinity              = prof_k%skin%salinity + Salinity_k
        tcelsius_k                        = tcelsius_k + Temperature_k
        wind10_sq_k                       = 0.0_JPRB
        wind10_direction_k              = wind10_direction_k + Rel_Azimuth_k * 180.0_JPRB / pi
        wind10_K = wind10_K + Wind_Speed_k
        transmission_aux_k%thermal_path1%tau_surf(0, i) = transmission_aux_k%thermal_path1%tau_surf(0, i) + Transmittance_k
      END SELECT
      IF (opts%rt_mw%fastem_version > 0) prof_k%skin%t = prof_k%skin%t + tcelsius_k
      windangle_k   = wind10_direction_k * quadcof(iquadrant, 2)
      windratio_k   = 0.0_JPRB
      windratio_k   = windangle_k / (1.0_JPRB + windratio * windratio)
      IF (abs(prof%s2m%v) >= windlimit) THEN
        prof_k%s2m%u = prof_k%s2m%u + windratio_k * prof%s2m%v / (prof%s2m%v * prof%s2m%v)
        prof_k%s2m%v = prof_k%s2m%v - windratio_k * prof%s2m%u / (prof%s2m%v * prof%s2m%v)
      ELSE
        IF (abs(prof%s2m%u) > windlimit) THEN
          prof_k%s2m%u = prof_k%s2m%u + windscale * windratio_k
        ENDIF
      ENDIF
!       wind10_k = wind10_k + wind10_sq_k * 2 * wind10
      IF (wind10 > 0._JPRB) THEN
        wind10_sq_k = wind10_sq_k + 0.5_JPRB * wind10_k / wind10
      ELSE
        wind10_sq_k = 0.0_JPRB
      ENDIF
      prof_k%s2m%u          = prof_k%s2m%u + 2 * wind10_sq_k * prof%s2m%u
      prof_k%s2m%v          = prof_k%s2m%v + 2 * wind10_sq_k * prof%s2m%v
      prof_k%skin%fastem(:) = 0._JPRB
    ELSE
!--------------------
!2. Land/ice surfaces
!--------------------
!Coherent surface scattering model coefficients (input with the profile)
      perm_static           = prof%skin%fastem(1)
      perm_infinite         = prof%skin%fastem(2)
      freqr = prof%skin%fastem(3)
      small_rough           = prof%skin%fastem(4)
      large_rough           = prof%skin%fastem(5)
      chan = chanprof(i)%chan
      freq_ghz              = coef%frequency_ghz(chan)
!Simple Debye + Fresnel model gives reflectivities
      fen = freq_ghz / freqr
      fen_sq                = fen * fen
      den1 = 1.0_JPRB + fen_sq
      perm_Real             = (perm_static + perm_infinite * fen_sq) / den1
      perm_imag             = fen * (perm_static - perm_infinite) / den1
      permittivity          = Cmplx(perm_Real, perm_imag, jprb)
      perm1 = sqrt(permittivity - geom%sinzen_sq)
      perm2 = permittivity * geom%coszen
      rhth = (geom%coszen - perm1) / (geom%coszen + perm1)
      rvth = (perm2 - perm1) / (perm2 + perm1)
!    fresnel_v_real = dble(rvth)
      fresnel_v_Real        = Real(rvth)
      fresnel_v_imag        = Aimag(rvth)
      fresnel(1)            = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
!    fresnel_h_real = dble(rhth)
      fresnel_h_Real        = Real(rhth)
      fresnel_h_imag        = Aimag(rhth)
      fresnel(2)            = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
!Small scale roughness correction
      delta = 4.0_JPRB * pi * coef%ff_cwn(chan) * 0.1_JPRB * small_rough
      delta2                = delta * delta
      small_rough_cor       = Exp( - delta2 * geom%coszen_sq)
!Large scale roughness correction
      qdepol                = 0.35_JPRB - 0.35_JPRB * Exp( - 0.60_JPRB * freq_ghz * large_rough * large_rough)
      emissfactor_v         = 1.0_JPRB - fresnel(1) * small_rough_cor
      emissfactor_h         = 1.0_JPRB - fresnel(2) * small_rough_cor
      emissfactor           = emissfactor_h - emissfactor_v
      emissstokes(i, 1)     = emissfactor_v + qdepol * emissfactor
      emissstokes(i, 2)     = emissfactor_h - qdepol * emissfactor
!reflect_v(i)  = 1.0_JPRB - emiss_v(i)
!reflect_h(i)  = 1.0_JPRB - emiss_h(i)
!.......end of forward part....................................
!
! * Now run K code of fastem
!
      emissstokes_k(i, 2)   = emissstokes_k(i, 2) - reflectstokes_k(i, 2)
      emissstokes_k(i, 1)   = emissstokes_k(i, 1) - reflectstokes_k(i, 1)
      emissfactor_h_k       = emissstokes_k(i, 2)
      qdepol_k              =  - emissstokes_k(i, 2) * emissfactor
      emissfactor_k         =  - emissstokes_k(i, 2) * qdepol
      emissfactor_v_k       = emissstokes_k(i, 1)
      qdepol_k              = qdepol_k + emissstokes_k(i, 1) * emissfactor
      emissfactor_k         = emissfactor_k + emissstokes_k(i, 1) * qdepol
      emissfactor_v_k       = emissfactor_v_k - emissfactor_k
      emissfactor_h_k       = emissfactor_h_k + emissfactor_k
      fresnel_h_k           =  - emissfactor_h_k * small_rough_cor
      small_rough_cor_k     =  - emissfactor_h_k * fresnel(2)
      fresnel_v_k           =  - emissfactor_v_k * small_rough_cor
      small_rough_cor_k     = small_rough_cor_k - emissfactor_v_k * fresnel(1)
!Large scale roughness correction
      large_rough_k         = qdepol_k * 0.35_JPRB * 0.60_JPRB * freq_ghz * 2 * large_rough *      &
        & Exp( - 0.60_JPRB * freq_ghz * large_rough * large_rough)
!Small scale roughness correction
      delta2_k              =  - small_rough_cor_k * geom%coszen_sq * small_rough_cor
      delta_k               = delta2_k * 2 * delta
      small_rough_k         = 4.0_JPRB * pi * coef%ff_cwn(chan) * 0.1_JPRB * delta_k
!1.3.1) Fresnel reflection coefficients
!Simple Debye + Fresnel model gives reflectivities
!------
      fresnel_h_real_k      = fresnel_h_k * 2 * fresnel_h_real
      fresnel_h_imag_k      = fresnel_h_k * 2 * fresnel_h_imag
      rhth_k                = CMPLX(fresnel_h_real_k,  - fresnel_h_imag_k, jprb)
      fresnel_v_real_k      = fresnel_v_k * 2 * fresnel_v_real
      fresnel_v_imag_k      = fresnel_v_k * 2 * fresnel_v_imag
      rvth_k                = CMPLX(fresnel_v_real_k,  - fresnel_v_imag_k, jprb)
      perm1_k               =  - rvth_k * 2 * perm2 / (perm2 + perm1) ** 2
      perm2_k               = rvth_k * 2 * perm1 / (perm2 + perm1) ** 2
      perm1_k               = perm1_k - rhth_k * 2 * geom%coszen / (geom%coszen + perm1) ** 2
      permittivity_k        = perm2_k * geom%coszen
      permittivity_k        = permittivity_k + perm1_k * 0.5_JPRB / perm1
      perm_Real_k           = Real(permittivity_k)
      perm_imag_k           =  - Aimag(permittivity_k)
      fen_k = perm_imag_k * (perm_static - perm_infinite) / den1
      perm_static_k         = perm_imag_k * fen / den1
      perm_infinite_k       =  - perm_imag_k * fen / den1
      den1_k                =  - perm_imag_k * fen * (perm_static - perm_infinite) / (den1 * den1)
      perm_static_k         = perm_static_k + perm_real_k / den1
      perm_infinite_k       = perm_infinite_k + perm_real_k * fen_sq / den1
      fen_sq_k              = perm_real_k * perm_infinite / den1
      den1_k                = den1_k - perm_real_k * (perm_static + perm_infinite * fen_sq) / (den1 * den1)
      fen_sq_k              = fen_sq_k + den1_k
      fen_k = fen_k + fen_sq_k * 2 * fen
      freqr_k               =  - fen_k * freq_ghz / freqr ** 2
      prof_k%skin%fastem(1) = prof_k%skin%fastem(1) + perm_static_k
      prof_k%skin%fastem(2) = prof_k%skin%fastem(2) + perm_infinite_k
      prof_k%skin%fastem(3) = prof_k%skin%fastem(3) + freqr_k
      prof_k%skin%fastem(4) = prof_k%skin%fastem(4) + small_rough_k
      prof_k%skin%fastem(5) = prof_k%skin%fastem(5) + large_rough_k
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_mw_k
