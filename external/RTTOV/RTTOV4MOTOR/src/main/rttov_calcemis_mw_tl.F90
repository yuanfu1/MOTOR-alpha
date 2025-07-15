! Description:
!> @file
!!   TL of MW surface emissivity calculation
!
!> @brief
!!   TL of MW surface emissivity calculation
!!
!! @param[in]     opts                options to configure the simulations
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in]     profiles_tl         profile perturbations
!! @param[in]     geometry            internal geometry structure
!! @param[in]     coef                optical depth coefficients structure
!! @param[in]     chanprof            specifies channels and profiles to simulate
!! @param[in]     transmission_aux    RTTOV internal auxiliary transmission structure
!! @param[in]     transmission_aux_tl auxiliary transmission perturbations
!! @param[in]     calcemis            flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity_tl       updated with emissivity TL on exit
!! @param[in,out] reflectivity_tl     updated with surface reflectance TL on exit
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
SUBROUTINE rttov_calcemis_mw_tl( &
            & opts,                &
            & profiles,            &
            & profiles_tl,         &
            & geometry,            &
            & coef,                &
            & chanprof,            &
            & transmission_aux,    &
            & transmission_aux_tl, &
            & calcemis,            &
            & emissivity_tl,       &
            & reflectivity_tl)

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
  USE rttov_tessem_mod, ONLY : rttov_tessem_tl
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options),          INTENT(IN)            :: opts
  TYPE(rttov_chanprof),         INTENT(IN)            :: chanprof(:)
  TYPE(rttov_profile),          INTENT(IN),    TARGET :: profiles(:)
  TYPE(rttov_geometry),         INTENT(IN),    TARGET :: geometry(size(profiles))
  TYPE(rttov_coef),             INTENT(IN)            :: coef
  TYPE(rttov_transmission_aux), INTENT(IN)            :: transmission_aux
  LOGICAL(KIND=jplm),           INTENT(IN)            :: calcemis   (size(chanprof))
  TYPE(rttov_profile),          INTENT(IN),    TARGET :: profiles_tl(size(profiles))
  TYPE(rttov_transmission_aux), INTENT(IN)            :: transmission_aux_tl
  REAL(KIND=jprb),              INTENT(INOUT)         :: emissivity_tl  (size(chanprof))
  REAL(KIND=jprb),              INTENT(INOUT)         :: reflectivity_tl(size(chanprof))
!INTF_END

#include "rttov_fastem5_tl.interface"

!local constants:
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
  REAL(KIND=jprb) :: del1           , del2
  REAL(KIND=jprb) :: einf
  REAL(KIND=jprb) :: fen, fen_sq
  REAL(KIND=jprb) :: den1           , den2
  REAL(KIND=jprb) :: perm_free
  REAL(KIND=jprb) :: sigma
  REAL(KIND=jprb) :: perm_real1     , perm_real2
  REAL(KIND=jprb) :: perm_imag1     , perm_imag2    , perm_imag3
  REAL(KIND=jprb) :: perm_Real      , perm_imag
  REAL(KIND=jprb) :: perm_static    , perm_infinite
  REAL(KIND=jprb) :: freq_ghz       , freq_ghz_sq
  REAL(KIND=jprb) :: fresnel_v_Real , fresnel_v_imag
  REAL(KIND=jprb) :: fresnel_h_Real , fresnel_h_imag
  REAL(KIND=jprb) :: fresnel_v      , fresnel_h
  REAL(KIND=jprb) :: small_rough_cor, foam_cor
  REAL(KIND=jprb) :: large_rough_cor(2)
  REAL(KIND=jprb) :: small_rough     , large_rough
  REAL(KIND=jprb) :: variance        , varm
  REAL(KIND=jprb) :: wind10
  REAL(KIND=jprb) :: wind10_sq       , windsec
  REAL(KIND=jprb) :: wind10_direction, windangle     , windratio ! Note wind azimuth is in radians
  REAL(KIND=jprb) :: opdpsfc         , freqr
  REAL(KIND=jprb) :: zrough_v        , zrough_h
  REAL(KIND=jprb) :: zreflmod_v      , zreflmod_h
  REAL(KIND=jprb) :: delta           , delta2
  REAL(KIND=jprb) :: qdepol          , emissfactor
  REAL(KIND=jprb) :: emissfactor_v   , emissfactor_h
  REAL(KIND=jprb) :: emissstokes     (size(chanprof), 4)
  REAL(KIND=jprb) :: emissstokes_tl  (size(chanprof), 4)
  REAL(KIND=jprb) :: reflectstokes_tl(size(chanprof), 4)
  REAL(KIND=jprb) :: zc(12), zx(9)
  REAL(KIND=jprb) :: azimuthal_emiss, u19, phi       , dfreq
  REAL(KIND=jprb) :: tbfixed      (4, 4, 3)                      ! Surface brightness temperature azimuthal
                                                                 !   variation terms for 37, 19, 10, 7 GHz
  REAL(KIND=jprb) :: efixed       (4, 4, 3)                      ! Emissivity azimuthal variation terms for
                                                                 !   7, 10, 19, 37 GHz
  REAL(KIND=jprb) :: einterpolated(4, 3   )                      ! Emissivity azimuthal variation terms for
                                                                 !   interpolated to required frequency
  REAL(KIND=jprb) :: a1e, a2e, a3e                               ! coefficients used in azimuthal emissivity model
  REAL(KIND=jprb), POINTER :: c(:)
  COMPLEX(KIND=jprb) :: perm1       , perm2
  COMPLEX(KIND=jprb) :: rhth        , rvth
  COMPLEX(KIND=jprb) :: permittivity
  INTEGER(KIND=jpim) :: i, j, chan    , istokes, m
  INTEGER(KIND=jpim) :: iquadrant                              ! Determines which quadrant (NE, SE, SW, NW) the
                                                               !   wind is blowing to
  INTEGER(KIND=jpim) :: pol_id                                 ! polarisation indice
  INTEGER(KIND=jpim) :: ifreq       , i_freq, j_stokes, ich    ! indices used in azimuthal emissivity model
  INTEGER(KIND=jpim) :: jcof        , jcofm1
  TYPE(rttov_profile),  POINTER :: prof
  TYPE(rttov_profile),  POINTER :: prof_tl
  TYPE(rttov_geometry), POINTER :: geom
! TL variables
  REAL   (KIND=jprb) :: tcelsius_tl
  REAL   (KIND=jprb) :: tcelsius_sq_tl
  REAL   (KIND=jprb) :: tcelsius_cu_tl
  REAL   (KIND=jprb) :: f1_tl             , f2_tl
  REAL   (KIND=jprb) :: del1_tl           , del2_tl
  REAL   (KIND=jprb) :: einf_tl
  REAL   (KIND=jprb) :: fen_tl            , fen_sq_tl
  REAL   (KIND=jprb) :: den1_tl           , den2_tl
  REAL   (KIND=jprb) :: sigma_tl
  REAL   (KIND=jprb) :: perm_real1_tl     , perm_real2_tl
  REAL   (KIND=jprb) :: perm_imag1_tl     , perm_imag2_tl    , perm_imag3_tl
  REAL   (KIND=jprb) :: perm_Real_tl      , perm_imag_tl
  REAL   (KIND=jprb) :: perm_static_tl    , perm_infinite_tl
  REAL   (KIND=jprb) :: fresnel_v_Real_tl , fresnel_v_imag_tl
  REAL   (KIND=jprb) :: fresnel_h_Real_tl , fresnel_h_imag_tl
  REAL   (KIND=jprb) :: fresnel_v_tl      , fresnel_h_tl
  REAL   (KIND=jprb) :: small_rough_cor_tl, foam_cor_tl
  REAL   (KIND=jprb) :: large_rough_cor_tl(2)
  REAL   (KIND=jprb) :: small_rough_tl     , large_rough_tl
  REAL   (KIND=jprb) :: variance_tl        , varm_tl
  REAL   (KIND=jprb) :: wind10_tl
  REAL   (KIND=jprb) :: wind10_sq_tl       , windsec_tl
  REAL   (KIND=jprb) :: wind10_direction_tl, windangle_tl     , windratio_tl  ! Note wind azimuth is in radians
  REAL   (KIND=jprb) :: opdpsfc_tl         , freqr_tl
  REAL   (KIND=jprb) :: zrough_v_tl        , zrough_h_tl
  REAL   (KIND=jprb) :: zreflmod_v_tl      , zreflmod_h_tl
  REAL   (KIND=jprb) :: delta_tl           , delta2_tl
  REAL   (KIND=jprb) :: qdepol_tl          , emissfactor_tl
  REAL   (KIND=jprb) :: emissfactor_v_tl   , emissfactor_h_tl
  REAL   (KIND=jprb) :: zx_tl(9)
  REAL   (KIND=jprb) :: azimuthal_emiss_tl, u19_tl, phi_tl
  REAL   (KIND=jprb) :: tbfixed_tl      (4, 4, 3)              ! Surface brightness temperature azimuthal
                                                               !   variation terms for 37, 19, 10, 7 GHz
  REAL   (KIND=jprb) :: efixed_tl       (4, 4, 3)              ! Emissivity azimuthal variation terms for
                                                               !   7, 10, 19, 37 GHz
  REAL   (KIND=jprb) :: einterpolated_tl(4, 3   )              ! Emissivity azimuthal variation terms for
                                                               !   interpolated to required frequency
  REAL   (KIND=jprb) :: a1e_tl           , a2e_tl, a3e_tl      ! coefficients used in azimuthal emissivity model
  COMPLEX(KIND=jprb) :: perm1_tl         , perm2_tl
  COMPLEX(KIND=jprb) :: rhth_tl          , rvth_tl
  COMPLEX(KIND=jprb) :: permittivity_tl
! rttov_fastem4_tl   internal variables
  REAL   (KIND=jprb) :: Zenith_Angle     , Salinity         , Transmittance, Rel_Azimuth, Wind_Speed
  REAL(KIND=jprb), DIMENSION(4) :: jemissivity, jreflectivity
  REAL(KIND=jprb) :: Salinity_tl   , Transmittance_tl, Rel_Azimuth_tl
  REAL(KIND=jprb) :: Temperature_tl, Wind_Speed_tl
  REAL(KIND=jprb), DIMENSION(4) :: jemissivity_tl, jreflectivity_tl
  INTEGER(KIND=jpim) :: nchannels   ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! The rttov_check_options subroutine guarantees the selected FASTEM version is valid
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels         = size(chanprof)

  DO i = 1, nchannels
    IF (.NOT. calcemis(i)) CYCLE
    chan = chanprof(i)%chan
    prof => profiles(chanprof(i)%prof)
    prof_tl => profiles_tl(chanprof(i)%prof)
    geom => geometry(chanprof(i)%prof)
    pol_id = coef%fastem_polar(chan) + 1
!-------------------------------
!0. Point to fastem coefficients
!-------------------------------
    c => fastem3_coef
!---------------
!1. Sea surfaces
!---------------
    IF (prof%skin%surftype == surftype_sea) THEN
!-------------------------------------------
!1.1 Calculate channel independent variables
!-------------------------------------------
! no TL on wind direction, but TL on wind speed
      wind10_sq    = prof%s2m%u * prof%s2m%u + prof%s2m%v * prof%s2m%v
      wind10       = Sqrt(wind10_sq)
      wind10_sq_tl = 2.0_JPRB * prof%s2m%u * prof_tl%s2m%u + 2.0_JPRB * prof%s2m%v * prof_tl%s2m%v
      freq_ghz     = coef%frequency_ghz(chan)
      freq_ghz_sq  = freq_ghz * freq_ghz
      IF (wind10 > 0._JPRB) THEN
        wind10_tl = 0.5_JPRB * wind10_sq_tl / wind10
      ELSE
        wind10_tl = 0.0_JPRB
      ENDIF
      IF (prof%s2m%u >= 0.0_JPRB .AND. prof%s2m%v >= 0.0_JPRB) iquadrant = 1
      IF (prof%s2m%u >= 0.0_JPRB .AND. prof%s2m%v < 0.0_JPRB ) iquadrant = 2
      IF (prof%s2m%u < 0.0_JPRB .AND. prof%s2m%v >= 0.0_JPRB ) iquadrant = 4
      IF (prof%s2m%u < 0.0_JPRB .AND. prof%s2m%v < 0.0_JPRB  ) iquadrant = 3
      IF (abs(prof%s2m%v) >= windlimit) THEN
        windratio = prof%s2m%u / prof%s2m%v
      ELSE
        windratio = 0.0_JPRB
        IF (abs(prof%s2m%u) > windlimit) THEN
          windratio = windscale * prof%s2m%u
        ENDIF
      ENDIF
      windangle        = atan(abs(windratio))
      wind10_direction = quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2)
      windratio_tl     = 0.0_JPRB
      IF (abs(prof%s2m%v) >= windlimit) THEN
        windratio_tl = (prof%s2m%v * prof_tl%s2m%u - prof_tl%s2m%v * prof%s2m%u) / (prof%s2m%v * prof%s2m%v)
      ELSE
        windratio_tl = 0.0_JPRB
        IF (abs(prof%s2m%u) > windlimit) THEN
          windratio_tl = windscale * prof_tl%s2m%u
        ENDIF
      ENDIF
      windangle_tl        = windratio_tl / (1.0_JPRB + windratio * windratio)
      wind10_direction_tl = windangle_tl * quadcof(iquadrant, 2)
      SELECT CASE (opts%rt_mw%fastem_version)
      CASE (0)
        CALL rttov_tessem_tl(freq_ghz, prof%zenangle, wind10, prof%skin%t, prof%skin%salinity, &
                             wind10_tl, prof_tl%skin%t, prof_tl%skin%salinity, &
                             emissstokes_tl(i,2), emissstokes_tl(i,1))
        emissstokes_tl(i,3:4) = 0._jprb
        reflectstokes_tl(i,:) = - emissstokes_tl(i,:)
      CASE (1:3)
        windsec           = wind10 * geom%seczen
        windsec_tl        = wind10_tl * geom%seczen
!Set values for temperature polynomials (convert from kelvin to celsius)
        tcelsius          = prof%skin%t - 273.15_JPRB
        tcelsius_sq       = tcelsius * tcelsius                                                                 !quadratic
        tcelsius_cu       = tcelsius_sq * tcelsius                                                              !cubic
        tcelsius_tl       = prof_tl%skin%t
        tcelsius_sq_tl    = 2 * tcelsius * tcelsius_tl
        tcelsius_cu_tl    = 3 * tcelsius_sq * tcelsius_tl
!Define two relaxation frequencies, f1 and f2
        f1 = c(1) + c(2) * tcelsius + c(3) * tcelsius_sq
        f2 = c(4) + c(5) * tcelsius + c(6) * tcelsius_sq + c(7) * tcelsius_cu
        f1_tl             = c(2) * tcelsius_tl + c(3) * tcelsius_sq_tl
        f2_tl             = c(5) * tcelsius_tl + c(6) * tcelsius_sq_tl + c(7) * tcelsius_cu_tl
!Static permittivity estatic = del1+del2+einf
        del1 = c(8) + c(9) * tcelsius + c(10) * tcelsius_sq + c(11) * tcelsius_cu
        del2 = c(12) + c(13) * tcelsius + c(14) * tcelsius_sq + c(15) * tcelsius_cu
        einf = c(18) + c(19) * tcelsius
        del1_tl           = c(9) * tcelsius_tl + c(10) * tcelsius_sq_tl + c(11) * tcelsius_cu_tl
        del2_tl           = c(13) * tcelsius_tl + c(14) * tcelsius_sq_tl + c(15) * tcelsius_cu_tl
        einf_tl           = c(19) * tcelsius_tl
!-----------------------------------------------------
!1.2 calculate permittivity using double-debye formula
!-----------------------------------------------------
        fen = 2.0_JPRB * c(20) * freq_ghz * 0.001_JPRB
        fen_sq            = fen * fen
        den1 = 1.0_JPRB + fen_sq * f1 * f1
        den2 = 1.0_JPRB + fen_sq * f2 * f2
        perm_real1        = del1 / den1
        perm_real2        = del2 / den2
        perm_imag1        = del1 * fen * f1 / den1
        perm_imag2        = del2 * fen * f2 / den2
        perm_free         = 8.854E-03_JPRB
        sigma             = 2.906_JPRB + 0.09437_JPRB * tcelsius
        perm_imag3        = sigma / (2.0_JPRB * c(20) * perm_free * freq_ghz)
        perm_Real         = perm_real1 + perm_real2 + einf
!        perm_imag    = perm_imag1 + perm_imag2 + perm_imag3 + perm_imag3
        perm_imag         = perm_imag1 + perm_imag2 + perm_imag3
        permittivity      = Cmplx(perm_Real, perm_imag, jprb)
        den1_tl           = 2 * fen_sq * f1 * f1_tl
        den2_tl           = 2 * fen_sq * f2 * f2_tl
        perm_real1_tl     = (den1 * del1_tl - del1 * den1_tl) / (den1 * den1)
        perm_real2_tl     = (den2 * del2_tl - del2 * den2_tl) / (den2 * den2)
        perm_imag1_tl     = fen * (den1 * (del1_tl * f1 + del1 * f1_tl) - (del1 * f1 * den1_tl)) / (den1 * den1)
        perm_imag2_tl     = fen * (den2 * (del2_tl * f2 + del2 * f2_tl) - (del2 * f2 * den2_tl)) / (den2 * den2)
        sigma_tl          = 0.09437_JPRB * tcelsius_tl
        perm_imag3_tl     = sigma_tl / (2.0_JPRB * c(20) * perm_free * freq_ghz)
        perm_Real_tl      = perm_real1_tl + perm_real2_tl + einf_tl
!        perm_imag_tl    = perm_imag1_tl + perm_imag2_tl + perm_imag3_tl + perm_imag3_tl
        perm_imag_tl      = perm_imag1_tl + perm_imag2_tl + perm_imag3_tl
        permittivity_tl   = Cmplx(perm_Real_tl, perm_imag_tl, jprb)
!-------------------------------------------------------------
!1.3 calculate complex reflection coefficients and corrections
!-------------------------------------------------------------
!1.3.1) Fresnel reflection coefficients
!------
        perm1             = sqrt(permittivity - geom%sinzen_sq)
        perm2             = permittivity * geom%coszen
        rhth = (geom%coszen - perm1) / (geom%coszen + perm1)
        rvth = (perm2 - perm1) / (perm2 + perm1)
!    fresnel_v_real = dble(rvth)
        fresnel_v_Real    = Real(rvth)
        fresnel_v_imag    = Aimag(rvth)
        fresnel_v         = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
!    fresnel_h_real = dble(rhth)
        fresnel_h_Real    = Real(rhth)
        fresnel_h_imag    = Aimag(rhth)
        fresnel_h         = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
        perm1_tl          = 0.5_JPRB * permittivity_tl / perm1
        perm2_tl          = permittivity_tl * geom%coszen
        rhth_tl           =  - 2 * geom%coszen * perm1_tl / (geom%coszen + perm1) ** 2
        rvth_tl           = 2 * (perm1 * perm2_tl - perm1_tl * perm2) / (perm2 + perm1) ** 2
!    fresnel_v_real_tl = dble(rvth_tl)
        fresnel_v_Real_tl = Real(rvth_tl)
        fresnel_v_imag_tl = Aimag(rvth_tl)
        fresnel_v_tl      = 2 * fresnel_v_Real * fresnel_v_Real_tl + 2 * fresnel_v_imag * fresnel_v_imag_tl
!    fresnel_h_real_tl = dble(rhth_tl)
        fresnel_h_Real_tl = Real(rhth_tl)
        fresnel_h_imag_tl = Aimag(rhth_tl)
        fresnel_h_tl      = 2 * fresnel_h_Real * fresnel_h_Real_tl + 2 * fresnel_h_imag * fresnel_h_imag_tl
!1.3.2) Small scale correction to reflection coefficients
!------
        IF (freq_ghz >= 15.0) THEN
          small_rough_cor    = Exp(c(21) * wind10 * geom%coszen_sq / (freq_ghz_sq))
          small_rough_cor_tl = small_rough_cor * c(21) * wind10_tl * geom%coszen_sq / (freq_ghz_sq)
        ELSE
          small_rough_cor    = 1.0
          small_rough_cor_tl = 0.0
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
        large_rough_cor(1)    =                                                                                              &
          & (zc(1) + zc(2) * geom%seczen + zc(3) * geom%seczen_sq + zc(4) * wind10 + zc(5) * wind10_sq + zc(6) * windsec) /  &
          & 100._JPRB
        large_rough_cor(2)    =                                                                                              &
          & (zc(7) + zc(8) * geom%seczen + zc(9) * geom%seczen_sq + zc(10) * wind10 + zc(11) * wind10_sq + zc(12) * windsec) &
          &  / 100._JPRB
!    large_rough_cor(:) = large_rough_cor(:) * 0.01
        large_rough_cor_tl(1) = (zc(4) * wind10_tl + zc(5) * wind10_sq_tl + zc(6) * windsec_tl) / 100._JPRB
        large_rough_cor_tl(2) = (zc(10) * wind10_tl + zc(11) * wind10_sq_tl + zc(12) * windsec_tl) / 100._JPRB
! For Fastem-3 do not compute rough surface effects if theta > 60 degrees
        IF (opts%rt_mw%fastem_version <= 2.0_JPRB .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) THEN
          emissstokes(i, 1)    = 1.0_JPRB - fresnel_v * small_rough_cor + large_rough_cor(1)
          emissstokes(i, 2)    = 1.0_JPRB - fresnel_h * small_rough_cor + large_rough_cor(2)
          emissstokes_tl(i, 1) =      &
            &  - fresnel_v_tl * small_rough_cor - fresnel_v * small_rough_cor_tl + large_rough_cor_tl(1)
          emissstokes_tl(i, 2) =      &
            &  - fresnel_h_tl * small_rough_cor - fresnel_h * small_rough_cor_tl + large_rough_cor_tl(2)
        ELSE
          emissstokes(i, 1)    = 1.0_JPRB - fresnel_v
          emissstokes(i, 2)    = 1.0_JPRB - fresnel_h
          emissstokes_tl(i, 1) =  - fresnel_v_tl
          emissstokes_tl(i, 2) =  - fresnel_h_tl
        ENDIF
        emissstokes(i, 3)    = 0.0_JPRB
        emissstokes(i, 4)    = 0.0_JPRB
        emissstokes_tl(i, 3) = 0.0_JPRB
        emissstokes_tl(i, 4) = 0.0_JPRB
!Apply foam correction
        IF (.NOT. opts%rt_mw%supply_foam_fraction) THEN
          foam_cor             = c(22) * (wind10 ** c(23))
          foam_cor_tl          = c(22) * c(23) * wind10_tl * (wind10 ** (c(23) - 1.0_JPRB))
        ELSE
          foam_cor             = prof%skin%foam_fraction
          foam_cor_tl          = prof_tl%skin%foam_fraction
        ENDIF
          
! Be careful do TL first because the next 2 lines of the direct model
! have variables in input/output of the statement
        emissstokes_tl(i, 1) =      &
          & emissstokes_tl(i, 1) - foam_cor_tl * emissstokes(i, 1) - foam_cor * emissstokes_tl(i, 1) + foam_cor_tl
        emissstokes_tl(i, 2) =      &
          & emissstokes_tl(i, 2) - foam_cor_tl * emissstokes(i, 2) - foam_cor * emissstokes_tl(i, 2) + foam_cor_tl
        emissstokes(i, 1)    = emissstokes(i, 1) - foam_cor * emissstokes(i, 1) + foam_cor
        emissstokes(i, 2)    = emissstokes(i, 2) - foam_cor * emissstokes(i, 2) + foam_cor
        IF (opts%rt_mw%fastem_version == 3) THEN
! Add azimuthal component from Fuzhong Weng (NOAA/NESDIS) based on work by Dr. Gene Poe (NRL)
! Assume 19m wind = 10m wind for now (fix later)
          u19 = wind10
! Angle between wind direction and satellite azimuthal view angle
          phi = pi - wind10_direction + prof%azangle * pi / 180.0_JPRB
          phi_tl            =  - 1.0_JPRB * wind10_direction_tl
          u19_tl            = wind10_tl
          tbfixed(:,:,:)    = 0.0_JPRB
          tbfixed_tl(:,:,:) = 0.0_JPRB
          DO ich = 0, 15
            a1e = c(141 + ich * 12) + u19 * (c(142 + ich * 12) + u19 * (c(143 + ich * 12) + u19 * c(144 + ich * 12)))
            a2e = c(145 + ich * 12) + u19 * (c(146 + ich * 12) + u19 * (c(147 + ich * 12) + u19 * c(148 + ich * 12)))
            a3e = c(149 + ich * 12) + u19 * (c(150 + ich * 12) + u19 * (c(151 + ich * 12) + u19 * c(152 + ich * 12)))
            a1e_tl = u19_tl * (c(142 + ich * 12) + u19 * (2.0 * c(143 + ich * 12) + 3.0 * u19 * c(144 + ich * 12)))
            a2e_tl = u19_tl * (c(146 + ich * 12) + u19 * (2.0 * c(147 + ich * 12) + 3.0 * u19 * c(148 + ich * 12)))
            a3e_tl = u19_tl * (c(150 + ich * 12) + u19 * (2.0 * c(151 + ich * 12) + 3.0 * u19 * c(152 + ich * 12)))
            i_freq = int(ich / 4) + 1! 37, 19, 10, 7 GHz
            j_stokes                        = mod(ich, 4_jpim) + 1
            tbfixed(j_stokes, i_freq, 1)    = a1e
            tbfixed(j_stokes, i_freq, 2)    = a2e
            tbfixed(j_stokes, i_freq, 3)    = a3e
            tbfixed_tl(j_stokes, i_freq, 1) = a1e_tl
            tbfixed_tl(j_stokes, i_freq, 2) = a2e_tl
            tbfixed_tl(j_stokes, i_freq, 3) = a3e_tl
          ENDDO
          efixed_tl(:,:,:)      = 0.0_JPRB
          einterpolated_tl(:,:) = 0.0_JPRB
          DO M = 1, 3
            DO istokes = 1, 4
              efixed(1, istokes, M)    = tbfixed(istokes, 4, M)   ! 7   GHz
              efixed(2, istokes, M)    = tbfixed(istokes, 3, M)   ! 10  GHz
              efixed(3, istokes, M)    = tbfixed(istokes, 2, M)   ! 19  GHz
              efixed(4, istokes, M)    = tbfixed(istokes, 1, M)   ! 37  GHz
              efixed_tl(1, istokes, M) = tbfixed_tl(istokes, 4, M)! 7  GHz
              efixed_tl(2, istokes, M) = tbfixed_tl(istokes, 3, M)! 10  GHz
              efixed_tl(3, istokes, M) = tbfixed_tl(istokes, 2, M)! 19  GHz
              efixed_tl(4, istokes, M) = tbfixed_tl(istokes, 1, M)! 37  GHz
            ENDDO
! Interpolate results to required frequency based on 7, 10, 19, 37 GHz
            IF (freq_ghz .LE. freqfixed(1)) THEN
              einterpolated(:, M)    = efixed(1, :, M)
              einterpolated_tl(:, M) = efixed_tl(1, :, M)
            ELSE IF (freq_ghz .GE. freqfixed(4)) THEN
              einterpolated(:, M)    = efixed(4, :, M)
              einterpolated_tl(:, M) = efixed_tl(4, :, M)
            ELSE
              IF (freq_ghz .LT. freqfixed(2)                                 ) ifreq = 2
              IF (freq_ghz .LT. freqfixed(3) .AND. freq_ghz .GE. freqfixed(2)) ifreq = 3
              IF (freq_ghz .GE. freqfixed(3)                                 ) ifreq = 4
              dfreq = (freq_ghz - freqfixed(ifreq - 1)) / (freqfixed(ifreq) - freqfixed(ifreq - 1))
              einterpolated(:, M)    = efixed(ifreq - 1, :, M) + dfreq * (efixed(ifreq, :, M) - efixed(ifreq - 1, :, M))
              einterpolated_tl(:, M) =      &
                & efixed_tl(ifreq - 1, :, M) + dfreq * (efixed_tl(ifreq, :, M) - efixed_tl(ifreq - 1, :, M))
            ENDIF
          ENDDO
          DO istokes = 1, 4
            azimuthal_emiss    = 0.0_JPRB
            azimuthal_emiss_tl = 0.0_JPRB
            DO M = 1, 3
              IF (istokes .LE. 2) THEN
                azimuthal_emiss    = azimuthal_emiss +      &
                  & einterpolated(istokes, M) * cos(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
                azimuthal_emiss_tl = azimuthal_emiss_tl +                                                                    &
                  & (einterpolated_tl(istokes, M) * cos(m * phi) - einterpolated(istokes, M) * m * sin(m * phi) * phi_tl) *  &
                  & (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
              ELSE
                azimuthal_emiss    = azimuthal_emiss +      &
                  & einterpolated(istokes, M) * sin(m * phi) * (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
                azimuthal_emiss_tl = azimuthal_emiss_tl +                                                                    &
                  & (einterpolated_tl(istokes, M) * sin(m * phi) + einterpolated(istokes, M) * m * cos(m * phi) * phi_tl) *  &
                  & (1.0_JPRB - geom%coszen) / (1.0_JPRB - 0.6018_JPRB)
              ENDIF
            ENDDO
            emissstokes(i, istokes)    = emissstokes(i, istokes) + azimuthal_emiss
            emissstokes_tl(i, istokes) = emissstokes_tl(i, istokes) + azimuthal_emiss_tl
          ENDDO
        ENDIF
! Only apply non-specular correction for Fastem-3 if theta < 60 degrees
        IF ((opts%rt_mw%fastem_version == 2 .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) < 0.9999_JPRB .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) > 0.00001_JPRB) THEN
!Convert windspeed to slope variance using the Cox and Munk model
          variance    = 0.00512_JPRB * wind10 + 0.0030_JPRB
          varm        = variance * c(138)
          variance    = varm * (c(139) * freq_ghz + c(140))
          variance_tl = 0.00512_JPRB * wind10_tl
          varm_tl     = variance_tl * c(138)
          variance_tl = varm_tl * (c(139) * freq_ghz + c(140))
          IF (variance > varm) THEN
            variance    = varm
            variance_tl = varm_tl
          ENDIF
          IF (variance < 0.0_JPRB) THEN
            variance    = 0.0_JPRB
            variance_tl = 0.0_JPRB
          ENDIF
!Compute surface to space optical depth
          opdpsfc     =  - log(transmission_aux%thermal_path1%tau_surf(0, i)) / geom%seczen
          opdpsfc_tl  =  - transmission_aux_tl%thermal_path1%tau_surf(0, i) / &
                      & (transmission_aux%thermal_path1%tau_surf(0, i) * geom%seczen)
!Define nine predictors for the effective angle calculation
          zx(1)       = 1.0_JPRB
          zx(2)       = variance
          zx(4)       = 1.0_JPRB / geom%coszen
          zx(3)       = zx(2) * zx(4)
          zx(5)       = zx(3) * zx(3)
          zx(6)       = zx(4) * zx(4)
          zx(7)       = zx(2) * zx(2)
          zx(8)       = log(opdpsfc)
          zx(9)       = zx(8) * zx(8)
          zx_tl(1)    = 0._JPRB
          zx_tl(2)    = variance_tl
          zx_tl(4)    = 0._JPRB
          zx_tl(3)    = zx_tl(2) * zx(4)
          zx_tl(5)    = 2 * zx_tl(3) * zx(3)
          zx_tl(6)    = 2 * zx_tl(4) * zx(4)
          zx_tl(7)    = 2 * zx_tl(2) * zx(2)
          zx_tl(8)    = opdpsfc_tl / opdpsfc
          zx_tl(9)    = 2 * zx_tl(8) * zx(8)
          zrough_v    = 1.0_JPRB
          zrough_h    = 1.0_JPRB
          zrough_v_tl = 0._JPRB
          zrough_h_tl = 0._JPRB
          DO jcof = 1, 7
            jcofm1      = jcof - 1
!Switched h to v Deblonde SSMIS june 7, 2001
            zrough_h    =      &
              & zrough_h + zx(jcof) * (c(96 + jcofm1 * 3) + zx(8) * c(97 + jcofm1 * 3) + zx(9) * c(98 + jcofm1 * 3))
            zrough_v    =      &
              & zrough_v + zx(jcof) * (c(117 + jcofm1 * 3) + zx(8) * c(118 + jcofm1 * 3) + zx(9) * c(119 + jcofm1 * 3))
            zrough_h_tl = zrough_h_tl + zx(jcof) * (zx_tl(8) * c(97 + jcofm1 * 3) + zx_tl(9) * c(98 + jcofm1 * 3)) +      &
              & zx_tl(jcof) * (c(96 + jcofm1 * 3) + zx(8) * c(97 + jcofm1 * 3) + zx(9) * c(98 + jcofm1 * 3))
            zrough_v_tl = zrough_v_tl + zx(jcof) * (zx_tl(8) * c(118 + jcofm1 * 3) + zx_tl(9) * c(119 + jcofm1 * 3)) +      &
              & zx_tl(jcof) * (c(117 + jcofm1 * 3) + zx(8) * c(118 + jcofm1 * 3) + zx(9) * c(119 + jcofm1 * 3))
          ENDDO
          zreflmod_v             =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v) /            &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          zreflmod_h             =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h) /            &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          zreflmod_v_tl          = transmission_aux_tl%thermal_path1%tau_surf(0, i) * ( -         &
            & zrough_v * transmission_aux%thermal_path1%tau_surf(0, i) ** (zrough_v - 1.0_JPRB) * &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) +                        &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v)) /           &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) ** 2
          zreflmod_v_tl          = zreflmod_v_tl -                                                &
            & (transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v *                        &
            & Log(transmission_aux%thermal_path1%tau_surf(0, i)) * zrough_v_tl) /                 &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          zreflmod_h_tl          = transmission_aux_tl%thermal_path1%tau_surf(0, i) * ( -         &
            & zrough_h * transmission_aux%thermal_path1%tau_surf(0, i) ** (zrough_h - 1.0_JPRB) * &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) +                        &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h)) /           &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i)) ** 2
          zreflmod_h_tl          = zreflmod_h_tl -                                                &
            & (transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h *                        &
            & Log(transmission_aux%thermal_path1%tau_surf(0, i)) * zrough_h_tl) /                 &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          reflectstokes_tl(i, 1) = zreflmod_v_tl * (1.0 - emissstokes(i, 1)) - zreflmod_v * emissstokes_tl(i, 1)
          reflectstokes_tl(i, 2) = zreflmod_h_tl * (1.0 - emissstokes(i, 2)) - zreflmod_h * emissstokes_tl(i, 2)
!          zreflmod_v_tl = 0.0
!          zreflmod_h_tl = 0.0
          reflectstokes_tl(i, 3) =  - 0.5_JPRB *      &
            & ((zreflmod_v_tl + zreflmod_h_tl) * emissstokes(i, 3) + (zreflmod_v + zreflmod_h) * emissstokes_tl(i, 3))
          reflectstokes_tl(i, 4) =  - 0.5_JPRB *      &
            & ((zreflmod_v_tl + zreflmod_h_tl) * emissstokes(i, 4) + (zreflmod_v + zreflmod_h) * emissstokes_tl(i, 4))
!           reflectstokes_tl(i,3)  =     -emissstokes_tl(i,3)
!           reflectstokes_tl(i,4)  =     -emissstokes_tl(i,4)
        ELSE
          reflectstokes_tl(i, :) =  - emissstokes_tl(i, :)
        ENDIF
      CASE (4:max_fastem_version)
        Zenith_Angle     = 1.0_JPRB / geom%seczen
        Zenith_Angle     = acos(Zenith_Angle) * 180.0_JPRB / pi
        Salinity         = prof%skin%salinity
        Salinity_tl      = prof_tl%skin%salinity
        Transmittance    = transmission_aux%thermal_path1%tau_surf(0, i)
        Transmittance_tl = transmission_aux_tl%thermal_path1%tau_surf(0, i)
        Temperature_tl   = prof_tl%skin%t
        Wind_Speed       = wind10
        Wind_Speed_tl    = wind10_tl
! relative azimuth angle in degree
        Rel_Azimuth      = (wind10_direction * 180.0_JPRB / pi - prof%azangle)
        Rel_Azimuth_tl   = wind10_direction_tl * 180.0_JPRB / pi
        CALL rttov_fastem5_TL( &
              & opts%rt_mw%fastem_version,&
              & freq_ghz,         &
              & Zenith_Angle,     &
              & prof%skin%t,      &
              & Salinity,         &
              & Wind_Speed,       &
              & Temperature_tl,   &
              & Salinity_tl,      &
              & Wind_Speed_tl,    &
              & jemissivity,      &
              & jreflectivity,    &
              & jemissivity_tl,   &
              & jreflectivity_tl, &
              & Transmittance,    &
              & Rel_Azimuth,      &
              & Transmittance_tl, &
              & Rel_Azimuth_tl,   & ! Input, may not be used
              & Supply_Foam_Fraction = opts%rt_mw%supply_foam_fraction, &
              & Foam_Fraction = prof%skin%foam_fraction, &
              & Foam_Fraction_tl = prof_tl%skin%foam_fraction)
        emissstokes_tl(i, 1:4)   = jemissivity_tl(1:4)
        reflectstokes_tl(i, 1:4) = jreflectivity_tl(1:4)
      END SELECT
!--------------------
!2. Land/ice surfaces
!--------------------
    ELSE
!Coherent surface scattering model coefficients (input with the profile)
      perm_static            = prof%skin%fastem(1)
      perm_infinite          = prof%skin%fastem(2)
      freqr = prof%skin%fastem(3)
      small_rough            = prof%skin%fastem(4)
      large_rough            = prof%skin%fastem(5)
      freq_ghz               = coef%frequency_ghz(chan)
      perm_static_tl         = prof_tl%skin%fastem(1)
      perm_infinite_tl       = prof_tl%skin%fastem(2)
      freqr_tl               = prof_tl%skin%fastem(3)
      small_rough_tl         = prof_tl%skin%fastem(4)
      large_rough_tl         = prof_tl%skin%fastem(5)
!Simple Debye + Fresnel model gives reflectivities
      fen = freq_ghz / freqr
      fen_sq                 = fen * fen
      den1 = 1.0_JPRB + fen_sq
      perm_Real              = (perm_static + perm_infinite * fen_sq) / den1
      perm_imag              = fen * (perm_static - perm_infinite) / den1
      permittivity           = Cmplx(perm_Real, perm_imag, jprb)
      perm1 = sqrt(permittivity - geom%sinzen_sq)
      perm2 = permittivity * geom%coszen
      rhth = (geom%coszen - perm1) / (geom%coszen + perm1)
      rvth = (perm2 - perm1) / (perm2 + perm1)
!    fresnel_v_real = dble(rvth)
      fresnel_v_Real         = Real(rvth)
      fresnel_v_imag         = Aimag(rvth)
      fresnel_v              = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
!    fresnel_h_real = dble(rhth)
      fresnel_h_Real         = Real(rhth)
      fresnel_h_imag         = Aimag(rhth)
      fresnel_h              = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
      fen_tl                 =  - freq_ghz * freqr_tl / freqr ** 2
      fen_sq_tl              = 2 * fen_tl * fen
      den1_tl                = fen_sq_tl
      perm_Real_tl           = (den1 * (perm_static_tl + perm_infinite_tl * fen_sq + perm_infinite * fen_sq_tl) -      &
        & den1_tl * (perm_static + perm_infinite * fen_sq)) / (den1 * den1)
      perm_imag_tl           = (                                                                         &
        & den1 * (fen_tl * (perm_static - perm_infinite) + fen * (perm_static_tl - perm_infinite_tl)) -  &
        & den1_tl * fen * (perm_static - perm_infinite)) / (den1 * den1)
      permittivity_tl        = Cmplx(perm_Real_tl, perm_imag_tl, jprb)
      perm1_tl               = 0.5_JPRB * permittivity_tl / perm1
      perm2_tl               = permittivity_tl * geom%coszen
      rhth_tl                =  - 2 * geom%coszen * perm1_tl / (geom%coszen + perm1) ** 2
      rvth_tl                = 2 * (perm1 * perm2_tl - perm1_tl * perm2) / (perm2 + perm1) ** 2
!    fresnel_v_real_tl = dble(rvth_tl)
      fresnel_v_Real_tl      = Real(rvth_tl)
      fresnel_v_imag_tl      = Aimag(rvth_tl)
      fresnel_v_tl           = 2 * fresnel_v_Real * fresnel_v_Real_tl + 2 * fresnel_v_imag * fresnel_v_imag_tl
!    fresnel_h_real_tl = dble(rhth_tl)
      fresnel_h_Real_tl      = Real(rhth_tl)
      fresnel_h_imag_tl      = Aimag(rhth_tl)
      fresnel_h_tl           = 2 * fresnel_h_Real * fresnel_h_Real_tl + 2 * fresnel_h_imag * fresnel_h_imag_tl
!Small scale roughness correction
      delta = 4.0_JPRB * pi * coef%ff_cwn(chan) * 0.1_JPRB * small_rough
      delta2                 = delta * delta
      small_rough_cor        = Exp( - delta2 * geom%coszen_sq)
      delta_tl               = 4.0_JPRB * pi * coef%ff_cwn(chan) * 0.1_JPRB * small_rough_tl
      delta2_tl              = 2 * delta * delta_tl
      small_rough_cor_tl     =  - delta2_tl * geom%coszen_sq * small_rough_cor
!Large scale roughness correction
      qdepol                 = 0.35_JPRB - 0.35_JPRB * Exp( - 0.60_JPRB * freq_ghz * large_rough * large_rough)
      qdepol_tl              =  - 0.35_JPRB * ( - 0.60_JPRB * freq_ghz * 2 * large_rough_tl * large_rough) *      &
        & Exp( - 0.60_JPRB * freq_ghz * large_rough * large_rough)
      emissfactor_v          = 1.0_JPRB - fresnel_v * small_rough_cor
      emissfactor_h          = 1.0_JPRB - fresnel_h * small_rough_cor
      emissfactor            = emissfactor_h - emissfactor_v
      emissstokes(i, 1)      = emissfactor_v + qdepol * emissfactor
      emissstokes(i, 2)      = emissfactor_h - qdepol * emissfactor
      emissstokes(i, 3)      = 0.0_JPRB
      emissstokes(i, 4)      = 0.0_JPRB
!reflect_v(i)  = 1.0_JPRB - emiss_v(i)
!reflect_h(i)  = 1.0_JPRB - emiss_h(i)
      emissfactor_v_tl       =  - fresnel_v_tl * small_rough_cor - fresnel_v * small_rough_cor_tl
      emissfactor_h_tl       =  - fresnel_h_tl * small_rough_cor - fresnel_h * small_rough_cor_tl
      emissfactor_tl         = emissfactor_h_tl - emissfactor_v_tl
      emissstokes_tl(i, 1)   = emissfactor_v_tl + qdepol_tl * emissfactor + qdepol * emissfactor_tl
      emissstokes_tl(i, 2)   = emissfactor_h_tl - qdepol_tl * emissfactor - qdepol * emissfactor_tl
      emissstokes_tl(i, 3)   = 0.0_JPRB
      emissstokes_tl(i, 4)   = 0.0_JPRB
      reflectstokes_tl(i, :) =  - emissstokes_tl(i, :)
    ENDIF
! Now calc channel emissivity after mixing v and h pol
!
    emissfactor_v_tl   = pol_v(1, pol_id) * coef%pol_fac_v(chan) + &
                         pol_v(2, pol_id) * geom%sinview_sq + pol_v(3, pol_id) * geom%cosview_sq
    emissfactor_h_tl   = pol_h(1, pol_id) * coef%pol_fac_h(chan) + &
                         pol_h(2, pol_id) * geom%sinview_sq + pol_h(3, pol_id) * geom%cosview_sq
    emissivity_tl(i)   = emissstokes_tl(i, 1) * emissfactor_v_tl + emissstokes_tl(i, 2) * emissfactor_h_tl + &
                         emissstokes_tl(i, 3) * pol_s3(0, pol_id) + emissstokes_tl(i, 4) * pol_s3(1, pol_id)
    reflectivity_tl(i) = reflectstokes_tl(i, 1) * emissfactor_v_tl + reflectstokes_tl(i, 2) * emissfactor_h_tl + &
                         reflectstokes_tl(i, 3) * pol_s3(0, pol_id) + reflectstokes_tl(i, 4) * pol_s3(1, pol_id)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_mw_tl
