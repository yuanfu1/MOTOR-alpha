! Description:
!> @file
!!   Computes MW surface emissivities.
!
!> @brief
!!   Computes MW surface emissivities.
!!
!! @details
!!   This subroutine provides emissivity values for each
!!   channel/profile where calcemis is set to TRUE.
!!
!!   For land and sea-ice surfaces the FASTEM parameterisation is
!!   used with the coefficients in profiles(:)%skin%fastem(1:5)
!!   (see the user guide).
!!
!!   For sea surfaces FASTEM versions 1-6 are available, selected
!!   in the options structure.
!!
!!   The TESSEM2 model is also available and is recommended for
!!   channels above 200GHz (in particular MetopSG ICI) although it
!!   will provide emissivities at all frequencies.
!!
!!   See the user guide for more information on the models and
!!   references.
!!
!! @param[in]     opts                options to configure the simulations
!! @param[in]     profiles            input atmospheric profiles and surface variables
!! @param[in]     geometry            internal geometry structure
!! @param[in]     coef                optical depth coefficients structure
!! @param[in]     chanprof            specifies channels and profiles to simulate
!! @param[in]     transmission_aux    RTTOV internal auxiliary transmission structure
!! @param[in]     calcemis            flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity          updated with surface emissivities on exit
!! @param[in,out] reflectivity        updated with surface reflectances for downwelling radiation on exit
!! @param[out]    err                 status on exit
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
SUBROUTINE rttov_calcemis_mw( &
            & opts,             &
            & profiles,         &
            & geometry,         &
            & coef,             &
            & chanprof,         &
            & transmission_aux, &
            & calcemis,         &
            & emissivity,       &
            & reflectivity,     &
            & err)
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_options,          &
       & rttov_chanprof,         &
       & rttov_coef,             &
       & rttov_profile,          &
       & rttov_transmission_aux, &
       & rttov_geometry
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : &
       & pi,              &
       & surftype_sea,    &
       & pol_v,           &
       & pol_h,           &
       & pol_s3,          &
       & max_fastem_version
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE mod_rttov_fastem3_coef, ONLY : fastem3_coef
  USE rttov_tessem_mod, ONLY : rttov_tessem
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options),          INTENT(IN)            :: opts
  TYPE(rttov_chanprof),         INTENT(IN)            :: chanprof(:)
  TYPE(rttov_profile),          INTENT(IN)   , TARGET :: profiles(:)
  TYPE(rttov_geometry),         INTENT(IN)   , TARGET :: geometry(size(profiles))
  TYPE(rttov_coef),             INTENT(IN)            :: coef
  TYPE(rttov_transmission_aux), INTENT(IN)            :: transmission_aux
  LOGICAL(KIND=jplm),           INTENT(IN)            :: calcemis    (size(chanprof))
  REAL   (KIND=jprb),           INTENT(INOUT)         :: emissivity  (size(chanprof))
  REAL   (KIND=jprb),           INTENT(INOUT)         :: reflectivity(size(chanprof))
  INTEGER(KIND=jpim),           INTENT(OUT)           :: err
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_fastem5.interface"
!local constants:
  REAL(KIND=jprb), PARAMETER :: windscale = 999999.0_JPRB
  REAL(KIND=jprb), PARAMETER :: windlimit = 0.0001_JPRB
  REAL(KIND=jprb), PARAMETER :: quadcof  (4, 2  ) =      &
    & Reshape((/0.0_JPRB, 1.0_JPRB, 1.0_JPRB, 2.0_JPRB, 1.0_JPRB,  - 1.0_JPRB, 1.0_JPRB,  - 1.0_JPRB/), (/4, 2/))
  REAL(KIND=jprb), PARAMETER :: freqfixed(4)      = Reshape((/7.0_JPRB, 10.0_JPRB, 19.0_JPRB, 37.0_JPRB/), (/4/))

!local variables:
  REAL(KIND=jprb)     :: tcelsius
  REAL(KIND=jprb)     :: tcelsius_sq
  REAL(KIND=jprb)     :: tcelsius_cu
  REAL(KIND=jprb)     :: einf                                        ! Debye parameter Epsilon infinity
  REAL(KIND=jprb)     :: fen, fen_sq                                 ! intermediate Debye variable
  REAL(KIND=jprb)     :: del1           , del2                       ! intermediate Debye variable
  REAL(KIND=jprb)     :: den1           , den2                       ! intermediate Debye variable
  REAL(KIND=jprb)     :: f1, f2                                      ! intermediate Debye variable
  REAL(KIND=jprb)     :: perm_free                                   ! permittivity (space)
  REAL(KIND=jprb)     :: sigma                                       ! saline water conductivity
  REAL(KIND=jprb)     :: perm_real1     , perm_real2                 ! permittivity (real part)
  REAL(KIND=jprb)     :: perm_imag1     , perm_imag2    , perm_imag3 !    .... imaginary part
  REAL(KIND=jprb)     :: perm_Real      , perm_imag                  ! permittivity (real, imaginary part)
  REAL(KIND=jprb)     :: perm_static                                 ! static land permittivity
  REAL(KIND=jprb)     :: perm_infinite                               ! infinite frequency land permittivity
  REAL(KIND=jprb)     :: freq_ghz       , freq_ghz_sq                ! frequency in GHz , and squared
  REAL(KIND=jprb)     :: fresnel_v_Real , fresnel_v_imag
  REAL(KIND=jprb)     :: fresnel_h_Real , fresnel_h_imag
  REAL(KIND=jprb)     :: fresnel_v      , fresnel_h
  REAL(KIND=jprb)     :: small_rough_cor, foam_cor
  REAL(KIND=jprb)     :: large_rough_cor(2)
  REAL(KIND=jprb)     :: small_rough, large_rough                    ! small and large scale roughness
  REAL(KIND=jprb)     :: emissstokes  (size(chanprof), 4)
  REAL(KIND=jprb)     :: reflectstokes(size(chanprof), 4)
  REAL(KIND=jprb)     :: variance        , varm
  REAL(KIND=jprb)     :: wind10
  REAL(KIND=jprb)     :: wind10_sq       , windsec       , windratio
  REAL(KIND=jprb)     :: wind10_direction, windangle                 ! Note wind azimuth is in radians
  REAL(KIND=jprb)     :: opdpsfc         , freqr
  REAL(KIND=jprb)     :: zrough_v        , zrough_h
  REAL(KIND=jprb)     :: zreflmod_v      , zreflmod_h
  REAL(KIND=jprb)     :: delta           , delta2
  REAL(KIND=jprb)     :: qdepol          , emissfactor
  REAL(KIND=jprb)     :: emissfactor_v   , emissfactor_h
  REAL(KIND=jprb)     :: zc(12)                                      ! large scale correction
  REAL(KIND=jprb)     :: zx(9 )                                      ! effective path coefficients
  REAL(KIND=jprb)     :: azimuthal_emiss, u19, phi       , dfreq
  REAL(KIND=jprb)     :: tbfixed      (4, 4, 3)  ! Surface brightness temperature azimuthal variation terms for 37, 19, 10, 7 GHz
  REAL(KIND=jprb)     :: efixed       (4, 4, 3)  ! Emissivity azimuthal variation terms for 7, 10, 19, 37 GHz
  REAL(KIND=jprb)     :: einterpolated(4, 3   )  ! Emissivity azimuthal variation terms for interpolated to required frequency
  REAL(KIND=jprb)     :: a1e, a2e, a3e           ! coefficients used in azimuthal emissivity model
  REAL(KIND=jprb), POINTER :: c(:)! pointer to FASTEM coefs
  COMPLEX(KIND=jprb) :: perm1       , perm2          ! permittivity
  COMPLEX(KIND=jprb) :: rhth        , rvth           ! Fresnel reflectivity complex variables
  COMPLEX(KIND=jprb) :: permittivity                 ! permittivity
  INTEGER(KIND=jpim) :: i, j, chan, istokes, ifreq, m
  INTEGER(KIND=jpim) :: iquadrant                    ! Determines which quadrant (NE, SE, SW, NW) the wind is blowing to
  INTEGER(KIND=jpim) :: pol_id                       ! polarisation indice
  INTEGER(KIND=jpim) :: i_freq      , j_stokes, ich  ! indices used in azimuthal emissivity model
! == pol_id +1
!   1 average of vertical and horizontal
!   2 nominal vertical at nadir, rotating
!      with view angle
!   3 nominal horizontal at nadir, rotating
!      with view angle
!   4 vertical
!   5 horizontal
!   6 + 45 minus -45 (3rd stokes vector)
!   7 left circular - right circular (4th stokes vector)
  INTEGER(KIND=jpim) :: jcof        , jcofm1
  TYPE(rttov_profile),  POINTER :: prof
  TYPE(rttov_geometry), POINTER :: geom
! rttov_fastem4 internal variables
  REAL   (KIND=jprb) :: Zenith_Angle     , Wind_Speed, Salinity, Transmittance, Rel_Azimuth
  REAL(KIND=jprb), DIMENSION(4) :: jemissivity, jreflectivity
  INTEGER(KIND=jpim) :: nchannels   ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  TRY
! The rttov_check_options subroutine guarantees the selected FASTEM version is valid
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW', 0_jpim, ZHOOK_HANDLE)
  nchannels         = size(chanprof)
!Loop over channels
  DO i = 1, nchannels
    IF (.NOT. calcemis(i)) CYCLE
    chan = chanprof(i)%chan
    prof => profiles(chanprof(i)%prof)
    geom => geometry(chanprof(i)%prof)
    pol_id = coef%fastem_polar(chan) + 1_jpim
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
      wind10_sq = prof%s2m%u * prof%s2m%u + prof%s2m%v * prof%s2m%v
      wind10    = Sqrt(wind10_sq)
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
      freq_ghz         = coef%frequency_ghz(chan)
      freq_ghz_sq      = freq_ghz * freq_ghz
      SELECT CASE (opts%rt_mw%fastem_version)
      CASE (0)
        CALL rttov_tessem(freq_ghz, prof%zenangle, wind10, prof%skin%t, prof%skin%salinity, &
                          emissstokes(i,2), emissstokes(i,1))
        emissstokes(i,3:4) = 0._jprb
        reflectstokes(i,:) = 1._jprb - emissstokes(i,:)
      CASE (1:3)
        windsec        = wind10 * geom%seczen
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
! perm_free = 8.854E-3_JPRB not 8.854E-12 as multiplied by 1E9 for GHz
        perm_free      = 8.854E-3_JPRB
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
        fresnel_v_Real = Dble(rvth)
        fresnel_v_imag = Aimag(rvth)
        fresnel_v      = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
        fresnel_h_Real = Dble(rhth)
        fresnel_h_imag = Aimag(rhth)
        fresnel_h      = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
!1.3.2) Small scale correction to reflection coefficients
!------
        IF (freq_ghz >= 15.0_jprb) THEN
          small_rough_cor = Exp(c(21) * wind10 * geom%coszen_sq / (freq_ghz_sq))
        ELSE
          small_rough_cor = 1.0_jprb
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
        large_rough_cor(1) =      &
          & zc(1) + zc(2) * geom%seczen + zc(3) * geom%seczen_sq + zc(4) * wind10 + zc(5) * wind10_sq + zc(6) * windsec
        large_rough_cor(2) =      &
          & zc(7) + zc(8) * geom%seczen + zc(9) * geom%seczen_sq + zc(10) * wind10 + zc(11) * wind10_sq + zc(12) * windsec
        large_rough_cor(:) = large_rough_cor(:) * 0.01_JPRB
! For Fastem-3 do not compute rough surface effects if theta > 60 degrees
        IF (opts%rt_mw%fastem_version <= 2.0_JPRB .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) THEN
          emissstokes(i, 1) = 1.0_JPRB - fresnel_v * small_rough_cor + large_rough_cor(1)
          emissstokes(i, 2) = 1.0_JPRB - fresnel_h * small_rough_cor + large_rough_cor(2)
        ELSE
          emissstokes(i, 1) = 1.0_JPRB - fresnel_v
          emissstokes(i, 2) = 1.0_JPRB - fresnel_h
        ENDIF
        emissstokes(i, 3) = 0.0_JPRB
        emissstokes(i, 4) = 0.0_JPRB
!Apply foam correction
        IF (.NOT. opts%rt_mw%supply_foam_fraction) THEN
          foam_cor          = c(22) * (wind10 ** c(23))
        ELSE
          foam_cor          = prof%skin%foam_fraction
        ENDIF
        emissstokes(i, 1) = emissstokes(i, 1) - foam_cor * emissstokes(i, 1) + foam_cor
        emissstokes(i, 2) = emissstokes(i, 2) - foam_cor * emissstokes(i, 2) + foam_cor
        IF (opts%rt_mw%fastem_version == 3) THEN
! Add azimuthal component from Fuzhong Weng (NOAA/NESDIS) based on work by Dr. Gene Poe (NRL)
! Angle between wind direction and satellite azimuthal view angle
          phi = pi - (wind10_direction - prof%azangle * pi / 180.0_JPRB)
! Assume 19m wind = 10m wind for now (fix later).
          u19 = wind10
          DO ich = 0, 15
            a1e = c(141 + ich * 12) + u19 * (c(142 + ich * 12) + u19 * (c(143 + ich * 12) + u19 * c(144 + ich * 12)))
            a2e = c(145 + ich * 12) + u19 * (c(146 + ich * 12) + u19 * (c(147 + ich * 12) + u19 * c(148 + ich * 12)))
            a3e = c(149 + ich * 12) + u19 * (c(150 + ich * 12) + u19 * (c(151 + ich * 12) + u19 * c(152 + ich * 12)))
            i_freq = int(ich / 4_jpim) + 1                      ! 37, 19, 10, 7 GHz
            j_stokes                     = mod(ich, 4_jpim) + 1
            tbfixed(j_stokes, i_freq, 1) = a1e                  !* prof%skin%t
            tbfixed(j_stokes, i_freq, 2) = a2e                  !* prof%skin%t
            tbfixed(j_stokes, i_freq, 3) = a3e                  !* prof%skin%t
          ENDDO
          DO M = 1, 3
            DO istokes = 1, 4
              efixed(1, istokes, M) = tbfixed(istokes, 4, M)!/prof%skin%t  ! 7  GHz
              efixed(2, istokes, M) = tbfixed(istokes, 3, M)!/prof%skin%t  ! 10  GHz
              efixed(3, istokes, M) = tbfixed(istokes, 2, M)!/prof%skin%t  ! 19  GHz
              efixed(4, istokes, M) = tbfixed(istokes, 1, M)!/prof%skin%t  ! 37  GHz
            ENDDO
! Interpolate results to required frequency based on 7, 10, 19, 37 GHz
            IF (freq_ghz .LE. freqfixed(1)) THEN
              DO istokes = 1, 4
                einterpolated(istokes, M) = efixed(1, istokes, M)
              ENDDO
            ELSE IF (freq_ghz .GE. freqfixed(4)) THEN
              DO istokes = 1, 4
                einterpolated(istokes, M) = efixed(4, istokes, M)
              ENDDO
            ELSE
              IF (freq_ghz .LT. freqfixed(2)                                 ) ifreq = 2
              IF (freq_ghz .LT. freqfixed(3) .AND. freq_ghz .GE. freqfixed(2)) ifreq = 3
              IF (freq_ghz .GE. freqfixed(3)                                 ) ifreq = 4
              dfreq = (freq_ghz - freqfixed(ifreq - 1)) / (freqfixed(ifreq) - freqfixed(ifreq - 1))
              DO istokes = 1, 4
                einterpolated(istokes, M) =      &
                  & efixed(ifreq - 1, istokes, M) + dfreq * (efixed(ifreq, istokes, M) - efixed(ifreq - 1, istokes, M))
              ENDDO
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
! Only apply non-specular correction for Fastem-3 if theta < 60 degrees
        IF ((opts%rt_mw%fastem_version == 2 .OR. (opts%rt_mw%fastem_version == 3 .AND. geom%seczen <= 2.0_JPRB)) .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) < 0.9999_JPRB .AND. &
          & transmission_aux%thermal_path1%tau_surf(0, i) > 0.00001_JPRB) THEN
!Convert windspeed to slope variance using the Cox and Munk model
          variance = 0.00512_JPRB * wind10 + 0.0030_JPRB
          varm     = variance * c(138)
          variance = varm * (c(139) * freq_ghz + c(140))
          IF (variance > varm    ) variance = varm
          IF (variance < 0.0_JPRB) variance = 0.0_JPRB
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
          zreflmod_v          =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_v) / &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          zreflmod_h          =      &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i) ** zrough_h) / &
            & (1.0_JPRB - transmission_aux%thermal_path1%tau_surf(0, i))
          reflectstokes(i, 1) = zreflmod_v * (1.0_JPRB - emissstokes(i, 1))
          reflectstokes(i, 2) = zreflmod_h * (1.0_JPRB - emissstokes(i, 2))
          reflectstokes(i, 3) =  - 0.5_JPRB * (zreflmod_v + zreflmod_h) * emissstokes(i, 3)
          reflectstokes(i, 4) =  - 0.5_JPRB * (zreflmod_v + zreflmod_h) * emissstokes(i, 4)
        ELSE
          reflectstokes(i, 1) = 1.0_JPRB - emissstokes(i, 1)
          reflectstokes(i, 2) = 1.0_JPRB - emissstokes(i, 2)
          reflectstokes(i, 3) = 0.0_JPRB
          reflectstokes(i, 4) = 0.0_JPRB
        ENDIF
      CASE (4:max_fastem_version)
        Zenith_Angle  = 1.d0 / geom%seczen
        Zenith_Angle  = acos(Zenith_Angle) * 180.0_JPRB / pi
        Salinity      = prof%skin%salinity
        Transmittance = transmission_aux%thermal_path1%tau_surf(0, i)
! relative azimuth angle in degree
        Rel_Azimuth   = (wind10_direction * 180.0_JPRB / pi - prof%azangle)
        Wind_Speed    = wind10
        CALL rttov_fastem5( &
              & opts%rt_mw%fastem_version, &
              & freq_ghz,      &
              & Zenith_Angle,  &
              & prof%skin%t,   &
              & Salinity,      &
              & Wind_Speed,    &
              & jemissivity,   &
              & jreflectivity, &
              & Transmittance, &
              & Rel_Azimuth,   & ! Input, may not be used
              & Supply_Foam_Fraction = opts%rt_mw%supply_foam_fraction, &
              & Foam_Fraction = prof%skin%foam_fraction)
        emissstokes(i, 1:4)   = jemissivity(1:4)
        reflectstokes(i, 1:4) = jreflectivity(1:4)
      END SELECT
!--------------------
!2. Land/ice surfaces
!--------------------
    ELSE
! Test input FASTEM land coefficients
! only coefs 1-3 are checked
      IF (Any(prof%skin%fastem(1:3) == 0.0_JPRB)) err = errorstatus_fatal
      THROWM( ERR .NE. 0 , "some profile fastem(1:3) values are 0.0 ")
!Coherent surface scattering model coefficients (input with the profile)
      perm_static         = prof%skin%fastem(1)
      perm_infinite       = prof%skin%fastem(2)
      freqr               = prof%skin%fastem(3)
      small_rough         = prof%skin%fastem(4)
      large_rough         = prof%skin%fastem(5)
      freq_ghz            = coef%frequency_ghz(chan)
!Simple Debye + Fresnel model gives reflectivities
      fen = freq_ghz / freqr
      fen_sq              = fen * fen
      den1 = 1.0_JPRB + fen_sq
      perm_Real           = (perm_static + perm_infinite * fen_sq) / den1
      perm_imag           = fen * (perm_static - perm_infinite) / den1
      permittivity        = Cmplx(perm_Real, perm_imag, jprb)
      perm1               = sqrt(permittivity - geom%sinzen_sq)
      perm2               = permittivity * geom%coszen
      rhth = (geom%coszen - perm1) / (geom%coszen + perm1)
      rvth = (perm2 - perm1) / (perm2 + perm1)
      fresnel_v_Real      = Dble(rvth)
      fresnel_v_imag      = Aimag(rvth)
      fresnel_v           = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
      fresnel_h_Real      = Dble(rhth)
      fresnel_h_imag      = Aimag(rhth)
      fresnel_h           = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag
!Small scale roughness correction
      delta               = 4.0_JPRB * pi * coef%ff_cwn(chan) * 0.1_JPRB * small_rough
      delta2              = delta * delta
      small_rough_cor     = Exp( - delta2 * geom%coszen_sq)
!Large scale roughness correction
      qdepol              = 0.35_JPRB - 0.35_JPRB * Exp( - 0.60_JPRB * freq_ghz * large_rough * large_rough)
      emissfactor_v       = 1.0_JPRB - fresnel_v * small_rough_cor
      emissfactor_h       = 1.0_JPRB - fresnel_h * small_rough_cor
      emissfactor         = emissfactor_h - emissfactor_v
      emissstokes(i, 1)   = emissfactor_v + qdepol * emissfactor
      emissstokes(i, 2)   = emissfactor_h - qdepol * emissfactor
      emissstokes(i, 3)   = 0.0_JPRB
      emissstokes(i, 4)   = 0.0_JPRB
      reflectstokes(i, :) = 1.0_JPRB - emissstokes(i, :)
! End of if sea else land if else endif loop
    ENDIF
!
! Now calc channel emissivity after mixing v and h pol
!
    emissfactor_v   = pol_v(1, pol_id) * coef%pol_fac_v(chan) + &
                      pol_v(2, pol_id) * geom%sinview_sq + pol_v(3, pol_id) * geom%cosview_sq
    emissfactor_h   = pol_h(1, pol_id) * coef%pol_fac_h(chan) + &
                      pol_h(2, pol_id) * geom%sinview_sq + pol_h(3, pol_id) * geom%cosview_sq
    emissivity(i)   = emissstokes(i, 1) * emissfactor_v + emissstokes(i, 2) * emissfactor_h + &
                      emissstokes(i, 3) * pol_s3(0, pol_id) + emissstokes(i, 4) * pol_s3(1, pol_id)
    reflectivity(i) = reflectstokes(i, 1) * emissfactor_v + reflectstokes(i, 2) * emissfactor_h + &
                      reflectstokes(i, 3) * pol_s3(0, pol_id) + reflectstokes(i, 4) * pol_s3(1, pol_id)
! End loop over channels
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW', 1_jpim, ZHOOK_HANDLE)
  CATCH
  IF (LHOOK) CALL DR_HOOK('RTTOV_CALCEMIS_MW', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_calcemis_mw
