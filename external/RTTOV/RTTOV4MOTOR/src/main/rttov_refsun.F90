! Description:
!> @file
!!   Compute the fraction of solar radiance that is reflected by a wind
!!   roughened water surface.
!
!> @brief
!!   Compute the fraction of solar radiance that is reflected by a wind
!!   roughened water surface.
!!
!! @details
!!   The water surface is regarded as a collection of small mirror-like facets
!!   each randomly tilted with respect to the local horizon. A time passes, the
!!   tilt of a facet varies under the influence of the wind. When the open
!!   ocean reflects the solar disk, these fluctuating facets produce a dancing
!!   pattern known as sun glint.
!!   The fraction of solar radiance reflected by a wind roughened water surface
!!   is computed assuming that the slope of the facets obey to a
!!   two-dimensional Gaussian random process whose spectrum is specified by
!!   a parameterisation (see below).
!!   The total variance of the slope of the facet is obtained from the
!!   frequency spectrum of the surface wave and the inverse function of the
!!   dispersive relation of the full-gravity-capillary wave.
!!   The shadowing of the surface of the facets on the backsides of the waves
!!   and deep in the throughs between waves is considered.
!!   A coordinate system is used where the average water surface lies in the
!!   X-Y plane. The coordinate system is right-handed and is located at the
!!   reflection point. The Z axis points toward the zenith and the X axis
!!   points in the direction formed by the projection of the reflected ray
!!   on the average water surface.
!!
!!   Yoshimori et al used the Joint North Sea Wave Project (JONSWAP) wave
!!   spectral model for the wave spectrum. An alternative parameterisation
!!   by Elfouhaily et al has been implemented which yields smaller biases
!!   compared to observations in glint-affected conditions.
!!
!!   K. Yoshimori, K. Itoh and Y. Ichioka: 'Optical charateristics of a
!!   wind-roughened water surface:a two dimensional theory'. Applied Optics,
!!   Vol.34, No.27, 20 September 1995, pp.6236-6247.
!!
!!   Hasselmann, K. et al. (1973) Measurements of wind-wave growth and swell
!!   during the Joint North Sea Wave Project (JONSWAP)
!!   Dtsch.Hydrogr.Z., 12, 95 pp.
!!
!!   Elfouhaily, T., B. Chapron, K. Katsaros, and D. Vandemark (1997),
!!   A unified directional spectrum for long and short wind-driven waves,
!!   J. Geophys. Res., 102(C7), 15781-15796, doi:10.1029/97JC00467.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     aux            internal structure containing auxiliary profile variables
!! @param[in,out] sunglint       internal structure for sea surface BRDF model variables
!! @param[in]     raytracing     raytracing structure
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
SUBROUTINE rttov_refsun( &
              opts,       &
              profiles,   &
              coef,       &
              aux,        &
              sunglint,   &
              raytracing)

  USE rttov_types, ONLY : &
      rttov_options,     &
      rttov_profile_aux, &
      rttov_profile,     &
      rttov_raytracing,  &
      rttov_sunglint,    &
      rttov_coef
!INTF_OFF
  USE rttov_const, ONLY : &
      realtol,          &
      deg2rad, pi,      &
      gravity,          &
      surftype_sea,     &
      min_windsp,       &
      max_sol_zen,      &
      max_exp_exponent, &
      min_exponent,     &
      nk_elf, cm_elf,   &
      km_elf, x0_elf,   &
      a0_elf, ap_elf
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_math_mod, ONLY : exponential
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_coef),        INTENT(IN)    :: coef
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_sunglint),    INTENT(INOUT) :: sunglint
  TYPE(rttov_raytracing),  INTENT(IN)    :: raytracing
!INTF_END

  REAL(jprb)    :: ff(coef%ws_nomega)        ! Working space

  INTEGER(jpim) :: j, i, ilaysur, nprofiles
  REAL(jprb)    :: theta_fi                  ! First order shadowing factor for the reflected ray
  REAL(jprb)    :: theta_csi                 ! First order shadowing factor for the incident ray surface
  REAL(jprb)    :: omega_1                   ! Frequency of the surface wave
  REAL(jprb)    :: k
  REAL(jprb)    :: sigma_a2_r, sigma_b2_r, sigma2_r
  REAL(jprb)    :: bb, aa, hinc, s2, s4, sa, sb, sqrt2pi

  REAL(jprb)    :: omega_m2_r
  REAL(jprb)    :: omega_1_array(coef%ws_nomega)
  REAL(jprb)    :: omega_1_4_r(coef%ws_nomega), omega_1_5_r(coef%ws_nomega)
  REAL(jprb)    :: psi_exp(coef%ws_nomega)
  REAL(jprb)    :: beta_exponent, log_min_exponent, max_exp_exponent_4rt

  REAL(jprb)    :: karr(0:nk_elf), sk2o(0:nk_elf), sk2p(0:nk_elf)
  REAL(jprb)    :: sqrt_k_kp
  REAL(jprb)    :: cos2wangl
  REAL(jprb)    :: windsp, k0, x, kp, gamma, sigma
  REAL(jprb)    :: cp, omega, alpha_p, ustar, alpha_m, am
  REAL(jprb)    :: gamma_exponent, fp, bl, bh

  REAL(jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  sqrt2pi = SQRT(2_jpim * pi)
  k = 3.3_jprb
  sigma_a2_r = 1._jprb / 0.07_jprb**2
  sigma_b2_r = 1._jprb / 0.09_jprb**2

  max_exp_exponent_4rt = SQRT(SQRT(max_exp_exponent))
  log_min_exponent = LOG(min_exponent)

  IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN
    DO i = 1, coef%ws_nomega
      omega_1_array(i) = (i - 1) * 0.1_jprb
    ENDDO
    omega_1_4_r(2:coef%ws_nomega) = 1._jprb / (omega_1_array(2:coef%ws_nomega)**4_jpim)
    omega_1_4_r(1) = 0._jprb
    omega_1_5_r(2:coef%ws_nomega) = omega_1_4_r(2:coef%ws_nomega) / omega_1_array(2:coef%ws_nomega)
    omega_1_5_r(1) = 0._jprb
  ELSE
    ! nk_elf = 318 : 0-1 @ 0.1 [10], 1-20 @ 0.5 [38], 20-100 @ 1 [80], 100-2000 @ 10 [190]
    DO i = 0, 10
      karr(i) = 0.1_jprb * i
    ENDDO
    DO i = 1, 38
      karr(i+10) = 1 + 0.5_jprb * i
    ENDDO
    DO i = 1, 80
      karr(i+48) = 20 + 1._jprb * i
    ENDDO
    DO i = 1, 190
      karr(i+128) = 100 + 10._jprb * i
    ENDDO
  ENDIF

  sunglint%s(:)%glint = 0._jprb
  sunglint%s(:)%omega = 0._jprb
  DO j = 1, nprofiles
    IF (profiles(j)%sunzenangle < 0.0 .OR. &
        profiles(j)%sunzenangle > max_sol_zen) CYCLE
    IF (profiles(j)%skin%surftype /= surftype_sea) CYCLE

    ! --------------------------------------------------------------------
    ! Calculate profile variables
    ! --------------------------------------------------------------------

    ! Compute the wind speed

    sunglint%s(j)%windsp = SQRT(profiles(j)%s2m%u**2_jpim + profiles(j)%s2m%v**2_jpim)

    ! Compute the angle between the wind direction and the U axis

    IF (sunglint%s(j)%windsp > min_windsp) THEN
      sunglint%s(j)%wangl = MOD(ATAN2(profiles(j)%s2m%v, profiles(j)%s2m%u), 2._jprb * pi)
    ELSE
      sunglint%s(j)%wangl = 0._jprb
    ENDIF

    ! Compute the difference between the sun azimuth angle and the
    ! azimuth angle of the direction formed by the projection on the
    ! mean water surface of the surface-to-sun direction
    ! Note defn of azangle in user guide.

    sunglint%s(j)%dazng = MOD(profiles(j)%azangle - profiles(j)%sunazangle, 360._jprb)
    IF (sunglint%s(j)%dazng < 0.) sunglint%s(j)%dazng = sunglint%s(j)%dazng + 360._jprb

    ! Satellite and solar zenith angles at surface
    ilaysur = aux%s(j)%nearestlev_surf - 1
    sunglint%s(j)%zensat = ACOS(1._jprb / raytracing%pathsat(ilaysur, j))
    sunglint%s(j)%zensun = ACOS(1._jprb / raytracing%pathsun(ilaysur, j))

    IF (sunglint%s(j)%windsp <= min_windsp) THEN

      ! In absence of wind the water surface is a mirror-like surface
      IF (ABS(sunglint%s(j)%zensat - sunglint%s(j)%zensun) < realtol .AND. &
          ABS(sunglint%s(j)%dazng - 180.0_jprb) < realtol) THEN
        sunglint%s(j)%glint = 1._jprb
      ELSE
        sunglint%s(j)%glint = 0._jprb
      ENDIF

    ELSE

      ! ------------------------------------------------------------------
      ! Wave spectrum parameterisation
      ! ------------------------------------------------------------------

      ! Compute the angle between the wind direction and the view direction
      ! Note defn of azangle in user guide
      sunglint%s(j)%wangl = &
          MOD(sunglint%s(j)%wangl - MOD(90.0_jprb - profiles(j)%azangle, 360.0_jprb) * deg2rad, 2 * pi)
      cos2wangl = COS(2 * sunglint%s(j)%wangl)

      IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN

        ! JONSWAP wave frequency spectrum parameterisation

        sunglint%s(j)%x_u     = profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp**2_jpim
        sunglint%s(j)%alfa1   = 0.076_jprb / sunglint%s(j)%x_u**0.22_jprb
        sunglint%s(j)%omega_m =      &
            2._jprb * pi * (3.5_jprb * (gravity / sunglint%s(j)%windsp) * (1._jprb / sunglint%s(j)%x_u**0.33_jprb))
        omega_m2_r = 1._jprb / sunglint%s(j)%omega_m**2_jpim

        !-----------Compute the total variance of the slope of the facet---------------
        !                                                                              |
        !           To compute the total variance of the slope of the facet (or the    |
        !           total variance of the displacement of the water surface) the       |
        !           knowledge of the frequency spectrum of the surface wave and the    |
        !           inverse function of the dispersive relation of the full-gravity-   |
        !           capillary wave is needed.                                          |
        !                                                                              |
        !           The dispersive relation of the full-gravity-capillary wave is      |
        !           written as:                                                        |
        !                                                                              |
        !           omega_1(k)=k*sqrt((gravity/k+gamma*k/ro)*tanh(h*k))                |
        !                                                                              |
        !           k    = wave-number                                                 |
        !           gamma= surface-tension constant                                    |
        !           ro   = density of water                                            |
        !           h    = water depth                                                 |
        !                                                                              |
        !           The inverse function of the dispersive relation of the             |
        !           full-gravity-capillary wave has been pre-computed numerically      |
        !           using the following parameters:                                    |
        !                                                                              |
        !           gamma = 0.0735   [N/m]                                             |
        !           ro    = 999.1    [kg/m^3]                                          |
        !           h     = 50       [m]                                               |
        !                                                                              |
        !           The frequency spectrum (PSI) of the surface wave is computed       |
        !           using the JONSWAP wave-spectral model.                             |
        !                                                                              |
        !           The total variance of the slope is then evaluated by integrating   |
        !           the quantity (psi*k_omega**2):                                     |
        !                                                                              |
        !                     -infinity                                                |
        !                    -                                                         |
        !           gamma_sq=- (k_omega(omega_1)**2*psi(omega_1))  d omega_1           |
        !                    -                                                         |
        !                   -zero                                                      |
        !                                                                              |
        !           The integral is computed applying Simpson's rule using 302         |
        !           quadrature points.                                                 |
        !                                                                              |
        !------------------------------------------------------------------------------

        ! Compute JONSWAP frequency spectrum of the surface wave

        ! Move expensive calculations outside omega loop where possible;
        ! enable vectorised exponential
        psi_exp(:) = -(5._jprb / 4._jprb) * & ! stop exp overflow
          MIN(sunglint%s(j)%omega_m**4_jpim * omega_1_4_r(:), max_exp_exponent)
        CALL exponential(psi_exp, sunglint%psi(:,j))
        sunglint%psi(:,j) = sunglint%psi(:,j) * sunglint%s(j)%alfa1 * gravity**2_jpim * omega_1_5_r(:)

        ff(1) = 0._jprb
        DO i = 2, coef%ws_nomega
          omega_1 = omega_1_array(i)

          IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

            IF (omega_1 <= sunglint%s(j)%omega_m) THEN
              sigma2_r = sigma_a2_r
            ELSE
              sigma2_r = sigma_b2_r
            ENDIF

            beta_exponent = MIN((omega_1 - sunglint%s(j)%omega_m)**2_jpim * &
                                0.5_jprb * sigma2_r * omega_m2_r, max_exp_exponent)

            ! stop beta underflow
            ! IF ([sunglint%beta(i, j) == EXP(-beta_exponent)] > min_exponent)
            IF (-beta_exponent > log_min_exponent) THEN
              sunglint%beta(i, j) = EXP(-beta_exponent) ! Moved inside this IF clause to avoid unnecessary exponentials
              sunglint%psi(i, j) = sunglint%psi(i, j) * (k**sunglint%beta(i, j))
            ELSE
              sunglint%beta(i, j) = 0._jprb
            ENDIF

          ELSE
            sunglint%psi(i, j) = 0._jprb
          ENDIF
          ff(i) = sunglint%psi(i, j) * coef%ws_k_omega(i)**2_jpim
        ENDDO

        ! Compute the Simpson's integral

        bb   = 301._jprb
        aa   = 0._jprb
        hinc = (bb - aa) / coef%ws_nomega
        sa   = ff(1)
        sb   = ff(coef%ws_nomega)
        s4 = SUM(ff(1::2))
        s2 = SUM(ff(2::2))
        sunglint%s(j)%gamma_sq = hinc / 3.0_jprb * (sa + sb + 4.0_jprb * s4 + 2.0_jprb * s2)

        ! Compute the rms of the slope of the facet along the X (gamma_o)and Y (gamma_p) axis

        sunglint%s(j)%gamma_o = SQRT(0.25_jprb * (2 + cos2wangl) * sunglint%s(j)%gamma_sq)
        sunglint%s(j)%gamma_p = SQRT(0.25_jprb * (2 - cos2wangl) * sunglint%s(j)%gamma_sq)

      ELSEIF (opts%rt_ir%solar_sea_brdf_model == 2) THEN

        ! Elfouhaily et al wave spectrum parameterisation
        ! Numbers in parentheses below refer to equations in the paper

        ! For wind speeds below ~0.6m/s the equations from
        ! the paper yield entirely negative spectra
        windsp = MAX(sunglint%s(j)%windsp, 0.6_jprb)

        sunglint%sk2(0,j) = 0._jprb
        sk2o(0) = 0._jprb
        sk2p(0) = 0._jprb

        ! Active values independent of wavenumber (k)
        k0 = gravity / windsp**2                                      ! (3)
        x = k0 * profiles(j)%s2m%wfetc                                ! paragraph after (4)
        sunglint%s(j)%omega_c = &                                     ! (37) 0.84 < omega_c < 5
          0.84_jprb * (TANH((x / x0_elf)**0.4))**(-0.75)
        kp = k0 * sunglint%s(j)%omega_c**2                            ! (3)

        IF (sunglint%s(j)%omega_c < 1) THEN                           ! (3)
          gamma = 1.7_jprb
        ELSE
          gamma = 1.7_jprb + 6 * LOG10(sunglint%s(j)%omega_c)
        ENDIF
        sigma = 0.08_jprb * (1 + 4 / sunglint%s(j)%omega_c**3)        ! (3)

        cp = SQRT(gravity * (1 + (kp / km_elf)**2) / kp)              ! phase speed cp = angular_freq/kp
        omega = windsp / cp                                           ! paragraph after (31)
        alpha_p = 6E-3_jprb * SQRT(omega)                             ! (34)

        ! Air friction velocity: see Elfouhaily Fig 7b. This simple expression is a reasonable approximation.
        ! See Plant (1982) A Relationship Between Wind Stress and Wave Slope, between eqns (17) and (18)
        ustar = windsp / 25
        ! Donelan et al (from Elfouhaily) below differs at higher wind speed, but the simple ustar above
        ! gave better obs-sim stats for SEVIRI
        ! z0 = 3.7E-5_jprb * windsp**2 * &                            ! (66)
        !     (windsp / cp)**0.9 / gravity
        ! kappa_vonkarman = 0.41_jprb
        ! ustar = windsp * kappa_vonkarman / LOG(10. / z0)            ! from (60)

        IF (ustar < cm_elf) THEN                                      ! (44)
          alpha_m = 1E-2_jprb * (1 + LOG(ustar / cm_elf))
        ELSE
          alpha_m = 1E-2_jprb * (1 + 3 * LOG(ustar / cm_elf))
        ENDIF

        ! Angular spread
        am = 0.13_jprb * ustar / cm_elf                               ! (59)

        DO i = 1, nk_elf
          k = karr(i)
          sqrt_k_kp = SQRT(k / kp)

          ! This gives very similar values to the expression in Yoshimori
          sunglint%c(i,j) = SQRT(gravity * (1 + (k / km_elf)**2) / k) ! (24) phase speed c = angular_freq/k

          ! Long-wave spectrum
          sunglint%lpm(i,j) = EXP(-5 * (kp / k)**2 / 4)               ! (2)
          gamma_exponent = &
            MIN((sqrt_k_kp - 1)**2 / (2 * sigma**2), max_exp_exponent)
          IF (-gamma_exponent > log_min_exponent) THEN
            ! More efficient to do these calculations only when required
            sunglint%gamma_exp(i,j) = EXP(-gamma_exponent)            ! (3)
            sunglint%jp(i,j) = gamma**sunglint%gamma_exp(i,j)         ! (3)
          ELSE
            sunglint%gamma_exp(i,j) = 0._jprb
            sunglint%jp(i,j) = 1._jprb
          ENDIF

          sunglint%fpexp(i,j) = &
            EXP(-MIN((sqrt_k_kp - 1) * omega / &                      ! (32)
                      SQRT(10._jprb), max_exp_exponent))
          fp = sunglint%lpm(i,j) * sunglint%jp(i,j) * &               ! (32)
               sunglint%fpexp(i,j)
          bl = 0.5_jprb * alpha_p * cp * fp / sunglint%c(i,j)         ! (31)

          ! Short-wave spectrum
          sunglint%fm(i,j) = EXP(-(k / km_elf - 1)**2 / 4)            ! (41)
          bh = 0.5_jprb * alpha_m * cm_elf * sunglint%fm(i,j) / &     ! (40)
               sunglint%c(i,j)

          ! Combined spectrum
          sunglint%sk2(i,j) = (bl + bh) / k                           ! (30) x k**2 (integrate psi * k**2)

          IF (sunglint%sk2(i,j) > 0._jprb) THEN
            ! Angular spreading factor
            sunglint%dk(i,j) = &                                      ! (57)
              TANH(a0_elf + ap_elf * (sunglint%c(i,j) / cp)**2.5 + &
                   am * (cm_elf / sunglint%c(i,j))**2.5)

            ! The following expressions represent the wave spectrum parallel to and
            ! perpendicular to the view direction respectively. The first is derived
            ! by evaluating the following integral over phi from -pi to +pi:
            !      (1 + dk * cos(2(phi-wangl))) * cos(phi)^2 dphi
            ! wangl is the angle between the wind and view directions
            ! The perpendicular expression comes from replacing wangl with wangl + pi/2
            ! The above expression comes from (67) in Elfouhaily
            sk2o(i) = sunglint%sk2(i,j) * 0.5_jprb * (1 + 0.5_jprb * sunglint%dk(i,j) * cos2wangl)
            sk2p(i) = sunglint%sk2(i,j) * 0.5_jprb * (1 - 0.5_jprb * sunglint%dk(i,j) * cos2wangl)
          ELSE
            sk2o(i) = 0._jprb
            sk2p(i) = 0._jprb
          ENDIF
        ENDDO

        ! Trapezium integration (k grid is not uniform)
        sunglint%s(j)%gamma_o = SQRT(0.5_jprb * SUM((sk2o(0:nk_elf-1) + sk2o(1:nk_elf)) * (karr(1:nk_elf) - karr(0:nk_elf-1))))
        sunglint%s(j)%gamma_p = SQRT(0.5_jprb * SUM((sk2p(0:nk_elf-1) + sk2p(1:nk_elf)) * (karr(1:nk_elf) - karr(0:nk_elf-1))))

      ENDIF ! wave spectrum option


      ! ------------------------------------------------------------------
      ! Yoshimori wave facet reflectance model
      ! ------------------------------------------------------------------

      ! Obtain angles csi and alfa

      sunglint%s(j)%csi     = ATAN(TAN(sunglint%s(j)%zensun) * COS(pi - sunglint%s(j)%dazng * deg2rad))
      sunglint%s(j)%alfa    = ASIN(SIN(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad))
      sunglint%s(j)%omega   = ABS(sunglint%s(j)%csi + sunglint%s(j)%zensat) / 2._jprb

      ! Compute the value of the function that represents the
      ! shadowing of the incident and reflected ray

      ! First order shadowing (the slope of the facet is negative)

      sunglint%s(j)%gammax  = TAN((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
      IF (ABS(sunglint%s(j)%zensat) > 0._jprb) THEN
        IF ((1._jprb / TAN(ABS(sunglint%s(j)%zensat)) - &
            sunglint%s(j)%gammax * sunglint%s(j)%zensat / ABS(sunglint%s(j)%zensat)) >= 0._jprb) THEN
          theta_fi = 1._jprb
        ELSE
          theta_fi = 0._jprb
        ENDIF
      ELSE
        theta_fi = 1._jprb
      ENDIF
      IF (ABS(sunglint%s(j)%csi) > 0._jprb) THEN
        IF ((1._jprb / TAN(ABS(sunglint%s(j)%csi)) + &
            sunglint%s(j)%gammax * sunglint%s(j)%csi / ABS(sunglint%s(j)%csi)) >= 0._jprb) THEN
          theta_csi = 1._jprb
        ELSE
          theta_csi = 0._jprb
        ENDIF
      ELSE
        theta_csi = 1._jprb
      ENDIF
      IF (ABS(sunglint%s(j)%csi) > pi / 2._jprb) theta_csi = 0._jprb

      ! Second order shadowing (the facet cannot be seen)

      IF ((sunglint%s(j)%zensat * sunglint%s(j)%csi) <= 0._jprb) THEN
        IF (ABS(sunglint%s(j)%zensat) < realtol .AND. ABS(sunglint%s(j)%csi) < realtol) THEN
          sunglint%s(j)%q_shad = 1._jprb
        ELSE IF (sunglint%s(j)%zensat >= ABS(sunglint%s(j)%csi)) THEN
          sunglint%s(j)%a_shad   = 1._jprb / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o)
          sunglint%s(j)%lambda_a = &
              EXP(-sunglint%s(j)%a_shad**2_jpim / 2._jprb) / (sqrt2pi * sunglint%s(j)%a_shad) - &
              0.5_jprb * rttov_erfcx(sunglint%s(j)%a_shad / SQRT(2._jprb))
          sunglint%s(j)%q_shad   = theta_fi / (1._jprb + sunglint%s(j)%lambda_a)
        ELSE IF (ABS(sunglint%s(j)%csi) > sunglint%s(j)%zensat) THEN
          sunglint%s(j)%b_shad   = 1._jprb / (TAN(ABS(sunglint%s(j)%csi)) * sunglint%s(j)%gamma_o)
          sunglint%s(j)%lambda_b = &
              EXP(-sunglint%s(j)%b_shad**2_jpim / 2._jprb) / (sqrt2pi * sunglint%s(j)%b_shad) - &
              0.5_jprb * rttov_erfcx(sunglint%s(j)%b_shad / SQRT(2._jprb))
          sunglint%s(j)%q_shad   = theta_csi / (1._jprb + sunglint%s(j)%lambda_b)
        ENDIF
      ELSE
        sunglint%s(j)%a_shad   = 1._jprb / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o)
        sunglint%s(j)%b_shad   = 1._jprb / (TAN(ABS(sunglint%s(j)%csi)) * sunglint%s(j)%gamma_o)
        sunglint%s(j)%lambda_a = &
            EXP(-sunglint%s(j)%a_shad**2_jpim / 2._jprb) / (sqrt2pi * sunglint%s(j)%a_shad) - &
            0.5_jprb * rttov_erfcx(sunglint%s(j)%a_shad / SQRT(2._jprb))
        sunglint%s(j)%lambda_b = &
            EXP(-sunglint%s(j)%b_shad**2_jpim / 2._jprb) / (sqrt2pi * sunglint%s(j)%b_shad) - &
            0.5_jprb * rttov_erfcx(sunglint%s(j)%b_shad / SQRT(2._jprb))
        sunglint%s(j)%q_shad   = 1._jprb / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b)
        sunglint%s(j)%q_shad   = sunglint%s(j)%q_shad * theta_fi * theta_csi
      ENDIF

      ! Compute the probability density that the slope of the facet at
      ! a certain point is gammax when the incident ray and the ray from
      ! the reflection of the incident ray on the local surface do not
      ! intersect with any other surface.

      sunglint%s(j)%c_shad  = ABS(COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat)) * sunglint%s(j)%gamma_p
      sunglint%s(j)%p_prime = 1._jprb / (sqrt2pi * sunglint%s(j)%c_shad) * &
          EXP(-0.5_jprb * (TAN(sunglint%s(j)%alfa) / sunglint%s(j)%c_shad)**2_jpim)
      sunglint%s(j)%g_shad  = &
          (1._jprb - TAN((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb) * TAN(sunglint%s(j)%zensat))
      sunglint%s(j)%fac1    = 2._jprb * (COS(sunglint%s(j)%alfa))**2_jpim * &
          (COS((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2))**2_jpim
      sunglint%s(j)%pxy_gammaxy = (1._jprb / (sqrt2pi * sunglint%s(j)%gamma_o)) * &
          EXP(-sunglint%s(j)%gammax**2_jpim / (2._jprb * sunglint%s(j)%gamma_o**2_jpim))

      ! Compute effective distribution function

      sunglint%s(j)%glint = &
        sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * &
        sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1

    ENDIF ! wind speed > minimum
  ENDDO ! nprofiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  FUNCTION RTTOV_ERFCX(X)

!     Description:
!     RTTOV_ERFCX - To compute complementary error function
!
!     Method:
!     Chebyshev fitting
!
!     Owner:
!     Argonne National Laboratory
!
!     History:
!     Version      Date        Comment
!     1            15/07/2003  Marco Matricardi. ECMWF.
!                              Based on the routine CALREF written by
!                              W.J.Cody, Mathematics and Computer Science
!                              Division Argonne National Laboratory
!                              Argonne, IL 60439
!
    USE parkind1, ONLY : jprb
!INTF_OFF
    USE parkind1, ONLY : jpim
!INTF_ON
    IMPLICIT NONE

    REAL (jprb), INTENT(IN) :: X            ! Interpolation argument.
    REAL (jprb)             :: RTTOV_ERFCX
!INTF_END

    REAL    (jprb) :: A(5),B(4),C(9),D(8),P(6),Q(5)
    REAL    (jprb) :: Y,THRESH,YSQ,XSMALL,XNUM,XDEN,XBIG,SQRPI,DEL
    INTEGER (jpim) :: I
!-----End of header-------------------------------------------------------------

    DATA A/3.16112374387056560E00_JPRB,1.13864154151050156E02_JPRB,          &
           3.77485237685302021E02_JPRB,3.20937758913846947E03_JPRB,          &
           1.85777706184603153E-1_JPRB/
    DATA B/2.36012909523441209E01_JPRB,2.44024637934444173E02_JPRB,          &
           1.28261652607737228E03_JPRB,2.84423683343917062E03_JPRB/
    DATA C/5.64188496988670089E-1_JPRB,8.88314979438837594E0_JPRB,           &
           6.61191906371416295E01_JPRB,2.98635138197400131E02_JPRB,          &
           8.81952221241769090E02_JPRB,1.71204761263407058E03_JPRB,          &
           2.05107837782607147E03_JPRB,1.23033935479799725E03_JPRB,          &
           2.15311535474403846E-8_JPRB/
    DATA D/1.57449261107098347E01_JPRB,1.17693950891312499E02_JPRB,          &
           5.37181101862009858E02_JPRB,1.62138957456669019E03_JPRB,          &
           3.29079923573345963E03_JPRB,4.36261909014324716E03_JPRB,          &
           3.43936767414372164E03_JPRB,1.23033935480374942E03_JPRB/
    DATA P/3.05326634961232344E-1_JPRB,3.60344899949804439E-1_JPRB,          &
           1.25781726111229246E-1_JPRB,1.60837851487422766E-2_JPRB,          &
           6.58749161529837803E-4_JPRB,1.63153871373020978E-2_JPRB/
    DATA Q/2.56852019228982242E00_JPRB,1.87295284992346047E00_JPRB,          &
           5.27905102951428412E-1_JPRB,6.05183413124413191E-2_JPRB,          &
           2.33520497626869185E-3_JPRB/
    DATA THRESH /0.46875_JPRB/
    DATA XSMALL /5.96E-8_JPRB/
    DATA XBIG   /9.194_JPRB  /
    DATA SQRPI  /5.6418958354775628695E-1_JPRB/

    Y = ABS(X)

    IF (Y .LE. THRESH) THEN

      YSQ = 0.

      IF (Y .GT. XSMALL)THEN
        YSQ = Y * Y
      ENDIF

      XNUM = A(5) * YSQ
      XDEN = YSQ

        DO I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
        ENDDO

      RTTOV_ERFCX = 1 - X * (XNUM + A(4)) / (XDEN + B(4))

    ELSE IF (Y .LE. 4.) THEN

      XNUM = C(9)*Y
      XDEN = Y

      DO I = 1, 7
        XNUM = (XNUM + C(I)) * Y
        XDEN = (XDEN + D(I)) * Y
      ENDDO

      RTTOV_ERFCX = (XNUM + C(8)) / (XDEN + D(8))
      YSQ = AINT(Y * 16.) / 16.
      DEL = (Y - YSQ) * (Y + YSQ)
      RTTOV_ERFCX = EXP(-YSQ * YSQ) * EXP(-DEL) * RTTOV_ERFCX

    ELSE

      IF (Y .GE. XBIG) THEN
        RTTOV_ERFCX = 0.
      ELSE
        YSQ = 1. / (Y * Y)
        XNUM = P(6) * YSQ
        XDEN = YSQ
          DO I = 1, 4
            XNUM = (XNUM + P(I)) * YSQ
            XDEN = (XDEN + Q(I)) * YSQ
          ENDDO
        RTTOV_ERFCX = YSQ * (XNUM + P(5)) / (XDEN + Q(5))
        RTTOV_ERFCX = (SQRPI -  RTTOV_ERFCX) / Y
        YSQ = AINT(Y * 16.) / 16.
        DEL = (Y - YSQ) * (Y + YSQ)
        RTTOV_ERFCX = EXP(-YSQ * YSQ) * EXP(-DEL) * RTTOV_ERFCX
      ENDIF
    ENDIF

    IF (X < 0)THEN
      RTTOV_ERFCX = 2. - RTTOV_ERFCX
    ENDIF
  END FUNCTION RTTOV_ERFCX

END SUBROUTINE rttov_refsun
