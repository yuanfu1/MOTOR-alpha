! Description:
!> @file
!!   TL of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!
!> @brief
!!   TL of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     profiles_tl    input atmospheric profile and surface variable perturbations
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     aux            internal structure containing auxiliary profile variables
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_tl    sea surface BRDF model variable perturbations
!! @param[in]     raytracing     raytracing structure
!! @param[in]     raytracing_tl  raytracing variable perturbations
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
SUBROUTINE rttov_refsun_tl( &
              opts,          &
              profiles,      &
              profiles_tl,   &
              coef,          &
              aux,           &
              sunglint,      &
              sunglint_tl,   &
              raytracing,    &
              raytracing_tl)

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
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_tl(SIZE(profiles))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_raytracing),  INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),  INTENT(IN)    :: raytracing_tl
  TYPE(rttov_sunglint),    INTENT(IN)    :: sunglint
  TYPE(rttov_sunglint),    INTENT(INOUT) :: sunglint_tl
  TYPE(rttov_coef),        INTENT(IN)    :: coef
!INTF_END

  REAL(jprb)    :: psi_tl(coef%ws_nomega)              ! The frequency spectrum of the surface wave
  REAL(jprb)    :: ff_tl (coef%ws_nomega)              ! Working space

  INTEGER(jpim) :: j, i, ilaysur, nprofiles
  REAL(jprb)    :: csi_tl                              ! Angle between the zenith and the projection of the
                                                            ! incident ray on the x-z plane
  REAL(jprb)    :: alfa_tl                             ! Angle between the incident ray and the x-z plane
  REAL(jprb)    :: c_shad_tl                           ! The average magnitude of tan(alfa)
  REAL(jprb)    :: p_prime_tl                          ! The probability density of tan(alfa)
  REAL(jprb)    :: pxy_gammaxy_tl                      ! The joint probability density of the along-view and
                                                            ! cross view slope
  REAL(jprb)    :: gamma_sq_tl                         ! Total variance of the slope of the facet
  REAL(jprb)    :: gamma_o_tl                          ! The mean square of the along-view (X axis) slope
  REAL(jprb)    :: gamma_p_tl                          ! The mean square of the cross-view (Y axis) slope
  REAL(jprb)    :: g_shad_tl                           ! Normalization function
  REAL(jprb)    :: gammax_tl                           ! The x-slope of the facet
  REAL(jprb)    :: theta_fi                            ! First order shadowing factor for the reflected ray
  REAL(jprb)    :: theta_csi                           ! First ordesr shadowing factor for the incident ray
  REAL(jprb)    :: q_shad_tl                           ! Second order shadowing factor
  REAL(jprb)    :: q_shad_a_tl                         ! Second order shadowing factor
  REAL(jprb)    :: q_shad_b_tl                         ! Second order shadowing factor
  REAL(jprb)    :: zensat_tl                           ! Zenith angle of satellite viewing angle at surface
  REAL(jprb)    :: zensun_tl                           ! Zenith angle of sun at surface
  REAL(jprb)    :: omega_1                             ! Frequency of the surface wave
  REAL(jprb)    :: fac1_tl
  REAL(jprb)    :: a_shad_tl
  REAL(jprb)    :: b_shad_tl
  REAL(jprb)    :: lambda_a_tl
  REAL(jprb)    :: lambda_b_tl
  REAL(jprb)    :: x_u_tl
  REAL(jprb)    :: alfa1_tl
  REAL(jprb)    :: k
  REAL(jprb)    :: sigma_a2_r
  REAL(jprb)    :: sigma_b2_r
  REAL(jprb)    :: omega_m_tl
  REAL(jprb)    :: sigma2_r
  REAL(jprb)    :: beta_tl
  REAL(jprb)    :: windsp_tl
  REAL(jprb)    :: wangl_tl
  REAL(jprb)    :: s2_tl, s4_tl, sa_tl, sb_tl
  REAL(jprb)    :: bb, aa, hinc, sqrt2pi

  REAL(jprb)    :: max_exp_exponent_4rt, alfa1_r, logk
  REAL(jprb)    :: exp_max_exp_exponent_r, omega_m_3, omega_m_3_r
  REAL(jprb)    :: omega_1_array(coef%ws_nomega)
  REAL(jprb)    :: omega_1_4_r(coef%ws_nomega)

  REAL(jprb)    :: karr(0:nk_elf), sk2_tl, sk2o_tl(0:nk_elf), sk2p_tl(0:nk_elf)
  REAL(jprb)    :: sqrt_k_kp, cm_c_pow25, c_cp_pow25
  REAL(jprb)    :: cos2wangl, cos2wangl_tl
  REAL(jprb)    :: windsp, k0, x, kp, gamma, sigma
  REAL(jprb)    :: k0_tl, x_tl, kp_tl, gamma_tl, sigma_tl
  REAL(jprb)    :: cp, omega, alpha_p, ustar, am
  REAL(jprb)    :: cp_tl, omega_tl, alpha_p_tl, ustar_tl, alpha_m_tl, am_tl
  REAL(jprb)    :: fp, fp_tl, bl_tl, bh_tl
  REAL(jprb)    :: omega_c_tl, lpm_tl, jp_tl, gamma_exp_tl, fpexp_tl, dk_tl
  REAL(jprb)    :: c_r, kp_r

  REAL(jprb)    :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)
  sqrt2pi = SQRT(2_jpim * pi)
  k = 3.3_jprb
  sigma_a2_r = 1._jprb / 0.07_jprb**2
  sigma_b2_r = 1._jprb / 0.09_jprb**2

  max_exp_exponent_4rt = SQRT(SQRT(max_exp_exponent))
  logk = LOG(k) ! log(3.3)
  exp_max_exp_exponent_r = EXP(-max_exp_exponent)

  IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN
    DO i = 1, coef%ws_nomega
      omega_1_array(i) = (i - 1) * 0.1_jprb
    ENDDO
    omega_1_4_r(2:coef%ws_nomega) = 1._jprb / (omega_1_array(2:coef%ws_nomega)**4_jpim)
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

  sunglint_tl%s(:)%glint = 0._jprb
  sunglint_tl%s(:)%omega = 0._jprb
  DO j = 1, nprofiles
    IF (profiles(j)%sunzenangle < 0.0 .OR. &
        profiles(j)%sunzenangle > max_sol_zen) CYCLE
    IF (profiles(j)%skin%surftype /= surftype_sea) CYCLE

    IF (sunglint%s(j)%windsp <= min_windsp) THEN

      sunglint_tl%s(j)%glint = 0._jprb

    ELSE

      ! --------------------------------------------------------------------
      ! Calculate profile variables
      ! --------------------------------------------------------------------

      ! Compute the wind speed

      IF (sunglint%s(j)%windsp > min_windsp) THEN
        windsp_tl = (profiles(j)%s2m%u * profiles_tl(j)%s2m%u + &
                     profiles(j)%s2m%v * profiles_tl(j)%s2m%v) / sunglint%s(j)%windsp
      ELSE
        windsp_tl = (profiles_tl(j)%s2m%u + profiles_tl(j)%s2m%v) / SQRT(2._jprb)
      ENDIF

      ! Compute the angle between the wind direction and the U axis

      IF (sunglint%s(j)%windsp > min_windsp) THEN
        wangl_tl = (1._jprb / (profiles(j)%s2m%u**2_jpim + profiles(j)%s2m%v**2_jpim)) * &
            (profiles(j)%s2m%u * profiles_tl(j)%s2m%v - profiles(j)%s2m%v * profiles_tl(j)%s2m%u)
      ELSE
        wangl_tl = 0._jprb
      ENDIF

      ! Satellite and solar zenith angles at surface
      ilaysur = aux%s(j)%nearestlev_surf - 1
      IF (ABS(raytracing%pathsat(ilaysur, j) - 1._jprb) > realtol) THEN
        zensat_tl = raytracing_tl%pathsat(ilaysur, j) / &
            (raytracing%pathsat(ilaysur, j) * SQRT(raytracing%pathsat(ilaysur, j)**2_jpim - 1._jprb))
      ELSE
        zensat_tl = 0._jprb
      ENDIF
      IF (ABS(raytracing%pathsun(ilaysur, j) - 1._jprb) > realtol) THEN
        zensun_tl = raytracing_tl%pathsun(ilaysur, j) / &
            (raytracing%pathsun(ilaysur, j) * SQRT(raytracing%pathsun(ilaysur, j)**2_jpim - 1._jprb))
      ELSE
        zensun_tl = 0._jprb
      ENDIF

      ! ------------------------------------------------------------------
      ! Wave spectrum parameterisation
      ! ------------------------------------------------------------------

      cos2wangl = COS(2 * sunglint%s(j)%wangl)
      cos2wangl_tl = -2 * wangl_tl * SIN(2 * sunglint%s(j)%wangl)

      IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN

        ! JONSWAP wave frequency spectrum parameterisation

        x_u_tl     = profiles_tl(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp**2_jpim - &
            2._jprb * windsp_tl * profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp**3_jpim
        alfa1_tl   =  - 0.22_jprb * X_U_tl * sunglint%s(j)%alfa1 / sunglint%s(j)%x_u
        omega_m_tl =  - windsp_tl * sunglint%s(j)%omega_m / sunglint%s(j)%windsp - &
            x_u_tl * 0.33_jprb * sunglint%s(j)%omega_m / sunglint%s(j)%x_u

        ! Compute JONSWAP frequency spectrum of the surface wave

        alfa1_r = 1._jprb / sunglint%s(j)%alfa1
        omega_m_3 = sunglint%s(j)%omega_m**3_jpim
        omega_m_3_r = 1._jprb / omega_m_3
        psi_tl(1) = 0._jprb
        ff_tl(1) = 0._jprb

        DO i = 2, coef%ws_nomega
          omega_1 = omega_1_array(i)

          IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

            IF (sunglint%beta(i, j) >= exp_max_exp_exponent_r) THEN
              IF (omega_1 <= sunglint%s(j)%omega_m) THEN
                sigma2_r = sigma_a2_r
              ELSE
                sigma2_r = sigma_b2_r
              ENDIF

              beta_tl = omega_m_tl * sunglint%beta(i, j) * omega_1 * &
                  (omega_1 - sunglint%s(j)%omega_m) * sigma2_r * omega_m_3_r
            ELSE
              beta_tl = 0._jprb
            ENDIF

            psi_tl(i) = alfa1_tl * sunglint%psi(i, j) * alfa1_r - &
                omega_m_tl * sunglint%psi(i, j) * 5._jprb * omega_m_3 * omega_1_4_r(i)

            IF (sunglint%beta(i, j) > min_exponent) psi_tl(i) = psi_tl(i) + beta_tl * sunglint%psi(i, j) * logk

          ELSE
            psi_tl(i) = 0._jprb
          ENDIF
          ff_tl(i) = psi_tl(i) * coef%ws_k_omega(i)**2_jpim
        ENDDO

        ! Compute the Simpson's integral

        bb    = 301._jprb
        aa    = 0._jprb
        hinc  = (bb - aa) / coef%ws_nomega
        sa_tl = ff_tl(1)
        sb_tl = ff_tl(coef%ws_nomega)
        s4_tl = SUM(ff_tl(1::2))
        s2_tl = SUM(ff_tl(2::2))
        gamma_sq_tl = hinc / 3.0_jprb * (sa_tl + sb_tl + 4.0_jprb * s4_tl + 2.0_jprb * s2_tl)

        ! Compute the rms of the slope of the facet along the X (gamma_o)and Y (gamma_p) axis

        gamma_o_tl = 0.125_jprb * ((2 + cos2wangl) * gamma_sq_tl + &
                     cos2wangl_tl * sunglint%s(j)%gamma_sq) / sunglint%s(j)%gamma_o

        gamma_p_tl = 0.125_jprb * ((2 - cos2wangl) * gamma_sq_tl - &
                     cos2wangl_tl * sunglint%s(j)%gamma_sq) / sunglint%s(j)%gamma_p

      ELSEIF (opts%rt_ir%solar_sea_brdf_model == 2) THEN

        ! Elfouhaily et al wave spectrum parameterisation

        windsp = MAX(sunglint%s(j)%windsp, 0.6_jprb)
        IF (sunglint%s(j)%windsp < 0.6_jprb) windsp_tl = 0._jprb

        sk2o_tl(0) = 0._jprb
        sk2p_tl(0) = 0._jprb

        ! Active values independent of wavenumber (k)
        k0 = gravity / windsp**2
        k0_tl = -2 * gravity * windsp_tl / windsp**3

        x = k0 * profiles(j)%s2m%wfetc
        x_tl = k0_tl * profiles(j)%s2m%wfetc + k0 * profiles_tl(j)%s2m%wfetc

        omega_c_tl = 0.84_jprb * (-0.75_jprb) * (TANH((x / x0_elf)**0.4))**(-1.75) * &
                     (1._jprb / COSH((x / x0_elf)**0.4)**2) * &
                     (0.4_jprb * x_tl / (x**0.6_jprb * x0_elf**0.4_jprb))

        kp = k0 * sunglint%s(j)%omega_c**2
        kp_r = 1._jprb / kp
        kp_tl = k0_tl * sunglint%s(j)%omega_c**2 + &
                2 * k0 * omega_c_tl * sunglint%s(j)%omega_c

        IF (sunglint%s(j)%omega_c < 1) THEN
          gamma = 1.7_jprb
          gamma_tl = 0._jprb
        ELSE
          ! log10(x) = log(x) / log(10)
          gamma = 1.7_jprb + 6 * LOG10(sunglint%s(j)%omega_c)
          gamma_tl = 6 * omega_c_tl / (sunglint%s(j)%omega_c * LOG(10._jprb))
        ENDIF
        sigma = 0.08_jprb * (1 + 4 / sunglint%s(j)%omega_c**3)
        sigma_tl = -0.08_jprb * 4 * 3 * omega_c_tl / sunglint%s(j)%omega_c**4

        cp = SQRT(gravity * (1 + (kp / km_elf)**2) * kp_r)
        cp_tl = 0.5_jprb  * gravity * kp_tl * (1 / km_elf**2 - kp_r**2) / cp

        omega = windsp / cp
        omega_tl = windsp_tl / cp - windsp * cp_tl / cp**2

        alpha_p = 6E-3_jprb * SQRT(omega)
        alpha_p_tl = 6E-3_jprb * 0.5_jprb * omega_tl / SQRT(omega)

        ! Air friction velocity
        ustar = windsp / 25
        ustar_tl = windsp_tl / 25

        IF (ustar < cm_elf) THEN
          alpha_m_tl = 1E-2_jprb * ustar_tl / ustar
        ELSE
          alpha_m_tl = 1E-2_jprb * 3 * ustar_tl / ustar
        ENDIF

        ! Angular spread
        am = 0.13_jprb * ustar / cm_elf
        am_tl = 0.13_jprb * ustar_tl / cm_elf

        DO i = 1, nk_elf
          k = karr(i)
          sqrt_k_kp = SQRT(k * kp_r)
          c_r = 1._jprb / sunglint%c(i,j)

          ! Long-wave spectrum
          lpm_tl = -5 * 2 * kp_tl * kp * sunglint%lpm(i,j) / (4 * k**2)

          IF (sunglint%gamma_exp(i,j) >= exp_max_exp_exponent_r) THEN
            gamma_exp_tl = sunglint%gamma_exp(i,j) * (sqrt_k_kp - 1) * &
              (0.5_jprb * kp_tl * sqrt_k_kp / (kp * sigma**2) + &
               (sqrt_k_kp - 1) * sigma_tl / sigma**3)
            IF (sunglint%gamma_exp(i,j) > min_exponent) THEN
              jp_tl = sunglint%jp(i,j) * (gamma_exp_tl * LOG(gamma) + &
                        sunglint%gamma_exp(i,j) * gamma_tl / gamma)
            ELSE
              jp_tl = 0._jprb
            ENDIF
          ELSE
            jp_tl = 0._jprb
          ENDIF

          fpexp_tl = -(sunglint%fpexp(i,j) / SQRT(10._jprb)) * &
              ((sqrt_k_kp - 1) * omega_tl - 0.5_jprb * omega * kp_tl * sqrt_k_kp * kp_r)

          fp = sunglint%lpm(i,j) * sunglint%jp(i,j) * sunglint%fpexp(i,j)
          fp_tl = lpm_tl * sunglint%jp(i,j) * sunglint%fpexp(i,j) + &
                  sunglint%lpm(i,j) * jp_tl * sunglint%fpexp(i,j) + &
                  sunglint%lpm(i,j) * sunglint%jp(i,j) * fpexp_tl

          bl_tl = 0.5_jprb * (alpha_p_tl * cp * fp + &
                              alpha_p * cp_tl * fp + &
                              alpha_p * cp * fp_tl) * c_r

          ! Short-wave spectrum
          bh_tl = 0.5_jprb * alpha_m_tl * cm_elf * sunglint%fm(i,j) * c_r

          ! Combined spectrum
          sk2_tl = (bl_tl + bh_tl) / k

          IF (sunglint%sk2(i,j) > 0._jprb) THEN
            ! Angular spreading factor
            cm_c_pow25 = (cm_elf * c_r)**2.5
            c_cp_pow25 = (sunglint%c(i,j) / cp)**2.5
            dk_tl = (am_tl * cm_c_pow25 - 2.5_jprb * ap_elf * cp_tl * c_cp_pow25 / cp) / &
                    (COSH(a0_elf + ap_elf * c_cp_pow25 + am * cm_c_pow25))**2

            sk2o_tl(i) = sk2_tl * 0.5_jprb * (1 + 0.5_jprb * sunglint%dk(i,j) * cos2wangl) + &
                         sunglint%sk2(i,j) * 0.25_jprb * (dk_tl * cos2wangl + sunglint%dk(i,j) * cos2wangl_tl)
            sk2p_tl(i) = sk2_tl * 0.5_jprb * (1 - 0.5_jprb * sunglint%dk(i,j) * cos2wangl) - &
                         sunglint%sk2(i,j) * 0.25_jprb * (dk_tl * cos2wangl + sunglint%dk(i,j) * cos2wangl_tl)
          ELSE
            sk2o_tl(i) = 0._jprb
            sk2p_tl(i) = 0._jprb
          ENDIF
        ENDDO

        ! Trapezium integration
        gamma_o_tl = 0.25_jprb * SUM((sk2o_tl(0:nk_elf-1) + sk2o_tl(1:nk_elf)) * &
                     (karr(1:nk_elf) - karr(0:nk_elf-1))) / sunglint%s(j)%gamma_o
        gamma_p_tl = 0.25_jprb * SUM((sk2p_tl(0:nk_elf-1) + sk2p_tl(1:nk_elf)) * &
                     (karr(1:nk_elf) - karr(0:nk_elf-1))) / sunglint%s(j)%gamma_p

      ENDIF ! wave spectrum option


      ! ------------------------------------------------------------------
      ! Yoshimori wave facet reflectance model
      ! ------------------------------------------------------------------

      ! Obtain angles csi and alfa

      csi_tl = &
          zensun_tl * COS(pi - sunglint%s(j)%dazng * deg2rad) / (COS(sunglint%s(j)%zensun)**2_jpim * &
          (1._jprb + (TAN(sunglint%s(j)%zensun) * COS(pi - sunglint%s(j)%dazng * deg2rad))**2_jpim))
      IF (ABS((SIN(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad))**2_jpim - 1._jprb) > realtol) THEN
        alfa_tl = zensun_tl * COS(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad) / &
            (SQRT(1._jprb - (SIN(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad))**2_jpim))
      ELSE
        alfa_tl = 0._jprb
      ENDIF
      IF ((sunglint%s(j)%csi + sunglint%s(j)%zensat) >= 0._jprb) THEN
        sunglint_tl%s(j)%omega = (csi_tl + zensat_tl) / 2._jprb
      ELSE
        sunglint_tl%s(j)%omega =  -(csi_tl + zensat_tl) / 2._jprb
      ENDIF

      ! Compute the value of the function that represents the
      ! shadowing of the incident and reflected ray

      ! First order shadowing (the slope of the facet is negative)

      gammax_tl = &
          (csi_tl - zensat_tl) / (2._jprb * (COS((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb))**2_jpim)
      IF (ABS(sunglint%s(j)%zensat) > 0._jprb) THEN
        IF (1._jprb / TAN(ABS(sunglint%s(j)%zensat)) - &
            sunglint%s(j)%gammax * sunglint%s(j)%zensat / ABS(sunglint%s(j)%zensat) >= 0._jprb) THEN
          theta_fi = 1._jprb
        ELSE
          theta_fi = 0._jprb
        ENDIF
      ELSE
        theta_fi = 1._jprb
      ENDIF
      IF (ABS(sunglint%s(j)%csi) > 0._jprb) THEN
        IF (1._jprb / TAN(ABS(sunglint%s(j)%csi)) + &
            sunglint%s(j)%gammax * sunglint%s(j)%csi / ABS(sunglint%s(j)%csi) >= 0._jprb) THEN
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
          q_shad_tl = 0._jprb
        ELSE IF (sunglint%s(j)%zensat >= ABS(sunglint%s(j)%csi)) THEN
          a_shad_tl   = &
              -zensat_tl / (SIN(sunglint%s(j)%zensat)**2_jpim * sunglint%s(j)%gamma_o) - &
              gamma_o_tl / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o**2_jpim)
          lambda_a_tl = -a_shad_tl * (EXP(-sunglint%s(j)%a_shad**2_jpim / 2._jprb)) / &
              (sqrt2pi * sunglint%s(j)%a_shad**2_jpim)
          q_shad_a_tl = -1._jprb / (1._jprb + sunglint%s(j)%lambda_a)**2_jpim * lambda_a_tl
          q_shad_tl   = q_shad_a_tl * theta_fi
        ELSE IF (ABS(sunglint%s(j)%csi) > sunglint%s(j)%zensat) THEN
          IF (sunglint%s(j)%csi >= 0._jprb) THEN
            b_shad_tl = -csi_tl / (SIN(sunglint%s(j)%csi)**2_jpim * sunglint%s(j)%gamma_o) - &
                 gamma_o_tl / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2_jpim)
          ELSE
            b_shad_tl = &
                csi_tl / (SIN(sunglint%s(j)%csi)**2_jpim * sunglint%s(j)%gamma_o) +  &
                gamma_o_tl / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2_jpim)
          ENDIF
          lambda_b_tl = -b_shad_tl * (EXP(-sunglint%s(j)%b_shad**2_jpim / 2._jprb)) / &
              (sqrt2pi * sunglint%s(j)%b_shad**2_jpim)
          q_shad_b_tl = -lambda_b_tl / (1._jprb + sunglint%s(j)%lambda_b)**2_jpim
          q_shad_tl   = q_shad_b_tl * theta_csi
        ENDIF
      ELSE
        a_shad_tl = &
            -zensat_tl / (SIN(sunglint%s(j)%zensat)**2_jpim * sunglint%s(j)%gamma_o) - &
            gamma_o_tl / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o**2_jpim)
        IF (sunglint%s(j)%csi >= 0._jprb) THEN
          b_shad_tl = -csi_tl / (SIN(sunglint%s(j)%csi)**2_jpim * sunglint%s(j)%gamma_o) - &
              gamma_o_tl / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2_jpim)
        ELSE
          b_shad_tl = csi_tl / (SIN(sunglint%s(j)%csi)**2_jpim * sunglint%s(j)%gamma_o) + &
              gamma_o_tl / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2_jpim)
        ENDIF
        lambda_a_tl = -a_shad_tl * (EXP(-sunglint%s(j)%a_shad**2_jpim / 2._jprb)) / &
            (sqrt2pi * sunglint%s(j)%a_shad**2_jpim)
        lambda_b_tl = -b_shad_tl * (EXP(-sunglint%s(j)%b_shad**2_jpim / 2._jprb)) / &
            (sqrt2pi * sunglint%s(j)%b_shad**2_jpim)
        q_shad_tl   = -(lambda_a_tl + lambda_b_tl) / &
            (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b)**2_jpim
        q_shad_tl   = q_shad_tl * theta_fi * theta_csi
      ENDIF

      ! Compute the probability density that the slope of the facet at
      ! a certain point is gammax when the incident ray and the ray from
      ! the reflection of the incident ray on the local surface do not
      ! intersect with any other surface.

      IF (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat) >= 0._jprb) THEN
        c_shad_tl = gamma_p_tl * (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat)) - &
            csi_tl * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%csi) - &
            zensat_tl * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%zensat)
      ELSE
        c_shad_tl = -gamma_p_tl * (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat)) + &
            csi_tl * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%csi) + &
            zensat_tl * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%zensat)
      ENDIF
      p_prime_tl             = -alfa_tl * sunglint%s(j)%p_prime * TAN(sunglint%s(j)%alfa) / &
          (COS(sunglint%s(j)%alfa) * sunglint%s(j)%c_shad)**2_jpim + c_shad_tl * sunglint%s(j)%p_prime * &
          ((TAN(sunglint%s(j)%alfa))**2_jpim / sunglint%s(j)%c_shad**3_jpim - 1._jprb / sunglint%s(j)%c_shad)
      g_shad_tl              = -0.5_jprb * (csi_tl - zensat_tl) * TAN(sunglint%s(j)%zensat) / &
          COS(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat))**2_jpim - &
          zensat_tl * TAN(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%ZENSAT)) / &
          COS(sunglint%s(j)%zensat)**2_jpim
      fac1_tl                = -alfa_tl * 2._jprb * sunglint%s(j)%fac1 * TAN(sunglint%s(j)%alfa) - &
          (csi_tl - zensat_tl) * sunglint%s(j)%fac1 * TAN((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
      pxy_gammaxy_tl         = gamma_o_tl * sunglint%s(j)%pxy_gammaxy *                                           &
          (sunglint%s(j)%gammax**2_jpim / sunglint%s(j)%gamma_o**3_jpim - 1._jprb / sunglint%s(j)%gamma_o) -  &
          gammax_tl * sunglint%s(j)%pxy_gammaxy * sunglint%s(j)%gammax / sunglint%s(j)%gamma_o**2_jpim

      ! Compute effective distribution function

      sunglint_tl%s(j)%glint = &
          q_shad_tl * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1 +  &
          g_shad_tl * sunglint%s(j)%q_shad * sunglint%s(j)%p_prime * sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1 +  &
          p_prime_tl * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1 +  &
          pxy_gammaxy_tl * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime / sunglint%s(j)%fac1 -  &
          fac1_tl * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * sunglint%s(j)%pxy_gammaxy /  &
          sunglint%s(j)%fac1**2_jpim

    ENDIF ! wind speed > minimum
  ENDDO ! nprofiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_refsun_tl
