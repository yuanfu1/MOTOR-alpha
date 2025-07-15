! Description:
!> @file
!!   AD of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!
!> @brief
!!   AD of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad    gradient wrt atmospheric profile and surface variables
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     aux            internal structure containing auxiliary profile variables
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_ad    gradient wrt sea surface BRDF model variables
!! @param[in]     raytracing     raytracing structure
!! @param[in,out] raytracing_ad  gradient wrt raytracing variables
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
SUBROUTINE rttov_refsun_ad( &
              opts,          &
              profiles,      &
              profiles_ad,   &
              coef,          &
              aux,           &
              sunglint,      &
              sunglint_ad,   &
              raytracing,    &
              raytracing_ad)

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
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_ad(SIZE(profiles))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_raytracing),  INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),  INTENT(INOUT) :: raytracing_ad
  TYPE(rttov_sunglint),    INTENT(IN)    :: sunglint
  TYPE(rttov_sunglint),    INTENT(INOUT) :: sunglint_ad
  TYPE(rttov_coef),        INTENT(IN)    :: coef
!INTF_END

  REAL(jprb)    :: psi_ad(coef%ws_nomega)              ! The frequency spectrum of the surface wave
  REAL(jprb)    :: ff_ad (coef%ws_nomega)              ! Working space

  INTEGER(jpim) :: j, i, ilaysur, nprofiles
  REAL(jprb)    :: csi_ad                              ! Angle between the zenith and the proJection of the
                                                            ! incident ray on the x-z plane
  REAL(jprb)    :: alfa_ad                             ! Angle between the incident ray and the x-z plane
  REAL(jprb)    :: c_shad_ad                           ! The average magnitude of tan(alfa)
  REAL(jprb)    :: p_prime_ad                          ! The probability density of tan(alfa)
  REAL(jprb)    :: pxy_gammaxy_ad                      ! The joint probability density of the along-view and
                                                            ! cross view slope
  REAL(jprb)    :: gamma_sq_ad                         ! Total variance of the slope of the facet
  REAL(jprb)    :: gamma_o_ad                          ! The mean square of the along-view (X axis) slope
  REAL(jprb)    :: gamma_p_ad                          ! The mean square of the cross-view (Y axis) slope
  REAL(jprb)    :: g_shad_ad                           ! Normalization function
  REAL(jprb)    :: gammax_ad                           ! The x-slope of the facet
  REAL(jprb)    :: theta_fi                            ! First order shadowing factor for the reflected ray
  REAL(jprb)    :: theta_csi                           ! First ordesr shadowing factor for the incident ray
  REAL(jprb)    :: q_shad_ad                           ! Second order shadowing factor
  REAL(jprb)    :: q_shad_a_ad                         ! Second order shadowing factor
  REAL(jprb)    :: q_shad_b_ad                         ! Second order shadowing factor
  REAL(jprb)    :: zensat_ad                           ! Zenith angle of satellite viewing angle at surface
  REAL(jprb)    :: zensun_ad                           ! Zenith angle of sun at surface
  REAL(jprb)    :: omega_1                             ! Frequency of the surface wave
  REAL(jprb)    :: fac1_ad
  REAL(jprb)    :: a_shad_ad
  REAL(jprb)    :: b_shad_ad
  REAL(jprb)    :: lambda_a_ad
  REAL(jprb)    :: lambda_b_ad
  REAL(jprb)    :: x_u_ad
  REAL(jprb)    :: alfa1_ad
  REAL(jprb)    :: k
  REAL(jprb)    :: sigma_a2_r
  REAL(jprb)    :: sigma_b2_r
  REAL(jprb)    :: omega_m_ad
  REAL(jprb)    :: sigma2_r
  REAL(jprb)    :: beta_ad
  REAL(jprb)    :: windsp_ad
  REAL(jprb)    :: wangl_ad
  REAL(jprb)    :: s2_ad, s4_ad, sa_ad, sb_ad
  REAL(jprb)    :: bb, aa, hinc, sqrt2pi

  REAL(jprb)    :: max_exp_exponent_4rt, alfa1_r, logk
  REAL(jprb)    :: exp_max_exp_exponent_r, omega_m_3, omega_m_3_r
  REAL(jprb)    :: omega_1_array(coef%ws_nomega)
  REAL(jprb)    :: omega_1_4_r(coef%ws_nomega)

  REAL(jprb)    :: karr(0:nk_elf), sk2_ad, sk2o_ad(0:nk_elf), sk2p_ad(0:nk_elf)
  REAL(jprb)    :: sqrt_k_kp, cm_c_pow25, c_cp_pow25
  REAL(jprb)    :: cos2wangl, cos2wangl_ad
  REAL(jprb)    :: windsp, k0, x, kp, gamma, sigma
  REAL(jprb)    :: k0_ad, x_ad, kp_ad, gamma_ad, sigma_ad
  REAL(jprb)    :: cp, omega, alpha_p, ustar, am
  REAL(jprb)    :: cp_ad, omega_ad, alpha_p_ad, ustar_ad, alpha_m_ad, am_ad
  REAL(jprb)    :: fp, fp_ad, bl_ad, bh_ad
  REAL(jprb)    :: omega_c_ad, lpm_ad, jp_ad, gamma_exp_ad, fpexp_ad, dk_ad
  REAL(jprb)    :: cosh_tmp, cp_r, c_r, kp_r, sigma_r2, sigma_r3

  REAL(jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)
  sqrt2pi = SQRT(2 * pi)
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

  DO j = 1, nprofiles
    IF (profiles(j)%sunzenangle < 0.0 .OR. &
        profiles(j)%sunzenangle > max_sol_zen) CYCLE
    IF (profiles(j)%skin%surftype /= surftype_sea) CYCLE

    windsp_ad      = 0._jprb
    wangl_ad       = 0._jprb
    q_shad_ad      = 0._jprb
    q_shad_a_ad    = 0._jprb
    q_shad_b_ad    = 0._jprb
    g_shad_ad      = 0._jprb
    p_prime_ad     = 0._jprb
    pxy_gammaxy_ad = 0._jprb
    fac1_ad        = 0._jprb
    gamma_o_ad     = 0._jprb
    gammax_ad      = 0._jprb
    alfa_ad        = 0._jprb
    csi_ad         = 0._jprb
    zensat_ad      = 0._jprb
    zensun_ad      = 0._jprb
    c_shad_ad      = 0._jprb
    gamma_p_ad     = 0._jprb
    gamma_sq_ad    = 0._jprb
    lambda_a_ad    = 0._jprb
    lambda_b_ad    = 0._jprb
    b_shad_ad      = 0._jprb
    a_shad_ad      = 0._jprb

    IF (sunglint%s(j)%windsp <= min_windsp) THEN

      sunglint_ad%s(j)%glint = 0._jprb

    ELSE

      ! ------------------------------------------------------------------
      ! Yoshimori wave facet reflectance model
      ! ------------------------------------------------------------------

      ! Compute effective distribution function

      q_shad_ad      = q_shad_ad + &
          sunglint_ad%s(j)%glint * sunglint%s(j)%g_shad * sunglint%s(j)%p_prime * &
          sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1
      g_shad_ad      = g_shad_ad + &
          sunglint_ad%s(j)%glint * sunglint%s(j)%q_shad * sunglint%s(j)%p_prime * &
          sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1
      p_prime_ad     = p_prime_ad + &
          sunglint_ad%s(j)%glint * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * &
          sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1
      pxy_gammaxy_ad = pxy_gammaxy_ad + &
          sunglint_ad%s(j)%glint * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * &
          sunglint%s(j)%p_prime / sunglint%s(j)%fac1
      fac1_ad        = fac1_ad - &
          sunglint_ad%s(j)%glint * sunglint%s(j)%q_shad * sunglint%s(j)%g_shad * &
          sunglint%s(j)%p_prime * sunglint%s(j)%pxy_gammaxy / sunglint%s(j)%fac1**2

      ! Compute the probability density that the slope of the facet at
      ! a certain point is gammax when the incident ray and the ray from
      ! the reflection of the incident ray on the local surface do not
      ! intersect with any other surface.

      gamma_o_ad     = gamma_o_ad + pxy_gammaxy_ad * sunglint%s(j)%pxy_gammaxy * &
          (sunglint%s(j)%gammax**2 / sunglint%s(j)%gamma_o**3 - 1._jprb / sunglint%s(j)%gamma_o)
      gammax_ad      = gammax_ad - &
          pxy_gammaxy_ad * sunglint%s(j)%pxy_gammaxy * sunglint%s(j)%gammax / sunglint%s(j)%gamma_o**2
      pxy_gammaxy_ad = 0._jprb
      alfa_ad        = alfa_ad - fac1_ad * 2._jprb * sunglint%s(j)%fac1 * TAN(sunglint%s(j)%alfa)
      csi_ad         = csi_ad - &
          fac1_ad * sunglint%s(j)%fac1 * TAN((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
      zensat_ad      = zensat_ad + &
          fac1_ad * sunglint%s(j)%fac1 * TAN((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb)
      fac1_ad        = 0._jprb

      csi_ad         = csi_ad - g_shad_ad * 0.5_jprb * TAN(sunglint%s(j)%zensat) / &
          COS(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat))**2
      zensat_ad      = zensat_ad + 0.5_jprb * g_shad_ad * TAN(sunglint%s(j)%zensat) / &
          COS(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat))**2
      zensat_ad      = zensat_ad - g_shad_ad * TAN(0.5_jprb * (sunglint%s(j)%csi - sunglint%s(j)%zensat)) / &
          COS(sunglint%s(j)%zensat)**2
      g_shad_ad      = 0._jprb

      alfa_ad        = alfa_ad - p_prime_ad * sunglint%s(j)%p_prime * TAN(sunglint%s(j)%alfa) / &
          (COS(sunglint%s(j)%alfa) * sunglint%s(j)%c_shad)**2
      c_shad_ad      = c_shad_ad + p_prime_ad * sunglint%s(j)%p_prime * &
          ((TAN(sunglint%s(j)%alfa))**2 / sunglint%s(j)%c_shad**3 - 1._jprb / sunglint%s(j)%c_shad)
      p_prime_ad     = 0._jprb

      IF (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat) >= 0._jprb) THEN
        zensat_ad  = zensat_ad - c_shad_ad * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%zensat)
        csi_ad     = csi_ad - c_shad_ad * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%csi)
        gamma_p_ad = gamma_p_ad + c_shad_ad * (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat))
        c_shad_ad  = 0._jprb
      ELSE
        zensat_ad  = zensat_ad + c_shad_ad * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%zensat)
        csi_ad     = csi_ad + c_shad_ad * sunglint%s(j)%gamma_p * SIN(sunglint%s(j)%csi)
        gamma_p_ad = gamma_p_ad - c_shad_ad * (COS(sunglint%s(j)%csi) + COS(sunglint%s(j)%zensat))
        c_shad_ad  = 0._jprb
      ENDIF

      ! Direct model - first order shadowing

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
          q_shad_ad = 0._jprb
        ELSE IF (sunglint%s(j)%zensat >= ABS(sunglint%s(j)%csi)) THEN
          q_shad_a_ad = q_shad_a_ad + q_shad_ad * theta_fi
          q_shad_ad   = 0._jprb
          lambda_a_ad = lambda_a_ad - q_shad_a_ad / (1._jprb + sunglint%s(j)%lambda_a)**2
          q_shad_a_ad = 0._jprb
          a_shad_ad   = a_shad_ad - lambda_a_ad * EXP(-sunglint%s(j)%a_shad**2 / 2._jprb) / &
              (sqrt2pi * sunglint%s(j)%a_shad**2)
          lambda_a_ad = 0._jprb
          zensat_ad   = zensat_ad - a_shad_ad / (SIN(sunglint%s(j)%zensat)**2 * sunglint%s(j)%gamma_o)
          gamma_o_ad  = gamma_o_ad - a_shad_ad / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o**2)
          a_shad_ad   = 0._jprb
        ELSE IF (ABS(sunglint%s(j)%csi) > sunglint%s(j)%zensat) THEN
          q_shad_b_ad = q_shad_b_ad + q_shad_ad * theta_csi
          q_shad_ad   = 0._jprb
          lambda_b_ad = lambda_b_ad - q_shad_b_ad / (1._jprb + sunglint%s(j)%lambda_b)**2
          q_shad_b_ad = 0._jprb
          b_shad_ad   = b_shad_ad - lambda_b_ad * EXP(-sunglint%s(j)%b_shad**2 / 2._jprb) / &
              (sqrt2pi * sunglint%s(j)%b_shad**2)
          lambda_b_ad = 0._jprb
          IF (sunglint%s(j)%csi >= 0._jprb) THEN
            csi_ad     = csi_ad - b_shad_ad / (SIN(sunglint%s(j)%csi)**2 * sunglint%s(j)%gamma_o)
            gamma_o_ad = gamma_o_ad - b_shad_ad / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2)
            b_shad_ad  = 0._jprb
          ELSE
            csi_ad     = csi_ad + b_shad_ad / (SIN(sunglint%s(j)%csi)**2 * sunglint%s(j)%gamma_o)
            gamma_o_ad = gamma_o_ad + b_shad_ad / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2)
            b_shad_ad  = 0._jprb
          ENDIF
        ENDIF
      ELSE
        q_shad_ad   = q_shad_ad * theta_fi * theta_csi
        lambda_a_ad = lambda_a_ad - q_shad_ad / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b)**2
        lambda_b_ad = lambda_b_ad - q_shad_ad / (1._jprb + sunglint%s(j)%lambda_a + sunglint%s(j)%lambda_b)**2
        q_shad_ad   = 0._jprb
        b_shad_ad   = b_shad_ad - lambda_b_ad * EXP(-sunglint%s(j)%b_shad**2 / 2._jprb) / &
            (sqrt2pi * sunglint%s(j)%b_shad**2)
        lambda_b_ad = 0._jprb
        a_shad_ad   = a_shad_ad - lambda_a_ad * EXP(-sunglint%s(j)%a_shad**2 / 2._jprb) / &
            (sqrt2pi * sunglint%s(j)%a_shad**2)
        lambda_a_ad = 0._jprb
        IF (sunglint%s(j)%csi >= 0._jprb) THEN
          csi_ad     = csi_ad - b_shad_ad / (SIN(sunglint%s(j)%csi)**2 * sunglint%s(j)%gamma_o)
          gamma_o_ad = gamma_o_ad - b_shad_ad / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2)
          b_shad_ad  = 0._jprb
        ELSE IF (sunglint%s(j)%csi < 0._jprb) THEN
          csi_ad     = csi_ad + b_shad_ad / (SIN(sunglint%s(j)%csi)**2 * sunglint%s(j)%gamma_o)
          gamma_o_ad = gamma_o_ad + b_shad_ad / (TAN(sunglint%s(j)%csi) * sunglint%s(j)%gamma_o**2)
          b_shad_ad  = 0._jprb
        ENDIF
        zensat_ad  = zensat_ad - a_shad_ad / (SIN(sunglint%s(j)%zensat)**2 * sunglint%s(j)%gamma_o)
        gamma_o_ad = gamma_o_ad - a_shad_ad / (TAN(sunglint%s(j)%zensat) * sunglint%s(j)%gamma_o**2)
        a_shad_ad  = 0._jprb
      ENDIF

      ! First order shadowing (the slope of the facet is negative)

      csi_ad    = csi_ad + &
          gammax_ad / (2._jprb * (COS((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb))**2)
      zensat_ad = zensat_ad - &
          gammax_ad / (2._jprb * (COS((sunglint%s(j)%csi - sunglint%s(j)%zensat) / 2._jprb))**2)
      gammax_ad = 0._jprb

      ! Obtain angles csi and alfa

      IF (sunglint%s(j)%csi + sunglint%s(j)%zensat >= 0._jprb) THEN
        csi_ad    = csi_ad + sunglint_ad%s(j)%omega / 2._jprb
        zensat_ad = zensat_ad + sunglint_ad%s(j)%omega / 2._jprb
      ELSE
        csi_ad    = csi_ad - sunglint_ad%s(j)%omega / 2._jprb
        zensat_ad = zensat_ad - sunglint_ad%s(j)%omega / 2._jprb
      ENDIF
      IF (ABS((SIN(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad))**2 - 1._jprb) > realtol) THEN
        zensun_ad = zensun_ad + alfa_ad * (COS(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad) / &
            SQRT(1._jprb - (SIN(sunglint%s(j)%zensun) * SIN(pi - sunglint%s(j)%dazng * deg2rad))**2))
      ENDIF
      alfa_ad = 0._jprb
      zensun_ad   = zensun_ad + &
          csi_ad * COS(pi - sunglint%s(j)%dazng * deg2rad) / (COS(sunglint%s(j)%zensun)**2 * &
            (1._jprb + (TAN(sunglint%s(j)%zensun) * COS(pi - sunglint%s(j)%dazng * deg2rad))**2))
      csi_ad = 0._jprb


      ! ------------------------------------------------------------------
      ! Wave spectrum parameterisation
      ! ------------------------------------------------------------------

      cos2wangl = COS(2 * sunglint%s(j)%wangl)
      cos2wangl_ad = 0._jprb

      IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN

        ! JONSWAP wave frequency spectrum parameterisation

        ! AD init

        sa_ad          = 0._jprb
        sb_ad          = 0._jprb
        s4_ad          = 0._jprb
        s2_ad          = 0._jprb
        ff_ad(:)       = 0._jprb
        beta_ad        = 0._jprb
        alfa1_ad       = 0._jprb
        omega_m_ad     = 0._jprb
        psi_ad(:)      = 0._jprb
        x_u_ad         = 0._jprb

        ! Compute the rms of the slope of the facet along the X (gamma_o)and Y (gamma_p) axis

        gamma_sq_ad = gamma_sq_ad + 0.125_jprb * (2 - cos2wangl) * gamma_p_ad / sunglint%s(j)%gamma_p
        cos2wangl_ad = cos2wangl_ad - 0.125_jprb * gamma_p_ad * sunglint%s(j)%gamma_sq / sunglint%s(j)%gamma_p

        gamma_sq_ad = gamma_sq_ad + 0.125_jprb * (2 + cos2wangl) * gamma_o_ad / sunglint%s(j)%gamma_o
        cos2wangl_ad = cos2wangl_ad + 0.125_jprb * gamma_o_ad * sunglint%s(j)%gamma_sq / sunglint%s(j)%gamma_o

        ! Compute the Simpson's integral

        bb    = 301._jprb
        aa    = 0._jprb
        hinc  = (bb - aa) / coef%ws_nomega
        sa_ad       = sa_ad + gamma_sq_ad * hinc / 3.0_jprb
        sb_ad       = sb_ad + gamma_sq_ad * hinc / 3.0_jprb
        s4_ad       = s4_ad + gamma_sq_ad * hinc / 3.0_jprb * 4.0_jprb
        s2_ad       = s2_ad + gamma_sq_ad * hinc / 3.0_jprb * 2.0_jprb
        gamma_sq_ad = 0._jprb
        ff_ad(2::2) = ff_ad(2::2) + s2_ad
        s2_ad = 0._jprb
        ff_ad(1::2) = ff_ad(1::2) + s4_ad
        s4_ad = 0._jprb
        ff_ad(1)              = ff_ad(1) + sa_ad
        ff_ad(coef%ws_nomega) = ff_ad(coef%ws_nomega) + sb_ad
        sa_ad = 0._jprb
        sb_ad = 0._jprb

        ! Compute JONSWAP frequency spectrum of the surface wave

        alfa1_r = 1._jprb / sunglint%s(j)%alfa1
        omega_m_3 = sunglint%s(j)%omega_m**3
        omega_m_3_r = 1._jprb / omega_m_3
        DO i = 2, coef%ws_nomega
          omega_1 = omega_1_array(i)

          psi_ad(i) = psi_ad(i) + ff_ad(i) * coef%ws_k_omega(i)**2

          IF (sunglint%s(j)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

            IF (sunglint%beta(i, j) > min_exponent) beta_ad = beta_ad + psi_ad(i) * sunglint%psi(i, j) * logk

            omega_m_ad = omega_m_ad - &
                psi_ad(i) * sunglint%psi(i, j) * 5._jprb * omega_m_3 * omega_1_4_r(i)
            alfa1_ad   = alfa1_ad + psi_ad(i) * sunglint%psi(i, j) * alfa1_r

            IF (sunglint%beta(i, j) >= exp_max_exp_exponent_r) THEN
              IF (omega_1 <= sunglint%s(j)%omega_m) THEN
                sigma2_r = sigma_a2_r
              ELSE
                sigma2_r = sigma_b2_r
              ENDIF

              omega_m_ad = omega_m_ad + beta_ad * sunglint%beta(i, j) * omega_1 * &
                  (omega_1 - sunglint%s(j)%omega_m) * sigma2_r * omega_m_3_r
            ENDIF
          ENDIF
          beta_ad = 0._jprb
        ENDDO
!         ff_ad(:) = 0._jprb
!         psi_ad(:) = 0._jprb

        windsp_ad                = windsp_ad - omega_m_ad * sunglint%s(j)%omega_m / sunglint%s(j)%windsp
        x_u_ad                   = x_u_ad - omega_m_ad * 0.33_jprb * sunglint%s(j)%omega_m / sunglint%s(j)%x_u
        omega_m_ad               = 0._jprb
        x_u_ad                   = x_u_ad - 0.22_jprb * alfa1_ad * sunglint%s(j)%alfa1 / sunglint%s(j)%x_u
        alfa1_ad                 = 0._jprb
        windsp_ad                =      &
            windsp_ad - x_u_ad * 2._jprb * profiles(j)%s2m%wfetc * gravity / sunglint%s(j)%windsp**3
        profiles_ad(j)%s2m%wfetc = profiles_ad(j)%s2m%wfetc + x_u_ad * gravity / sunglint%s(j)%windsp**2
        x_u_ad                   = 0._jprb

      ELSEIF (opts%rt_ir%solar_sea_brdf_model == 2) THEN

        ! Elfouhaily et al wave spectrum parameterisation

        ! AD init

        am_ad = 0._jprb
        cp_ad = 0._jprb
        alpha_p_ad = 0._jprb
        alpha_m_ad = 0._jprb
        omega_ad = 0._jprb
        gamma_ad = 0._jprb
        sigma_ad = 0._jprb
        kp_ad = 0._jprb

        ! Direct model

        windsp = MAX(sunglint%s(j)%windsp, 0.6_jprb)

        ! Active values independent of wavenumber (k)
        k0 = gravity / windsp**2
        x = k0 * profiles(j)%s2m%wfetc
        kp = k0 * sunglint%s(j)%omega_c**2
        kp_r = 1._jprb / kp

        IF (sunglint%s(j)%omega_c < 1) THEN
          gamma = 1.7_jprb
        ELSE
          gamma = 1.7_jprb + 6 * LOG10(sunglint%s(j)%omega_c)
        ENDIF
        sigma = 0.08_jprb * (1 + 4 / sunglint%s(j)%omega_c**3)
        sigma_r2 = 1._jprb / sigma**2
        sigma_r3 = 1._jprb / sigma**3

        cp = SQRT(gravity * (1 + (kp / km_elf)**2) * kp_r)
        cp_r = 1._jprb / cp
        omega = windsp * cp_r
        alpha_p = 6E-3_jprb * SQRT(omega)

        ! Air friction velocity
        ustar = windsp / 25

        ! Angular spread
        am = 0.13_jprb * ustar / cm_elf

        ! AD model

        ! Trapezium integration
        sk2o_ad(1:nk_elf-1) = 0.25_jprb * gamma_o_ad * (karr(2:nk_elf) - karr(0:nk_elf-2)) / sunglint%s(j)%gamma_o
        sk2o_ad(0)          = 0.25_jprb * gamma_o_ad * (karr(1) - karr(0)) / sunglint%s(j)%gamma_o
        sk2o_ad(nk_elf)     = 0.25_jprb * gamma_o_ad * (karr(nk_elf) - karr(nk_elf-1)) / sunglint%s(j)%gamma_o

        sk2p_ad(1:nk_elf-1) = 0.25_jprb * gamma_p_ad * (karr(2:nk_elf) - karr(0:nk_elf-2)) / sunglint%s(j)%gamma_p
        sk2p_ad(0)          = 0.25_jprb * gamma_p_ad * (karr(1) - karr(0)) / sunglint%s(j)%gamma_p
        sk2p_ad(nk_elf)     = 0.25_jprb * gamma_p_ad * (karr(nk_elf) - karr(nk_elf-1)) / sunglint%s(j)%gamma_p

        DO i = 1, nk_elf
          k = karr(i)
          sqrt_k_kp = SQRT(k * kp_r)
          c_r = 1._jprb / sunglint%c(i,j)

          sk2_ad = 0._jprb
          dk_ad = 0._jprb

          IF (sunglint%sk2(i,j) > 0._jprb) THEN
            ! Angular spreading factor
            sk2_ad = sk2_ad + sk2p_ad(i) * 0.5_jprb * (1 - 0.5_jprb * sunglint%dk(i,j) * cos2wangl)
            dk_ad = dk_ad - sk2p_ad(i) * 0.25_jprb * sunglint%sk2(i,j) * cos2wangl
            cos2wangl_ad = cos2wangl_ad - sk2p_ad(i) * 0.25_jprb * sunglint%sk2(i,j) * sunglint%dk(i,j)

            sk2_ad = sk2_ad + sk2o_ad(i) * 0.5_jprb * (1 + 0.5_jprb * sunglint%dk(i,j) * cos2wangl)
            dk_ad = dk_ad + sk2o_ad(i) * 0.25_jprb * sunglint%sk2(i,j) * cos2wangl
            cos2wangl_ad = cos2wangl_ad + sk2o_ad(i) * 0.25_jprb * sunglint%sk2(i,j) * sunglint%dk(i,j)

            cm_c_pow25 = (cm_elf * c_r)**2.5
            c_cp_pow25 = (sunglint%c(i,j) * cp_r)**2.5
            cosh_tmp = 1._jprb / COSH(a0_elf + ap_elf * c_cp_pow25 + am * cm_c_pow25)**2
            am_ad = am_ad + dk_ad * cm_c_pow25 * cosh_tmp
            cp_ad = cp_ad - dk_ad *  2.5_jprb * ap_elf * c_cp_pow25 * cp_r * cosh_tmp
          ELSE
            sk2o_ad(i) = 0._jprb
            sk2p_ad(i) = 0._jprb
          ENDIF

          ! Combined spectrum
          bl_ad = sk2_ad / k
          bh_ad = sk2_ad / k

          ! Short-wave spectrum
          alpha_m_ad = alpha_m_ad + 0.5_jprb * bh_ad * cm_elf * sunglint%fm(i,j) * c_r

          ! Long-wave spectrum
          fp = sunglint%lpm(i,j) * sunglint%jp(i,j) * sunglint%fpexp(i,j)

          alpha_p_ad = alpha_p_ad + 0.5_jprb * bl_ad * cp * fp * c_r
          cp_ad      = cp_ad + 0.5_jprb * bl_ad * alpha_p * fp * c_r
          fp_ad      = 0.5_jprb * bl_ad * alpha_p * cp * c_r

          lpm_ad   = fp_ad * sunglint%jp(i,j) * sunglint%fpexp(i,j)
          jp_ad    = fp_ad * sunglint%lpm(i,j) * sunglint%fpexp(i,j)
          fpexp_ad = fp_ad * sunglint%lpm(i,j) * sunglint%jp(i,j)

          omega_ad = omega_ad - (sunglint%fpexp(i,j) / SQRT(10._jprb)) * &
                     (sqrt_k_kp - 1) * fpexp_ad
          kp_ad    = kp_ad + (sunglint%fpexp(i,j) / SQRT(10._jprb)) * &
                     0.5_jprb * omega * sqrt_k_kp * fpexp_ad * kp_r

          IF (sunglint%gamma_exp(i,j) >= exp_max_exp_exponent_r) THEN
            IF (sunglint%gamma_exp(i,j) > min_exponent) THEN
              gamma_exp_ad = jp_ad * sunglint%jp(i,j) * LOG(gamma)
              gamma_ad = gamma_ad + jp_ad * sunglint%jp(i,j) * &
                         sunglint%gamma_exp(i,j) / gamma
            ELSE
              jp_ad = 0._jprb
            ENDIF

            kp_ad    = kp_ad + gamma_exp_ad * sunglint%gamma_exp(i,j) * &
                        (sqrt_k_kp - 1) * 0.5_jprb * sqrt_k_kp * kp_r * sigma_r2
            sigma_ad = sigma_ad + gamma_exp_ad * sunglint%gamma_exp(i,j) * &
                        (sqrt_k_kp - 1)**2 * sigma_r3
          ELSE
            jp_ad = 0._jprb
          ENDIF

          kp_ad = kp_ad - 5 * 2 * lpm_ad * kp * sunglint%lpm(i,j) / (4 * k**2)
        ENDDO

        ! Angular spread
        ustar_ad = 0.13_jprb * am_ad / cm_elf

        ! Air friction velocity
        IF (ustar < cm_elf) THEN
          ustar_ad = ustar_ad + 1E-2_jprb * alpha_m_ad / ustar
        ELSE
          ustar_ad = ustar_ad + 1E-2_jprb * 3 * alpha_m_ad / ustar
        ENDIF
        windsp_ad = ustar_ad / 25

        omega_ad = omega_ad + 6E-3_jprb * 0.5_jprb * alpha_p_ad / SQRT(omega)

        windsp_ad = windsp_ad + omega_ad * cp_r
        cp_ad = cp_ad - windsp * omega_ad * cp_r**2

        kp_ad = kp_ad + 0.5_jprb  * gravity * cp_ad * (1 / km_elf**2 - kp_r**2) * cp_r

        omega_c_ad = -0.08_jprb * 4 * 3 * sigma_ad / sunglint%s(j)%omega_c**4
        IF (sunglint%s(j)%omega_c < 1) THEN
          gamma_ad = 0._jprb
        ELSE
          ! log10(x) = log(x) / log(10)
          omega_c_ad = omega_c_ad + 6 * gamma_ad / (sunglint%s(j)%omega_c * LOG(10._jprb))
        ENDIF

        k0_ad = kp_ad * sunglint%s(j)%omega_c**2
        omega_c_ad = omega_c_ad + 2 * k0 * kp_ad * sunglint%s(j)%omega_c

        x_ad = 0.84_jprb * (-0.75_jprb) * (TANH((x / x0_elf)**0.4))**(-1.75) * &
               (1._jprb / COSH((x / x0_elf)**0.4)**2) * &
               (0.4_jprb * omega_c_ad / (x**0.6_jprb * x0_elf**0.4_jprb))

        k0_ad = k0_ad + x_ad * profiles(j)%s2m%wfetc
        profiles_ad(j)%s2m%wfetc = profiles_ad(j)%s2m%wfetc + x_ad * k0

        windsp_ad = windsp_ad - 2 * gravity * k0_ad / windsp**3

      ENDIF ! wave spectrum option

      wangl_ad = wangl_ad - 2 * cos2wangl_ad * SIN(2 * sunglint%s(j)%wangl)


      ! --------------------------------------------------------------------
      ! Calculate profile variables
      ! --------------------------------------------------------------------

      ! Satellite and solar zenith angles at surface
      ilaysur = aux%s(j)%nearestlev_surf - 1
      IF (ABS(raytracing%pathsun(ilaysur, j) - 1._jprb) > realtol) THEN
        raytracing_ad%pathsun(ilaysur, j) = raytracing_ad%pathsun(ilaysur, j) + zensun_ad / &
            (raytracing%pathsun(ilaysur, j) * SQRT(raytracing%pathsun(ilaysur, j)**2 - 1._jprb))
      ENDIF
      zensun_ad = 0._jprb
      IF (ABS(raytracing%pathsat(ilaysur, j) - 1._jprb) > realtol) THEN
        raytracing_ad%pathsat(ilaysur, j) = raytracing_ad%pathsat(ilaysur, j) + zensat_ad / &
            (raytracing%pathsat(ilaysur, j) * SQRT(raytracing%pathsat(ilaysur, j)**2 - 1._jprb))
      ENDIF
      zensat_ad = 0._jprb

      ! Compute the angle between the wind direction and the U axis

      IF (sunglint%s(j)%windsp > min_windsp) THEN
        profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v +      &
            profiles(j)%s2m%u * wangl_ad / (profiles(j)%s2m%u**2 + profiles(j)%s2m%v**2)
        profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u -      &
            profiles(j)%s2m%v * wangl_ad / (profiles(j)%s2m%u**2 + profiles(j)%s2m%v**2)
      ELSE
        wangl_ad = 0._jprb
      ENDIF

      ! Compute the wind speed

      IF (sunglint%s(j)%windsp > min_windsp) THEN
        profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u + windsp_ad * profiles(j)%s2m%u / sunglint%s(j)%windsp
        profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v + windsp_ad * profiles(j)%s2m%v / sunglint%s(j)%windsp
        windsp_ad            = 0._jprb
      ELSE
        profiles_ad(j)%s2m%u = profiles_ad(j)%s2m%u + windsp_ad / SQRT(2._jprb)
        profiles_ad(j)%s2m%v = profiles_ad(j)%s2m%v + windsp_ad / SQRT(2._jprb)
        windsp_ad            = 0._jprb
      ENDIF

    ENDIF ! wind speed > minimum
  ENDDO ! nprofiles
  sunglint_ad%s(:)%glint = 0._jprb
  sunglint_ad%s(:)%omega = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_refsun_ad
