! Description:
!> @file
!!   K of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!
!> @brief
!!   K of calculation of the fraction of solar radiance that is
!!   reflected by a wind roughened water surface.
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation
!! @param[in]     solar          flag to indicate channels with solar contribution
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in,out] profiles_k     Jacobian wrt atmospheric profile and surface variables
!! @param[in]     coef           optical depth coefficient structure
!! @param[in]     aux            internal structure containing auxiliary profile variables
!! @param[in]     sunglint       internal structure for sea surface BRDF model variables
!! @param[in,out] sunglint_k     Jacobian wrt sea surface BRDF model variables
!! @param[in]     raytracing     raytracing structure
!! @param[in,out] raytracing_k   Jacobian wrt raytracing variables
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
SUBROUTINE rttov_refsun_k( &
              opts,         &
              chanprof,     &
              calcrefl,     &
              solar,        &
              profiles,     &
              profiles_k,   &
              coef,         &
              aux,          &
              sunglint,     &
              sunglint_k,   &
              raytracing,   &
              raytracing_k)

  USE rttov_types, ONLY : &
      rttov_options,     &
      rttov_chanprof,    &
      rttov_profile_aux, &
      rttov_profile,     &
      rttov_raytracing,  &
      rttov_sunglint,    &
      rttov_coef
  USE parkind1, ONLY : jplm
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
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),           INTENT(IN)    :: calcrefl(SIZE(chanprof))
  LOGICAL(jplm),           INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile_aux), INTENT(IN)    :: aux
  TYPE(rttov_raytracing),  INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),  INTENT(INOUT) :: raytracing_k
  TYPE(rttov_sunglint),    INTENT(IN)    :: sunglint
  TYPE(rttov_sunglint),    INTENT(INOUT) :: sunglint_k
  TYPE(rttov_coef),        INTENT(IN)    :: coef
!INTF_END

  REAL(jprb)    :: psi_k(coef%ws_nomega)                            ! The frequency spectrum of the surface wave
  REAL(jprb)    :: ff_k (coef%ws_nomega)                            ! Working space

  INTEGER(jpim) :: j, i, ilaysur, prof, nchannels
  REAL(jprb)    :: csi_k                                            ! Angle between the zenith and the proJection of the
                                                                         ! incident ray on the x-z plane
  REAL(jprb)    :: alfa_k                                           ! Angle between the incident ray and the x-z plane
  REAL(jprb)    :: c_shad_k                                         ! The average magnitude of tan(alfa)
  REAL(jprb)    :: p_prime_k                                        ! The probability density of tan(alfa)
  REAL(jprb)    :: pxy_gammaxy_k                                    ! The joint probability density of the along-view and
                                                                         ! cross view slope
  REAL(jprb)    :: gamma_sq_k                                       ! Total variance of the slope of the facet
  REAL(jprb)    :: gamma_o_k                                        ! The mean square of the along-view (X axis) slope
  REAL(jprb)    :: gamma_p_k                                        ! The mean square of the cross-view (Y axis) slope
  REAL(jprb)    :: g_shad_k                                         ! Normalization function
  REAL(jprb)    :: gammax_k                                         ! The x-slope of the facet
  REAL(jprb)    :: theta_fi                                         ! First order shadowing factor for the reflected ray
  REAL(jprb)    :: theta_csi                                        ! First ordesr shadowing factor for the incident ray
  REAL(jprb)    :: q_shad_k                                         ! Second order shadowing factor
  REAL(jprb)    :: q_shad_a_k                                       ! Second order shadowing factor
  REAL(jprb)    :: q_shad_b_k                                       ! Second order shadowing factor
  REAL(jprb)    :: zensat_k                                         ! Zenith angle of satellite viewing angle at surface
  REAL(jprb)    :: zensun_k                                         ! Zenith angle of sun at surface
  REAL(jprb)    :: omega_1                                          ! Frequency of the surface wave
  REAL(jprb)    :: fac1_k
  REAL(jprb)    :: a_shad_k
  REAL(jprb)    :: b_shad_k
  REAL(jprb)    :: lambda_a_k
  REAL(jprb)    :: lambda_b_k
  REAL(jprb)    :: x_u_k
  REAL(jprb)    :: alfa1_k
  REAL(jprb)    :: k
  REAL(jprb)    :: sigma_a2_r
  REAL(jprb)    :: sigma_b2_r
  REAL(jprb)    :: omega_m_k
  REAL(jprb)    :: sigma2_r
  REAL(jprb)    :: beta_k
  REAL(jprb)    :: windsp_k
  REAL(jprb)    :: wangl_k
  REAL(jprb)    :: s2_k, s4_k, sa_k, sb_k
  REAL(jprb)    :: bb, aa, hinc, sqrt2pi

  REAL(jprb)    :: max_exp_exponent_4rt, alfa1_r, logk
  REAL(jprb)    :: exp_max_exp_exponent_r, omega_m_3, omega_m_3_r
  REAL(jprb)    :: omega_1_array(coef%ws_nomega)
  REAL(jprb)    :: omega_1_4_r(coef%ws_nomega)

  REAL(jprb)    :: karr(0:nk_elf), sk2_k, sk2o_k(0:nk_elf), sk2p_k(0:nk_elf)
  REAL(jprb)    :: sqrt_k_kp, cm_c_pow25, c_cp_pow25
  REAL(jprb)    :: cos2wangl, cos2wangl_k
  REAL(jprb)    :: windsp, k0, x, kp, gamma, sigma
  REAL(jprb)    :: k0_k, x_k, kp_k, gamma_k, sigma_k
  REAL(jprb)    :: cp, omega, alpha_p, ustar, am
  REAL(jprb)    :: cp_k, omega_k, alpha_p_k, ustar_k, alpha_m_k, am_k
  REAL(jprb)    :: fp, fp_k, bl_k, bh_k
  REAL(jprb)    :: omega_c_k, lpm_k, jp_k, gamma_exp_k, fpexp_k, dk_k
  REAL(jprb)    :: cosh_tmp, cp_r, c_r, kp_r, sigma_r2, sigma_r3

  REAL(jprb)    :: ZHOOK_HANDLE
!-----End of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_K', 0_jpim, ZHOOK_HANDLE)
  nchannels = SIZE(chanprof)
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

  DO j = 1, nchannels
    prof = chanprof(j)%prof
    IF (.NOT. (solar(j) .AND. calcrefl(j))) CYCLE
    IF (profiles(prof)%sunzenangle < 0.0 .OR. &
        profiles(prof)%sunzenangle > max_sol_zen) CYCLE
    IF (profiles(prof)%skin%surftype /= surftype_sea) CYCLE

    windsp_k      = 0._jprb
    wangl_k       = 0._jprb
    q_shad_k      = 0._jprb
    q_shad_a_k    = 0._jprb
    q_shad_b_k    = 0._jprb
    g_shad_k      = 0._jprb
    p_prime_k     = 0._jprb
    pxy_gammaxy_k = 0._jprb
    fac1_k        = 0._jprb
    gamma_o_k     = 0._jprb
    gammax_k      = 0._jprb
    alfa_k        = 0._jprb
    csi_k         = 0._jprb
    zensat_k      = 0._jprb
    zensun_k      = 0._jprb
    c_shad_k      = 0._jprb
    gamma_p_k     = 0._jprb
    gamma_sq_k    = 0._jprb
    lambda_a_k    = 0._jprb
    lambda_b_k    = 0._jprb
    b_shad_k      = 0._jprb
    a_shad_k      = 0._jprb

    IF (sunglint%s(prof)%windsp <= min_windsp) THEN

      sunglint_k%s(j)%glint = 0._jprb

    ELSE

      ! ------------------------------------------------------------------
      ! Yoshimori wave facet reflectance model
      ! ------------------------------------------------------------------

      ! Compute effective distribution function

      q_shad_k      = q_shad_k + &
          sunglint_k%s(j)%glint * sunglint%s(prof)%g_shad * sunglint%s(prof)%p_prime * &
          sunglint%s(prof)%pxy_gammaxy / sunglint%s(prof)%fac1
      g_shad_k      = g_shad_k + &
          sunglint_k%s(j)%glint * sunglint%s(prof)%q_shad * sunglint%s(prof)%p_prime * &
          sunglint%s(prof)%pxy_gammaxy / sunglint%s(prof)%fac1
      p_prime_k     = p_prime_k + &
          sunglint_k%s(j)%glint * sunglint%s(prof)%q_shad * sunglint%s(prof)%g_shad * &
          sunglint%s(prof)%pxy_gammaxy / sunglint%s(prof)%fac1
      pxy_gammaxy_k = pxy_gammaxy_k + &
          sunglint_k%s(j)%glint * sunglint%s(prof)%q_shad * sunglint%s(prof)%g_shad * &
          sunglint%s(prof)%p_prime / sunglint%s(prof)%fac1
      fac1_k        = fac1_k - &
          sunglint_k%s(j)%glint * sunglint%s(prof)%q_shad * sunglint%s(prof)%g_shad * &
          sunglint%s(prof)%p_prime * sunglint%s(prof)%pxy_gammaxy / sunglint%s(prof)%fac1**2

      ! Compute the probability density that the slope of the facet at
      ! a certain point is gammax when the incident ray and the ray from
      ! the reflection of the incident ray on the local surface do not
      ! intersect with any other surface.

      gamma_o_k     = gamma_o_k + pxy_gammaxy_k * sunglint%s(prof)%pxy_gammaxy * &
          (sunglint%s(prof)%gammax**2 / sunglint%s(prof)%gamma_o**3 - 1._jprb / sunglint%s(prof)%gamma_o)
      gammax_k      = gammax_k - &
          pxy_gammaxy_k * sunglint%s(prof)%pxy_gammaxy * sunglint%s(prof)%gammax / sunglint%s(prof)%gamma_o**2
      pxy_gammaxy_k = 0._jprb
      alfa_k        = alfa_k - fac1_k * 2._jprb * sunglint%s(prof)%fac1 * TAN(sunglint%s(prof)%alfa)
      csi_k         = csi_k - &
          fac1_k * sunglint%s(prof)%fac1 * TAN((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)
      zensat_k      = zensat_k + &
          fac1_k * sunglint%s(prof)%fac1 * TAN((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb)
      fac1_k        = 0._jprb

      csi_k         = csi_k - g_shad_k * 0.5_jprb * TAN(sunglint%s(prof)%zensat) / &
          COS(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%zensat))**2
      zensat_k      = zensat_k + 0.5_jprb * g_shad_k * TAN(sunglint%s(prof)%zensat) / &
          COS(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%zensat))**2
      zensat_k      = zensat_k - g_shad_k * TAN(0.5_jprb * (sunglint%s(prof)%csi - sunglint%s(prof)%zensat)) / &
          COS(sunglint%s(prof)%zensat)**2
      g_shad_k      = 0._jprb

      alfa_k        = alfa_k - p_prime_k * sunglint%s(prof)%p_prime * TAN(sunglint%s(prof)%alfa) / &
          (COS(sunglint%s(prof)%alfa) * sunglint%s(prof)%c_shad)**2
      c_shad_k      = c_shad_k + p_prime_k * sunglint%s(prof)%p_prime * &
          ((TAN(sunglint%s(prof)%alfa))**2 / sunglint%s(prof)%c_shad**3 - 1._jprb / sunglint%s(prof)%c_shad)
      p_prime_k     = 0._jprb

      IF (COS(sunglint%s(prof)%csi) + COS(sunglint%s(prof)%zensat) >= 0._jprb) THEN
        zensat_k  = zensat_k - c_shad_k * sunglint%s(prof)%gamma_p * SIN(sunglint%s(prof)%zensat)
        csi_k     = csi_k - c_shad_k * sunglint%s(prof)%gamma_p * SIN(sunglint%s(prof)%csi)
        gamma_p_k = gamma_p_k + c_shad_k * (COS(sunglint%s(prof)%csi) + COS(sunglint%s(prof)%zensat))
        c_shad_k  = 0._jprb
      ELSE
        zensat_k  = zensat_k + c_shad_k * sunglint%s(prof)%gamma_p * SIN(sunglint%s(prof)%zensat)
        csi_k     = csi_k + c_shad_k * sunglint%s(prof)%gamma_p * SIN(sunglint%s(prof)%csi)
        gamma_p_k = gamma_p_k - c_shad_k * (COS(sunglint%s(prof)%csi) + COS(sunglint%s(prof)%zensat))
        c_shad_k  = 0._jprb
      ENDIF

      ! Direct model - first order shadowing

      IF (ABS(sunglint%s(prof)%zensat) > 0._jprb) THEN
        IF (1._jprb / TAN(ABS(sunglint%s(prof)%zensat)) - &
            sunglint%s(prof)%gammax * sunglint%s(prof)%zensat / ABS(sunglint%s(prof)%zensat) >= 0._jprb) THEN
          theta_fi = 1._jprb
        ELSE
          theta_fi = 0._jprb
        ENDIF
      ELSE
        theta_fi = 1._jprb
      ENDIF
      IF (ABS(sunglint%s(prof)%csi) > 0._jprb) THEN
        IF (1._jprb / TAN(ABS(sunglint%s(prof)%csi)) + &
            sunglint%s(prof)%gammax * sunglint%s(prof)%csi / ABS(sunglint%s(prof)%csi) >= 0._jprb) THEN
          theta_csi = 1._jprb
        ELSE
          theta_csi = 0._jprb
        ENDIF
      ELSE
        theta_csi = 1._jprb
      ENDIF
      IF (ABS(sunglint%s(prof)%csi) > pi / 2._jprb) theta_csi = 0._jprb

      ! Second order shadowing (the facet cannot be seen)

      IF ((sunglint%s(prof)%zensat * sunglint%s(prof)%csi) <= 0._jprb) THEN
        IF (ABS(sunglint%s(prof)%zensat) < realtol .AND. ABS(sunglint%s(prof)%csi) < realtol) THEN
          q_shad_k = 0._jprb
        ELSE IF (sunglint%s(prof)%zensat >= ABS(sunglint%s(prof)%csi)) THEN
          q_shad_a_k = q_shad_a_k + q_shad_k * theta_fi
          q_shad_k   = 0._jprb
          lambda_a_k = lambda_a_k - q_shad_a_k / (1._jprb + sunglint%s(prof)%lambda_a)**2
          q_shad_a_k = 0._jprb
          a_shad_k   = a_shad_k - lambda_a_k * EXP(-sunglint%s(prof)%a_shad**2 / 2._jprb) / &
              (sqrt2pi * sunglint%s(prof)%a_shad**2)
          lambda_a_k = 0._jprb
          zensat_k   = zensat_k - a_shad_k / (SIN(sunglint%s(prof)%zensat)**2 * sunglint%s(prof)%gamma_o)
          gamma_o_k  = gamma_o_k - a_shad_k / (TAN(sunglint%s(prof)%zensat) * sunglint%s(prof)%gamma_o**2)
          a_shad_k   = 0._jprb
        ELSE IF (ABS(sunglint%s(prof)%csi) > sunglint%s(prof)%zensat) THEN
          q_shad_b_k = q_shad_b_k + q_shad_k * theta_csi
          q_shad_k   = 0._jprb
          lambda_b_k = lambda_b_k - q_shad_b_k / (1._jprb + sunglint%s(prof)%lambda_b)**2
          q_shad_b_k = 0._jprb
          b_shad_k   = b_shad_k - lambda_b_k * EXP(-sunglint%s(prof)%b_shad**2 / 2._jprb) / &
              (sqrt2pi * sunglint%s(prof)%b_shad**2)
          lambda_b_k = 0._jprb
          IF (sunglint%s(prof)%csi >= 0._jprb) THEN
            csi_k     = csi_k - b_shad_k / (SIN(sunglint%s(prof)%csi)**2 * sunglint%s(prof)%gamma_o)
            gamma_o_k = gamma_o_k - b_shad_k / (TAN(sunglint%s(prof)%csi) * sunglint%s(prof)%gamma_o**2)
            b_shad_k  = 0._jprb
          ELSE
            csi_k     = csi_k + b_shad_k  / (SIN(sunglint%s(prof)%csi)**2 * sunglint%s(prof)%gamma_o)
            gamma_o_k = gamma_o_k + b_shad_k / (TAN(sunglint%s(prof)%csi) * sunglint%s(prof)%gamma_o**2)
            b_shad_k  = 0._jprb
          ENDIF
        ENDIF
      ELSE
        q_shad_k   = q_shad_k * theta_fi * theta_csi
        lambda_a_k = lambda_a_k - q_shad_k / (1._jprb + sunglint%s(prof)%lambda_a + sunglint%s(prof)%lambda_b)**2
        lambda_b_k = lambda_b_k - q_shad_k / (1._jprb + sunglint%s(prof)%lambda_a + sunglint%s(prof)%lambda_b)**2
        q_shad_k   = 0._jprb
        b_shad_k   = b_shad_k - lambda_b_k * EXP(-sunglint%s(prof)%b_shad**2 / 2._jprb) / &
            (sqrt2pi * sunglint%s(prof)%b_shad**2)
        lambda_b_k = 0._jprb
        a_shad_k   = a_shad_k - lambda_a_k * EXP(-sunglint%s(prof)%a_shad**2 / 2._jprb) / &
            (sqrt2pi * sunglint%s(prof)%a_shad**2)
        lambda_a_k = 0._jprb
        IF (sunglint%s(prof)%csi >= 0._jprb) THEN
          csi_k     = csi_k - b_shad_k / (SIN(sunglint%s(prof)%csi)**2 * sunglint%s(prof)%gamma_o)
          gamma_o_k = gamma_o_k - b_shad_k / (TAN(sunglint%s(prof)%csi) * sunglint%s(prof)%gamma_o**2)
          b_shad_k  = 0._jprb
        ELSE IF (sunglint%s(prof)%csi < 0._jprb) THEN
          csi_k     = csi_k + b_shad_k / (SIN(sunglint%s(prof)%csi)**2 * sunglint%s(prof)%gamma_o)
          gamma_o_k = gamma_o_k + b_shad_k / (TAN(sunglint%s(prof)%csi) * sunglint%s(prof)%gamma_o**2)
          b_shad_k  = 0._jprb
        ENDIF
        zensat_k  = zensat_k - a_shad_k / (SIN(sunglint%s(prof)%zensat)**2 * sunglint%s(prof)%gamma_o)
        gamma_o_k = gamma_o_k - a_shad_k / (TAN(sunglint%s(prof)%zensat) * sunglint%s(prof)%gamma_o**2)
        a_shad_k  = 0._jprb
      ENDIF

      ! First order shadowing (the slope of the facet is negative)

      csi_k    = csi_k + &
          gammax_k / (2._jprb * (COS((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb))**2)
      zensat_k = zensat_k - &
          gammax_k / (2._jprb * (COS((sunglint%s(prof)%csi - sunglint%s(prof)%zensat) / 2._jprb))**2)
      gammax_k = 0._jprb

      ! Obtain angles csi and alfa

      IF (sunglint%s(prof)%csi + sunglint%s(prof)%zensat >= 0._jprb) THEN
        csi_k    = csi_k + sunglint_k%s(j)%omega / 2._jprb
        zensat_k = zensat_k + sunglint_k%s(j)%omega / 2._jprb
      ELSE
        csi_k    = csi_k - sunglint_k%s(j)%omega / 2._jprb
        zensat_k = zensat_k - sunglint_k%s(j)%omega / 2._jprb
      ENDIF
      IF (ABS((SIN(sunglint%s(prof)%zensun) * SIN(pi - sunglint%s(prof)%dazng * deg2rad))**2 - 1._jprb) > realtol) THEN
        zensun_k = zensun_k + alfa_k * (COS(sunglint%s(prof)%zensun) * SIN(pi - sunglint%s(prof)%dazng * deg2rad) / &
            SQRT(1._jprb - (SIN(sunglint%s(prof)%zensun) * SIN(pi - sunglint%s(prof)%dazng * deg2rad))**2))
      ENDIF
      alfa_k = 0._jprb
      zensun_k   = zensun_k + &
          csi_k * COS(pi - sunglint%s(prof)%dazng * deg2rad) / (COS(sunglint%s(prof)%zensun)**2 * &
           (1._jprb + (TAN(sunglint%s(prof)%zensun) * COS(pi - sunglint%s(prof)%dazng * deg2rad))**2))
      csi_k = 0._jprb


      ! ------------------------------------------------------------------
      ! Wave spectrum parameterisation
      ! ------------------------------------------------------------------

      cos2wangl = COS(2 * sunglint%s(prof)%wangl)
      cos2wangl_k = 0._jprb

      IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN

        ! JONSWAP wave frequency spectrum parameterisation

        ! AD init

        sa_k          = 0._jprb
        sb_k          = 0._jprb
        s4_k          = 0._jprb
        s2_k          = 0._jprb
        ff_k(:)       = 0._jprb
        beta_k        = 0._jprb
        alfa1_k       = 0._jprb
        omega_m_k     = 0._jprb
        psi_k(:)      = 0._jprb
        x_u_k         = 0._jprb

        ! Compute the rms of the slope of the facet along the X (gamma_o)and Y (gamma_p) axis

        gamma_sq_k = gamma_sq_k + 0.125_jprb * (2 - cos2wangl) * gamma_p_k / sunglint%s(prof)%gamma_p
        cos2wangl_k = cos2wangl_k - 0.125_jprb * gamma_p_k * sunglint%s(prof)%gamma_sq / sunglint%s(prof)%gamma_p

        gamma_sq_k = gamma_sq_k + 0.125_jprb * (2 + cos2wangl) * gamma_o_k / sunglint%s(prof)%gamma_o
        cos2wangl_k = cos2wangl_k + 0.125_jprb * gamma_o_k * sunglint%s(prof)%gamma_sq / sunglint%s(prof)%gamma_o

        ! Compute the Simpson's integral

        bb    = 301._jprb
        aa    = 0._jprb
        hinc  = (bb - aa) / coef%ws_nomega
        sa_k       = sa_k + gamma_sq_k * hinc / 3.0_jprb
        sb_k       = sb_k + gamma_sq_k * hinc / 3.0_jprb
        s4_k       = s4_k + gamma_sq_k * hinc / 3._jprb * 4._jprb
        s2_k       = s2_k + gamma_sq_k * hinc / 3._jprb * 2._jprb
        gamma_sq_k = 0._jprb
        ff_k(2::2) = ff_k(2::2) + s2_k
        s2_k = 0._jprb
        ff_k(1::2) = ff_k(1::2) + s4_k
        s4_k = 0._jprb
        ff_k(1)              = ff_k(1) + sa_k
        ff_k(coef%ws_nomega) = ff_k(coef%ws_nomega) + sb_k
        sa_k = 0._jprb
        sb_k = 0._jprb

        ! Compute JONSWAP frequency spectrum of the surface wave

        alfa1_r = 1._jprb / sunglint%s(prof)%alfa1
        omega_m_3 = sunglint%s(prof)%omega_m**3
        omega_m_3_r = 1._jprb / omega_m_3
        DO i = 2, coef%ws_nomega
          omega_1 = omega_1_array(i)

          psi_k(i) = psi_k(i) + ff_k(i) * coef%ws_k_omega(i)**2

          IF (sunglint%s(prof)%omega_m < max_exp_exponent_4rt * omega_1) THEN ! stop exp overflow

            IF (sunglint%beta(i, prof) > min_exponent) beta_k = beta_k + psi_k(i) * sunglint%psi(i, prof) * logk

            omega_m_k = omega_m_k - &
                psi_k(i) * sunglint%psi(i, prof) * 5._jprb * omega_m_3 * omega_1_4_r(i)
            alfa1_k   = alfa1_k + psi_k(i) * sunglint%psi(i, prof) * alfa1_r

            IF (sunglint%beta(i, prof) >= exp_max_exp_exponent_r) THEN
              IF (omega_1 <= sunglint%s(prof)%omega_m) THEN
                sigma2_r = sigma_a2_r
              ELSE
                sigma2_r = sigma_b2_r
              ENDIF

              omega_m_k = omega_m_k + beta_k * sunglint%beta(i, prof) * omega_1 * &
                  (omega_1 - sunglint%s(prof)%omega_m) * sigma2_r * omega_m_3_r
            ENDIF
          ENDIF
          beta_k = 0._jprb
        ENDDO
!         ff_k(:) = 0._jprb
!         psi_k(:) = 0._jprb

        windsp_k                = windsp_k - omega_m_k * sunglint%s(prof)%omega_m / sunglint%s(prof)%windsp
        x_u_k                   = x_u_k - omega_m_k * 0.33_jprb * sunglint%s(prof)%omega_m / sunglint%s(prof)%x_u
        omega_m_k               = 0._jprb
        x_u_k                   = x_u_k - 0.22_jprb * alfa1_k * sunglint%s(prof)%alfa1 / sunglint%s(prof)%x_u
        alfa1_k                 = 0._jprb
        windsp_k                =      &
            windsp_k - x_u_k * 2._jprb * profiles(prof)%s2m%wfetc * gravity / sunglint%s(prof)%windsp**3
        profiles_k(j)%s2m%wfetc = profiles_k(j)%s2m%wfetc + x_u_k * gravity / sunglint%s(prof)%windsp**2
        x_u_k                   = 0._jprb

      ELSEIF (opts%rt_ir%solar_sea_brdf_model == 2) THEN

        ! Elfouhaily et al wave spectrum parameterisation

        ! AD init

        am_k = 0._jprb
        cp_k = 0._jprb
        alpha_p_k = 0._jprb
        alpha_m_k = 0._jprb
        omega_k = 0._jprb
        gamma_k = 0._jprb
        sigma_k = 0._jprb
        kp_k = 0._jprb

        ! Direct model

        windsp = MAX(sunglint%s(prof)%windsp, 0.6_jprb)

        ! Active values independent of wavenumber (k)
        k0 = gravity / windsp**2
        x = k0 * profiles(prof)%s2m%wfetc
        kp = k0 * sunglint%s(prof)%omega_c**2
        kp_r = 1._jprb / kp

        IF (sunglint%s(prof)%omega_c < 1) THEN
          gamma = 1.7_jprb
        ELSE
          gamma = 1.7_jprb + 6 * LOG10(sunglint%s(prof)%omega_c)
        ENDIF
        sigma = 0.08_jprb * (1 + 4 / sunglint%s(prof)%omega_c**3)
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
        sk2o_k(1:nk_elf-1) = 0.25_jprb * gamma_o_k * (karr(2:nk_elf) - karr(0:nk_elf-2)) / sunglint%s(prof)%gamma_o
        sk2o_k(0)          = 0.25_jprb * gamma_o_k * (karr(1) - karr(0)) / sunglint%s(prof)%gamma_o
        sk2o_k(nk_elf)     = 0.25_jprb * gamma_o_k * (karr(nk_elf) - karr(nk_elf-1)) / sunglint%s(prof)%gamma_o

        sk2p_k(1:nk_elf-1) = 0.25_jprb * gamma_p_k * (karr(2:nk_elf) - karr(0:nk_elf-2)) / sunglint%s(prof)%gamma_p
        sk2p_k(0)          = 0.25_jprb * gamma_p_k * (karr(1) - karr(0)) / sunglint%s(prof)%gamma_p
        sk2p_k(nk_elf)     = 0.25_jprb * gamma_p_k * (karr(nk_elf) - karr(nk_elf-1)) / sunglint%s(prof)%gamma_p

        DO i = 1, nk_elf
          k = karr(i)
          sqrt_k_kp = SQRT(k * kp_r)
          c_r = 1._jprb / sunglint%c(i,prof)

          sk2_k = 0._jprb
          dk_k = 0._jprb

          IF (sunglint%sk2(i,prof) > 0._jprb) THEN
            ! Angular spreading factor
            sk2_k = sk2_k + sk2p_k(i) * 0.5_jprb * (1 - 0.5_jprb * sunglint%dk(i,prof) * cos2wangl)
            dk_k = dk_k - sk2p_k(i) * 0.25_jprb * sunglint%sk2(i,prof) * cos2wangl
            cos2wangl_k = cos2wangl_k - sk2p_k(i) * 0.25_jprb * sunglint%sk2(i,prof) * sunglint%dk(i,prof)

            sk2_k = sk2_k + sk2o_k(i) * 0.5_jprb * (1 + 0.5_jprb * sunglint%dk(i,prof) * cos2wangl)
            dk_k = dk_k + sk2o_k(i) * 0.25_jprb * sunglint%sk2(i,prof) * cos2wangl
            cos2wangl_k = cos2wangl_k + sk2o_k(i) * 0.25_jprb * sunglint%sk2(i,prof) * sunglint%dk(i,prof)

            cm_c_pow25 = (cm_elf * c_r)**2.5
            c_cp_pow25 = (sunglint%c(i,prof) * cp_r)**2.5
            cosh_tmp = 1._jprb / COSH(a0_elf + ap_elf * c_cp_pow25 + am * cm_c_pow25)**2
            am_k = am_k + dk_k * cm_c_pow25 * cosh_tmp
            cp_k = cp_k - dk_k *  2.5_jprb * ap_elf * c_cp_pow25 * cp_r * cosh_tmp
          ELSE
            sk2o_k(i) = 0._jprb
            sk2p_k(i) = 0._jprb
          ENDIF

          ! Combined spectrum
          bl_k = sk2_k / k
          bh_k = sk2_k / k

          ! Short-wave spectrum
          alpha_m_k = alpha_m_k + 0.5_jprb * bh_k * cm_elf * sunglint%fm(i,prof) * c_r

          ! Long-wave spectrum
          fp = sunglint%lpm(i,prof) * sunglint%jp(i,prof) * sunglint%fpexp(i,prof)

          alpha_p_k = alpha_p_k + 0.5_jprb * bl_k * cp * fp * c_r
          cp_k      = cp_k + 0.5_jprb * bl_k * alpha_p * fp * c_r
          fp_k      = 0.5_jprb * bl_k * alpha_p * cp * c_r

          lpm_k   = fp_k * sunglint%jp(i,prof) * sunglint%fpexp(i,prof)
          jp_k    = fp_k * sunglint%lpm(i,prof) * sunglint%fpexp(i,prof)
          fpexp_k = fp_k * sunglint%lpm(i,prof) * sunglint%jp(i,prof)

          omega_k = omega_k - (sunglint%fpexp(i,prof) / SQRT(10._jprb)) * &
                     (sqrt_k_kp - 1) * fpexp_k
          kp_k    = kp_k + (sunglint%fpexp(i,prof) / SQRT(10._jprb)) * &
                     0.5_jprb * omega * sqrt_k_kp * fpexp_k * kp_r

          IF (sunglint%gamma_exp(i,prof) >= exp_max_exp_exponent_r) THEN
            IF (sunglint%gamma_exp(i,prof) > min_exponent) THEN
              gamma_exp_k = jp_k * sunglint%jp(i,prof) * LOG(gamma)
              gamma_k = gamma_k + jp_k * sunglint%jp(i,prof) * &
                         sunglint%gamma_exp(i,prof) / gamma
            ELSE
              jp_k = 0._jprb
            ENDIF

            kp_k    = kp_k + gamma_exp_k * sunglint%gamma_exp(i,prof) * &
                        (sqrt_k_kp - 1) * 0.5_jprb * sqrt_k_kp * kp_r * sigma_r2
            sigma_k = sigma_k + gamma_exp_k * sunglint%gamma_exp(i,prof) * &
                        (sqrt_k_kp - 1)**2 * sigma_r3
          ELSE
            jp_k = 0._jprb
          ENDIF

          kp_k = kp_k - 5 * 2 * lpm_k * kp * sunglint%lpm(i,prof) / (4 * k**2)
        ENDDO

        ! Angular spread
        ustar_k = 0.13_jprb * am_k / cm_elf

        ! Air friction velocity
        IF (ustar < cm_elf) THEN
          ustar_k = ustar_k + 1E-2_jprb * alpha_m_k / ustar
        ELSE
          ustar_k = ustar_k + 1E-2_jprb * 3 * alpha_m_k / ustar
        ENDIF
        windsp_k = ustar_k / 25

        omega_k = omega_k + 6E-3_jprb * 0.5_jprb * alpha_p_k / SQRT(omega)

        windsp_k = windsp_k + omega_k * cp_r
        cp_k = cp_k - windsp * omega_k * cp_r**2

        kp_k = kp_k + 0.5_jprb  * gravity * cp_k * (1 / km_elf**2 - kp_r**2) * cp_r

        omega_c_k = -0.08_jprb * 4 * 3 * sigma_k / sunglint%s(prof)%omega_c**4
        IF (sunglint%s(prof)%omega_c < 1) THEN
          gamma_k = 0._jprb
        ELSE
          ! log10(x) = log(x) / log(10)
          omega_c_k = omega_c_k + 6 * gamma_k / (sunglint%s(prof)%omega_c * LOG(10._jprb))
        ENDIF

        k0_k = kp_k * sunglint%s(prof)%omega_c**2
        omega_c_k = omega_c_k + 2 * k0 * kp_k * sunglint%s(prof)%omega_c

        x_k = 0.84_jprb * (-0.75_jprb) * (TANH((x / x0_elf)**0.4))**(-1.75) * &
               (1._jprb / COSH((x / x0_elf)**0.4)**2) * &
               (0.4_jprb * omega_c_k / (x**0.6_jprb * x0_elf**0.4_jprb))

        k0_k = k0_k + x_k * profiles(prof)%s2m%wfetc
        profiles_k(j)%s2m%wfetc = profiles_k(j)%s2m%wfetc + x_k * k0

        windsp_k = windsp_k - 2 * gravity * k0_k / windsp**3

      ENDIF ! wave spectrum option

      wangl_k = wangl_k - 2 * cos2wangl_k * SIN(2 * sunglint%s(prof)%wangl)


      ! --------------------------------------------------------------------
      ! Calculate profile variables
      ! --------------------------------------------------------------------

      ! Satellite and solar zenith angles at surface
      ilaysur = aux%s(prof)%nearestlev_surf - 1
      IF (ABS(raytracing%pathsun(ilaysur, prof) - 1._jprb) > realtol) THEN
        raytracing_k%pathsun(ilaysur, j) = raytracing_k%pathsun(ilaysur, j) + zensun_k / &
            (raytracing%pathsun(ilaysur, prof) * SQRT(raytracing%pathsun(ilaysur, prof)**2 - 1._jprb))
      ENDIF
      zensun_k = 0._jprb
      IF (ABS(raytracing%pathsat(ilaysur, prof) - 1._jprb) > realtol) THEN
        raytracing_k%pathsat(ilaysur, j) = raytracing_k%pathsat(ilaysur, j) + zensat_k / &
            (raytracing%pathsat(ilaysur, prof) * SQRT(raytracing%pathsat(ilaysur, prof)**2 - 1._jprb))
      ENDIF
      zensat_k = 0._jprb

      ! Compute the angle between the wind direction and the U axis

      IF (sunglint%s(prof)%windsp > min_windsp) THEN
        profiles_k(j)%s2m%v = profiles_k(j)%s2m%v +      &
            profiles(prof)%s2m%u * wangl_k / (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2)
        profiles_k(j)%s2m%u = profiles_k(j)%s2m%u -      &
            profiles(prof)%s2m%v * wangl_k / (profiles(prof)%s2m%u**2 + profiles(prof)%s2m%v**2)
      ELSE
        wangl_k = 0._jprb
      ENDIF

      ! Compute the wind speed

      IF (sunglint%s(prof)%windsp > min_windsp) THEN
        profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + windsp_k * profiles(prof)%s2m%u / sunglint%s(prof)%windsp
        profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + windsp_k * profiles(prof)%s2m%v / sunglint%s(prof)%windsp
        windsp_k            = 0._jprb
      ELSE
        profiles_k(j)%s2m%u = profiles_k(j)%s2m%u + windsp_k / SQRT(2._jprb)
        profiles_k(j)%s2m%v = profiles_k(j)%s2m%v + windsp_k / SQRT(2._jprb)
        windsp_k            = 0._jprb
      ENDIF

    ENDIF ! wind speed > minimum
  ENDDO ! nchannels
  sunglint_k%s(:)%glint = 0._jprb
  sunglint_k%s(:)%omega = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_REFSUN_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_refsun_k
