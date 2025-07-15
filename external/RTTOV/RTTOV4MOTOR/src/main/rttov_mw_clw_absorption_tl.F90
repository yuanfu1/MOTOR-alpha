! Description:
!> @file
!!   TL of MW CLW absorption optical depth calculation.
!!
!> @brief
!!   TL of MW CLW absorption optical depth calculation.
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     coef            rttov_coef structure
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     raytracing      raytracing structure
!! @param[in]     raytracing_tl   raytracing perturbations
!! @param[in]     profiles        input atmospheric profiles on user levels
!! @param[in]     profiles_tl     atmospheric profile perturbations
!! @param[in]     aux             profile_aux structure
!! @param[in,out] aux_tl          perturbations in profile_aux structure
!! @param[in,out] opdp_path_tl    optical depth profile perturbations
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_mw_clw_absorption_tl( &
            opts,          &
            coef,          &
            chanprof,      &
            raytracing,    &
            raytracing_tl, &
            profiles,      &
            profiles_tl,   &
            aux,           &
            aux_tl,        &
            opdp_path_tl)

  USE rttov_types, ONLY : &
      rttov_options,     &
      rttov_coef,        &
      rttov_chanprof,    &
      rttov_raytracing,  &
      rttov_profile,     &
      rttov_profile_aux, &
      rttov_opdp_path
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : &
      pi, gravity, dcoeff,                &
      mw_clw_scheme_liebe,                &
      mw_clw_scheme_rosenkranz,           &
      mw_clw_scheme_tkc, tkc_t_c,         &
      tkc_b_1, tkc_d_1, tkc_b_2, tkc_d_2, &
      tkc_s_1, tkc_s_2, tkc_s_3
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),      INTENT(IN)     :: opts
  TYPE(rttov_coef),         INTENT(IN)     :: coef
  TYPE(rttov_chanprof),     INTENT(IN)     :: chanprof(:)
  TYPE(rttov_raytracing),   INTENT(IN)     :: raytracing
  TYPE(rttov_raytracing),   INTENT(IN)     :: raytracing_tl
  TYPE(rttov_profile),      INTENT(IN)     :: profiles(:)
  TYPE(rttov_profile),      INTENT(IN)     :: profiles_tl(:)
  TYPE(rttov_profile_aux),  INTENT(IN)     :: aux
  TYPE(rttov_profile_aux),  INTENT(INOUT)  :: aux_tl
  TYPE(rttov_opdp_path),    INTENT(INOUT)  :: opdp_path_tl
!INTF_END

  INTEGER(jpim) :: prof, chan, lay, lev, j
  INTEGER(jpim) :: nlayers, nlevels, nprofiles, nchanprof
  REAL(jprb)    :: ztemp
  REAL(jprb)    :: tc, theta
  REAL(jprb)    :: tc_tl, theta_tl
  REAL(jprb)    :: od_tl(profiles(1)%nlayers)

  ! Liebe variables
  REAL(jprb)    :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL(jprb)    :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL(jprb)    :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL(jprb)    :: z1f_sq_z1_sq, z2f_sq_z2_sq
  REAL(jprb)    :: z34_dif_tl, z45_dif_tl, z1_sq_tl, z2_sq_tl, z1_div_tl, z2_div_tl
  REAL(jprb)    :: z1_den_tl, z2_den_tl, zastar_tl, z1_prod_tl, z2_prod_tl, z3_prod_tl, z4_prod_tl
  REAL(jprb)    :: zbstar_tl, zbstar_sq_tl, za2star_tl, za2star_sq_tl, zdiv_tl, zgstar_tl
  REAL(jprb)    :: z1f_sq_z1_sq_tl, z2f_sq_z2_sq_tl
  ! Rosenkranz variables
  REAL(jprb)    :: nu, abs_eps, zabs_eps, abs_eps_tl, zabs_eps_tl
  COMPLEX(jprb) :: inu, eps, fr, fb, eps_tl, fr_tl, fb_tl
  COMPLEX(jprb) :: log1, log2, log1_tl, log2_tl
  REAL(jprb)    :: log_abs_z1_m_inu_sq, log_abs_z1star_m_inu_sq
  REAL(jprb)    :: log_abs_z1_m_inu_sq_tl, log_abs_z1star_m_inu_sq_tl
  ! TKC variables
  REAL(jprb)    :: freq, perm_re, perm_im
  REAL(jprb)    :: term1_1, term1_2, term2_1, term2_2
  REAL(jprb)    :: perm_re_tl, perm_im_tl
  REAL(jprb)    :: term1_1_tl, term1_2_tl, term2_1_tl, term2_2_tl

  nchanprof = SIZE(chanprof)
  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  ! -----------------------------------------------------------
  ! Calculate layer CLW amounts (g.m^-3.km)
  ! -----------------------------------------------------------

  ztemp = 1.0_jprb / (4.3429_jprb * gravity)
  DO prof = 1, nprofiles
    DO lay = 1, nlayers
      lev = lay + 1
      aux_tl%clw(lay,prof) = &
          0.1820_jprb * 100.0_jprb * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * ztemp * 0.5_jprb * &
          (raytracing%pathsat(lay,prof) * (profiles_tl(prof)%clw(lev-1) + profiles_tl(prof)%clw(lev)) + &
           raytracing_tl%pathsat(lay,prof) * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev)))

      IF (opts%interpolation%lgradp) THEN
        aux_tl%clw(lay,prof) = aux_tl%clw(lay,prof) + &
          0.1820_jprb * 100.0_jprb * (profiles_tl(prof)%p(lev) - profiles_tl(prof)%p(lev-1)) * ztemp * &
          raytracing%pathsat(lay,prof) * 0.5_jprb * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev))
      ENDIF
    ENDDO
  ENDDO

  ! -----------------------------------------------------------
  ! Calculate per-profile permittivity-related variables
  ! -----------------------------------------------------------

  DO prof = 1, nprofiles
    IF (opts%rt_mw%clw_scheme == mw_clw_scheme_liebe) THEN

      ! Liebe 1989

      ! Debye terms for permittivity calculations:
      !   calculate individual debye terms for temperature at each level
      DO lev = 1, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        theta = 300.0_jprb / profiles(prof)%t(lev) - 1.0_jprb
        theta_tl = -300.0_jprb * profiles_tl(prof)%t(lev) / profiles(prof)%t(lev) ** 2
        aux_tl%debye_prof(1,lev,prof) = -dcoeff(2) * theta_tl + 2 * dcoeff(3) * theta_tl * theta
        aux_tl%debye_prof(2,lev,prof) = dcoeff(4) * aux_tl%debye_prof(1,lev,prof)
        aux_tl%debye_prof(3,lev,prof) = dcoeff(5) * theta_tl
        aux_tl%debye_prof(4,lev,prof) = dcoeff(7) * aux_tl%debye_prof(3,lev,prof)
        aux_tl%debye_prof(5,lev,prof) = 0._jprb
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      ! Rosenkranz 2015

      DO lev = 2, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb
        theta = 300._jprb / profiles(prof)%t(lev)

        tc_tl = profiles_tl(prof)%t(lev)
        theta_tl = -300._jprb * profiles_tl(prof)%t(lev) / profiles(prof)%t(lev)**2

        ! These parameters are directly from the Rosenkranz paper
        aux_tl%ros_eps_s(lay,prof) = (-43.7527_jprb * 0.05_jprb * theta**(-0.95_jprb) + &
                                       299.504_jprb * 1.47_jprb * theta**0.47_jprb - &
                                       399.364_jprb * 2.11_jprb * theta**1.11_jprb + &
                                       221.327_jprb * 2.31_jprb * theta**1.31_jprb) * theta_tl

        aux_tl%ros_dr(lay,prof) = -tc_tl * aux%ros_dr(lay,prof) / 226.45_jprb
        aux_tl%ros_gammar(lay,prof) = 651.4728_jprb * tc_tl * aux%ros_gammar(lay,prof) / (tc + 133.07_jprb)**2

        aux_tl%ros_db(lay,prof) = -tc_tl * aux%ros_db(lay,prof) / 103.05_jprb
        aux_tl%ros_nu1(lay,prof) = (0.1454962_jprb + 2._jprb * 0.063267156_jprb * tc + &
                                     3._jprb * 0.00093786645_jprb * tc**2) * tc_tl
        aux_tl%ros_z1(lay,prof) = CMPLX(-0.75_jprb * aux_tl%ros_nu1(lay,prof), &
                                         aux_tl%ros_nu1(lay,prof), KIND=jprb)

        ! Precalculated per-profile values for efficiency
        aux_tl%ros_log_abs_z1_sq(lay,prof) = 2._jprb * aux_tl%ros_nu1(lay,prof) / aux%ros_nu1(lay,prof)
        aux_tl%ros_log_z1star_z2(lay,prof) = CONJG(aux_tl%ros_z1(lay,prof)) / CONJG(aux%ros_z1(lay,prof))
        aux_tl%ros_div1(lay,prof) = &
          -(aux_tl%ros_log_z1star_z2(lay,prof) - aux_tl%ros_log_abs_z1_sq(lay,prof)) * aux%ros_div1(lay,prof)**2
        aux_tl%ros_div2(lay,prof) = &
          -(CONJG(aux_tl%ros_log_z1star_z2(lay,prof)) - aux_tl%ros_log_abs_z1_sq(lay,prof)) * aux%ros_div2(lay,prof)**2
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_tkc) THEN

      ! Turner, Kneifel, Cadeddu (2016)

      DO lev = 2, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb
        tc_tl = profiles_tl(prof)%t(lev)

        ! Static dielectric permittivity
        aux_tl%tkc_eps_s(lay,prof) = tc_tl * (tkc_s_1 + 2 * tkc_s_2 * tc + 3 * tkc_s_3 * tc**2)

        aux_tl%tkc_delta_1(lay,prof) = -tkc_b_1 * tc_tl * aux%tkc_delta_1(lay,prof)
        aux_tl%tkc_delta_2(lay,prof) = -tkc_b_2 * tc_tl * aux%tkc_delta_2(lay,prof)

        aux_tl%tkc_tau_1(lay,prof) = -aux%tkc_tau_1(lay,prof) * tkc_d_1 * tc_tl / (tc + tkc_t_c)**2
        aux_tl%tkc_tau_2(lay,prof) = -aux%tkc_tau_2(lay,prof) * tkc_d_2 * tc_tl / (tc + tkc_t_c)**2
      ENDDO

    ENDIF
  ENDDO

  ! -----------------------------------------------------------
  ! Calculate the optical depths
  ! -----------------------------------------------------------

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    od_tl = 0._jprb
    IF (opts%rt_mw%clw_scheme == mw_clw_scheme_liebe) THEN

      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        zf = coef%frequency_ghz(chan)
        ! Direct
        zf_sq        = zf * zf
        z1_sq        = aux%debye_prof(1,lev,prof) * aux%debye_prof(1,lev,prof)
        z2_sq        = aux%debye_prof(2,lev,prof) * aux%debye_prof(2,lev,prof)
        z34_dif      = aux%debye_prof(3,lev,prof) - aux%debye_prof(4,lev,prof)
        z45_dif      = aux%debye_prof(4,lev,prof) - aux%debye_prof(5,lev,prof)
        z1f_sq_z1_sq = zf_sq + z1_sq
        z2f_sq_z2_sq = zf_sq + z2_sq
        z1_div       = 1.0_jprb / z1f_sq_z1_sq
        z2_div       = 1.0_jprb / z2f_sq_z2_sq
        z1_den       = z34_dif * z1_div
        z2_den       = z45_dif * z2_div
        zastar       = aux%debye_prof(3,lev,prof) - zf_sq * (z1_den + z2_den)
        z1_prod      = z34_dif * aux%debye_prof(1,lev,prof)
        z2_prod      = z1_prod * z1_div
        z3_prod      = z45_dif * aux%debye_prof(2,lev,prof)
        z4_prod      = z3_prod * z2_div
        zbstar       =  - zf * (z2_prod + z4_prod)
        zbstar_sq    = zbstar * zbstar
        za2star      = zastar + 2.0_jprb
        za2star_sq   = za2star * za2star
        zdiv         = za2star_sq + zbstar_sq
        zgstar       =  - 3.0_jprb * zbstar / zdiv
        ! Tangent-linear
        !zf_tl    = 0
        !zf_sq_tl = 0
        z1_sq_tl        = 2.0_jprb * aux%debye_prof(1,lev,prof) * aux_tl%debye_prof(1,lev,prof)
        z2_sq_tl        = 2.0_jprb * aux%debye_prof(2,lev,prof) * aux_tl%debye_prof(2,lev,prof)
        z34_dif_tl      = aux_tl%debye_prof(3,lev,prof) - aux_tl%debye_prof(4,lev,prof)
        z45_dif_tl      = aux_tl%debye_prof(4,lev,prof) - aux_tl%debye_prof(5,lev,prof)
        z1f_sq_z1_sq_tl = z1_sq_tl
        z2f_sq_z2_sq_tl = z2_sq_tl
        z1_div_tl       =  - z1f_sq_z1_sq_tl / (z1f_sq_z1_sq * z1f_sq_z1_sq)
        z2_div_tl       =  - z2f_sq_z2_sq_tl / (z2f_sq_z2_sq * z2f_sq_z2_sq)
        z1_den_tl       = z34_dif * z1_div_tl + z34_dif_tl * z1_div
        z2_den_tl       = z45_dif * z2_div_tl + z45_dif_tl * z2_div
        zastar_tl       = aux_tl%debye_prof(3,lev,prof) - zf_sq * (z1_den_tl + z2_den_tl)
        z1_prod_tl      = z34_dif_tl * aux%debye_prof(1,lev,prof) + z34_dif * aux_tl%debye_prof(1,lev,prof)
        z2_prod_tl      = z1_prod_tl * z1_div + z1_prod * z1_div_tl
        z3_prod_tl      = z45_dif_tl * aux%debye_prof(2,lev,prof) + z45_dif * aux_tl%debye_prof(2,lev,prof)
        z4_prod_tl      = z3_prod_tl * z2_div + z3_prod * z2_div_tl
        zbstar_tl       =  - zf * (z2_prod_tl + z4_prod_tl)
        zbstar_sq_tl    = 2.0_jprb * zbstar * zbstar_tl
        za2star_tl      = zastar_tl
        za2star_sq_tl   = 2.0_jprb * za2star * za2star_tl
        zdiv_tl         = za2star_sq_tl + zbstar_sq_tl
        zgstar_tl       = -3.0_jprb * (zbstar_tl * zdiv - zbstar * zdiv_tl) / (zdiv * zdiv)

        od_tl(lay)      = -1.5_jprb * zf * (zgstar_tl * aux%clw(lay,prof) + zgstar * aux_tl%clw(lay,prof))
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      nu = coef%frequency_ghz(chan)
      inu = CMPLX(0._jprb, nu, KIND=jprb)
      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE

        ! NB nu, inu and ros_z2 are constants

        fr = aux%ros_gammar(lay,prof) / (aux%ros_gammar(lay,prof) + inu)
        fr_tl = aux_tl%ros_gammar(lay,prof) * inu / (aux%ros_gammar(lay,prof) + inu)**2

        log_abs_z1_m_inu_sq     = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) - nu)**2)
        log_abs_z1star_m_inu_sq = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) + nu)**2)

        log_abs_z1_m_inu_sq_tl     = 2._jprb * aux_tl%ros_nu1(lay,prof) * &
                                     (0.5625_jprb * aux%ros_nu1(lay,prof) + (aux%ros_nu1(lay,prof) - nu)) / &
                                     (0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) - nu)**2)
        log_abs_z1star_m_inu_sq_tl = 2._jprb * aux_tl%ros_nu1(lay,prof) * &
                                     (0.5625_jprb * aux%ros_nu1(lay,prof) + (aux%ros_nu1(lay,prof) + nu)) / &
                                     (0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) + nu)**2)

        log1 = LOG((aux%ros_z2 - inu) * (CONJG(aux%ros_z1(lay,prof)) + inu))
        log1_tl = CONJG(aux_tl%ros_z1(lay,prof)) / (CONJG(aux%ros_z1(lay,prof)) + inu)
        log2 = LOG((CONJG(aux%ros_z2) - inu) * (aux%ros_z1(lay,prof) + inu))
        log2_tl = aux_tl%ros_z1(lay,prof) / (aux%ros_z1(lay,prof) + inu)

        fb = 0.5_jprb * (log1 - log_abs_z1_m_inu_sq) * aux%ros_div1(lay,prof) + &
             0.5_jprb * (log2 - log_abs_z1star_m_inu_sq) * aux%ros_div2(lay,prof)

        fb_tl = 0.5_jprb * (log1_tl - log_abs_z1_m_inu_sq_tl) * aux%ros_div1(lay,prof) + &
                0.5_jprb * (log1 - log_abs_z1_m_inu_sq) * aux_tl%ros_div1(lay,prof) + &
                0.5_jprb * (log2_tl - log_abs_z1star_m_inu_sq_tl) * aux%ros_div2(lay,prof) + &
                0.5_jprb * (log2 - log_abs_z1star_m_inu_sq) * aux_tl%ros_div2(lay,prof)

        ! Compute the complex relative permittivity
        eps = aux%ros_eps_s(lay,prof) - aux%ros_dr(lay,prof) * (1 - fr) - aux%ros_db(lay,prof) * (1 - fb)

        eps_tl = aux_tl%ros_eps_s(lay,prof) - &
                 aux_tl%ros_dr(lay,prof) * (1 - fr) + &
                 aux%ros_dr(lay,prof) * fr_tl - &
                 aux_tl%ros_db(lay,prof) * (1 - fb) + &
                 aux%ros_db(lay,prof) * fb_tl

        ! Compute optical depth as for Liebe (see notes above)
        abs_eps = REAL(eps + 2)**2 + AIMAG(eps)**2
        abs_eps_tl = 2._jprb * (REAL(eps_tl) * REAL(eps + 2) + AIMAG(eps_tl) * AIMAG(eps))
        zabs_eps = 1._jprb / abs_eps
        zabs_eps_tl = -abs_eps_tl / abs_eps**2
        zgstar = -3._jprb * AIMAG(eps) * zabs_eps
        zgstar_tl = -3._jprb * (AIMAG(eps_tl) * zabs_eps + AIMAG(eps) * zabs_eps_tl)

        od_tl(lay) = -1.5_jprb * nu * (zgstar_tl * aux%clw(lay,prof) + zgstar * aux_tl%clw(lay,prof))
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_tkc) THEN

      ! Turner, Kneifel, Cadeddu (2016)

      freq = coef%frequency_ghz(chan) * 1.E9_jprb

      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE

        ! Direct

        ! Relaxation terms
        term1_1 = aux%tkc_tau_1(lay,prof)**2 * aux%tkc_delta_1(lay,prof) / &
                  (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)
        term1_2 = aux%tkc_tau_2(lay,prof)**2 * aux%tkc_delta_2(lay,prof) / &
                  (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        term2_1 = aux%tkc_tau_1(lay,prof) * aux%tkc_delta_1(lay,prof) / &
                  (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)
        term2_2 = aux%tkc_tau_2(lay,prof) * aux%tkc_delta_2(lay,prof) / &
                  (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        ! Permittivity coefficients
        perm_re = aux%tkc_eps_s(lay,prof) - (2._jprb * pi * freq)**2 * (term1_1 + term1_2)
        perm_im = -2._jprb * pi * freq * (term2_1 + term2_2)

        ! Compute optical depth as for Liebe (see notes above)
        zdiv    = (perm_re + 2)**2 + perm_im**2
        zgstar  = -3.0_jprb * perm_im / zdiv


        ! TL

        ! Relaxation terms
        term1_1_tl = &
          (2 * aux_tl%tkc_tau_1(lay,prof) * aux%tkc_tau_1(lay,prof) * aux%tkc_delta_1(lay,prof) + &
           aux%tkc_tau_1(lay,prof)**2 * aux_tl%tkc_delta_1(lay,prof) - &
           term1_1 * (2._jprb * pi * freq)**2 * 2._jprb * aux_tl%tkc_tau_1(lay,prof) * aux%tkc_tau_1(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)

        term1_2_tl = &
          (2 * aux_tl%tkc_tau_2(lay,prof) * aux%tkc_tau_2(lay,prof) * aux%tkc_delta_2(lay,prof) + &
           aux%tkc_tau_2(lay,prof)**2 * aux_tl%tkc_delta_2(lay,prof) - &
           term1_2 * (2._jprb * pi * freq)**2 * 2._jprb * aux_tl%tkc_tau_2(lay,prof) * aux%tkc_tau_2(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        term2_1_tl = &
          (aux_tl%tkc_tau_1(lay,prof) * aux%tkc_delta_1(lay,prof) + &
           aux%tkc_tau_1(lay,prof) * aux_tl%tkc_delta_1(lay,prof) - &
           term2_1 * (2._jprb * pi * freq)**2 * 2._jprb * aux_tl%tkc_tau_1(lay,prof) * aux%tkc_tau_1(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)

        term2_2_tl = &
          (aux_tl%tkc_tau_2(lay,prof) * aux%tkc_delta_2(lay,prof) + &
           aux%tkc_tau_2(lay,prof) * aux_tl%tkc_delta_2(lay,prof) - &
           term2_2 * (2._jprb * pi * freq)**2 * 2._jprb * aux_tl%tkc_tau_2(lay,prof) * aux%tkc_tau_2(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        ! Permittivity coefficients
        perm_re_tl = aux_tl%tkc_eps_s(lay,prof) - (2._jprb * pi * freq)**2 * (term1_1_tl + term1_2_tl)
        perm_im_tl = -2._jprb * pi * freq * (term2_1_tl + term2_2_tl)

        ! Compute optical depth as for Liebe (see notes above)
        zdiv_tl    = 2 * perm_re_tl * (perm_re + 2) + 2 * perm_im_tl * perm_im
        zgstar_tl  = -3.0_jprb * (perm_im_tl * zdiv - perm_im * zdiv_tl) / zdiv**2
        od_tl(lay) = -1.5_jprb * coef%frequency_ghz(chan) * &
                     (zgstar_tl * aux%clw(lay,prof) + zgstar * aux_tl%clw(lay,prof))
      ENDDO

    ENDIF

    ! Add CLW optical depths to totals
    ztemp = 0.0_jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ! Note that optical depth in the calculations is negative
      ztemp = ztemp + od_tl(lay)
      opdp_path_tl%atm_level(lev,j) = opdp_path_tl%atm_level(lev,j) + ztemp
    ENDDO
  ENDDO ! chanprof

END SUBROUTINE rttov_mw_clw_absorption_tl
