! Description:
!> @file
!!   K of MW CLW absorption optical depth calculation.
!!
!> @brief
!!   K of MW CLW absorption optical depth calculation.
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     coef            rttov_coef structure
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     raytracing      raytracing structure
!! @param[in,out] raytracing_k    raytracing increments
!! @param[in]     profiles        input atmospheric profiles on user levels
!! @param[in,out] profiles_k      atmospheric profile increments
!! @param[in]     aux             profile_aux structure
!! @param[in,out] aux_k           increments in profile_aux structure
!! @param[in]     opdp_path_k     optical depth profile increments
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
SUBROUTINE rttov_mw_clw_absorption_k( &
            opts,         &
            coef,         &
            chanprof,     &
            raytracing,   &
            raytracing_k, &
            profiles,     &
            profiles_k,   &
            aux,          &
            aux_k,        &
            opdp_path_k)

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
  TYPE(rttov_raytracing),   INTENT(INOUT)  :: raytracing_k
  TYPE(rttov_profile),      INTENT(IN)     :: profiles(:)
  TYPE(rttov_profile),      INTENT(INOUT)  :: profiles_k(:)
  TYPE(rttov_profile_aux),  INTENT(IN)     :: aux
  TYPE(rttov_profile_aux),  INTENT(INOUT)  :: aux_k
  TYPE(rttov_opdp_path),    INTENT(IN)     :: opdp_path_k
!INTF_END

  INTEGER(jpim) :: prof, chan, lay, lev, j
  INTEGER(jpim) :: nlayers, nlevels, nchanprof
  REAL(jprb)    :: ztemp
  REAL(jprb)    :: tc, theta
  REAL(jprb)    :: tc_k, theta_k
  REAL(jprb)    :: od_k(profiles(1)%nlayers)

  ! Liebe variables
  REAL(jprb)    :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL(jprb)    :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL(jprb)    :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL(jprb)    :: z1f_sq_z1_sq, z2f_sq_z2_sq
  REAL(jprb)    :: z34_dif_k, z45_dif_k, z1_sq_k, z2_sq_k, z1_div_k, z2_div_k
  REAL(jprb)    :: z1_den_k, z2_den_k, zastar_k, z1_prod_k, z2_prod_k, z3_prod_k, z4_prod_k
  REAL(jprb)    :: zbstar_k, zbstar_sq_k, za2star_k, za2star_sq_k, zdiv_k, zgstar_k
  REAL(jprb)    :: z1f_sq_z1_sq_k, z2f_sq_z2_sq_k
  ! Rosenkranz variables
  REAL(jprb)    :: nu, abs_eps, zabs_eps, abs_eps_k, zabs_eps_k
  REAL(jprb)    :: eps_r_k, eps_i_k
  COMPLEX(jprb) :: inu, eps, fr, fb, eps_k, fr_k, fb_k
  COMPLEX(jprb) :: log1, log2, log1_k, log2_k
  REAL(jprb)    :: log_abs_z1_m_inu_sq, log_abs_z1star_m_inu_sq
  REAL(jprb)    :: log_abs_z1_m_inu_sq_k, log_abs_z1star_m_inu_sq_k
  ! TKC variables
  REAL(jprb)    :: freq, perm_re, perm_im
  REAL(jprb)    :: term1_1, term1_2, term2_1, term2_2
  REAL(jprb)    :: perm_re_k, perm_im_k
  REAL(jprb)    :: term1_1_k, term1_2_k, term2_1_k, term2_2_k

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  ! -----------------------------------------------------------
  ! Calculate the optical depths
  ! -----------------------------------------------------------

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    ! Add CLW optical depths to totals
    ztemp = 0._jprb
    od_k = 0._jprb
    DO lev = nlevels, 2, -1
      lay = lev - 1
      ! Note that optical depth in the calculations is negative
      ztemp = ztemp + opdp_path_k%atm_level(lev,j)
      od_k(lay) = od_k(lay) + ztemp
    ENDDO

    IF (opts%rt_mw%clw_scheme == mw_clw_scheme_liebe) THEN

      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        ! Direct
        zf = coef%frequency_ghz(chan)
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

        ! Adjoint
        zgstar_k                  = od_k(lay) * (-1.5_jprb * zf * aux%clw(lay,prof))
        aux_k%clw(lay,j)          = aux_k%clw(lay,j) + od_k(lay) * (-1.5_jprb * zf * zgstar)
        zbstar_k                  =  - 3.0_jprb * zgstar_k / zdiv
        zdiv_k                    = 3.0_jprb * zgstar_k * zbstar / (zdiv * zdiv)
        za2star_sq_k              = zdiv_k
        zbstar_sq_k               = zdiv_k
        za2star_k                 = 2.0_jprb * za2star * za2star_sq_k
        zastar_k                  = za2star_k
        zbstar_k                  = zbstar_k + 2.0_jprb * zbstar * zbstar_sq_k
        z2_prod_k                 =  - zf * zbstar_k
        z4_prod_k                 =  - zf * zbstar_k
        z3_prod_k                 = z2_div * z4_prod_k
        z2_div_k                  = z3_prod * z4_prod_k
        z45_dif_k                 = aux%debye_prof(2,lev,prof) * z3_prod_k
        aux_k%debye_prof(2,lev,j) = aux_k%debye_prof(2,lev,j) + z45_dif * z3_prod_k
        z1_prod_k                 = z1_div * z2_prod_k
        z1_div_k                  = z1_prod * z2_prod_k
        z34_dif_k                 = aux%debye_prof(1,lev,prof) * z1_prod_k
        aux_k%debye_prof(1,lev,j) = aux_k%debye_prof(1,lev,j) + z34_dif * z1_prod_k
        aux_k%debye_prof(3,lev,j) = aux_k%debye_prof(3,lev,j) + zastar_k
        z1_den_k                  =  - zf_sq * zastar_k
        z2_den_k                  =  - zf_sq * zastar_k
        z2_div_k                  = z2_div_k + z45_dif * z2_den_k
        z45_dif_k                 = z45_dif_k + z2_div * z2_den_k
        z1_div_k                  = z1_div_k + z34_dif * z1_den_k
        z34_dif_k                 = z34_dif_k + z1_div * z1_den_k
        z2f_sq_z2_sq_k            =  - z2_div_k / (z2f_sq_z2_sq * z2f_sq_z2_sq)
        z1f_sq_z1_sq_k            =  - z1_div_k / (z1f_sq_z1_sq * z1f_sq_z1_sq)
        z2_sq_k                   = z2f_sq_z2_sq_k
        z1_sq_k                   = z1f_sq_z1_sq_k
        aux_k%debye_prof(4,lev,j) = aux_k%debye_prof(4,lev,j) + z45_dif_k
        aux_k%debye_prof(5,lev,j) = aux_k%debye_prof(5,lev,j) - z45_dif_k
        aux_k%debye_prof(3,lev,j) = aux_k%debye_prof(3,lev,j) + z34_dif_k
        aux_k%debye_prof(4,lev,j) = aux_k%debye_prof(4,lev,j) - z34_dif_k
        aux_k%debye_prof(2,lev,j) = aux_k%debye_prof(2,lev,j) + z2_sq_k * 2.0_jprb * aux%debye_prof(2,lev,prof)
        aux_k%debye_prof(1,lev,j) = aux_k%debye_prof(1,lev,j) + z1_sq_k * 2.0_jprb * aux%debye_prof(1,lev,prof)
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      ! Rosenkranz
      ! Various per-profile quantities have been precalculated.
      ! Complex arithmetic is avoided in some places, but tests suggested that
      ! replacing it real calculations did not always speed the code up and it
      ! often reduces readability.

      nu = coef%frequency_ghz(chan)
      inu = CMPLX(0._jprb, nu, KIND=jprb)
      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE

        ! Direct
        fr = aux%ros_gammar(lay,prof) / (aux%ros_gammar(lay,prof) + inu)

        log_abs_z1_m_inu_sq     = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) - nu)**2)
        log_abs_z1star_m_inu_sq = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) + nu)**2)

        log1 = LOG((aux%ros_z2 - inu) * (CONJG(aux%ros_z1(lay,prof)) + inu))
        log2 = LOG((CONJG(aux%ros_z2) - inu) * (aux%ros_z1(lay,prof) + inu))

        fb = 0.5_jprb * (log1 - log_abs_z1_m_inu_sq) * aux%ros_div1(lay,prof) + &
             0.5_jprb * (log2 - log_abs_z1star_m_inu_sq) * aux%ros_div2(lay,prof)

        ! Compute the complex relative permittivity
        eps = aux%ros_eps_s(lay,prof) - aux%ros_dr(lay,prof) * (1 - fr) - aux%ros_db(lay,prof) * (1 - fb)

        ! Compute optical depth as for Liebe (see notes above)
        abs_eps = REAL(eps + 2)**2 + AIMAG(eps)**2
        zabs_eps = 1._jprb / abs_eps
        zgstar = -3._jprb * AIMAG(eps) * zabs_eps

        ! Adjoint - NB nu, inu and ros_z2 are constants

        ! Compute optical depth as for Liebe (see notes above)
        zgstar_k = -1.5_jprb * nu * aux%clw(lay,prof) * od_k(lay)
        aux_k%clw(lay,j) = aux_k%clw(lay,j) - 1.5_jprb * nu * zgstar * od_k(lay)

        eps_i_k = -3._jprb * zgstar_k * zabs_eps
        zabs_eps_k = -3._jprb * zgstar_k * AIMAG(eps)
        abs_eps_k = -zabs_eps_k / abs_eps**2
        eps_r_k = 2._jprb * abs_eps_k * REAL(eps + 2)
        eps_i_k = eps_i_k + 2._jprb * abs_eps_k * AIMAG(eps)
        eps_k = CMPLX(eps_r_k, eps_i_k, KIND=jprb)

        ! Compute the complex relative permittivity
        aux_k%ros_eps_s(lay,j) = aux_k%ros_eps_s(lay,j) + eps_r_k
        aux_k%ros_dr(lay,j) = aux_k%ros_dr(lay,j) - REAL(CONJG(1 - fr) * eps_k)
        fr_k = aux%ros_dr(lay,prof) * eps_k
        aux_k%ros_db(lay,j) = aux_k%ros_db(lay,j) - REAL(CONJG(1 - fb) * eps_k)
        fb_k = aux%ros_db(lay,prof) * eps_k

        log1_k = 0.5_jprb * CONJG(aux%ros_div1(lay,prof)) * fb_k
        log_abs_z1_m_inu_sq_k = -0.5_jprb * REAL(CONJG(aux%ros_div1(lay,prof)) * fb_k)
        aux_k%ros_div1(lay,j) = aux_k%ros_div1(lay,j) + &
            0.5_jprb * CONJG(log1 - log_abs_z1_m_inu_sq) * fb_k
        log2_k = 0.5_jprb * CONJG(aux%ros_div2(lay,prof)) * fb_k
        log_abs_z1star_m_inu_sq_k = -0.5_jprb * REAL(CONJG(aux%ros_div2(lay,prof)) * fb_k)
        aux_k%ros_div2(lay,j) = aux_k%ros_div2(lay,j) + &
            0.5_jprb * CONJG(log2 - log_abs_z1star_m_inu_sq) * fb_k

        aux_k%ros_z1(lay,j) = aux_k%ros_z1(lay,j) + &
            log2_k / CONJG(aux%ros_z1(lay,prof) + inu)
        aux_k%ros_z1(lay,j) = aux_k%ros_z1(lay,j) + &
            CONJG(log1_k) / (CONJG(aux%ros_z1(lay,prof)) + inu)

        aux_k%ros_nu1(lay,j) = aux_k%ros_nu1(lay,j) + &
            2._jprb * log_abs_z1star_m_inu_sq_k * &
            (0.5625_jprb * aux%ros_nu1(lay,prof) + (aux%ros_nu1(lay,prof) + nu)) / &
            (0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) + nu)**2)
        aux_k%ros_nu1(lay,j) = aux_k%ros_nu1(lay,j) + &
            2._jprb * log_abs_z1_m_inu_sq_k * &
            (0.5625_jprb * aux%ros_nu1(lay,prof) + (aux%ros_nu1(lay,prof) - nu)) / &
            (0.5625_jprb * aux%ros_nu1(lay,prof)**2 + (aux%ros_nu1(lay,prof) - nu)**2)

        aux_k%ros_gammar(lay,j) = aux_k%ros_gammar(lay,j) + &
            REAL(fr_k * CONJG(inu / (aux%ros_gammar(lay,prof) + inu)**2))
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


        ! Adjoint

        ! Compute optical depth as for Liebe (see notes above)
        zgstar_k = -1.5_jprb * coef%frequency_ghz(chan) * aux%clw(lay,prof) * od_k(lay)
        aux_k%clw(lay,j) = aux_k%clw(lay,j) - 1.5_jprb * coef%frequency_ghz(chan) * zgstar * od_k(lay)

        perm_im_k = -3.0_jprb * zdiv * zgstar_k / zdiv**2
        zdiv_k = 3.0_jprb * perm_im * zgstar_k / zdiv**2

        perm_re_k = 2 * (perm_re + 2) * zdiv_k
        perm_im_k = perm_im_k + 2 * perm_im * zdiv_k

        ! Permittivity coefficients
        term2_1_k = -2._jprb * pi * freq * perm_im_k
        term2_2_k = -2._jprb * pi * freq * perm_im_k
        term1_1_k = -(2._jprb * pi * freq)**2 * perm_re_k
        term1_2_k = -(2._jprb * pi * freq)**2 * perm_re_k
        aux_k%tkc_eps_s(lay,j) = aux_k%tkc_eps_s(lay,j) + perm_re_k

        ! Relaxation terms
        aux_k%tkc_tau_2(lay,j) = aux_k%tkc_tau_2(lay,j) + &
          term2_2_k * (aux%tkc_delta_2(lay,prof) - &
          term2_2 * (2._jprb * pi * freq)**2 * 2._jprb * aux%tkc_tau_2(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)
        aux_k%tkc_delta_2(lay,j) = aux_k%tkc_delta_2(lay,j) + &
          term2_2_k * aux%tkc_tau_2(lay,prof) / (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        aux_k%tkc_tau_1(lay,j) = aux_k%tkc_tau_1(lay,j) + &
          term2_1_k * (aux%tkc_delta_1(lay,prof) - &
          term2_1 * (2._jprb * pi * freq)**2 * 2._jprb * aux%tkc_tau_1(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)
        aux_k%tkc_delta_1(lay,j) = aux_k%tkc_delta_1(lay,j) + &
          term2_1_k * aux%tkc_tau_1(lay,prof) / (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)

        aux_k%tkc_tau_2(lay,j) = aux_k%tkc_tau_2(lay,j) + &
          term1_2_k * (2 * aux%tkc_tau_2(lay,prof) * aux%tkc_delta_2(lay,prof) - &
          term1_2 * (2._jprb * pi * freq)**2 * 2._jprb * aux%tkc_tau_2(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)
        aux_k%tkc_delta_2(lay,j) = aux_k%tkc_delta_2(lay,j) + &
          term1_2_k * aux%tkc_tau_2(lay,prof)**2 / (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_2(lay,prof))**2)

        aux_k%tkc_tau_1(lay,j) = aux_k%tkc_tau_1(lay,j) + &
          term1_1_k * (2 * aux%tkc_tau_1(lay,prof) * aux%tkc_delta_1(lay,prof) - &
          term1_1 * (2._jprb * pi * freq)**2 * 2._jprb * aux%tkc_tau_1(lay,prof)) / &
          (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)
        aux_k%tkc_delta_1(lay,j) = aux_k%tkc_delta_1(lay,j) + &
          term1_1_k * aux%tkc_tau_1(lay,prof)**2 / (1._jprb + (2._jprb * pi * freq * aux%tkc_tau_1(lay,prof))**2)
      ENDDO

    ENDIF
  ENDDO ! chanprof

  ! -----------------------------------------------------------
  ! Calculate per-profile permittivity-related variables
  ! -----------------------------------------------------------

  DO j = 1, nchanprof
    prof = chanprof(j)%prof

    IF (opts%rt_mw%clw_scheme == mw_clw_scheme_liebe) THEN

      ! Liebe 1989

      DO lev = nlevels, 1, -1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        theta = 300.0_jprb / profiles(prof)%t(lev) - 1.0_jprb
        aux_k%debye_prof(3,lev,j) = &
            aux_k%debye_prof(3,lev,j) + aux_k%debye_prof(4,lev,j) * dcoeff(7)
        theta_k = aux_k%debye_prof(3,lev,j) * dcoeff(5)
        aux_k%debye_prof(1,lev,j) = &
            aux_k%debye_prof(1,lev,j) + aux_k%debye_prof(2,lev,j) * dcoeff(4)
        theta_k = theta_k + aux_k%debye_prof(1,lev,j) * (2 * dcoeff(3) * theta - dcoeff(2))
        profiles_k(j)%t(lev) = profiles_k(j)%t(lev) - theta_k * 300.0_jprb / profiles(prof)%t(lev)**2
      ENDDO
      !aux_k%debye_prof(:,:) = 0.

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      ! Rosenkranz 2015

      DO lev = nlevels, 2, -1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb
        theta = 300._jprb / profiles(prof)%t(lev)

        ! Precalculated per-profile values for efficiency
        aux_k%ros_log_z1star_z2(lay,j) = aux_k%ros_log_z1star_z2(lay,j) - &
            CONJG(aux_k%ros_div2(lay,j)) * aux%ros_div2(lay,prof)**2
        aux_k%ros_log_abs_z1_sq(lay,j) = aux_k%ros_log_abs_z1_sq(lay,j) + &
            REAL(aux_k%ros_div2(lay,j) * CONJG(aux%ros_div2(lay,prof)**2))

        aux_k%ros_log_z1star_z2(lay,j) = aux_k%ros_log_z1star_z2(lay,j) - &
            aux_k%ros_div1(lay,j) * CONJG(aux%ros_div1(lay,prof)**2)
        aux_k%ros_log_abs_z1_sq(lay,j) = aux_k%ros_log_abs_z1_sq(lay,j) + &
            REAL(aux_k%ros_div1(lay,j) * CONJG(aux%ros_div1(lay,prof)**2))

        aux_k%ros_z1(lay,j) = aux_k%ros_z1(lay,j) + &
            CONJG(aux_k%ros_log_z1star_z2(lay,j)) / CONJG(aux%ros_z1(lay,prof))
        aux_k%ros_nu1(lay,j) = aux_k%ros_nu1(lay,j) + &
            2._jprb * aux_k%ros_log_abs_z1_sq(lay,j) / aux%ros_nu1(lay,prof)

        ! These parameters are directly from the Rosenkranz paper
        aux_k%ros_nu1(lay,j) = aux_k%ros_nu1(lay,j) - &
            0.75_jprb * REAL(aux_k%ros_z1(lay,j)) + AIMAG(aux_k%ros_z1(lay,j))
        tc_k = (0.1454962_jprb + 2._jprb * 0.063267156_jprb * tc + &
                3._jprb * 0.00093786645_jprb * tc**2) * aux_k%ros_nu1(lay,j)
        tc_k = tc_k - aux_k%ros_db(lay,j) * aux%ros_db(lay,prof) / 103.05_jprb

        tc_k = tc_k + 651.4728_jprb * aux_k%ros_gammar(lay,j) * &
                      aux%ros_gammar(lay,prof) / (tc + 133.07_jprb)**2
        tc_k = tc_k - aux_k%ros_dr(lay,j) * aux%ros_dr(lay,prof) / 226.45_jprb

        theta_k = aux_k%ros_eps_s(lay,j) * &
            (-43.7527_jprb * 0.05_jprb * theta**(-0.95_jprb) + &
              299.504_jprb * 1.47_jprb * theta**0.47_jprb - &
              399.364_jprb * 2.11_jprb * theta**1.11_jprb + &
              221.327_jprb * 2.31_jprb * theta**1.31_jprb)
        profiles_k(j)%t(lev) = profiles_k(j)%t(lev) - &
            300._jprb * theta_k / profiles(prof)%t(lev)**2
        profiles_k(j)%t(lev) = profiles_k(j)%t(lev) + tc_k
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_tkc) THEN

      ! Turner, Kneifel, Cadeddu (2016)

      DO lev = 2, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb

        tc_k = 0._jprb
        tc_k = tc_k - aux_k%tkc_tau_2(lay,j) * aux%tkc_tau_2(lay,prof) * tkc_d_2 / (tc + tkc_t_c)**2
        tc_k = tc_k - aux_k%tkc_tau_1(lay,j) * aux%tkc_tau_1(lay,prof) * tkc_d_1 / (tc + tkc_t_c)**2

        tc_k = tc_k - aux_k%tkc_delta_2(lay,j) * tkc_b_2 * aux%tkc_delta_2(lay,prof)
        tc_k = tc_k - aux_k%tkc_delta_1(lay,j) * tkc_b_1 * aux%tkc_delta_1(lay,prof)

        tc_k = tc_k + aux_k%tkc_eps_s(lay,j) * (tkc_s_1 + 2 * tkc_s_2 * tc + 3 * tkc_s_3 * tc**2)

        profiles_k(j)%t(lev) = profiles_k(j)%t(lev) + tc_k
      ENDDO

    ENDIF
  ENDDO

  ! -----------------------------------------------------------
  ! Calculate layer CLW amounts (g.m^-3.km)
  ! -----------------------------------------------------------

  ztemp = 1.0_jprb / (4.3429_jprb * gravity)
  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    DO lay = 1, nlayers
      lev = lay + 1

      IF (opts%interpolation%lgradp) THEN
        profiles_k(j)%p(lev) = profiles_k(j)%p(lev) + &
          0.1820_jprb * 100.0_jprb * aux_k%clw(lay,j) * ztemp * raytracing%pathsat(lay,prof) * &
          0.5_jprb * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev))
        profiles_k(j)%p(lev-1) = profiles_k(j)%p(lev-1) - &
          0.1820_jprb * 100.0_jprb * aux_k%clw(lay,j) * ztemp * raytracing%pathsat(lay,prof) * &
          0.5_jprb * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev))
      ENDIF

      raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) + &
          0.1820_jprb * 100.0_jprb * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * ztemp * &
          0.5_jprb * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev)) * aux_k%clw(lay,j)
      profiles_k(j)%clw(lev-1) = profiles_k(j)%clw(lev-1) + &
          0.1820_jprb * 100.0_jprb * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * ztemp * &
          raytracing%pathsat(lay,prof) * 0.5_jprb * aux_k%clw(lay,j)
      profiles_k(j)%clw(lev) = profiles_k(j)%clw(lev) + &
          0.1820_jprb * 100.0_jprb * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * ztemp * &
          raytracing%pathsat(lay,prof) * 0.5_jprb * aux_k%clw(lay,j)
    ENDDO
  ENDDO

END SUBROUTINE rttov_mw_clw_absorption_k
