! Description:
!> @file
!!   Calculates optical depths due to CLW absorption for MW channels.
!!
!> @brief
!!   Calculates optical depths due to CLW absorption for MW channels.
!!
!! @details
!!   This implements the MW CLW absorption models. This subroutine
!!   is only called for "clear-sky" MW simulations: it is not used
!!   by RTTOV-SCATT.
!!
!!   The comments in the code give more details about the method.
!!   There are three permittivity parameterisations: Liebe (1989),
!!   Rosenkranz (2015), and Turner, Kneifel, Cadeddu (2016)
!!
!!   NB The output optical depths are NEGATIVE numbers.
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     coef            rttov_coef structure
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     raytracing      raytracing structure
!! @param[in]     profiles        input atmospheric profiles on user levels
!! @param[in,out] aux             profile_aux structure
!! @param[in,out] opdp_path       calculated optical depth profiles
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
SUBROUTINE rttov_mw_clw_absorption( &
            opts,       &
            coef,       &
            chanprof,   &
            raytracing, &
            profiles,   &
            aux,        &
            opdp_path)

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
      tkc_a_1, tkc_b_1, tkc_c_1, tkc_d_1, &
      tkc_a_2, tkc_b_2, tkc_c_2, tkc_d_2, &
      tkc_s_0, tkc_s_1, tkc_s_2, tkc_s_3
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),      INTENT(IN)     :: opts
  TYPE(rttov_coef),         INTENT(IN)     :: coef
  TYPE(rttov_chanprof),     INTENT(IN)     :: chanprof(:)
  TYPE(rttov_raytracing),   INTENT(IN)     :: raytracing
  TYPE(rttov_profile),      INTENT(IN)     :: profiles(:)
  TYPE(rttov_profile_aux),  INTENT(INOUT)  :: aux
  TYPE(rttov_opdp_path),    INTENT(INOUT)  :: opdp_path
!INTF_END

  INTEGER(jpim) :: prof, chan, lay, lev, j
  INTEGER(jpim) :: nlayers, nlevels, nprofiles, nchanprof
  REAL(jprb)    :: ztemp
  REAL(jprb)    :: tc, theta
  REAL(jprb)    :: od(profiles(1)%nlayers)

  ! Liebe variables
  REAL(jprb)    :: zf, zf_sq, z34_dif, z45_dif, z1_sq, z2_sq, z1_div, z2_div
  REAL(jprb)    :: z1_den, z2_den, zastar, z1_prod, z2_prod, z3_prod, z4_prod
  REAL(jprb)    :: zbstar, zbstar_sq, za2star, za2star_sq, zdiv, zgstar
  REAL(jprb)    :: z1f_sq_z1_sq, z2f_sq_z2_sq
  ! Rosenkranz variables
  REAL(jprb)    :: nu, abs_eps, zabs_eps
  COMPLEX(jprb) :: inu, eps, fr, fb, log1, log2
  REAL(jprb)    :: log_abs_z1_m_inu_sq, log_abs_z1star_m_inu_sq
  ! TKC variables
  REAL(jprb)    :: freq, perm_re, perm_im
  REAL(jprb)    :: term1_1, term1_2, term2_1, term2_2

  nchanprof = SIZE(chanprof)
  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  ! -----------------------------------------------------------
  ! Calculate layer CLW amounts (g.m^-3.km)
  ! -----------------------------------------------------------

  ! 0.1820 is constant in relationship between refractivity and
  ! absorption - see Liebe MPM 1989 Eq. 5a.
  ! 4.3429 converts dB to Np (Nepers - also dimensionless).
  ! 100*dp is the depth of the layer in Pascals, g is gravity
  ! dp/g = rho.dz where dz is the depth of the layer in m and rho
  ! is air density in kg.m^-3.
  ! The cloud absorption coefficient is Np.(g.m^-3)^-1.km^-1 so we
  ! must convert CLW input profile from kg.kg-1 to g.m^-3 and
  ! multiply by layer thickness in km.

  ztemp = 1._jprb / (4.3429_jprb * gravity)
  DO prof = 1, nprofiles
    DO lay = 1, nlayers
      lev = lay + 1
      aux%clw(lay,prof) = &
          0.1820_jprb * 100._jprb * (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * ztemp * &
          raytracing%pathsat(lay,prof) * 0.5_jprb * (profiles(prof)%clw(lev-1) + profiles(prof)%clw(lev))
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
      !   The model for the Debeye calculations is a hybrid of LIEBE MPM 1993
      !   and the PIOM laboratory measurements reported by ELLISON et al. 1999.
      DO lev = 1, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        theta = 300.0_jprb / profiles(prof)%t(lev) - 1.0_jprb
        aux%debye_prof(1,lev,prof) = dcoeff(1) - dcoeff(2) * theta + dcoeff(3) * theta * theta
        aux%debye_prof(2,lev,prof) = dcoeff(4) * aux%debye_prof(1,lev,prof)
        aux%debye_prof(3,lev,prof) = dcoeff(5) * theta + dcoeff(6)
        aux%debye_prof(4,lev,prof) = dcoeff(7) * aux%debye_prof(3,lev,prof)
        aux%debye_prof(5,lev,prof) = dcoeff(8)
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      ! Rosenkranz 2015

      aux%ros_z2 = CMPLX(-4500, 2000, KIND=jprb)  ! This is a constant
      DO lev = 2, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb
        theta = 300._jprb / profiles(prof)%t(lev)

        ! These parameters are directly from the Rosenkranz paper
        aux%ros_eps_s(lay,prof) = -43.7527_jprb * theta**0.05_jprb + &
                                   299.504_jprb * theta**1.47_jprb - &
                                   399.364_jprb * theta**2.11_jprb + &
                                   221.327_jprb * theta**2.31_jprb

        aux%ros_dr(lay,prof) = 80.69715_jprb * EXP(-tc / 226.45_jprb)
        aux%ros_gammar(lay,prof) = 1164.023_jprb * EXP(-651.4728_jprb / (tc + 133.07_jprb))

        aux%ros_db(lay,prof) = 4.008724_jprb * EXP(-tc / 103.05_jprb)
        aux%ros_nu1(lay,prof) = 10.46012_jprb + 0.1454962_jprb * tc + &
                                0.063267156_jprb * tc**2 + 0.00093786645_jprb * tc**3
        aux%ros_z1(lay,prof) = CMPLX(-0.75_jprb * aux%ros_nu1(lay,prof), &
                                     aux%ros_nu1(lay,prof), KIND=jprb)

        ! Precalculated per-profile values
        aux%ros_log_abs_z1_sq(lay,prof) = LOG(1.5625_jprb * aux%ros_nu1(lay,prof)**2)   ! LOG(ABS(z1)**2)
        aux%ros_log_z1star_z2(lay,prof) = LOG(CONJG(aux%ros_z1(lay,prof)) * aux%ros_z2) ! LOG(CONJG(z1)*z2)
        aux%ros_div1(lay,prof) = &
          1._jprb  / (aux%ros_log_z1star_z2(lay,prof) - aux%ros_log_abs_z1_sq(lay,prof))
        aux%ros_div2(lay,prof) = &
          1._jprb  / (CONJG(aux%ros_log_z1star_z2(lay,prof)) - aux%ros_log_abs_z1_sq(lay,prof))
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_tkc) THEN

      ! Turner, Kneifel, Cadeddu (2016)

      DO lev = 2, nlevels
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
        lay = lev - 1

        tc = profiles(prof)%t(lev) - 273.15_jprb

        ! Static dielectric permittivity
        aux%tkc_eps_s(lay,prof) = tkc_s_0 + tkc_s_1 * tc + tkc_s_2 * tc**2 + tkc_s_3 * tc**3

        aux%tkc_delta_1(lay,prof) = tkc_a_1 * EXP(-tkc_b_1 * tc)
        aux%tkc_delta_2(lay,prof) = tkc_a_2 * EXP(-tkc_b_2 * tc)

        aux%tkc_tau_1(lay,prof) = tkc_c_1 * EXP(tkc_d_1 / (tc + tkc_t_c))
        aux%tkc_tau_2(lay,prof) = tkc_c_2 * EXP(tkc_d_2 / (tc + tkc_t_c))
      ENDDO

    ENDIF
  ENDDO

  ! -----------------------------------------------------------
  ! Calculate the optical depths
  ! -----------------------------------------------------------

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan

    od = 0._jprb
    IF (opts%rt_mw%clw_scheme == mw_clw_scheme_liebe) THEN

      ! Liebe 1989

      ! Calculate real (zastar) and imaginary (zbstar) parts of complex refractive
      ! index without using complex Fortran as this is computationally expensive.
      !
      ! The refractive index is based on Liebe, Hufford and Manabe (1991) Int. J.
      ! Infrared and Millimeter Waves, 12, 659-671 except that for cloud
      ! with T>273.15 K the coefficients for this model have been recalculated to
      ! fit the measurements made at PIOM laboratory (Lamkaouchi, Balana and Ellison
      ! (1997) New permittivity data for sea water (30-100 GHz), Extension of ESA
      ! report 11197/94/NL/CN. The new model is described in English, Poulsen and
      ! Smith (1999), Proceedings of ATOVS workshop at ECMWF 2-5/9/99.
      !
      ! The optical depth calculation combines equations 5a and 16 from Liebe MPM
      ! 1989 along with unit conversion and path length calculation - see comments
      ! alongside the layer CLW amount calculation above.
      !
      ! The code is written in a step by step format to ease interpretation of
      ! tangent-linear and adjoint code.

      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE
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
        zgstar       = -3.0_jprb * zbstar / zdiv
        od(lay)      = -1.5_jprb * zf * zgstar * aux%clw(lay,prof)
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_rosenkranz) THEN

      ! Rosenkranz 2015

      ! Various per-profile quantities have been precalculated.
      ! Complex arithmetic is avoided in some places, but tests suggested that
      ! replacing it real calculations did not always speed the code up and it
      ! often reduces readability.

      nu = coef%frequency_ghz(chan)
      inu = CMPLX(0._jprb, nu, KIND=jprb)
      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE

        fr = aux%ros_gammar(lay,prof) / (aux%ros_gammar(lay,prof) + inu)

        ! This is the original Fb from Rosenkranz. This requires the principal value of each LOG
        ! (i.e. the complex argument must lie between -pi and pi) so care must be taken when
        ! simplifying the expression to aid TL/AD coding and/or to make it more efficient.

        ! fb = 0.5_jprb * LOG((z2 - inu) / (z1 - inu)) / LOG(z2 / z1) + &
        !      0.5_jprb * LOG((CONJG(z2) - inu) / (CONJG(z1) - inu)) / LOG(CONJG(z2) / CONJG(z1))

        ! This is the above expression written using the aux variables

        ! fb = 0.5_jprb * LOG((aux%ros_z2 - inu) / (aux%ros_z1(lay,prof) - inu)) / &
        !                 LOG(aux%ros_z2 / aux%ros_z1(lay,prof)) + &
        !      0.5_jprb * LOG((CONJG(aux%ros_z2) - inu) / (CONJG(aux%ros_z1(lay,prof)) - inu)) / &
        !                 LOG(CONJG(aux%ros_z2) / CONJG(aux%ros_z1(lay,prof)))

        ! In the following version of the expression the complex divisions within each LOG have been
        ! re-expressed with a complex numerator and strictly real denominator so that the LOGs of
        ! the numerator and denominator can be separated without causing problems. This simplifies
        ! the expression for the TL/AD and eliminates the complex divisions.
        ! Note that LOG(CONJG(z)) == CONJG(LOG(z))

        log_abs_z1_m_inu_sq     = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + &  ! LOG(ABS(z1-inu)**2)
                                      (aux%ros_nu1(lay,prof) - nu)**2)
        log_abs_z1star_m_inu_sq = LOG(0.5625_jprb * aux%ros_nu1(lay,prof)**2 + &  ! LOG(ABS(CONJG(z1)-inu)**2)
                                      (aux%ros_nu1(lay,prof) + nu)**2)

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
        od(lay) = -1.5_jprb * nu * zgstar * aux%clw(lay,prof)
      ENDDO

    ELSEIF (opts%rt_mw%clw_scheme == mw_clw_scheme_tkc) THEN

      ! Turner, Kneifel, Cadeddu (2016)

      ! Adapted from mw_scatt_coef/perm_water_TKC_16.
      ! Above 500GHz this extrapolates the frequency-dependence.
      ! The fitting formula operates on deg Celsius and Hz.

      freq = coef%frequency_ghz(chan) * 1.E9_jprb

      DO lay = 1, nlayers
        lev = lay + 1
        IF (profiles(prof)%p(MIN(lev+1, nlevels)) < opts%rt_mw%clw_cloud_top) CYCLE

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
        od(lay) = -1.5_jprb * coef%frequency_ghz(chan) * zgstar * aux%clw(lay,prof)
      ENDDO

    ENDIF

    ! Add CLW optical depths to totals
    ztemp = 0._jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ! Note that optical depth in the calculations is negative
      ztemp = ztemp + od(lay)
      opdp_path%atm_level(lev,j) = opdp_path%atm_level(lev,j) + ztemp
    ENDDO
  ENDDO ! chanprof

END SUBROUTINE rttov_mw_clw_absorption
