! Description:
!> @file
!!   AD of atmospheric radiation path calculations.
!!
!> @brief
!!   AD of atmospheric radiation path calculations.
!!
!! @param[in]     opts            RTTOV options structure
!! @param[in]     dosolar         flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel  flag for strict plane parallel geometry
!! @param[in]     profiles        profiles structure (on user levels or coefficient levels)
!! @param[in,out] profiles_ad     profiles structure containing increments
!! @param[in]     profiles_dry    profiles structure containing gas profiles in units of ppmv dry
!! @param[in,out] profiles_dry_ad profiles structure containing gas increments in units of ppmv dry
!! @param[in]     auxs            RTTOV profile_aux_s structure
!! @param[in]     angles          geometry structure
!! @param[in]     raytracing      raytracing structure
!! @param[in,out] raytracing_ad   raytracing structure containing increments
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
SUBROUTINE rttov_locpat_ad(opts, dosolar, plane_parallel, profiles, profiles_ad, &
      profiles_dry, profiles_dry_ad, auxs, angles, raytracing, raytracing_ad)

  USE rttov_types, ONLY : rttov_options, rttov_profile_aux_s, rttov_profile, &
                          rttov_geometry, rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4, &
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt, &
                          eqrad, omega, grave, t0, max_sol_zen, co2_conc

  USE parkind1, ONLY : jpim, jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_math_mod, ONLY : divide, invsqrt, sintosec_ad, reciprocal, reciprocal_ad
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),       INTENT(IN)    :: opts
  LOGICAL(jplm),             INTENT(IN)    :: dosolar
  LOGICAL(jplm),             INTENT(IN)    :: plane_parallel
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_ad(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_dry(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_dry_ad(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_geometry),      INTENT(IN)    :: angles(:)
  TYPE(rttov_raytracing),    INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),    INTENT(INOUT) :: raytracing_ad
!INTF_END

  INTEGER(jpim) :: iprof, lay, nlevels, nlayers, nprofiles
  REAL   (jprb) :: rearth(SIZE(profiles))
  REAL   (jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (jprb) :: gravh_r(SIZE(profiles))
  REAL   (jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (jprb) :: rlh(SIZE(profiles))
  REAL   (jprb) :: dflat, fac, c, scale_eps
  REAL   (jprb) :: pres_ad(profiles(1)%nlevels,SIZE(profiles))
  REAL   (jprb) :: qwet(profiles(1)%nlevels)
  REAL   (jprb) :: qwet_ad(profiles(1)%nlevels)
  REAL   (jprb) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_AD', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

  ! Scale water profile in air density calculation
  scale_eps = 1._jprb - mh2o/mair

  dflat = (1._jprb - flatt) ** 2_jpim
  fac = (omega ** 2_jpim * (eqrad * 1000._jprb)) / (grave)
  beta = 5._jprb * fac / 2._jprb - flatt - 17._jprb * fac * flatt / 14._jprb
  eta = flatt * (5._jprb * fac - flatt) / 8._jprb


  ! Calculate the earth's radius at latitude lat assuming the
  ! earth is an ellipsoid of revolution
  rearth(:) = SQRT(eqrad ** 2_jpim * dflat / &
               (angles(:)%sinlat ** 2_jpim + &
                dflat * angles(:)%coslat ** 2_jpim))

  rlh(:) = rearth(:) * 5.E2_jprb ! DAR: used to be gravl / gravh


  ! The value of earth's gravity at surface at latitude lat is
  ! computed using the international gravity formula.
  gravl(:) = grave * &
    (1._jprb + beta * angles(:)%sinlat ** 2_jpim +      &
    eta * (2._jprb * angles(:)%sinlat * angles(:)%coslat) ** 2_jpim)

  gravh_r(:) = rearth(:) / (2.E-3_jprb * gravl(:)) ! computational useful form

  c = 1000._jprb * rgc / (100._jprb * mair)


! Start AD

  IF (opts%interpolation%lgradp) pres_ad = 0._jprb


  ! Layer thickness
  IF (opts%rt_ir%addsolar .OR. opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO iprof = 1, nprofiles
      DO lay = nlevels - 1, 1,  - 1
        raytracing_ad%hgpl(lay,iprof)   = raytracing_ad%hgpl(lay,iprof) + raytracing_ad%ltick(lay,iprof)
        raytracing_ad%hgpl(lay+1,iprof) = raytracing_ad%hgpl(lay+1,iprof) - raytracing_ad%ltick(lay,iprof)
      ENDDO
    ENDDO
  ENDIF


  IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
    raytracing_ad%r = 0._jprb ! moved from init_raytracing
  ENDIF
  raytracing_ad%ppw = 0._jprb


  ! Compute local path angles
  IF (plane_parallel) THEN
    raytracing_ad%zasat         = 0._jprb
    raytracing_ad%pathsat       = 0._jprb

    IF (dosolar) THEN
      raytracing_ad%zasun   = 0._jprb
      raytracing_ad%pathsun = 0._jprb
      raytracing_ad%patheff = 0._jprb
    ENDIF
  ELSE
    ! Compute secant of the zenith angle at the lower boundary of each -
    ! layer for the satellite-to-surface line of view and of each
    ! layer for the sun-to-surface line of view (if possible and required)

    IF (dosolar) THEN
      raytracing_ad%pathsat = raytracing_ad%pathsat + raytracing_ad%patheff
      raytracing_ad%pathsun = raytracing_ad%pathsun + raytracing_ad%patheff

      CALL calc_seczen_ad(&
        raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
        raytracing_ad%ratoesun, raytracing_ad%zasun, raytracing_ad%pathsun,&
        do_sat = .FALSE._jplm) ! solar call
  !  ELSE
  !    raytracing_ad%pathsun(:,:)  = 0._jprb
    ENDIF

    CALL calc_seczen_ad(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
      raytracing_ad%ratoesat, raytracing_ad%zasat, raytracing_ad%pathsat, &
      do_sat = .TRUE._jplm) ! satellite call

    ! Compute atmospheric refractive index (if required) for dry air given temperature
    ! and pressure as a function of wavenumber using an updated version of edlen equation
    ! not yet corrected for CO2 contribution - used in calc_seczen
    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) CALL calc_refractivity_ad()
  ENDIF


  ! Compute height of pressure levels
  CALL calc_hgpl_ad()


  ! Calculate moist air density (for later calculation).
  ! Density of dry air is adjusted to account for the presence of water vapour
  ! by replacing temp with virtual temp
  DO iprof = 1, nprofiles

    ! Calculate water vapour in units of ppmv wet
    qwet(:) = profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1.E-6_jprb)

    CALL divide(-c * raytracing%dmair(:,iprof) * raytracing_ad%dmair(:,iprof), &
      c * profiles(iprof)%t(:), &
      profiles_ad(iprof)%t(:), acc = .TRUE._jplm)

    CALL divide(-raytracing_ad%dmair(:,iprof) * scale_eps, &
      c * profiles(iprof)%t(:), &
      raytracing_ad%ppw(:,iprof), acc = .TRUE._jplm)

      qwet_ad(:) = &!qwet_ad(:) + &
        1.E-6_jprb * profiles(iprof)%p(:) * raytracing_ad%ppw(:,iprof)

    IF (opts%interpolation%lgradp) THEN
      CALL divide(raytracing_ad%dmair(:,iprof), &
        c * profiles(iprof)%t(:), &
        pres_ad(:,iprof), acc = .TRUE._jplm)

      pres_ad(:,iprof) = pres_ad(:,iprof) + &
        1.E-6_jprb * raytracing_ad%ppw(:,iprof) * qwet(:)

      profiles_ad(iprof)%p(:) = profiles_ad(iprof)%p(:) + &
        pres_ad(:,iprof)
    ENDIF

    profiles_dry_ad(iprof)%q(:) = profiles_dry_ad(iprof)%q(:) + &
      qwet_ad(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1.E-6_jprb)**2
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_AD', 1_jpim, ZHOOK_HANDLE)
CONTAINS

  SUBROUTINE calc_refractivity_ad()

    REAL(jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(jprb) :: t_norm((profiles(1)%nlevels))
    REAL(jprb) :: dry_pp_h_ad((profiles(1)%nlevels))
    REAL(jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(jprb) :: DU_AD((profiles(1)%nlevels)), DV_AD((profiles(1)%nlevels))

    REAL(jprb) :: disp ! the value of the refractive index given by the dispersion equation.
    REAL(jprb) :: disp_moist
    REAL(jprb) :: dispco2
    INTEGER(jpim) :: iprof

    ! Dry air dispersion equation constant in updated Edlen eqn
    disp = 1.E-8_jprb * ( &
      d1 + d2 / (d3 - (waver * ctom) ** 2_jpim) + &
           d4 / (d5 - (waver * ctom) ** 2_jpim))

    ! Moist air dispersion constant
    disp_moist = htop * (ew1 - ew2 * (waver * ctom) ** 2_jpim) * 1.E-10_jprb

    DO iprof = 1, nprofiles

      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:,iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0

      U = 1._jprb + 1.E-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
      CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R)

      raytracing_ad%ppw(:,iprof) = -&! raytracing_ad%ppw(:,iprof) - &!
        (raytracing_ad%r(:,iprof) * disp_moist)

      ! apply CO2 correction
      IF (opts%rt_all%co2_data) THEN
        raytracing_ad%dispco2(:,iprof) = &! raytracing_ad%dispco2(:,iprof) + &
             raytracing%refractivity(:,iprof) * raytracing_ad%r(:,iprof)

        raytracing_ad%refractivity(:,iprof) = &!raytracing_ad%refractivity(:,iprof) + &
          raytracing%dispco2(:,iprof) * raytracing_ad%r(:,iprof)
      ELSE
        dispco2 = disp * (1._jprb + dco2 * (co2_conc * 1.E-6_jprb - 0.0003_jprb))
        raytracing_ad%refractivity(:,iprof) = &! raytracing_ad%refractivity(:,iprof) +
          raytracing_ad%r(:,iprof) * dispco2 ! AD
      ENDIF

      dry_pp_h_ad(:) = &!dry_pp_h_ad(:) + &
        raytracing_ad%refractivity(:,iprof) * (U * V_R)/ed1

      DU_ad(:) = &!DU_ad(:) + &
        dry_pp_h(:) * raytracing_ad%refractivity(:,iprof) * V_R(:) / ed1

      DV_ad(:) = &!DV_ad(:) &
        -raytracing%refractivity(:,iprof) * raytracing_ad%refractivity(:,iprof) * V_R(:)

      profiles_ad(iprof)%t(:) = profiles_ad(iprof)%t(:) + &
        ed4 * DV_ad(:) + & !combine two lines from TL
        DU_ad(:) * 1.E-8_jprb * dry_pp_h(:) * (-ed3)

      dry_pp_h_ad(:) = dry_pp_h_ad(:) + 1.E-8_jprb * (ed2 - ed3 * t_norm(:)) * DU_ad(:)

      raytracing_ad%ppw(:,iprof) = raytracing_ad%ppw(:,iprof) - htop * dry_pp_h_ad(:)

      IF (opts%interpolation%lgradp) THEN
        pres_ad(:,iprof) = pres_ad(:,iprof) + htop * dry_pp_h_ad(:)
      ENDIF
    ENDDO

    IF (opts%rt_all%co2_data) THEN
      DO iprof = 1, nprofiles
        profiles_dry_ad(iprof)%co2(:) = profiles_dry_ad(iprof)%co2(:) + &
          raytracing_ad%dispco2(:,iprof) * disp * &
          dco2 * 1.E-6_jprb
      ENDDO
    ENDIF

  END SUBROUTINE calc_refractivity_ad


  SUBROUTINE calc_hgpl_ad()
    REAL(jprb)    :: dp(profiles(1)%nlevels)
    REAL(jprb)    :: dp_ad(profiles(1)%nlevels)
    REAL(jprb)    :: ztemp(profiles(1)%nlevels,SIZE(profiles))
    INTEGER(jpim) :: i_above, i_below
    REAL(jprb)    :: dp_above, dp_below, t_surf, q_surf, ppw_surf, dmair_surf
    REAL(jprb)    :: dp_above_ad, dp_below_ad, t_surf_ad, q_surf_ad, ppw_surf_ad, dmair_surf_ad
    INTEGER(jpim) :: iprof, lev

    CALL invsqrt(raytracing%ztemp(:,:), ztemp(:,:))

    DO iprof = 1, nprofiles
      dp(1:nlevels-1) = profiles(iprof)%p(1:nlevels-1) - profiles(iprof)%p(2:nlevels)

      DO lev = nlevels, auxs(iprof)%nearestlev_surf + 1, -1

        raytracing_ad%ztemp(lev,iprof) = &!raytracing_ad%ztemp(lev,iprof) -
          -0.5_jprb * ztemp(lev,iprof) * 1.E-3_jprb * raytracing_ad%hgpl(lev,iprof)

        raytracing_ad%int(lev,iprof) = &!raytracing_ad%int(lev,iprof) - &
          -2._jprb * 100._jprb * gravh_r(iprof) * raytracing_ad%ztemp(lev,iprof)

        raytracing_ad%hgpl(lev - 1,iprof) = raytracing_ad%hgpl(lev - 1,iprof) + 1.E3_jprb * &
          2._jprb * raytracing_ad%ztemp(lev,iprof) * &
          (1.E3_jprb * raytracing%hgpl(lev - 1,iprof) - rlh(iprof))

        IF (raytracing%dmair(lev,iprof) == raytracing%dmair(lev-1,iprof)) THEN
          raytracing_ad%dmair(lev,iprof) = raytracing_ad%dmair(lev,iprof) - &
            raytracing_ad%int(lev,iprof) * &
            raytracing%int(lev,iprof) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing_ad%dmair(lev-1,iprof) = raytracing_ad%dmair(lev-1,iprof) + &
            (dp(lev-1) * raytracing_ad%int(lev,iprof) / raytracing%dmair(lev-1,iprof) + &
             raytracing%int(lev,iprof) * raytracing_ad%int(lev,iprof)) / &
            (raytracing%dmair(lev,iprof) - raytracing%dmair(lev-1,iprof))

          raytracing_ad%dmair(lev,iprof) = raytracing_ad%dmair(lev,iprof) - &
            (dp(lev-1) * raytracing_ad%int(lev,iprof) / raytracing%dmair(lev,iprof) + &
             raytracing%int(lev,iprof) * raytracing_ad%int(lev,iprof)) / &
            (raytracing%dmair(lev,iprof) - raytracing%dmair(lev-1,iprof))
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          dp_ad(lev-1) = &!dp_ad(lev-1) + &
            raytracing_ad%int(lev,iprof) * raytracing%int(lev,iprof) / dp(lev-1)
        ENDIF
      ENDDO

      DO lev = 1, auxs(iprof)%nearestlev_surf - 1
        raytracing_ad%ztemp(lev,iprof) = &!raytracing_ad%ztemp(lev,iprof) -
          -0.5_jprb * ztemp(lev,iprof) * 1.E-3_jprb * raytracing_ad%hgpl(lev,iprof)

        raytracing_ad%int(lev,iprof) = &!raytracing_ad%int(lev,iprof) - &
          -2._jprb * 100._jprb * gravh_r(iprof) * raytracing_ad%ztemp(lev,iprof)

        raytracing_ad%hgpl(lev + 1,iprof) = raytracing_ad%hgpl(lev + 1,iprof) + 1.E3_jprb * &
          2._jprb * raytracing_ad%ztemp(lev,iprof) * &
          (1.E3_jprb * raytracing%hgpl(lev + 1,iprof) - rlh(iprof))

        IF (raytracing%dmair(lev+1,iprof) == raytracing%dmair(lev,iprof)) THEN
          raytracing_ad%dmair(lev,iprof) = raytracing_ad%dmair(lev,iprof) - &
            raytracing_ad%int(lev,iprof) * &
            raytracing%int(lev,iprof) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing_ad%dmair(lev,iprof) = raytracing_ad%dmair(lev,iprof) + &
            (dp(lev) * raytracing_ad%int(lev,iprof) / raytracing%dmair(lev,iprof) + &
             raytracing%int(lev,iprof) * raytracing_ad%int(lev,iprof)) / &
            (raytracing%dmair(lev+1,iprof) - raytracing%dmair(lev,iprof))

          raytracing_ad%dmair(lev+1,iprof) = raytracing_ad%dmair(lev+1,iprof) - &
            (dp(lev) * raytracing_ad%int(lev,iprof) / raytracing%dmair(lev+1,iprof) + &
             raytracing%int(lev,iprof) * raytracing_ad%int(lev,iprof)) / &
            (raytracing%dmair(lev+1,iprof) - raytracing%dmair(lev,iprof))
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          dp_ad(lev) = &!dp_ad(lev) + &
            raytracing_ad%int(lev,iprof) * raytracing%int(lev,iprof) / dp(lev)
        ENDIF
      ENDDO

      IF (opts%interpolation%lgradp) THEN
        DO lev = nlevels - 1, 1, -1
          pres_ad(lev,iprof) = pres_ad(lev,iprof) + dp_ad(lev)
          pres_ad(lev+1,iprof) = pres_ad(lev+1,iprof) - dp_ad(lev)
        ENDDO
      ENDIF

      ! Set the height of the surface
      i_below     = auxs(iprof)%nearestlev_surf
      i_above     = auxs(iprof)%nearestlev_surf - 1
      dp_above    = profiles(iprof)%p(i_above) - profiles(iprof)%s2m%p
      dp_below    = profiles(iprof)%s2m%p - profiles(iprof)%p(i_below)

      IF (dp_below /= 0._jprb .AND. opts%config%fix_hgpl) THEN

        q_surf        = (profiles_dry(iprof)%q(i_below) * dp_above + &
                         profiles_dry(iprof)%q(i_above) * dp_below) / (dp_above + dp_below)
        t_surf        = (profiles(iprof)%t(i_below) * dp_above + &
                         profiles(iprof)%t(i_above) * dp_below) / (dp_above + dp_below)
        ppw_surf      = profiles(iprof)%s2m%p * 1.E-6_jprb * q_surf / (1._jprb + q_surf * 1.E-6_jprb)
        dmair_surf    = (profiles(iprof)%s2m%p - ppw_surf * scale_eps) / (c * t_surf)


        dp_above_ad   = 0._jprb
        dp_below_ad   = 0._jprb
        q_surf_ad     = 0._jprb
        t_surf_ad     = 0._jprb
        ppw_surf_ad   = 0._jprb
        dmair_surf_ad = 0._jprb

        raytracing_ad%ztemp(i_below,iprof) = &!&raytracing_ad%ztemp(i_below,iprof) - &
          -0.5_jprb * ztemp(i_below,iprof) * 1.E-3_jprb * raytracing_ad%hgpl(i_below,iprof)

        raytracing_ad%int(i_below,iprof) = &!raytracing_ad%int(i_below,iprof) - &
          -2._jprb * raytracing_ad%ztemp(i_below,iprof) * 1.E2_jprb * gravh_r(iprof)

        IF (raytracing%dmair(i_below,iprof) == dmair_surf) THEN
          dp_below_ad = dp_below_ad + &
            raytracing_ad%int(i_below,iprof) / dmair_surf

          dmair_surf_ad = dmair_surf_ad - &
            raytracing_ad%int(i_below,iprof) * raytracing%int(i_below,iprof) / dmair_surf
        ELSE
          dp_below_ad = dp_below_ad + &
            raytracing_ad%int(i_below,iprof) * raytracing%int(i_below,iprof) / dp_below

          raytracing_ad%dmair(i_below,iprof) = raytracing_ad%dmair(i_below,iprof) + &
            (dp_below / raytracing%dmair(i_below,iprof) - raytracing%int(i_below,iprof)) * &
            raytracing_ad%int(i_below,iprof) / (raytracing%dmair(i_below,iprof) - dmair_surf)

          dmair_surf_ad = dmair_surf_ad + &
            (raytracing%int(i_below,iprof) - dp_below / dmair_surf) * &
            raytracing_ad%int(i_below,iprof) / (raytracing%dmair(i_below,iprof) - dmair_surf)
        ENDIF

        profiles_ad(iprof)%s2m%p = profiles_ad(iprof)%s2m%p + dmair_surf_ad / (c * t_surf)
        ppw_surf_ad = ppw_surf_ad - dmair_surf_ad * scale_eps / (c * t_surf)
        t_surf_ad = t_surf_ad - dmair_surf_ad * dmair_surf / t_surf

        profiles_ad(iprof)%s2m%p = profiles_ad(iprof)%s2m%p + &
          ppw_surf_ad * q_surf * 1.E-6_jprb / (1._jprb + q_surf * 1.E-6_jprb)
        q_surf_ad = q_surf_ad + ppw_surf_ad * (profiles(iprof)%s2m%p - ppw_surf) * &
          1.E-6_jprb / (1._jprb + q_surf * 1.E-6_jprb)

        profiles_ad(iprof)%t(i_below) = profiles_ad(iprof)%t(i_below) + &
          t_surf_ad * dp_above / (dp_above + dp_below)
        profiles_ad(iprof)%t(i_above) = profiles_ad(iprof)%t(i_above) + &
          t_surf_ad * dp_below / (dp_above + dp_below)
        dp_above_ad = dp_above_ad + &
          t_surf_ad * (profiles(iprof)%t(i_below) - t_surf) / (dp_above + dp_below)
        dp_below_ad = dp_below_ad + &
          t_surf_ad * (profiles(iprof)%t(i_above) - t_surf) / (dp_above + dp_below)

        profiles_dry_ad(iprof)%q(i_below) = profiles_dry_ad(iprof)%q(i_below) + &
          q_surf_ad * dp_above / (dp_above + dp_below)
        profiles_dry_ad(iprof)%q(i_above) = profiles_dry_ad(iprof)%q(i_above) + &
          q_surf_ad * dp_below / (dp_above + dp_below)
        dp_above_ad = dp_above_ad + &
          q_surf_ad * (profiles_dry(iprof)%q(i_below) - q_surf) / (dp_above + dp_below)
        dp_below_ad = dp_below_ad + &
          q_surf_ad * (profiles_dry(iprof)%q(i_above) - q_surf) / (dp_above + dp_below)

        IF (opts%interpolation%lgradp) THEN
          pres_ad(i_below,iprof) = pres_ad(i_below,iprof) - dp_below_ad
          pres_ad(i_above,iprof) = pres_ad(i_above,iprof) + dp_above_ad
        ENDIF
        profiles_ad(iprof)%s2m%p = profiles_ad(iprof)%s2m%p + dp_below_ad - dp_above_ad

      ELSE
        raytracing_ad%hgpl(i_below,iprof) = 0._jprb
      ENDIF

    ENDDO
  END SUBROUTINE calc_hgpl_ad


  SUBROUTINE calc_seczen_ad(ra,za,path,ra_ad,za_ad,path_ad,do_sat)
    LOGICAL(jplm), INTENT(IN)    :: do_sat ! do satellite secant or solar secant
    REAL(jprb),    INTENT(IN)    :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb),    INTENT(INOUT) :: ra_ad(:,:), za_ad(:,:), path_ad(:,:)

    REAL(jprb)    :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    INTEGER(jpim) :: n, iprof

    IF (opts%rt_ir%addsolar) THEN
      CALL sintosec_ad(za(1:nlayers,:), za_ad(1:nlayers,:), path, path_ad, acc = .TRUE._jplm)
    ELSE
      CALL sintosec_ad(za(1:nlayers,:), za_ad(1:nlayers,:), path, path_ad, acc = .FALSE._jplm) ! first use of za
    ENDIF

!   path_ad = 0._jprb;

    DO iprof = 1, nprofiles
      IF (do_sat) THEN
        ra_ad(:,iprof) = &!ra_ad(:,iprof) +
                         za_ad(1:nlayers,iprof) * angles(iprof)%sinzen
      ELSE ! do_sun
        IF (profiles(iprof)%sunzenangle >= 0._jprb .AND. &
            profiles(iprof)%sunzenangle < max_sol_zen) THEN
          ra_ad(:,iprof) = &!ra_ad(:,iprof) +
                           za_ad(1:nlayers,iprof) * angles(iprof)%sinzen_sun
        ELSE
          ra_ad(:,iprof) = 0._jprb
        ENDIF
      ENDIF
    ENDDO
!  za_ad = 0._jprb;

    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
      IF (do_sat) THEN
        n = 2
      ELSE ! do sun
        n = nlevels
      ENDIF
      DO iprof = 1, nprofiles
        z(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
        raytracing_ad%r(2:nlevels,iprof) = raytracing_ad%r(2:nlevels,iprof) - & !first use
          z(:,iprof) * raytracing%r(n,iprof) * raytracing%r_r(2:nlevels,iprof)**2_jpim * ra_ad(:,iprof)

        raytracing_ad%r(n,iprof) = raytracing_ad%r(n,iprof) + &
          SUM(raytracing%r_r(2:nlevels,iprof) * z(:,iprof) * ra_ad(:,iprof))

        ra_ad(:,iprof) = &
          raytracing%r_r(2:nlevels,iprof) * raytracing%r(n,iprof) * ra_ad(:,iprof)
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      IF (do_sat .OR. (profiles(iprof)%sunzenangle >= 0._jprb .AND. &
                       profiles(iprof)%sunzenangle < max_sol_zen)) THEN
        raytracing_ad%z_r(:,iprof) = &!raytracing_ad%z_r(:,iprof) + &
          rearth(iprof) * ra_ad(:,iprof)
      ELSE
        raytracing_ad%z_r(:,iprof) = 0._jprb
      ENDIF
    ENDDO
!   ra_ad = 0._jprb;

    IF (opts%config%fix_hgpl) THEN
      DO iprof = 1, nprofiles
        raytracing_ad%z_r(auxs(iprof)%nearestlev_surf-1,iprof) = 0._jprb
      ENDDO
    ENDIF
    CALL reciprocal_ad(raytracing%z_r,raytracing_ad%hgpl(2:nlayers+1, :), &
         raytracing_ad%z_r, acc = .TRUE._jplm)
  END SUBROUTINE calc_seczen_ad

END SUBROUTINE rttov_locpat_ad
