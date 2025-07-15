! Description:
!> @file
!!   Jacobian of atmospheric radiation path calculations.
!!
!> @brief
!!   Jacobian of atmospheric radiation path calculations.
!!
!! @param[in]     opts            RTTOV options structure
!! @param[in]     dosolar         flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel  flag for strict plane parallel geometry
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     profiles        profiles structure (on user levels or coefficient levels)
!! @param[in,out] profiles_k      profiles structure containing increments
!! @param[in]     profiles_dry    profiles structure containing gas profiles in units of ppmv dry
!! @param[in,out] profiles_dry_k  profiles structure containing gas increments in units of ppmv dry
!! @param[in]     auxs            RTTOV profile_aux_s structure
!! @param[in]     angles          geometry structure
!! @param[in]     raytracing      raytracing structure
!! @param[in,out] raytracing_k    raytracing structure containing increments
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
SUBROUTINE rttov_locpat_k(opts, dosolar, plane_parallel, chanprof, profiles, profiles_k, &
      profiles_dry, profiles_dry_k, auxs, angles, raytracing, raytracing_k)

  USE rttov_types, ONLY : rttov_options, rttov_profile_aux_s, rttov_profile, &
                          rttov_geometry, rttov_raytracing, rttov_chanprof
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4, &
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt, &
                          eqrad, omega, grave, t0, max_sol_zen, co2_conc

  USE parkind1, ONLY : jpim, jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_math_mod, ONLY : invsqrt, sintosec_k, reciprocal, reciprocal_k
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),       INTENT(IN)    :: opts
  LOGICAL(jplm),             INTENT(IN)    :: dosolar
  LOGICAL(jplm),             INTENT(IN)    :: plane_parallel
  TYPE(rttov_chanprof),      INTENT(IN)    :: chanprof(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_dry(:)
  TYPE(rttov_profile),       INTENT(INOUT) :: profiles_dry_k(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_geometry),      INTENT(IN)    :: angles(:)
  TYPE(rttov_raytracing),    INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),    INTENT(INOUT) :: raytracing_k
!INTF_END

  INTEGER(jpim) :: lay, prof, i, nlevels, nlayers, nchannels
  REAL   (jprb) :: rearth(SIZE(profiles))
  REAL   (jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (jprb) :: gravh_r(SIZE(profiles))
  REAL   (jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (jprb) :: rlh(SIZE(profiles))
  REAL   (jprb) :: dflat, fac, c, scale_eps
  REAL   (jprb) :: pres_k(profiles(1)%nlevels,SIZE(chanprof))
  REAL   (jprb) :: qwet(profiles(1)%nlevels)
  REAL   (jprb) :: qwet_k(profiles(1)%nlevels)
  INTEGER(jpim) :: map(SIZE(chanprof),2), prof_stat, last_prof
  REAL   (jprb) :: ztemp(profiles(1)%nlevels, 4)
  REAL   (jprb) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_K', 0_jpim, ZHOOK_HANDLE)

  nchannels = SIZE(chanprof)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

  map(1,1) = chanprof(1)%prof
  map(1,2) = chanprof(1)%chan

  prof_stat = 1 ! assume profs are contiguous and monotonic
  DO i = 2, nchannels
    map(i,1) = chanprof(i)%prof
    map(i,2) = chanprof(i)%chan

    IF (map(i,1) < map(i-1,1)) THEN ! they are not.
      prof_stat = -1
    ENDIF
  ENDDO

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


! Start K

  IF (opts%interpolation%lgradp) pres_k = 0._jprb


  ! Layer thickness
  IF (opts%rt_ir%addsolar .OR. opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO i = 1, nchannels
      DO lay = nlevels - 1, 1,  - 1
        raytracing_k%hgpl(lay,i)   = raytracing_k%hgpl(lay,i) + raytracing_k%ltick(lay,i)
        raytracing_k%hgpl(lay+1,i) = raytracing_k%hgpl(lay+1,i) - raytracing_k%ltick(lay,i)
      ENDDO
    ENDDO
  ENDIF


  IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
    raytracing_k%r = 0._jprb ! moved from init_raytracing
  ENDIF
  raytracing_k%ppw = 0._jprb


  ! Compute local path angles
  IF (plane_parallel) THEN
    raytracing_k%zasat         = 0._jprb
    raytracing_k%pathsat       = 0._jprb

    IF (dosolar) THEN
      raytracing_k%zasun   = 0._jprb
      raytracing_k%pathsun = 0._jprb
      raytracing_k%patheff = 0._jprb
    ENDIF
  ELSE
    ! Compute secant of the zenith angle at the lower boundary of each -
    ! layer for the satellite-to-surface line of view and of each
    ! layer for the sun-to-surface line of view (if possible and required)

    IF (dosolar) THEN
      raytracing_k%pathsat = raytracing_k%pathsat + raytracing_k%patheff
      raytracing_k%pathsun = raytracing_k%pathsun + raytracing_k%patheff

      CALL calc_seczen_k(&
        raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
        raytracing_k%ratoesun, raytracing_k%zasun, raytracing_k%pathsun,&
        do_sat = .FALSE._jplm) ! solar call
  !  ELSE
  !    raytracing_k%pathsun(:,:)  = 0._jprb
    ENDIF

    CALL calc_seczen_k(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
      raytracing_k%ratoesat, raytracing_k%zasat, raytracing_k%pathsat, &
      do_sat = .TRUE._jplm) ! satellite call

    ! Compute atmospheric refractive index (if required) for dry air given temperature
    ! and pressure as a function of wavenumber using an updated version of edlen equation
    ! not yet corrected for CO2 contribution - used in calc_seczen
    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) CALL calc_refractivity_k()
  ENDIF


  ! Compute height of pressure levels
  CALL calc_hgpl_k()


  ! Calculate moist air density (for later calculation).
  ! Density of dry air is adjusted to account for the presence of water vapour
  ! by replacing temp with virtual temp
  last_prof = -1
  DO i = 1, nchannels
    prof = chanprof(i)%prof
    IF (prof .NE. last_prof) THEN
      CALL reciprocal(&
        c * profiles(prof)%t(:), &
        ztemp(:,1))
      ztemp(:,2) = ztemp(:,1) * raytracing%dmair(:,prof)

    ! Calculate water vapour in units of ppmv wet
      qwet(:) = profiles_dry(prof)%q(:) / (1._jprb + profiles_dry(prof)%q(:) * 1.E-6_jprb)

      last_prof = prof
    ENDIF

    profiles_k(i)%t(:) = profiles_k(i)%t(:) - &
      c * raytracing_k%dmair(:,i) * ztemp(:,2)

    raytracing_k%ppw(:,i) = raytracing_k%ppw(:,i) - &
      raytracing_k%dmair(:,i) * scale_eps * ztemp(:,1)


    qwet_k(:) = &!qwet_k(:) + &
      1.E-6_jprb * profiles(prof)%p(:) * raytracing_k%ppw(:,i)

    IF (opts%interpolation%lgradp) THEN
      pres_k(:,i) = pres_k(:,i) + raytracing_k%dmair(:,i) * ztemp(:,1)

      pres_k(:,i) = pres_k(:,i) + &
        1.E-6_jprb * raytracing_k%ppw(:,i) * qwet(:)

      profiles_k(i)%p(:) = profiles_k(i)%p(:) + pres_k(:,i)
    ENDIF

    profiles_dry_k(i)%q(:) = profiles_dry_k(i)%q(:) + &
      qwet_k(:) / (1._jprb + profiles_dry(prof)%q(:) * 1.E-6_jprb)**2
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_K', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_refractivity_k()

    REAL(jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(jprb) :: t_norm((profiles(1)%nlevels))
    REAL(jprb) :: dry_pp_h_k((profiles(1)%nlevels))
    REAL(jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(jprb) :: DU_K((profiles(1)%nlevels)), DV_K((profiles(1)%nlevels))

    REAL(jprb) :: disp ! the value of the refractive index given by the dispersion equation.
    REAL(jprb) :: disp_moist
    REAL(jprb) :: dispco2
    INTEGER(jpim) :: prof, i, last_prof

    ! Dry air dispersion equation constant in updated Edlen eqn
    disp = 1.E-8_jprb * ( &
      d1 + d2 / (d3 - (waver * ctom) ** 2_jpim) + &
           d4 / (d5 - (waver * ctom) ** 2_jpim))

    ! Moist air dispersion constant
    disp_moist = htop * (ew1 - ew2 * (waver * ctom) ** 2_jpim) * 1.E-10_jprb

    last_prof = -1
    DO i = 1, nchannels
      ! Store temporary FWD variables
      prof = chanprof(i)%prof
      IF (prof .NE. last_prof) THEN
        dry_pp_h(:) = htop * (profiles(prof)%p(:) - raytracing%ppw(:,prof))
        t_norm(:) = profiles(prof)%t(:) - t0

        U = 1._jprb + 1.E-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
        CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R(:))

        ztemp(:,1) = U(:) * V_R / ed1
        ztemp(:,2) = dry_pp_h(:) * V_R(:) / ed1
        ztemp(:,3) = -raytracing%refractivity(:,prof) * V_R(:)
        ztemp(:,4) = 1.E-8_jprb * (ed2 - ed3 * t_norm(:))
        last_prof = prof
      ENDIF

      raytracing_k%ppw(:,i) = -&! raytracing_k%ppw(:,i) - &!
        (raytracing_k%r(:,i) * disp_moist)

      ! apply CO2 correction
      IF (opts%rt_all%co2_data) THEN
        raytracing_k%dispco2(:,i) = &! raytracing_k%dispco2(:,i) + &
             raytracing%refractivity(:,prof) * raytracing_k%r(:,i)

        raytracing_k%refractivity(:,i) = &!raytracing_k%refractivity(:,i) + &
          raytracing%dispco2(:,prof) * raytracing_k%r(:,i)
      ELSE
        dispco2 = disp * (1._jprb + dco2 * (co2_conc * 1.E-6_jprb - 0.0003_jprb))
        raytracing_k%refractivity(:,i) = &! raytracing_k%refractivity(:,i) +
          raytracing_k%r(:,i) * dispco2 ! AD
      ENDIF

      dry_pp_h_k(:) = &!dry_pp_h_k(:) + &
        raytracing_k%refractivity(:,i) * ztemp(:,1)

      DU_k(:) = &!DU_k(:) + &
        raytracing_k%refractivity(:,i) * ztemp(:,2)

      DV_k(:) = &!DV_k(:) &
        raytracing_k%refractivity(:,i) * ztemp(:,3)

      profiles_k(i)%t(:) = profiles_k(i)%t(:) + &
        DV_k(:) * ed4 + & !combine two lines from TL
        DU_k(:) * 1.E-8_jprb * dry_pp_h(:) * (-ed3)

      dry_pp_h_k(:) = dry_pp_h_k(:) + DU_k(:) * ztemp(:,4)

      raytracing_k%ppw(:,i) = raytracing_k%ppw(:,i) - htop * dry_pp_h_k(:)

      IF (opts%interpolation%lgradp) THEN
        profiles_k(i)%p(:) = profiles_k(i)%p(:) + htop * dry_pp_h_k(:)
      ENDIF
    ENDDO

    IF (opts%rt_all%co2_data) THEN
      DO i = 1, nchannels
        profiles_dry_k(i)%co2(:) = profiles_dry_k(i)%co2(:) + &
          raytracing_k%dispco2(:,i) * disp * &
          dco2 * 1.E-6_jprb
      ENDDO
    ENDIF

  END SUBROUTINE calc_refractivity_k


  SUBROUTINE calc_hgpl_k()
    REAL(jprb)    :: dp(profiles(1)%nlevels)
    REAL(jprb)    :: dp_k(profiles(1)%nlevels)
    REAL(jprb)    :: ztemp(profiles(1)%nlevels,SIZE(profiles))
    REAL(jprb)    :: dmair_r(profiles(1)%nlevels)
    REAL(jprb)    :: ddmair_r(profiles(1)%nlayers)
    INTEGER(jpim) :: i_above, i_below
    REAL(jprb)    :: dp_above, dp_below, t_surf, q_surf, ppw_surf, dmair_surf
    REAL(jprb)    :: dp_above_k, dp_below_k, t_surf_k, q_surf_k, ppw_surf_k, dmair_surf_k
    INTEGER(jpim) :: lev, prof, i, last_prof

    CALL invsqrt(raytracing%ztemp(:, :), ztemp(:, :))

    last_prof = -1
    DO i = 1, nchannels
      prof = chanprof(i)%prof

      IF (prof .NE. last_prof) THEN
        DO lev = 1, nlevels - 1
          dp(lev) = profiles(prof)%p(lev) - profiles(prof)%p(lev+1)
        ENDDO

        CALL reciprocal(raytracing%dmair(:,prof), dmair_r)
        WHERE (raytracing%dmair(2:nlevels,prof) /= raytracing%dmair(1:nlevels-1,prof))
          ddmair_r = 1._jprb / (raytracing%dmair(2:nlevels,prof) - raytracing%dmair(1:nlevels-1,prof))
        ENDWHERE

        last_prof = prof
      ENDIF

      DO lev = nlevels, auxs(prof)%nearestlev_surf + 1, -1

        raytracing_k%ztemp(lev,i) = &!raytracing_k%ztemp(lev,i) -
          -0.5_jprb * ztemp(lev,prof) * 1.E-3_jprb * raytracing_k%hgpl(lev,i)

        raytracing_k%int(lev,i) = &!raytracing_k%int(lev,i) - &
          -2._jprb * 100._jprb * gravh_r(prof) * raytracing_k%ztemp(lev,i)

        raytracing_k%hgpl(lev - 1,i) = raytracing_k%hgpl(lev - 1,i) + 1.E3_jprb * &
          2._jprb * raytracing_k%ztemp(lev,i) * &
          (1.E3_jprb * raytracing%hgpl(lev - 1,prof) - rlh(prof))

        IF (raytracing%dmair(lev,prof) == raytracing%dmair(lev-1,prof)) THEN
          raytracing_k%dmair(lev,i) = raytracing_k%dmair(lev,i) - &
            raytracing_k%int(lev,i) * &
            raytracing%int(lev,prof) / raytracing%dmair(lev,prof)
        ELSE
          raytracing_k%dmair(lev-1,i) = raytracing_k%dmair(lev-1,i) + &
            (dp(lev-1) * raytracing_k%int(lev,i) * dmair_r(lev-1) + &
             raytracing%int(lev,prof) * raytracing_k%int(lev,i)) * ddmair_r(lev-1)

          raytracing_k%dmair(lev,i) = raytracing_k%dmair(lev,i) - &
            (dp(lev-1) * raytracing_k%int(lev,i) * dmair_r(lev) + &
             raytracing%int(lev,prof) * raytracing_k%int(lev,i)) * ddmair_r(lev-1)
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          dp_k(lev-1) = &!dp_k(lev-1) + &
            raytracing_k%int(lev,i) * raytracing%int(lev,prof) / dp(lev-1)
        ENDIF
      ENDDO

      DO lev = 1, auxs(prof)%nearestlev_surf - 1
        raytracing_k%ztemp(lev,i) = &!raytracing_k%ztemp(lev,i) -
          -0.5_jprb * ztemp(lev,prof) * 1.E-3_jprb * raytracing_k%hgpl(lev,i)

        raytracing_k%int(lev,i) = &!raytracing_k%int(lev,i) - &
          -2._jprb * 100._jprb * gravh_r(prof) * raytracing_k%ztemp(lev,i)

        raytracing_k%hgpl(lev + 1,i) = raytracing_k%hgpl(lev + 1,i) + 1.E3_jprb * &
          2._jprb * raytracing_k%ztemp(lev,i) * &
          (1.E3_jprb * raytracing%hgpl(lev + 1,prof) - rlh(prof))

        IF (raytracing%dmair(lev+1,prof) == raytracing%dmair(lev,prof)) THEN
          raytracing_k%dmair(lev,i) = raytracing_k%dmair(lev,i) - &
            raytracing_k%int(lev,i) * &
            raytracing%int(lev,prof) / raytracing%dmair(lev,prof)
        ELSE
          raytracing_k%dmair(lev,i) = raytracing_k%dmair(lev,i) + &
            (dp(lev) * raytracing_k%int(lev,i) * dmair_r(lev) + &
             raytracing%int(lev,prof) * raytracing_k%int(lev,i)) * ddmair_r(lev)

          raytracing_k%dmair(lev+1,i) = raytracing_k%dmair(lev+1,i) - &
            (dp(lev) * raytracing_k%int(lev,i) * dmair_r(lev+1) + &
             raytracing%int(lev,prof) * raytracing_k%int(lev,i)) * ddmair_r(lev)
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          dp_k(lev) = &!dp_k(lev) + &
            raytracing_k%int(lev,i) * raytracing%int(lev,prof) / dp(lev)
        ENDIF
      ENDDO

      IF (opts%interpolation%lgradp) THEN
        DO lev = nlevels - 1, 1, -1
          pres_k(lev,i) = pres_k(lev,i) + dp_k(lev)
          pres_k(lev+1,i) = pres_k(lev+1,i) - dp_k(lev)
        ENDDO
      ENDIF

      ! Set the height of the surface
      i_below     = auxs(prof)%nearestlev_surf
      i_above     = auxs(prof)%nearestlev_surf - 1
      dp_above    = profiles(prof)%p(i_above) - profiles(prof)%s2m%p
      dp_below    = profiles(prof)%s2m%p - profiles(prof)%p(i_below)

      IF (dp_below /= 0._jprb .AND. opts%config%fix_hgpl) THEN

        q_surf        = (profiles_dry(prof)%q(i_below) * dp_above + &
                         profiles_dry(prof)%q(i_above) * dp_below) / (dp_above + dp_below)
        t_surf        = (profiles(prof)%t(i_below) * dp_above + &
                         profiles(prof)%t(i_above) * dp_below) / (dp_above + dp_below)
        ppw_surf      = profiles(prof)%s2m%p * 1.E-6_jprb * q_surf / (1._jprb + q_surf * 1.E-6_jprb)
        dmair_surf    = (profiles(prof)%s2m%p - ppw_surf * scale_eps) / (c * t_surf)


        dp_above_k   = 0._jprb
        dp_below_k   = 0._jprb
        q_surf_k     = 0._jprb
        t_surf_k     = 0._jprb
        ppw_surf_k   = 0._jprb
        dmair_surf_k = 0._jprb

        raytracing_k%ztemp(i_below,i) = &!raytracing_k%ztemp(i_below,i) - &
          -0.5_jprb * ztemp(i_below,prof) * 1.E-3_jprb * raytracing_k%hgpl(i_below,i)

        raytracing_k%int(i_below,i) = &!raytracing_k%int(i_below,i) - &
          -2._jprb * raytracing_k%ztemp(i_below,i) * 1.E2_jprb * gravh_r(prof)

        IF (raytracing%dmair(i_below,prof) == dmair_surf) THEN
          dp_below_k = dp_below_k + &
            raytracing_k%int(i_below,i) / dmair_surf

          dmair_surf_k = dmair_surf_k - &
            raytracing_k%int(i_below,i) * raytracing%int(i_below,prof) / dmair_surf
        ELSE
          dp_below_k = dp_below_k + &
            raytracing_k%int(i_below,i) * raytracing%int(i_below,prof) / dp_below

          raytracing_k%dmair(i_below,i) = raytracing_k%dmair(i_below,i) + &
            (dp_below / raytracing%dmair(i_below,prof) - raytracing%int(i_below,prof)) * &
            raytracing_k%int(i_below,i) / (raytracing%dmair(i_below,prof) - dmair_surf)

          dmair_surf_k = dmair_surf_k + &
            (raytracing%int(i_below,prof) - dp_below / dmair_surf) * &
            raytracing_k%int(i_below,i) / (raytracing%dmair(i_below,prof) - dmair_surf)
        ENDIF

        profiles_k(i)%s2m%p = profiles_k(i)%s2m%p + dmair_surf_k / (c * t_surf)
        ppw_surf_k = ppw_surf_k - dmair_surf_k * scale_eps / (c * t_surf)
        t_surf_k = t_surf_k - dmair_surf_k * dmair_surf / t_surf

        profiles_k(i)%s2m%p = profiles_k(i)%s2m%p + &
          ppw_surf_k * q_surf * 1.E-6_jprb / (1._jprb + q_surf * 1.E-6_jprb)
        q_surf_k = q_surf_k + ppw_surf_k * (profiles(prof)%s2m%p - ppw_surf) * &
          1.E-6_jprb / (1._jprb + q_surf * 1.E-6_jprb)

        profiles_k(i)%t(i_below) = profiles_k(i)%t(i_below) + &
          t_surf_k * dp_above / (dp_above + dp_below)
        profiles_k(i)%t(i_above) = profiles_k(i)%t(i_above) + &
          t_surf_k * dp_below / (dp_above + dp_below)
        dp_above_k = dp_above_k + &
          t_surf_k * (profiles(prof)%t(i_below) - t_surf) / (dp_above + dp_below)
        dp_below_k = dp_below_k + &
          t_surf_k * (profiles(prof)%t(i_above) - t_surf) / (dp_above + dp_below)

        profiles_dry_k(i)%q(i_below) = profiles_dry_k(i)%q(i_below) + &
          q_surf_k * dp_above / (dp_above + dp_below)
        profiles_dry_k(i)%q(i_above) = profiles_dry_k(i)%q(i_above) + &
          q_surf_k * dp_below / (dp_above + dp_below)
        dp_above_k = dp_above_k + &
          q_surf_k * (profiles_dry(prof)%q(i_below) - q_surf) / (dp_above + dp_below)
        dp_below_k = dp_below_k + &
          q_surf_k * (profiles_dry(prof)%q(i_above) - q_surf) / (dp_above + dp_below)

        IF (opts%interpolation%lgradp) THEN
          pres_k(i_below,i) = pres_k(i_below,i) - dp_below_k
          pres_k(i_above,i) = pres_k(i_above,i) + dp_above_k
        ENDIF
        profiles_k(i)%s2m%p = profiles_k(i)%s2m%p + dp_below_k - dp_above_k

      ELSE
        raytracing_k%hgpl(i_below,i) = 0._jprb
      ENDIF

    ENDDO
  END SUBROUTINE calc_hgpl_k


  SUBROUTINE calc_seczen_k(ra,za,path,ra_k,za_k,path_k,do_sat)

    LOGICAL(jplm), INTENT(IN)    :: do_sat ! do satellite secant or solar secant
    REAL(jprb),    INTENT(IN)    :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb),    INTENT(INOUT) :: ra_k(:,:), za_k(:,:), path_k(:,:)

    REAL(jprb)    :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    INTEGER(jpim) :: n, i, prof, last_prof

    IF (opts%rt_ir%addsolar) THEN
      CALL sintosec_k(za(1:nlayers,:), za_k(1:nlayers,:), path, path_k, acc = .TRUE._jplm, map = map, &
        prof_stat = prof_stat)
    ELSE
      CALL sintosec_k(za(1:nlayers,:), za_k(1:nlayers,:), path, path_k, acc = .FALSE._jplm, map = map, &
        prof_stat = prof_stat) ! first use of za
    ENDIF
!   path_k = 0._jprb;

    IF (do_sat) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        ra_k(:,i) = &!ra_k(:,i) +
                    za_k(1:nlayers,i) * angles(prof)%sinzen
      ENDDO
    ELSE ! do_sun
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF (profiles(prof)%sunzenangle >= 0._jprb .AND. &
            profiles(prof)%sunzenangle < max_sol_zen) THEN
          ra_k(:,i) = &!ra_k(:,i) +
                      za_k(1:nlayers,i) * angles(prof)%sinzen_sun
        ELSE
          ra_k(:,i) = 0._jprb
        ENDIF
      ENDDO
    ENDIF
!  za_k = 0._jprb;

    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
      IF (do_sat) THEN
        n = 2
      ELSE ! do sun
        n = nlevels
      ENDIF
      last_prof = -1
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        IF (prof .NE. last_prof) THEN
          z(:,prof) = rearth(prof) * raytracing%z_r(:,prof)
          ztemp(1:nlayers,1) = z(1:nlayers,prof) * raytracing%r(n,prof) * &
                       raytracing%r_r(2:nlevels,prof)**2_jpim
          ztemp(1:nlayers,2) = raytracing%r_r(2:nlevels,prof) * z(1:nlayers,prof)
          ztemp(1:nlayers,3) = raytracing%r_r(2:nlevels,prof) * raytracing%r(n,prof)
          last_prof = prof
        ENDIF

        ! First use but called twice so keep accumulation
        raytracing_k%r(2:nlevels,i) = raytracing_k%r(2:nlevels,i) - &
          ra_k(:,i) * ztemp(1:nlayers,1)

        raytracing_k%r(n,i) = raytracing_k%r(n,i) + SUM(ztemp(1:nlayers,2) * ra_k(:,i))

        ra_k(:,i) = ra_k(:,i) * ztemp(1:nlayers,3)
      ENDDO
    ENDIF

    DO i = 1, nchannels
      prof = chanprof(i)%prof
      IF (do_sat .OR. (profiles(prof)%sunzenangle >= 0._jprb .AND. &
                       profiles(prof)%sunzenangle < max_sol_zen)) THEN
        ! Note accumulation as first call is above
        raytracing_k%z_r(:,i) = &!raytracing_k%z_r(:,i) + &
          rearth(prof) * ra_k(:,i)
      ELSE
        raytracing_k%z_r(:,i) = 0._jprb
      ENDIF
    ENDDO
!   ra_k = 0._jprb;

    IF (opts%config%fix_hgpl) THEN
      DO i = 1, nchannels
        prof = chanprof(i)%prof
        raytracing_k%z_r(auxs(prof)%nearestlev_surf-1,i) = 0._jprb
      ENDDO
    ENDIF
    CALL reciprocal_k(raytracing%z_r, raytracing_k%hgpl(2:nlayers+1, :), &
      raytracing_k%z_r, acc = .TRUE._jplm, map = map, prof_stat = prof_stat)

  END SUBROUTINE calc_seczen_k

END SUBROUTINE rttov_locpat_k
