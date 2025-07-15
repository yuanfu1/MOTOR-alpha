! Description:
!> @file
!!   TL of atmospheric radiation path calculations.
!!
!> @brief
!!   TL of atmospheric radiation path calculations.
!!
!! @param[in]     opts            RTTOV options structure
!! @param[in]     dosolar         flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel  flag for strict plane parallel geometry
!! @param[in]     profiles        profiles structure (on user levels or coefficient levels)
!! @param[in]     profiles_tl     profiles structure containing perturbations
!! @param[in]     profiles_dry    profiles structure containing gas profiles in units of ppmv dry
!! @param[in]     profiles_dry_tl profiles structure containing gas perturbations in units of ppmv dry
!! @param[in]     auxs            RTTOV profile_aux_s structure
!! @param[in]     angles          geometry structure
!! @param[in]     raytracing      raytracing structure
!! @param[in,out] raytracing_tl   raytracing structure containing perturbations
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
SUBROUTINE rttov_locpat_tl(opts, dosolar, plane_parallel, profiles, profiles_tl, &
      profiles_dry, profiles_dry_tl, auxs, angles, raytracing, raytracing_tl)

  USE rttov_types, ONLY : rttov_options, rttov_profile_aux_s, rttov_profile, &
                          rttov_geometry, rttov_raytracing
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4,&
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt,  &
                          eqrad, omega, grave, t0, max_sol_zen, co2_conc

  USE parkind1, ONLY : jpim, jprb
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_math_mod, ONLY : divide, invsqrt, sintosec_tl, reciprocal, reciprocal_tl
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),       INTENT(IN)    :: opts
  LOGICAL(jplm),             INTENT(IN)    :: dosolar
  LOGICAL(jplm),             INTENT(IN)    :: plane_parallel
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_dry(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles_dry_tl(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)    :: auxs(:)
  TYPE(rttov_geometry),      INTENT(IN)    :: angles(:)
  TYPE(rttov_raytracing),    INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),    INTENT(INOUT) :: raytracing_tl
!INTF_END

  INTEGER(jpim) :: iprof, lay, nlevels , nlayers, nprofiles
  REAL   (jprb) :: rearth(SIZE(profiles))
  REAL   (jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (jprb) :: gravh_r(SIZE(profiles))
  REAL   (jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (jprb) :: rlh(SIZE(profiles))
  REAL   (jprb) :: dflat, fac, c, scale_eps
  REAL   (jprb) :: pres_tl(profiles(1)%nlevels,SIZE(profiles))
  REAL   (jprb) :: qwet(profiles(1)%nlevels)
  REAL   (jprb) :: qwet_tl(profiles(1)%nlevels)
  REAL   (jprb) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_TL', 0_jpim, ZHOOK_HANDLE)

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


  ! Calculate moist air density (for later calculation).
  ! Density of dry air is adjusted to account for the presence of water vapour
  ! by replacing temp with virtual temp
  c = 1000._jprb * rgc / (100._jprb * mair)
  DO iprof = 1, nprofiles

    ! Calculate water vapour in units of ppmv wet
    qwet(:) = profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1.E-6_jprb)
    qwet_tl(:) = profiles_dry_tl(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1.E-6_jprb)**2

    IF (opts%interpolation%lgradp) THEN
      pres_tl(:,iprof) = profiles_tl(iprof)%p(:)
      raytracing_tl%ppw(:,iprof) = 1.E-6_jprb * &
           (pres_tl(:,iprof) * qwet(:) + &
            profiles(iprof)%p(:) * qwet_tl(:))

      CALL divide(pres_tl(:,iprof) - raytracing_tl%ppw(:,iprof) * scale_eps - &
           c * profiles_tl(iprof)%t(:) * raytracing%dmair(:,iprof), &
           c * profiles(iprof)%t(:), &
           raytracing_tl%dmair(:,iprof), acc = .FALSE._jplm)

    ELSE
      raytracing_tl%ppw(:,iprof) = 1.E-6_jprb * &
           profiles(iprof)%p(:) * qwet_tl(:)

      CALL divide(-raytracing_tl%ppw(:,iprof) * scale_eps - &
           c * profiles_tl(iprof)%t(:) * raytracing%dmair(:,iprof), &
           c * profiles(iprof)%t(:), &
           raytracing_tl%dmair(:,iprof), acc = .FALSE._jplm)
    ENDIF
  ENDDO


  ! Compute height of pressure levels
  CALL calc_hgpl_tl()


  ! Compute local path angles
  IF (plane_parallel) THEN
    raytracing_tl%zasat         = 0._jprb
    raytracing_tl%pathsat       = 0._jprb

    IF (dosolar) THEN
      raytracing_tl%zasun         = 0._jprb
      raytracing_tl%pathsun       = 0._jprb
      raytracing_tl%patheff       = 0._jprb
    ENDIF
  ELSE
    ! Compute atmospheric refractive index (if required) for dry air given temperature
    ! and pressure as a function of wavenumber using an updated version of edlen equation
    ! not yet corrected for CO2 contribution - used in calc_seczen
    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) CALL calc_refractivity_tl()

    ! Compute secant of the zenith angle at the lower boundary of each -
    ! layer for the satellite-to-surface line of view and of each
    ! layer for the sun-to-surface line of view (if possible and required)
    CALL calc_seczen_tl(&
      raytracing%ratoesat, raytracing%zasat, raytracing%pathsat, &
      raytracing_tl%ratoesat, raytracing_tl%zasat, raytracing_tl%pathsat, &
      do_sat=.TRUE._jplm) ! satellite call

    IF (dosolar) THEN
      CALL calc_seczen_tl(&
        raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
        raytracing_tl%ratoesun, raytracing_tl%zasun, raytracing_tl%pathsun,&
        do_sat=.FALSE._jplm) ! solar call

      raytracing_tl%patheff = raytracing_tl%pathsat + raytracing_tl%pathsun
    ENDIF
  ENDIF ! plane-parallel


  ! Layer thickness
  IF (opts%rt_ir%addsolar .OR. opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO iprof = 1, nprofiles
      DO lay = 1, nlayers
        raytracing_tl%ltick(lay,iprof) = raytracing_tl%hgpl(lay,iprof) - raytracing_tl%hgpl(lay + 1,iprof)
      ENDDO
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT_TL', 1_jpim, ZHOOK_HANDLE)
CONTAINS

  SUBROUTINE calc_refractivity_tl()

    REAL(jprb) :: dry_pp_h((profiles(1)%nlevels))
    REAL(jprb) :: t_norm((profiles(1)%nlevels))
    REAL(jprb) :: dry_pp_h_tl((profiles(1)%nlevels))
    REAL(jprb) :: U((profiles(1)%nlevels)), V_R((profiles(1)%nlevels))
    REAL(jprb) :: DU((profiles(1)%nlevels)), DV((profiles(1)%nlevels))

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

    ! Compute and store (multiplicative) correction factor to account for presence of CO2 in the air
    ! Assume constant CO2 using figure from 2005 (376ppmv) unless co2 data supplied
    ! Calculation assumes ppmv over dry air
    IF (opts%rt_all%co2_data) THEN
      DO iprof = 1, nprofiles
        raytracing_tl%dispco2(:,iprof) = disp * &
          dco2 * profiles_dry_tl(iprof)%co2(:) * 1.E-6_jprb
      ENDDO
    ELSE
      dispco2 = disp * (1._jprb + dco2 * (co2_conc * 1.E-6_jprb - 0.0003_jprb))
    ENDIF

    DO iprof = 1, nprofiles
      ! Calculate useful temporary array to break refractivity calculation up
      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:,iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0
      U = 1._jprb + 1.E-8_jprb * (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)
      CALL reciprocal(1._jprb + ed4 * t_norm(:), V_R)

      IF (opts%interpolation%lgradp) THEN
        dry_pp_h_tl(:) = htop * (pres_tl(:,iprof) - raytracing_tl%ppw(:,iprof))
      ELSE
        dry_pp_h_tl(:) = -htop * raytracing_tl%ppw(:,iprof)
      ENDIF

      DU = 1.E-8_jprb * &
           (dry_pp_h(:) * (-ed3) * profiles_tl(iprof)%t(:) + &
           (ed2 - ed3 * t_norm(:)) * dry_pp_h_tl(:))

      DV = ed4 * profiles_tl(iprof)%t(:)

      raytracing_tl%refractivity(:,iprof) = dry_pp_h_tl(:) * (U * V_R)/ed1 + &
        (dry_pp_h(:) * DU/ed1 - raytracing%refractivity(:,iprof) * DV) * V_R

      ! Apply CO2 correction then compute moist contribution to refractivity
      ! and store final refractive index for atmospheric profile
      IF (opts%rt_all%co2_data) THEN
        raytracing_tl%r(:,iprof) = & !refractivity back to refractance
          raytracing_tl%refractivity(:,iprof) * raytracing%dispco2(:,iprof) +&
          (raytracing%refractivity(:,iprof)) * raytracing_tl%dispco2(:,iprof) - &
          raytracing_tl%ppw(:,iprof) * disp_moist
      ELSE
        raytracing_tl%r(:,iprof) = raytracing_tl%refractivity(:,iprof) * dispco2 - &
          raytracing_tl%ppw(:,iprof) * disp_moist
      ENDIF
    ENDDO

  END SUBROUTINE calc_refractivity_tl


  SUBROUTINE calc_hgpl_tl()
    REAL(jprb)    :: dp(profiles(1)%nlevels)
    REAL(jprb)    :: dp_tl(profiles(1)%nlevels)
    REAL(jprb)    :: ztemp(profiles(1)%nlevels,SIZE(profiles))
    INTEGER(jpim) :: i_above, i_below
    REAL(jprb)    :: dp_above, dp_below, t_surf, q_surf, ppw_surf, dmair_surf
    REAL(jprb)    :: dp_above_tl, dp_below_tl, t_surf_tl, q_surf_tl, ppw_surf_tl, dmair_surf_tl
    INTEGER(jpim) :: iprof, lev

    ! Store ztemp from direct temporarily.
    CALL invsqrt(raytracing%ztemp(:,:), ztemp(:,:))

    DO iprof = 1, nprofiles

      ! Set the height of the surface
      i_below     = auxs(iprof)%nearestlev_surf
      i_above     = auxs(iprof)%nearestlev_surf - 1
      dp_above    = profiles(iprof)%p(i_above) - profiles(iprof)%s2m%p
      dp_below    = profiles(iprof)%s2m%p - profiles(iprof)%p(i_below)

      IF (dp_below /= 0._jprb .AND. opts%config%fix_hgpl) THEN

        IF (opts%interpolation%lgradp) THEN
          dp_above_tl = profiles_tl(iprof)%p(i_above) - profiles_tl(iprof)%s2m%p
          dp_below_tl = profiles_tl(iprof)%s2m%p - profiles_tl(iprof)%p(i_below)
        ELSE
          dp_above_tl = -profiles_tl(iprof)%s2m%p
          dp_below_tl = profiles_tl(iprof)%s2m%p
        ENDIF

        q_surf        = (profiles_dry(iprof)%q(i_below) * dp_above + &
                         profiles_dry(iprof)%q(i_above) * dp_below) / (dp_above + dp_below)
        t_surf        = (profiles(iprof)%t(i_below) * dp_above + &
                         profiles(iprof)%t(i_above) * dp_below) / (dp_above + dp_below)
        ppw_surf      = profiles(iprof)%s2m%p * 1.E-6_jprb * q_surf / (1._jprb + q_surf * 1.E-6_jprb)
        dmair_surf    = (profiles(iprof)%s2m%p - ppw_surf * scale_eps) / (c * t_surf)

        q_surf_tl     = (profiles_dry_tl(iprof)%q(i_below) * dp_above + &
                         profiles_dry(iprof)%q(i_below) * dp_above_tl + &
                         profiles_dry_tl(iprof)%q(i_above) * dp_below + &
                         profiles_dry(iprof)%q(i_above) * dp_below_tl - &
                         (dp_above_tl + dp_below_tl) * q_surf) / (dp_above + dp_below)
        t_surf_tl     = (profiles_tl(iprof)%t(i_below) * dp_above + &
                         profiles(iprof)%t(i_below) * dp_above_tl + &
                         profiles_tl(iprof)%t(i_above) * dp_below + &
                         profiles(iprof)%t(i_above) * dp_below_tl - &
                         (dp_above_tl + dp_below_tl) * t_surf) / (dp_above + dp_below)
        ppw_surf_tl   = (profiles_tl(iprof)%s2m%p * q_surf + &
                         q_surf_tl * (profiles(iprof)%s2m%p - ppw_surf)) * &
                        1.E-6_jprb / (1._jprb + q_surf * 1.E-6_jprb)
        dmair_surf_tl = (profiles_tl(iprof)%s2m%p - ppw_surf_tl * scale_eps) / (c * t_surf) - &
                        t_surf_tl * dmair_surf / t_surf

        IF (raytracing%dmair(i_below,iprof) == dmair_surf) THEN
          raytracing_tl%int(i_below,iprof) = &
            dp_below_tl / dmair_surf - &
            dmair_surf_tl * raytracing%int(i_below,iprof) / dmair_surf
        ELSE
          raytracing_tl%int(i_below,iprof) = &
            dp_below_tl * raytracing%int(i_below,iprof) / dp_below + &
            (dp_below * &
            (raytracing_tl%dmair(i_below,iprof) / raytracing%dmair(i_below,iprof) - &
             dmair_surf_tl / dmair_surf) - &
            raytracing%int(i_below,iprof) * &
            (raytracing_tl%dmair(i_below,iprof) - dmair_surf_tl)) / &
            (raytracing%dmair(i_below,iprof) - dmair_surf)
        ENDIF

        raytracing_tl%ztemp(i_below,iprof) = &
          -2._jprb * raytracing_tl%int(i_below,iprof) * 1.E2_jprb * gravh_r(iprof)

        raytracing_tl%hgpl(i_below,iprof) = &
          -0.5_jprb * ztemp(i_below,iprof) * 1.E-3_jprb * raytracing_tl%ztemp(i_below,iprof)
      ELSE
        raytracing_tl%hgpl(i_below,iprof) = 0._jprb
      ENDIF

      dp(1:nlevels-1) = profiles(iprof)%p(1:nlevels-1) - profiles(iprof)%p(2:nlevels)

      IF (opts%interpolation%lgradp) THEN
        dp_tl(1:nlevels-1) = profiles_tl(iprof)%p(1:nlevels-1) - profiles_tl(iprof)%p(2:nlevels)
      ENDIF

      DO lev = auxs(iprof)%nearestlev_surf - 1, 1, -1

        IF (raytracing%dmair(lev+1,iprof) == raytracing%dmair(lev,iprof)) THEN
          raytracing_tl%int(lev,iprof) = &
            -raytracing_tl%dmair(lev,iprof) * &
            raytracing%int(lev,iprof) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing_tl%int(lev,iprof) = &
            (dp(lev) * &
            (raytracing_tl%dmair(lev,iprof) / raytracing%dmair(lev,iprof) - &
             raytracing_tl%dmair(lev+1,iprof) / raytracing%dmair(lev+1,iprof)) - &
            raytracing%int(lev,iprof) * &
            (raytracing_tl%dmair(lev+1,iprof) - raytracing_tl%dmair(lev,iprof))) / &
            (raytracing%dmair(lev+1,iprof) - raytracing%dmair(lev,iprof))
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          raytracing_tl%int(lev,iprof) = &
            raytracing_tl%int(lev,iprof) + &
            dp_tl(lev) * raytracing%int(lev,iprof) / dp(lev)
        ENDIF

        raytracing_tl%ztemp(lev,iprof) = &
          2._jprb * 1.E3_jprb * raytracing_tl%hgpl(lev + 1,iprof) * &
          (1.E3_jprb * raytracing%hgpl(lev + 1,iprof) - rlh(iprof)) - &
          2._jprb * 100._jprb * gravh_r(iprof) * raytracing_tl%int(lev,iprof)

        raytracing_tl%hgpl(lev,iprof) = &
          -0.5_jprb * ztemp(lev,iprof) * 1.E-3_jprb * raytracing_tl%ztemp(lev,iprof)
      ENDDO

      DO lev = auxs(iprof)%nearestlev_surf + 1, nlevels

        IF (raytracing%dmair(lev,iprof) == raytracing%dmair(lev-1,iprof)) THEN
          raytracing_tl%int(lev,iprof) = &
            -raytracing_tl%dmair(lev,iprof) * &
            raytracing%int(lev,iprof) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing_tl%int(lev,iprof) = &
            (dp(lev-1) * &
            (raytracing_tl%dmair(lev-1,iprof) / raytracing%dmair(lev-1,iprof) - &
             raytracing_tl%dmair(lev,iprof) / raytracing%dmair(lev,iprof)) - &
            raytracing%int(lev,iprof) * &
            (raytracing_tl%dmair(lev,iprof) - raytracing_tl%dmair(lev-1,iprof))) / &
            (raytracing%dmair(lev,iprof) - raytracing%dmair(lev-1,iprof))
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          raytracing_tl%int(lev,iprof) = &
            raytracing_tl%int(lev,iprof) + &
            dp_tl(lev-1) * raytracing%int(lev,iprof) / dp(lev-1)
        ENDIF

        raytracing_tl%ztemp(lev,iprof) = &
          2._jprb * 1.E3_jprb * raytracing_tl%hgpl(lev - 1,iprof) * &
          (1.E3_jprb * raytracing%hgpl(lev - 1,iprof) - rlh(iprof)) - &
          2._jprb * 100._jprb * gravh_r(iprof) * raytracing_tl%int(lev,iprof)

        raytracing_tl%hgpl(lev,iprof) = &
          -0.5_jprb * ztemp(lev,iprof) * 1.E-3_jprb * raytracing_tl%ztemp(lev,iprof)
      ENDDO
    ENDDO

  END SUBROUTINE calc_hgpl_tl


  SUBROUTINE calc_seczen_tl(ra,za,path,ra_tl,za_tl,path_tl,do_sat)

    LOGICAL(jplm), INTENT(IN)  :: do_sat ! do satellite secant or solar secant
    REAL(jprb),    INTENT(IN)  :: ra(:,:), za(:,:), path(:,:)
    REAL(jprb),    INTENT(OUT) :: ra_tl(:,:), za_tl(:,:), path_tl(:,:)

    REAL(jprb)    :: z(SIZE(ra(:,1)),SIZE(ra(1,:)))
    INTEGER(jpim) :: n, iprof

    CALL reciprocal_tl(raytracing%z_r,raytracing_tl%hgpl(2:nlayers+1, :),&
         raytracing_tl%z_r)
    IF (opts%config%fix_hgpl) THEN
      DO iprof = 1, nprofiles
        raytracing_tl%z_r(auxs(iprof)%nearestlev_surf-1,iprof) = 0._jprb
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      ! z is temporary array because we overwrote ra
      z(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
      ra_tl(:,iprof) = rearth(iprof) * raytracing_tl%z_r(:,iprof)
    ENDDO

    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
      IF (do_sat) THEN
        n = 2
      ELSE ! do sun
        n = nlevels
      ENDIF
      DO iprof = 1, nprofiles
        ra_tl(:,iprof) = &
          (z(:,iprof) * raytracing_tl%r(n,iprof)  - &
          z(:,iprof) * raytracing%r(n,iprof) * raytracing_tl%r(2:nlevels,iprof) * raytracing%r_r(2:nlevels,iprof) + &
          ra_tl(:,iprof) * raytracing%r(n,iprof)) * raytracing%r_r(2:nlevels,iprof)
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      IF (do_sat) THEN
        za_tl(1:nlayers,iprof) = ra_tl(:,iprof) * angles(iprof)%sinzen
      ELSE ! do_sun
        IF (profiles(iprof)%sunzenangle >= 0._jprb .AND. &
          profiles(iprof)%sunzenangle < max_sol_zen) THEN
          za_tl(1:nlayers,iprof) = ra_tl(:,iprof) * angles(iprof)%sinzen_sun
        ELSE
          za_tl(1:nlayers,iprof) = 0._jprb
        ENDIF
      ENDIF
    ENDDO

    ! Calculate path from zenith angle using path = 1/sqrt(1-za**2)
    CALL sintosec_tl(za(1:nlayers,:), za_tl(1:nlayers,:), path, path_tl)

  END SUBROUTINE calc_seczen_tl

END SUBROUTINE rttov_locpat_tl
