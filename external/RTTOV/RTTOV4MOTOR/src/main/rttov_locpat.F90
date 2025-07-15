! Description:
!> @file
!!   Calculates information related to the path of the radiation
!!   through the atmosphere.
!!
!> @brief
!!   Calculates information related to the path of the radiation
!!   through the atmosphere.
!!
!! @details
!!   This subroutine is always called on user levels and is also
!!   called on coefficient levels if the interpolator is enabled,
!!   otherwise the contents of raytracing on user levels is copied
!!   into the raytracing structure on coefficient levels in rttov_direct.
!!   Only quantities required on both user and coefficient levels are
!!   calculated here. A small number of additional raytracing path quantities
!!   are computed in rttov_predictor_precalc because they are required on
!!   coefficient levels only.
!!
!!   This subroutine carries out a number of calculations:
!!   - geometric altitude of each pressure level
!!   - secant of the local zenith angle at the lower boundary
!!     of each layer
!!   - secant of local solar zenith angle at lower boundary of each
!!     layer (if solar calculations are enabled)
!!   - layer thickness (only for VIS/IR scattering calculations)
!!
!!   The calculation of the height of the pressure levels is based
!!   on the hydrostatic equation. To account for the presence of
!!   water vapour, virtual temperatures are used. The variation of
!!   gravity with latitude is introduced using the international
!!   gravity formula. The bending of rays as they traverse the atmosphere
!!   is fully accounted for applying the Snell's law. For the computation
!!   of the refractive index of air an updated version of Edlen's formula
!!   is used.
!!
!!     K.P. Birch and M.J.Downs:'An Updated Edlen Equation for the
!!     refractive index of air'. Metrologia, 1993, 30, 155-162.
!!
!!     K.P. Birch and M.J.Downs:'Correction to the Updated Edlen Equation
!!     for the refractive index of air'. Metrologia, 1994, 31, 315-316.
!!
!! @param[in]     opts              RTTOV options structure
!! @param[in]     dosolar           flag indicating whether solar computations are being performed
!! @param[in]     plane_parallel    flag for strict plane parallel geometry
!! @param[in]     profiles          profiles structure (on user levels or coefficient levels)
!! @param[in]     profiles_dry      profiles structure containing gas profiles in units of ppmv dry
!! @param[in]     auxs              RTTOV profile_aux_s structure
!! @param[in]     angles            geometry structure
!! @param[in,out] raytracing        raytracing structure
!! @param[in]     chanprof          channels/profiles to simulate, optional (required for geometric_height)
!! @param[in,out] geometric_height  array for output geometric heights, optional
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
SUBROUTINE rttov_locpat(opts, dosolar, plane_parallel, profiles, profiles_dry, &
                        auxs, angles, raytracing, chanprof, geometric_height)

  USE rttov_types, ONLY : rttov_options, rttov_profile_aux_s, rttov_profile, &
                          rttov_geometry, rttov_raytracing, rttov_chanprof
  USE parkind1, ONLY : jplm, jprb
!INTF_OFF
  USE rttov_const, ONLY : d1, d2, d3, d4, d5, dco2, ed1, ed2, ed3, ed4,&
                          ew1, ew2, htop, ctom, waver, rgc, mair, mh2o, flatt, &
                          eqrad, omega, grave, t0, max_sol_zen, co2_conc

  USE parkind1, ONLY : jpim
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_math_mod, ONLY : divide, sintosec, reciprocal
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),       INTENT(IN)              :: opts
  LOGICAL(jplm),             INTENT(IN)              :: dosolar
  LOGICAL(jplm),             INTENT(IN)              :: plane_parallel
  TYPE(rttov_profile),       INTENT(IN)              :: profiles(:)
  TYPE(rttov_profile),       INTENT(IN)              :: profiles_dry(:)
  TYPE(rttov_profile_aux_s), INTENT(IN)              :: auxs(:)
  TYPE(rttov_geometry),      INTENT(IN)              :: angles(:)
  TYPE(rttov_raytracing),    INTENT(INOUT)           :: raytracing
  TYPE(rttov_chanprof),      INTENT(IN),    OPTIONAL :: chanprof(:)
  REAL(jprb),                INTENT(INOUT), OPTIONAL :: geometric_height(:,:)
!INTF_END

  INTEGER(jpim) :: iprof, lay, i, nlevels, nlayers, nprofiles
  REAL   (jprb) :: rearth(SIZE(profiles))
  REAL   (jprb) :: gravl(SIZE(profiles))   ! gravity at latitude lat [m/s^2]
  REAL   (jprb) :: gravh_r(SIZE(profiles))
  REAL   (jprb) :: eta, beta ! coefficients of the international gravity formula
  REAL   (jprb) :: rlh(SIZE(profiles))
  REAL   (jprb) :: dflat, fac, c, scale_eps, small_val
  REAL   (jprb) :: ZHOOK_HANDLE
!-----end of header-------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT', 0_jpim, ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlevels     = profiles(1)%nlevels
  nlayers     = nlevels - 1

  small_val = (TINY(1._jprb)) ** (0.333333_jprb)

  ! Scale water profile in air density calculation
  scale_eps = (1._jprb - mh2o/mair)

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

  gravh_r(:) = rearth(:) / (2.E-3_jprb * gravl(:)) ! computationally useful form


  ! Calculate moist air density (for later calculation).
  ! Density of dry air is adjusted to account for the presence of water vapour
  ! by replacing temp with virtual temp
  c = 1000._jprb * rgc / (100._jprb * mair)
  DO iprof = 1, nprofiles
    ! Calculate partial pressure of water vapour (convert ppmv dry to ppmv wet)
    raytracing%ppw(:,iprof) = profiles(iprof)%p(:) * 1.E-6_jprb * &
                              profiles_dry(iprof)%q(:) / (1._jprb + profiles_dry(iprof)%q(:) * 1.E-6_jprb)

    CALL divide(profiles(iprof)%p(:) - raytracing%ppw(:,iprof) * scale_eps, &
                c * profiles(iprof)%t(:), &
                raytracing%dmair(:,iprof), acc = .FALSE._jplm)
  ENDDO


  ! Compute height of pressure levels
  CALL calc_hgpl()


  ! Compute local path angles
  IF (plane_parallel) THEN
    DO iprof = 1, nprofiles
      raytracing%zasat(:,iprof)   = angles(iprof)%sinview
      raytracing%pathsat(:,iprof) = angles(iprof)%seczen
    ENDDO

    IF (dosolar) THEN
      DO iprof = 1, nprofiles
        IF (profiles(iprof)%sunzenangle >= 0._jprb .AND. &
            profiles(iprof)%sunzenangle < max_sol_zen) THEN
          raytracing%zasun(:,iprof)   = angles(iprof)%sinzen_sun
          raytracing%pathsun(:,iprof) = 1._jprb / angles(iprof)%coszen_sun
        ELSE
          ! These values are not used for radiance calculations, but should be
          ! assigned to avoid numerical problems on whole-array operations later
          raytracing%zasun(:,iprof)   = 0._jprb
          raytracing%pathsun(:,iprof) = 1._jprb
        ENDIF
      ENDDO
    ENDIF
  ELSE
    ! Compute atmospheric refractive index (if required) for dry air given temperature
    ! and pressure as a function of wavenumber using an updated version of edlen equation
    ! not yet corrected for CO2 contribution - used in calc_seczen
    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) CALL calc_refractivity()

    ! Compute secant of the zenith angle at the lower boundary of each -
    ! layer for the satellite-to-surface line of view and of each
    ! layer for the sun-to-surface line of view (if possible and required)
    CALL calc_seczen(raytracing%ratoesat, raytracing%zasat, raytracing%pathsat,&
                     do_sat=.TRUE._jplm) ! satellite call

    IF (dosolar) THEN
      CALL calc_seczen(raytracing%ratoesun, raytracing%zasun, raytracing%pathsun,&
                       do_sat=.FALSE._jplm) ! solar call
    ENDIF
  ENDIF ! plane-parallel

  IF (dosolar) THEN
    raytracing%patheff = raytracing%pathsat + raytracing%pathsun
  ENDIF
    
  ! Layer thickness
  IF (opts%rt_ir%addsolar .OR. opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    DO iprof = 1, nprofiles
      DO lay = 1, nlayers
        raytracing%ltick(lay,iprof) = raytracing%hgpl(lay,iprof) - raytracing%hgpl(lay + 1,iprof)
      ENDDO
    ENDDO
  ENDIF

  ! Copy altitudes into output geometric_height array (convert km -> m)
  ! Direct model only
  IF (PRESENT(geometric_height)) THEN
    DO i = 1, SIZE(chanprof)
      geometric_height(:,i) = raytracing%hgpl(:,chanprof(i)%prof) * 1.E3_jprb
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_LOCPAT', 1_jpim, ZHOOK_HANDLE)
CONTAINS

  SUBROUTINE calc_refractivity()

    REAL(jprb) :: dry_pp_h(profiles(1)%nlevels)
    REAL(jprb) :: t_norm(profiles(1)%nlevels)

    REAL(jprb) :: disp ! the value of the refractive index given by the dispersion equation.
    REAL(jprb) :: disp_moist
    REAL(jprb) :: dispco2
    INTEGER(jpim) :: iprof

    ! Dry air dispersion equation constant in updated Edlen eqn
    disp = 1.0e-8_jprb * ( &
      d1 + d2 / (d3 - (waver * ctom) ** 2_jpim) + &
           d4 / (d5 - (waver * ctom) ** 2_jpim))

    ! Moist air dispersion constant
    disp_moist = htop * (ew1 - ew2 * (waver * ctom) ** 2_jpim) * 1.E-10_jprb

    ! Compute and store (multiplicative) correction factor to account for presence of CO2 in the air
    ! Assume constant CO2 using figure from 2005 (376ppmv) unless co2 data supplied
    ! Calculation assumes ppmv over dry air
    IF (opts%rt_all%co2_data) THEN
      DO iprof = 1, nprofiles
        raytracing%dispco2(:,iprof) = disp * &
          (1._jprb + dco2 * (profiles_dry(iprof)%co2(:) * 1.E-6_jprb - 0.0003_jprb))
      ENDDO
    ELSE
      dispco2 = disp * (1._jprb + dco2 * (co2_conc * 1.E-6_jprb - 0.0003_jprb))
    ENDIF

    DO iprof = 1, nprofiles
      ! Calculate useful temporary array to break refractivity calculation up

      dry_pp_h(:) = htop * (profiles(iprof)%p(:) - raytracing%ppw(:,iprof))
      t_norm(:) = profiles(iprof)%t(:) - t0

      raytracing%refractivity(:,iprof) = (dry_pp_h(:) / ed1) *         &
        (1._jprb + 1.E-8_jprb * &
        (ed2 - ed3 * t_norm(:)) * dry_pp_h(:)) /  &
        (1._jprb + ed4 * t_norm(:))

      ! Apply CO2 correction then compute moist contribution to refractivity
      ! and store final refractive index for atmospheric profile
      IF (opts%rt_all%co2_data) THEN
        raytracing%r(:,iprof) = (raytracing%refractivity(:,iprof) * &
          raytracing%dispco2(:,iprof)) - &
          (raytracing%ppw(:,iprof) * disp_moist) + 1._jprb
      ELSE
        raytracing%r(:,iprof) = raytracing%refractivity(:,iprof) * dispco2 - &
          (raytracing%ppw(:,iprof) * disp_moist) + 1._jprb
      ENDIF
    ENDDO

    ! Store reciprocal of refractive index for performance
    CALL reciprocal(raytracing%r, raytracing%r_r)
  END SUBROUTINE calc_refractivity


  SUBROUTINE calc_hgpl()

    REAL(jprb)    :: dp(profiles(1)%nlevels)
    INTEGER(jpim) :: i_above, i_below
    REAL(jprb)    :: dp_above, dp_below, t_surf, q_surf, ppw_surf, dmair_surf
    INTEGER(jpim) :: iprof, lev

    !-------Compute height of pressure levels:levels above surface. ------------
    !       The height of pressure levels H is obtained by integrating             |
    !       the hydrostatic equation dPRES=-GRAVL(H)*DMAIR*dH between              |
    !       two adjacent pressure levels                                           |
    !                                                                              |
    !           -P2         - H2                                                   |
    !          |           |                                                       |
    !          | dP/DMAIR= | GRAVL(H)*dH                                           |
    !          |           |                                                       |
    !         -  P1       -   H1                                                   |
    !                                                                              |
    !       The integration of 1/DMAIR is carried out assuming DMAIR               |
    !       varies linearly with pressure. The integral is then                    |
    !       computed analytically.                                                 |
    !                                                                              |
    !       The value of the gravity as a function of altitude H can be            |
    !       expressed using the inverse-square law of gravitation:                 |
    !                                                                              |
    !       GRAVL(H)=GRAVL*(REARTH/(H+REARTH))**2=                                 |
    !                GRAVL*(1-2*H/REARTH+3*H**2/REARTH**2+terms of higher order)   |
    !                                                                              |
    !       If we eliminate the second and higher order terms we can write:        |
    !                                                                              |
    !       GRAVL(H)=GRAVL-2*GRAVL*H/REARTH=GRAVL-GRAVH*H                          |
    !                                                                              |
    !       Note that RLH = GRAVL / GRAVH in the equations below                   |
    !------------------------------------------------------------------------------

    ! Compute height of pressure levels above surface
    ! Lower boundary of first layer is nearest level below surface
    ! Upper boundary of last layer is TOA
    DO iprof = 1, nprofiles

      ! Set the height of the surface
      i_below     = auxs(iprof)%nearestlev_surf
      i_above     = auxs(iprof)%nearestlev_surf - 1
      dp_above    = profiles(iprof)%p(i_above) - profiles(iprof)%s2m%p
      dp_below    = profiles(iprof)%s2m%p - profiles(iprof)%p(i_below)

      IF (dp_below /= 0._jprb .AND. opts%config%fix_hgpl) THEN
        q_surf      = (profiles_dry(iprof)%q(i_below) * dp_above + &
                       profiles_dry(iprof)%q(i_above) * dp_below) / (dp_above + dp_below)
        t_surf      = (profiles(iprof)%t(i_below) * dp_above + &
                       profiles(iprof)%t(i_above) * dp_below) / (dp_above + dp_below)
        ppw_surf    = profiles(iprof)%s2m%p * 1.E-6_jprb * q_surf / (1._jprb + q_surf * 1.E-6_jprb)
        dmair_surf  = (profiles(iprof)%s2m%p - ppw_surf * scale_eps) / (c * t_surf)

        IF (raytracing%dmair(i_below,iprof) == dmair_surf) THEN
          ! If dmair constant across layer, by L'Hopital's rule: int = dp / dmair
          raytracing%int(i_below,iprof) = dp_below / dmair_surf
        ELSE
          raytracing%int(i_below,iprof) = &
            dp_below / (raytracing%dmair(i_below,iprof) - dmair_surf) * &
            LOG(raytracing%dmair(i_below,iprof) / dmair_surf)
        ENDIF

        raytracing%ztemp(i_below,iprof) = &
          rlh(iprof)**2_jpim - 1.E3_jprb * profiles(iprof)%elevation * &
           (2._jprb * rlh(iprof) - 1.E3_jprb * profiles(iprof)%elevation) - &
            2._jprb * raytracing%int(i_below,iprof) * 1.E2_jprb * gravh_r(iprof)

        raytracing%hgpl(i_below,iprof) = &
          (rlh(iprof) - SQRT(raytracing%ztemp(i_below,iprof))) * 1.E-3_jprb
      ELSE
        raytracing%hgpl(i_below,iprof) = profiles(iprof)%elevation

        ! Not used, but must be non-zero for TL/AD/K:
        raytracing%ztemp(i_below,iprof) = 1._jprb
      ENDIF

      dp(1:nlevels-1) = profiles(iprof)%p(1:nlevels-1) - profiles(iprof)%p(2:nlevels)

      DO lev = auxs(iprof)%nearestlev_surf - 1, 1, -1

        IF (raytracing%dmair(lev+1,iprof) == raytracing%dmair(lev,iprof)) THEN
          ! If dmair constant across layer, by L'Hopital's rule: int = dp / dmair
          raytracing%int(lev,iprof) = & ! integrated layer values for dmair
            -dp(lev) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing%int(lev,iprof) = & ! integrated layer values for dmair
            -dp(lev) / (raytracing%dmair(lev+1,iprof) - raytracing%dmair(lev,iprof)) * &
            LOG(raytracing%dmair(lev+1,iprof) / raytracing%dmair(lev,iprof))
        ENDIF

        raytracing%ztemp(lev,iprof) = &
          rlh(iprof)**2_jpim - 1.E3_jprb * raytracing%hgpl(lev + 1,iprof) * &
          (2._jprb * rlh(iprof) - 1.E3_jprb * raytracing%hgpl(lev + 1,iprof)) - &
           2._jprb * raytracing%int(lev,iprof) * 1.E2_jprb * gravh_r(iprof)

        IF (raytracing%ztemp(lev,iprof) > 0._jprb) THEN
          raytracing%hgpl(lev,iprof) = &
            (rlh(iprof) - SQRT(raytracing%ztemp(lev,iprof))) * 1.E-3_jprb
        ELSE
          raytracing%ztemp(lev,iprof) = small_val
          raytracing%hgpl(lev,iprof) = rlh(iprof) * 1.E-3_jprb
        ENDIF
      ENDDO

      ! Compute height of pressure levels below surface
      ! Upper boundary of first layer is nearest level below surface
      ! Lower boundary of last layer is bottom of profile grid
      DO lev = auxs(iprof)%nearestlev_surf + 1, nlevels

        IF (raytracing%dmair(lev,iprof) == raytracing%dmair(lev-1,iprof)) THEN
          ! If dmair constant across layer, by L'Hopital's rule: int = dp / dmair
          raytracing%int(lev,iprof) = & ! integrated layer values for dmair
            dp(lev-1) / raytracing%dmair(lev,iprof)
        ELSE
          raytracing%int(lev,iprof) = & ! integrated layer values for dmair
            dp(lev-1) / (raytracing%dmair(lev,iprof) - raytracing%dmair(lev-1,iprof)) * &
            LOG(raytracing%dmair(lev,iprof) / raytracing%dmair(lev-1,iprof))
        ENDIF

        raytracing%ztemp(lev,iprof) = &
          rlh(iprof) ** 2_jpim - 1.E3_jprb * raytracing%hgpl(lev - 1,iprof) * &
          (2._jprb * rlh(iprof) - 1.E3_jprb * raytracing%hgpl(lev - 1,iprof)) - &
          2._jprb * raytracing%int(lev,iprof) * 1.E2_jprb * gravh_r(iprof)

        raytracing%hgpl(lev,iprof) = &
          (rlh(iprof) - SQRT(raytracing%ztemp(lev,iprof))) * 1.E-3_jprb
      ENDDO

    ENDDO

  END SUBROUTINE calc_hgpl


  SUBROUTINE calc_seczen(ra,za,path,do_sat)
    ! Compute secant of the zenith angle at the lower boundary of each -
    ! layer for the sat/sun-to-surface line of view.
    ! The atmospheric layers are considered as concentric rings. If we trace
    ! a ray across these rings at any angle other than nadir, the local
    ! angle relative to the outward radial direction at the point of
    ! intersection will be different at each ring because due to the
    ! curvature of the Earth and to atmospheric refraction. The secant of
    ! the local/solar zenith angle PATHSAT/SUN is thus computed taking into
    ! account the geometry of the situation and the bending of rays as they
    ! traverse the atmosphere (by application of Snell's law).

    ! In the layer containing the surface the local path angle matches the
    ! input zenith angle if the elevation is zero and addrefrac is false.

    REAL(jprb),    INTENT(OUT) :: ra(:,:), za(:,:), path(:,:)
    LOGICAL(jplm), INTENT(IN)  :: do_sat ! do satellite secant or solar secant

    INTEGER(jpim) :: n, iprof

    ! Calculate and store pressure level reciprocal distance from earth centre
    DO iprof = 1, nprofiles
      CALL reciprocal(rearth(iprof) + raytracing%hgpl(2:nlayers+1,iprof), raytracing%z_r(:,iprof))
    ENDDO
    IF (opts%config%fix_hgpl) THEN
      ! Ensure that the local path angle is consistent with the user-
      ! specified zenangle in the layer containing the surface
      DO iprof = 1, nprofiles
        raytracing%z_r(auxs(iprof)%nearestlev_surf-1,iprof) = &
          1._jprb / (rearth(iprof) + profiles(iprof)%elevation)
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      ra(:,iprof) = rearth(iprof) * raytracing%z_r(:,iprof)
    ENDDO

    IF (opts%rt_all%addrefrac .OR. opts%rt_ir%pc%addpc) THEN
      IF (do_sat) THEN
        n = 2
      ELSE ! do sun
        n = nlevels
      ENDIF
      DO iprof = 1, nprofiles
        ra(:,iprof) = ra(:,iprof) * raytracing%r(n,iprof) * raytracing%r_r(2:nlevels,iprof)
      ENDDO
    ENDIF

    DO iprof = 1, nprofiles
      IF (do_sat) THEN
        za(1:nlayers,iprof) = ra(:,iprof) * angles(iprof)%sinzen
      ELSE ! do_sun
        IF (profiles(iprof)%sunzenangle >= 0._jprb .AND. &
            profiles(iprof)%sunzenangle < max_sol_zen) THEN
          za(1:nlayers,iprof) = ra(:,iprof) * angles(iprof)%sinzen_sun
        ELSE
          ! This value is not used for radiance calculations, but should be
          ! assigned to avoid numerical problems on whole-array operations later
          za(1:nlayers,iprof) = 0._jprb
        ENDIF
      ENDIF
    ENDDO

    ! Calculate path from zenith angle using path = 1/sqrt(1-za**2)
    CALL sintosec(za(1:nlayers,:), path)

  END SUBROUTINE calc_seczen

END SUBROUTINE rttov_locpat
