! Description:
!> @file
!!   Calculate various optical depth and transmittance quantities for solar radiation
!
!> @brief
!!   Calculate various optical depth and transmittance quantities for solar radiation
!!
!! @details
!!   Transmittances for solar radiation are required on "path2" which is the
!!   combined sun-surface-satellite path. These are calculated using the
!!   effective path length (pathsat+pathsun). However some variables are
!!   actually calculated along the sun-surface path in particular
!!   solar_path2\%od_singlelayer and solar_path2\%od_sfrac.
!!
!!   All solar_path1 variables are along the surface-satellite path, but note
!!   that these are derived explicitly from the solar predictor calculations
!!   (on path2) so that they are consistent with the solar_path2 variables.
!!   Therefore these will always be slightly different to equivalent the
!!   thermal_path1 quantities.
!!
!!
!! @param[in]     opts                        RTTOV options structure
!! @param[in]     nlayers                     number of layers in input profile
!! @param[in]     chanprof                    specifies channels and profiles to simulate
!! @param[in]     solar                       flags indicating which channels have solar contribution
!! @param[in]     aux                         auxiliary profile variables
!! @param[in]     coef                        optical depth coefficients structure
!! @param[in]     raytracing                  raytracing structure
!! @param[in]     ircld                       computed cloud column data
!! @param[in]     opdp_path                   optical depths calculated by rttov_opdep
!! @param[in,out] path2                       optical depths and transmittances on "path2"
!! @param[in,out] path1                       optical depths and transmittances on "path1"
!! @param[in,out] transmission                the output transmission structure
!! @param[in,out] transmission_aux            top-level auxiliary transmission structure
!! @param[in]     transmission_scatt_ir       cloud/aerosol layer optical depths
!! @param[in]     transmission_scatt_ir_dyn   cloud/aerosol level-to-space optical depths
!!
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
SUBROUTINE rttov_transmit_solar( &
              opts,                         &
              nlayers,                      &
              chanprof,                     &
              solar,                        &
              aux,                          &
              coef,                         &
              raytracing,                   &
              ircld,                        &
              opdp_path,                    &
              path2,                        &
              path1,                        &
              transmission,                 &
              transmission_aux,             &
              transmission_scatt_ir,        &
              transmission_scatt_ir_dyn)

  USE rttov_types, ONLY :  &
       rttov_options,               &
       rttov_chanprof,              &
       rttov_coef,                  &
       rttov_opdp_path,             &
       rttov_path_traj_trans,       &
       rttov_transmission,          &
       rttov_transmission_aux,      &
       rttov_transmission_scatt_ir, &
       rttov_profile_aux,           &
       rttov_ircld,                 &
       rttov_raytracing
  USE parkind1, ONLY : jpim,jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : max_optical_depth, vis_scatt_single, min_tau
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_opdp_path),             INTENT(IN)    :: opdp_path
  TYPE(rttov_path_traj_trans),       INTENT(INOUT) :: path2
  TYPE(rttov_path_traj_trans),       INTENT(INOUT) :: path1
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission
  TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
!INTF_END

  REAL   (KIND=jprb) :: od_surf_ac(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))                 ! sat to surface optical depth (path2)
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, col
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: levsurf
  INTEGER(KIND=jpim) :: nchanprof
  LOGICAL(KIND=jplm) :: do_scatt, do_single_scatt
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_single_scatt = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single

!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  DO j = 1, nchanprof
    IF (solar(j)) THEN
      chan = chanprof(j)%chan
      DO lay = 1, nlayers
        path2%od_singlelayer(lay,j) = -coef%ff_gam(chan) * &
          (opdp_path%sun_level_path2(lay+1,j) - opdp_path%sun_level_path2(lay,j))
      ENDDO
      path2%od_level(:,j) = MAX(coef%ff_gam(chan) * opdp_path%sun_level_path2(:,j), -max_optical_depth)
      path2%tau_level(:,j) = EXP(path2%od_level(:,j))
    ENDIF
  ENDDO

  ! solar_path1 transmittances: the top-level-to-space transmittance is set to 1.0
  DO j = 1, nchanprof
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      path1%tau_level(1,j) = 1.0_jprb
      DO lev = 2, nlevels
        lay = lev - 1
        path1%tau_level(lev,j) = EXP(path2%od_level(lev,j) * &
            raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof))
      ENDDO
    ENDIF
  ENDDO

!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchanprof
    chan    = chanprof(j)%chan
    prof    = chanprof(j)%prof
! as defined in rttov_calc_nearest_lev
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (solar(j)) THEN
      od_surf(j) = path2%od_level(levsurf,j) + &
          aux%s(prof)%pfraction_surf * (path2%od_level(levsurf-1,j) - path2%od_level(levsurf,j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        path2%od_frac(j) = od_surf(j) - path2%od_level(levsurf-1,j)
      ELSE
        path2%od_frac(j) = od_surf(j) - path2%od_level(levsurf,j)
      ENDIF
      path2%tau_surf(j)     = EXP(od_surf(j))
      path1%tau_surf(j)     = EXP(od_surf(j) * &
          raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof))

      IF (do_scatt) THEN
        DO col = 0, ircld%ncolumn(prof)
          od_surf_ac(j) = transmission_scatt_ir_dyn%opdpacsun(levsurf, col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j) - &
               transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j))

          transmission_aux%solar_path2%tau_surf_ac(col,j) = EXP(-od_surf_ac(j))

          transmission_aux%solar_path1%tau_surf_ac(col,j) = &
              EXP(-od_surf_ac(j) * (raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof)))

          transmission_aux%solar_path2%tau_surf(col,j) = path2%tau_surf(j) * &
              transmission_aux%solar_path2%tau_surf_ac(col,j)

          transmission_aux%solar_path1%tau_surf(col,j) = path1%tau_surf(j) * &
              transmission_aux%solar_path1%tau_surf_ac(col,j)

          IF (do_single_scatt) THEN
            IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
              transmission_aux%solar_path2%od_frac_ac(col,j) = &
                  od_surf_ac(j) - transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j)
            ELSE
              transmission_aux%solar_path2%od_frac_ac(col,j) = &
                  od_surf_ac(j) - transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j)
            ENDIF

            ! solar_path2%od_sfrac is actually along the sun-surface path (not "path2")
            transmission_aux%solar_path2%od_sfrac(col,j) = &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                raytracing%pathsun(levsurf-1, prof) / raytracing%patheff(levsurf-1, prof)

            transmission_aux%solar_path1%od_sfrac(col,j) = &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                raytracing%pathsat(levsurf-1, prof) / raytracing%patheff(levsurf-1, prof)

            IF (transmission_aux%solar_path1%tau_surf(col,j) > min_tau) THEN
              transmission_aux%solar_path1%tau_surf_r(col,j) = 1._jprb / transmission_aux%solar_path1%tau_surf(col,j)
            ELSE
              transmission_aux%solar_path1%tau_surf_r(col,j) = 0._jprb
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO

!---------------------------------------------------------------------------------------
!3. Store transmittances for other columns
!---------------------------------------------------------------------------------------

  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    IF (solar(j)) THEN
      IF (do_scatt) THEN
        DO col = 0, ircld%ncolumn(prof)
          DO lev = 1, nlevels
            transmission_aux%solar_path2%tau_level(lev,col,j) = &
                path2%tau_level(lev,j) * EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j))

            IF (lev > 1) THEN
              lay = lev - 1
              transmission_aux%solar_path1%tau_level(lev, col,j) = &
                  path1%tau_level(lev,j) * EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j) * &
                  raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof))
            ELSE
              ! for lev == 1, opdpacsun(lev, col,j) == 0.
              transmission_aux%solar_path1%tau_level(lev,col,j) = path1%tau_level(lev,j)
            ENDIF

            IF (do_single_scatt) THEN
              IF (transmission_aux%solar_path1%tau_level(lev,col,j) > min_tau) THEN
                transmission_aux%solar_path1%tau_level_r(lev,col,j) = &
                    1._jprb / transmission_aux%solar_path1%tau_level(lev,col,j)
              ELSE
                transmission_aux%solar_path1%tau_level_r(lev,col,j) = 0._jprb
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        IF (do_single_scatt) THEN
          DO lay = 1, nlayers
            ! solar_path2%od_singlelayer is actually along the sun-surface path (not "path2")
            transmission_aux%solar_path2%od_singlelayer(:,lay,j) = &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(:,lay,j)) * &
                raytracing%pathsun(lay,prof) / raytracing%patheff(lay,prof)

            transmission_aux%solar_path1%od_singlelayer(:,lay,j) = &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(:,lay,j)) * &
                raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof)
          ENDDO
        ENDIF

      ELSE
        col = 0
        transmission_aux%solar_path2%tau_level(:,col,j) = path2%tau_level(:,j)
        transmission_aux%solar_path2%tau_surf(col,j)    = path2%tau_surf(j)
        transmission_aux%solar_path1%tau_level(:,col,j) = path1%tau_level(:,j)
        transmission_aux%solar_path1%tau_surf(col,j)    = path1%tau_surf(j)
      ENDIF

      IF (do_single_scatt) THEN
        DO lay = 1, nlayers
          transmission_scatt_ir%opdpext(:,lay,j) = &
              path2%od_singlelayer(lay,j) + raytracing%patheff(lay,prof) * &
              (transmission_scatt_ir%opdpabs(:,lay,j) + transmission_scatt_ir%opdpsca(:,lay,j))
          WHERE (transmission_scatt_ir%opdpext(:,lay,j) > 0._jprb)
            transmission_scatt_ir%ssa_solar(:,lay,j) = raytracing%patheff(lay,prof) * &
                transmission_scatt_ir%opdpsca(:,lay,j) / transmission_scatt_ir%opdpext(:,lay,j)
          ENDWHERE
        ENDDO
      ENDIF

      ! Store clear-stream transmittances in transmission structure (includes aerosols if present)
      transmission%tausun_total_path2(j) = transmission_aux%solar_path2%tau_surf(0,j)
      transmission%tausun_total_path1(j) = transmission_aux%solar_path1%tau_surf(0,j)
      transmission%tausun_levels_path2(:,j) = transmission_aux%solar_path2%tau_level(:,0,j)
      transmission%tausun_levels_path1(:,j) = transmission_aux%solar_path1%tau_level(:,0,j)
    ENDIF ! solar channel
  ENDDO

! See notes in rttov_transmit.F90 regarding fac2
  IF (do_single_scatt) THEN
    DO j = 1, nchanprof
      IF (solar(j)) THEN
        prof = chanprof(j)%prof
        WHERE (transmission_aux%solar_path1%tau_level(:,0:ircld%ncolumn(prof),j) < min_tau)
          transmission_aux%solar_path1%fac2(:,0:ircld%ncolumn(prof),j) = 0.0_jprb
        ELSEWHERE
          transmission_aux%solar_path1%fac2(:,0:ircld%ncolumn(prof),j) = 1.0_jprb
        ENDWHERE
      ENDIF
    ENDDO
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_solar
