! Description:
!> @file
!!   TL of optical depth and transmittance calculations for solar radiation
!
!> @brief
!!   TL of optical depth and transmittance calculations for solar radiation
!!
!!
!! @param[in]     opts                           RTTOV options structure
!! @param[in]     nlayers                        number of layers in input profile
!! @param[in]     nprofiles                      number of profiles being simulated
!! @param[in]     chanprof                       specifies channels and profiles to simulate
!! @param[in]     solar                          flags indicating which channels have solar contribution
!! @param[in]     aux                            auxiliary profile variables
!! @param[in]     aux_tl                         auxiliary profile variable perturbations
!! @param[in]     coef                           optical depth coefficients structure
!! @param[in]     raytracing                     raytracing structure
!! @param[in]     raytracing_tl                  raytracing perturbations
!! @param[in]     ircld                          computed cloud column data
!! @param[in]     opdp_path_tl                   optical depth perturbations
!! @param[in]     path2                          optical depths and transmittances on "path2"
!! @param[in]     path1                          optical depths and transmittances on "path1"
!! @param[in,out] transmission_tl                output transmission perturbations
!! @param[in]     transmission_aux               top-level auxiliary transmission structure
!! @param[in,out] transmission_aux_tl            auxiliary transmission perturbations
!! @param[in]     transmission_scatt_ir          cloud/aerosol layer optical depths
!! @param[in,out] transmission_scatt_ir_tl       cloud/aerosol layer optical depth perturbations
!! @param[in]     transmission_scatt_ir_dyn      cloud/aerosol level-to-space optical depths
!! @param[in]     transmission_scatt_ir_dyn_tl   cloud/aerosol level-to-space optical depth perturbations
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
SUBROUTINE rttov_transmit_solar_tl( &
              opts,                            &
              nlayers,                         &
              nprofiles,                       &
              chanprof,                        &
              solar,                           &
              aux,                             &
              aux_tl,                          &
              coef,                            &
              raytracing,                      &
              raytracing_tl,                   &
              ircld,                           &
              opdp_path_tl,                    &
              path2,                           &
              path1,                           &
              transmission_tl,                 &
              transmission_aux,                &
              transmission_aux_tl,             &
              transmission_scatt_ir,           &
              transmission_scatt_ir_tl,        &
              transmission_scatt_ir_dyn,       &
              transmission_scatt_ir_dyn_tl)

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
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_const, ONLY : vis_scatt_single
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  INTEGER(KIND=jpim),                INTENT(IN)    :: nprofiles
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux_tl
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing_tl
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_opdp_path),             INTENT(IN)    :: opdp_path_tl
  TYPE(rttov_path_traj_trans),       INTENT(IN)    :: path2
  TYPE(rttov_path_traj_trans),       INTENT(IN)    :: path1
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission_tl
  TYPE(rttov_transmission_aux),      INTENT(IN)    :: transmission_aux
  TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn_tl
!INTF_END

  REAL   (KIND=jprb) :: od_singlelayer_tl(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_tl(nlayers+1,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_level_tl(nlayers+1,SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_tl(nlayers+1, SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsat_patheff_tl(nlayers,nprofiles)
  REAL   (KIND=jprb) :: pathsun_patheff_tl(nlayers,nprofiles)
  REAL   (KIND=jprb) :: od_surf, od_surf_ac
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, col, levsurf, nlevels
  INTEGER(KIND=jpim) :: nchanprof
  LOGICAL(KIND=jplm) :: do_scatt, do_single_scatt
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR_TL', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)
  nlevels   = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_single_scatt = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single

!----------------------------------------
!2. Compute layer to space optical depths
!----------------------------------------
  DO j = 1, nchanprof
    IF (solar(j)) THEN
      chan = chanprof(j)%chan
      DO lay = 1, nlayers
        od_singlelayer_tl(lay,j) = -coef%ff_gam(chan) * &
          (opdp_path_tl%sun_level_path2(lay+1,j) - opdp_path_tl%sun_level_path2(lay,j))
      ENDDO
      od_level_tl(:,j) = coef%ff_gam(chan) * opdp_path_tl%sun_level_path2(:,j)
      tausun_level_tl(:,j) = od_level_tl(:,j) * path2%tau_level(:,j)
    ENDIF
  ENDDO

  DO j = 1, nchanprof
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      tau_level_tl(1,j) = 0._jprb
      DO lev = 2, nlevels
        lay = lev - 1

        pathsat_patheff_tl(lay,prof) = &
            raytracing_tl%pathsat(lay,prof) / raytracing%patheff(lay,prof) - &
            raytracing%pathsat(lay,prof) * raytracing_tl%patheff(lay,prof) / &
            raytracing%patheff(lay,prof) ** 2_jpim

        pathsun_patheff_tl(lay,prof) = &
            raytracing_tl%pathsun(lay,prof) / raytracing%patheff(lay,prof) - &
            raytracing%pathsun(lay,prof) * raytracing_tl%patheff(lay,prof) / &
            raytracing%patheff(lay,prof) ** 2_jpim

        tau_level_tl(lev,j) = path1%tau_level(lev,j) * &
            (path2%od_level(lev,j) * pathsat_patheff_tl(lay,prof) + od_level_tl(lev,j) * &
            raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof))
      ENDDO
    ENDIF
  ENDDO

!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchanprof
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
    levsurf = aux%s(prof)%nearestlev_surf
    IF (solar(j)) THEN
      od_surf_tl(j) = od_level_tl(levsurf,j) +                                                      &
          aux_tl%s(prof)%pfraction_surf * (path2%od_level(levsurf-1,j) - path2%od_level(levsurf,j)) +  &
          aux%s(prof)%pfraction_surf * (od_level_tl(levsurf-1,j) - od_level_tl(levsurf,j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf-1,j)
      ELSE
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf,j)
      ENDIF
      tausun_surf_tl(j) = od_surf_tl(j) * path2%tau_surf(j)
      od_surf = path2%od_level(levsurf,j) + &
          aux%s(prof)%pfraction_surf * (path2%od_level(levsurf-1,j) - path2%od_level(levsurf,j))
      tau_surf_tl(j) = path1%tau_surf(j) * &
          (od_surf * pathsat_patheff_tl(levsurf-1,prof) + od_surf_tl(j) * &
          raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof))

      IF (do_scatt) THEN
        DO col = 0, ircld%ncolumn(prof)
          od_surf_ac_tl(j) = transmission_scatt_ir_dyn_tl%opdpacsun(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn_tl%opdpacsun(levsurf-1,col,j) - &
               transmission_scatt_ir_dyn_tl%opdpacsun(levsurf,col,j)) + aux_tl%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j) - &
               transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j))

          transmission_aux_tl%solar_path2%tau_surf_ac(col,j) = &
               - od_surf_ac_tl(j) * transmission_aux%solar_path2%tau_surf_ac(col,j)

          od_surf_ac = transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j) - &
               transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j))

          transmission_aux_tl%solar_path1%tau_surf_ac(col,j) = &
              transmission_aux%solar_path1%tau_surf_ac(col,j) * &
              (- od_surf_ac * pathsat_patheff_tl(levsurf-1,prof) - od_surf_ac_tl(j) * &
              raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof))

          transmission_aux_tl%solar_path2%tau_surf(col,j) = &
              tausun_surf_tl(j) * transmission_aux%solar_path2%tau_surf_ac(col,j) + &
              transmission_aux_tl%solar_path2%tau_surf_ac(col,j) * path2%tau_surf(j)

          transmission_aux_tl%solar_path1%tau_surf(col,j) = &
              tau_surf_tl(j) * transmission_aux%solar_path1%tau_surf_ac(col,j) + &
              transmission_aux_tl%solar_path1%tau_surf_ac(col,j) * path1%tau_surf(j)

          IF (do_single_scatt) THEN
            IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
              transmission_aux_tl%solar_path2%od_frac_ac(col,j) = &
                  od_surf_ac_tl(j) - transmission_scatt_ir_dyn_tl%opdpacsun(levsurf-1,col,j)
            ELSE
              transmission_aux_tl%solar_path2%od_frac_ac(col,j) = &
                  od_surf_ac_tl(j) - transmission_scatt_ir_dyn_tl%opdpacsun(levsurf,col,j)
            ENDIF

            transmission_aux_tl%solar_path2%od_sfrac(col,j) = &
                (-od_frac_tl(j) + transmission_aux_tl%solar_path2%od_frac_ac(col,j)) * &
                  raytracing%pathsun(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof) + &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                pathsun_patheff_tl(levsurf-1,prof)

            transmission_aux_tl%solar_path1%od_sfrac(col,j) = &
                (-od_frac_tl(j) + transmission_aux_tl%solar_path2%od_frac_ac(col,j)) * &
                  raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof) + &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                pathsat_patheff_tl(levsurf-1,prof)
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
    chan = chanprof(j)%chan
    IF (solar(j)) THEN
      IF (do_scatt) THEN
        DO col = 0, ircld%ncolumn(prof)
          DO lev = 1, nlevels
            transmission_aux_tl%solar_path2%tau_level(lev,col,j) = &
                (tausun_level_tl(lev,j) - &
                transmission_scatt_ir_dyn_tl%opdpacsun(lev,col,j) * path2%tau_level(lev,j)) * &
                EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j))

            IF (lev > 1) THEN
              lay = lev - 1
              transmission_aux_tl%solar_path1%tau_level(lev,col,j) =      &
                  EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j) * &
                  raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof)) * &
                  (tau_level_tl(lev,j) + path1%tau_level(lev,j) * &
                  (-transmission_scatt_ir_dyn_tl%opdpacsun(lev,col,j) * &
                  raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) + &
                  (-transmission_scatt_ir_dyn%opdpacsun(lev,col,j) * pathsat_patheff_tl(lay,prof))))
            ELSE
              ! for lev == 1, opdpacsun(lev,col,j) == 0.
              transmission_aux_tl%solar_path1%tau_level(lev,col,j) = tau_level_tl(lev,j)
            ENDIF
          ENDDO
        ENDDO

        IF (do_single_scatt) THEN
          DO lay = 1, nlayers
            lev = lay + 1
            transmission_aux_tl%solar_path2%od_singlelayer(:,lay,j) = &
                (od_singlelayer_tl(lay,j) + transmission_scatt_ir_tl%opdpaclsun(:,lay,j)) * &
                raytracing%pathsun(lay,prof) / raytracing%patheff(lay,prof) + &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(:,lay,j)) * &
                pathsun_patheff_tl(lay,prof)

            transmission_aux_tl%solar_path1%od_singlelayer(:,lay,j) = &
                (od_singlelayer_tl(lay,j) + transmission_scatt_ir_tl%opdpaclsun(:,lay,j)) * &
                raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) + &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(:,lay,j)) * &
                pathsat_patheff_tl(lay,prof)
          ENDDO
        ENDIF
      ELSE
        col = 0
        transmission_tl%tausun_levels_path2(:,j) = tausun_level_tl(:,j)
        transmission_tl%tausun_total_path2(j)    = tausun_surf_tl(j)
        transmission_tl%tausun_levels_path1(:,j) = tau_level_tl(:,j)
        transmission_tl%tausun_total_path1(j)    = tau_surf_tl(j)
        transmission_aux_tl%solar_path2%tau_level(:,col,j) = tausun_level_tl(:,j)
        transmission_aux_tl%solar_path2%tau_surf(col,j)    = tausun_surf_tl(j)
        transmission_aux_tl%solar_path1%tau_level(:,col,j) = tau_level_tl(:,j)
        transmission_aux_tl%solar_path1%tau_surf(col,j)    = tau_surf_tl(j)
      ENDIF

      IF (do_single_scatt) THEN
        DO lay = 1, nlayers
          transmission_scatt_ir_tl%opdpext(:,lay,j) = &
              od_singlelayer_tl(lay,j) + raytracing_tl%patheff(lay,prof) * &
              (transmission_scatt_ir%opdpabs(:,lay,j) + transmission_scatt_ir%opdpsca(:,lay,j)) + &
              raytracing%patheff(lay,prof) * &
              (transmission_scatt_ir_tl%opdpabs(:,lay,j) + transmission_scatt_ir_tl%opdpsca(:,lay,j))
          WHERE (transmission_scatt_ir%opdpext(:,lay,j) > 0._jprb)
            transmission_scatt_ir_tl%ssa_solar(:,lay,j) = &
                (raytracing_tl%patheff(lay,prof) * transmission_scatt_ir%opdpsca(:,lay,j) + &
                raytracing%patheff(lay,prof) * transmission_scatt_ir_tl%opdpsca(:,lay,j)) / &
                transmission_scatt_ir%opdpext(:,lay,j) - &
                transmission_scatt_ir_tl%opdpext(:,lay,j) * raytracing%patheff(lay,prof) * &
                transmission_scatt_ir%opdpsca(:,lay,j) / transmission_scatt_ir%opdpext(:,lay,j)**2
          ENDWHERE
        ENDDO
      ENDIF
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_solar_tl
