! Description:
!> @file
!!   K of optical depth and transmittance calculations for solar radiation
!
!> @brief
!!   K of optical depth and transmittance calculations for solar radiation
!!
!!
!! @param[in]     opts                           RTTOV options structure
!! @param[in]     nlayers                        number of layers in input profile
!! @param[in]     chanprof                       specifies channels and profiles to simulate
!! @param[in]     solar                          flags indicating which channels have solar contribution
!! @param[in]     aux                            auxiliary profile variables
!! @param[in,out] aux_k                          auxiliary profile variable increments
!! @param[in]     coef                           optical depth coefficients structure
!! @param[in]     raytracing                     raytracing structure
!! @param[in,out] raytracing_k                   raytracing increments
!! @param[in]     ircld                          computed cloud column data
!! @param[in,out] opdp_path_k                    optical depth increments
!! @param[in]     path2                          optical depths and transmittances on "path2"
!! @param[in]     path1                          optical depths and transmittances on "path1"
!! @param[in,out] transmission_k                 output transmission increments
!! @param[in]     transmission_aux               top-level auxiliary transmission structure
!! @param[in,out] transmission_aux_k             auxiliary transmission increments
!! @param[in]     transmission_scatt_ir          cloud/aerosol layer optical depths
!! @param[in,out] transmission_scatt_ir_k        cloud/aerosol layer optical depth increments
!! @param[in]     transmission_scatt_ir_dyn      cloud/aerosol level-to-space optical depths
!! @param[in,out] transmission_scatt_ir_dyn_k    cloud/aerosol level-to-space optical depth increments
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
SUBROUTINE rttov_transmit_solar_k( &
              opts,                           &
              nlayers,                        &
              chanprof,                       &
              solar,                          &
              aux,                            &
              aux_k,                          &
              coef,                           &
              raytracing,                     &
              raytracing_k,                   &
              ircld,                          &
              opdp_path_k,                    &
              path2,                          &
              path1,                          &
              transmission_k,                 &
              transmission_aux,               &
              transmission_aux_k,             &
              transmission_scatt_ir,          &
              transmission_scatt_ir_k,        &
              transmission_scatt_ir_dyn,      &
              transmission_scatt_ir_dyn_k)

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
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: solar(SIZE(chanprof))
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux_k
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
  TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_k
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_opdp_path),             INTENT(INOUT) :: opdp_path_k
  TYPE(rttov_path_traj_trans),       INTENT(IN)    :: path2
  TYPE(rttov_path_traj_trans),       INTENT(IN)    :: path1
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission_k
  TYPE(rttov_transmission_aux),      INTENT(IN)    :: transmission_aux
  TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux_k
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_dyn_k
!INTF_END

  REAL   (KIND=jprb) :: od_singlelayer_k(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_level_k(nlayers+1,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_level_k(nlayers+1,SIZE(chanprof))
  REAL   (KIND=jprb) :: tausun_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_k(nlayers+1,SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsat_patheff_k(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: pathsun_patheff_k(nlayers,SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf, od_surf_ac, zpatheff, tauacsunpath1
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, col, coli, levsurf, nlevels
  INTEGER(KIND=jpim) :: nchanprof
  LOGICAL(KIND=jplm) :: do_scatt, do_single_scatt
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR_K', 0_jpim, ZHOOK_HANDLE)
  nchanprof             = SIZE(chanprof)
  nlevels               = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_single_scatt = do_scatt .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single

!---------------------------------------------------------------------------------------
!K of store transmittances for other polarisations
!---------------------------------------------------------------------------------------
  od_level_k(:,:)       = 0._jprb
  tausun_level_k(:,:)   = 0._jprb
  tau_level_k(:,:)      = 0._jprb
  od_singlelayer_k(:,:) = 0._jprb
  tausun_surf_k(:)      = 0._jprb
  tau_surf_k(:)         = 0._jprb
  od_frac_k(:)          = 0._jprb
  od_surf_k(:)          = 0._jprb
  od_surf_ac_k(:)       = 0._jprb
  pathsat_patheff_k(:,:) = 0._jprb
  pathsun_patheff_k(:,:) = 0._jprb
  DO j = 1, nchanprof
    prof = chanprof(j)%prof
    chan = chanprof(j)%chan
    IF (solar(j)) THEN

      IF (do_single_scatt) THEN
        DO lay = nlayers, 1, -1
          WHERE (transmission_scatt_ir%opdpext(:,lay,j) > 0._jprb)
            transmission_scatt_ir_k%opdpext(:,lay,j) = transmission_scatt_ir_k%opdpext(:,lay,j) - &
                transmission_scatt_ir_k%ssa_solar(:,lay,j) * raytracing%patheff(lay,prof) * &
                transmission_scatt_ir%opdpsca(:,lay,j) / transmission_scatt_ir%opdpext(:,lay,j)**2
            transmission_scatt_ir_k%opdpsca(:,lay,j) = transmission_scatt_ir_k%opdpsca(:,lay,j) + &
                transmission_scatt_ir_k%ssa_solar(:,lay,j) * raytracing%patheff(lay,prof) / &
                transmission_scatt_ir%opdpext(:,lay,j)
          ENDWHERE
          DO coli = 0, 1
            IF (transmission_scatt_ir%opdpext(coli,lay,j) > 0._jprb) &
              raytracing_k%patheff(lay,j) = &
                raytracing_k%patheff(lay,j) + &
                  transmission_scatt_ir_k%ssa_solar(coli,lay,j) * transmission_scatt_ir%opdpsca(coli,lay,j) / &
                  transmission_scatt_ir%opdpext(coli,lay,j)
            IF (.NOT. opts%rt_ir%addclouds) EXIT
          ENDDO

          od_singlelayer_k(lay,j) = &
              od_singlelayer_k(lay,j) + SUM(transmission_scatt_ir_k%opdpext(:,lay,j))
          raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) + &
              SUM(transmission_scatt_ir_k%opdpext(:,lay,j) * &
                  (transmission_scatt_ir%opdpabs(:,lay,j) + transmission_scatt_ir%opdpsca(:,lay,j)))
          transmission_scatt_ir_k%opdpabs(:,lay,j) = &
              transmission_scatt_ir_k%opdpabs(:,lay,j) + &
              raytracing%patheff(lay,prof) * transmission_scatt_ir_k%opdpext(:,lay,j)
          transmission_scatt_ir_k%opdpsca(:,lay,j) = &
              transmission_scatt_ir_k%opdpsca(:,lay,j) + &
              raytracing%patheff(lay,prof) * transmission_scatt_ir_k%opdpext(:,lay,j)
        ENDDO

        DO coli = 0, 1
          DO lay = nlayers, 1,  - 1
            lev = lay + 1

            od_singlelayer_k(lay,j) = od_singlelayer_k(lay,j) + &
                transmission_aux_k%solar_path1%od_singlelayer(coli,lay,j) * &
                raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof)
            transmission_scatt_ir_k%opdpaclsun(coli,lay,j) = &
                transmission_scatt_ir_k%opdpaclsun(coli,lay,j) + &
                transmission_aux_k%solar_path1%od_singlelayer(coli,lay,j) * &
                raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof)
            pathsat_patheff_k(lay,j) = pathsat_patheff_k(lay,j) + &
                transmission_aux_k%solar_path1%od_singlelayer(coli,lay,j) * &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(coli,lay,j))

            od_singlelayer_k(lay,j) = od_singlelayer_k(lay,j) + &
                transmission_aux_k%solar_path2%od_singlelayer(coli,lay,j) * &
                raytracing%pathsun(lay,prof) / raytracing%patheff(lay,prof)
            transmission_scatt_ir_k%opdpaclsun(coli,lay,j) = &
                transmission_scatt_ir_k%opdpaclsun(coli,lay,j) + &
                transmission_aux_k%solar_path2%od_singlelayer(coli,lay,j) * &
                raytracing%pathsun(lay,prof) / raytracing%patheff(lay,prof)
            pathsun_patheff_k(lay,j) = pathsun_patheff_k(lay,j) + &
                transmission_aux_k%solar_path2%od_singlelayer(coli,lay,j) * &
                (path2%od_singlelayer(lay,j) + transmission_scatt_ir%opdpaclsun(coli,lay,j))
          ENDDO
          IF (.NOT. opts%rt_ir%addclouds) EXIT
        ENDDO
      ENDIF

      IF (do_scatt) THEN
        DO col = ircld%ncolumn(prof), 0, -1
          DO lev = nlevels, 1, -1
            IF (lev > 1) THEN
              lay = lev - 1

              tauacsunpath1 = EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j) * &
                  raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof))

              tau_level_k(lev,j) = tau_level_k(lev,j) +      &
                  tauacsunpath1 * transmission_aux_k%solar_path1%tau_level(lev,col,j)

              transmission_scatt_ir_dyn_k%opdpacsun(lev,col,j) = &
                  transmission_scatt_ir_dyn_k%opdpacsun(lev,col,j) - &
                  tauacsunpath1 * path1%tau_level(lev,j) * &
                  transmission_aux_k%solar_path1%tau_level(lev,col,j) * &
                  raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)

              pathsat_patheff_k(lay,j) = pathsat_patheff_k(lay,j) - &
                  tauacsunpath1 * path1%tau_level(lev,j) * &
                  transmission_aux_k%solar_path1%tau_level(lev,col,j) * &
                  transmission_scatt_ir_dyn%opdpacsun(lev,col,j)
            ELSE
              ! for lev == 1, opdpacsun(lev,col,j) == 0.
              tau_level_k(lev,j) = tau_level_k(lev,j) + &
                  transmission_aux_k%solar_path1%tau_level(lev,col,j)
            ENDIF

            tausun_level_k(lev,j) = tausun_level_k(lev,j) + &
                transmission_aux_k%solar_path2%tau_level(lev,col,j) * &
                EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j))
            transmission_scatt_ir_dyn_k%opdpacsun(lev,col,j) = &
                transmission_scatt_ir_dyn_k%opdpacsun(lev,col,j) - &
                transmission_aux_k%solar_path2%tau_level(lev,col,j) * path2%tau_level(lev,j) *  &
                EXP(-transmission_scatt_ir_dyn%opdpacsun(lev,col,j))
          ENDDO
        ENDDO
      ELSE
        col = 0
        transmission_k%tausun_levels_path2(:,j) = 0._jprb
        transmission_k%tausun_total_path2(j) = 0._jprb
        transmission_k%tausun_levels_path1(:,j) = 0._jprb
        transmission_k%tausun_total_path1(j) = 0._jprb

        tausun_level_k(:,j) = tausun_level_k(:,j) + transmission_aux_k%solar_path2%tau_level(:,col,j)
        tausun_surf_k(j) = tausun_surf_k(j) + transmission_aux_k%solar_path2%tau_surf(col,j)
        tau_level_k(:,j) = tau_level_k(:,j) + transmission_aux_k%solar_path1%tau_level(:,col,j)
        tau_surf_k(j) = tau_surf_k(j) + transmission_aux_k%solar_path1%tau_surf(col,j)
      ENDIF

      transmission_aux_k%solar_path2%tau_level(:,0:ircld%ncolumn(prof),j) = 0._jprb
    ENDIF
  ENDDO

!---------------------------------------------------------------------------------------
!K of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchanprof
    prof    = chanprof(j)%prof
    chan    = chanprof(j)%chan
    levsurf = aux%s(prof)%nearestlev_surf
    IF (solar(j)) THEN
      IF (do_scatt) THEN
        DO col = ircld%ncolumn(prof), 0, -1

          IF (do_single_scatt) THEN
            od_frac_k(j)  = od_frac_k(j) - &
                transmission_aux_k%solar_path1%od_sfrac(col,j) * &
                raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof)
            transmission_aux_k%solar_path2%od_frac_ac(col,j) = transmission_aux_k%solar_path2%od_frac_ac(col,j) + &
                transmission_aux_k%solar_path1%od_sfrac(col,j) * &
                raytracing%pathsat(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof)
            pathsat_patheff_k(levsurf-1,j) = pathsat_patheff_k(levsurf-1,j) + &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                transmission_aux_k%solar_path1%od_sfrac(col,j)

            od_frac_k(j) = od_frac_k(j) - &
                transmission_aux_k%solar_path2%od_sfrac(col,j) * &
                raytracing%pathsun(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof)
            transmission_aux_k%solar_path2%od_frac_ac(col,j) = transmission_aux_k%solar_path2%od_frac_ac(col,j) + &
                transmission_aux_k%solar_path2%od_sfrac(col,j) * &
                raytracing%pathsun(levsurf-1,prof) / raytracing%patheff(levsurf-1,prof)
            pathsun_patheff_k(levsurf-1,j) = pathsun_patheff_k(levsurf-1,j) + &
                (-path2%od_frac(j) + transmission_aux%solar_path2%od_frac_ac(col,j)) * &
                transmission_aux_k%solar_path2%od_sfrac(col,j)

            IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
              od_surf_ac_k(j) = &
                  od_surf_ac_k(j) + transmission_aux_k%solar_path2%od_frac_ac(col,j)
              transmission_scatt_ir_dyn_k%opdpacsun(levsurf-1,col,j) = &
                  transmission_scatt_ir_dyn_k%opdpacsun(levsurf-1,col,j) - &
                  transmission_aux_k%solar_path2%od_frac_ac(col,j)
            ELSE
              od_surf_ac_k(j) = &
                  od_surf_ac_k(j) + transmission_aux_k%solar_path2%od_frac_ac(col,j)
              transmission_scatt_ir_dyn_k%opdpacsun(levsurf,col,j) = &
                  transmission_scatt_ir_dyn_k%opdpacsun(levsurf,col,j) - &
                  transmission_aux_k%solar_path2%od_frac_ac(col,j)
            ENDIF
          ENDIF

          tau_surf_k(j) = &
              tau_surf_k(j) + transmission_aux_k%solar_path1%tau_surf(col,j) * &
              transmission_aux%solar_path1%tau_surf_ac(col,j)
          transmission_aux_k%solar_path1%tau_surf_ac(col,j) =     &
              transmission_aux_k%solar_path1%tau_surf_ac(col,j) + &
              transmission_aux_k%solar_path1%tau_surf(col,j) * path1%tau_surf(j)

          tausun_surf_k(j)                          =      &
              tausun_surf_k(j) + transmission_aux_k%solar_path2%tau_surf(col,j) * &
              transmission_aux%solar_path2%tau_surf_ac(col,j)
          transmission_aux_k%solar_path2%tau_surf_ac(col,j) =     &
              transmission_aux_k%solar_path2%tau_surf_ac(col,j) + &
              transmission_aux_k%solar_path2%tau_surf(col,j) * path2%tau_surf(j)

          od_surf_ac = transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j) - &
              transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j))

          pathsat_patheff_k(levsurf-1,j) = pathsat_patheff_k(levsurf-1,j) - &
              od_surf_ac * transmission_aux%solar_path1%tau_surf_ac(col,j) * &
              transmission_aux_k%solar_path1%tau_surf_ac(col,j)

          od_surf_ac_k(j) = od_surf_ac_k(j) - &
              transmission_aux%solar_path1%tau_surf_ac(col,j) * &
              transmission_aux_k%solar_path1%tau_surf_ac(col,j) * &
              raytracing%pathsat(levsurf-1, prof) / raytracing%patheff(levsurf-1, prof)

          od_surf_ac_k(j)                          =      &
              od_surf_ac_k(j) - transmission_aux_k%solar_path2%tau_surf_ac(col,j) * &
              transmission_aux%solar_path2%tau_surf_ac(col,j)

          transmission_scatt_ir_dyn_k%opdpacsun(levsurf,col,j)     = &
              transmission_scatt_ir_dyn_k%opdpacsun(levsurf,col,j) + &
              od_surf_ac_k(j) * (1._jprb - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_dyn_k%opdpacsun(levsurf-1,col,j) =      &
              transmission_scatt_ir_dyn_k%opdpacsun(levsurf-1,col,j) +  &
              od_surf_ac_k(j) *  aux%s(prof)%pfraction_surf
          aux_k%s(j)%pfraction_surf = aux_k%s(j)%pfraction_surf + od_surf_ac_k(j) &
               * (transmission_scatt_ir_dyn%opdpacsun(levsurf-1,col,j) - &
              transmission_scatt_ir_dyn%opdpacsun(levsurf,col,j))
          od_surf_ac_k(j) = 0._jprb
        ENDDO
      ENDIF

      od_surf = path2%od_level(levsurf,j) + &
          aux%s(prof)%pfraction_surf * (path2%od_level(levsurf-1,j) - path2%od_level(levsurf,j))
      pathsat_patheff_k(levsurf-1,j) = pathsat_patheff_k(levsurf-1,j) + &
          path1%tau_surf(j) * od_surf * tau_surf_k(j)
      od_surf_k(j) = od_surf_k(j) + &
          path1%tau_surf(j) * tau_surf_k(j) * &
          raytracing%pathsat(levsurf-1, prof) / raytracing%patheff(levsurf-1, prof)
      od_surf_k(j) = od_surf_k(j) + tausun_surf_k(j) * path2%tau_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_k(j) = od_surf_k(j) + od_frac_k(j)
        od_level_k(levsurf-1,j) = od_level_k(levsurf-1,j) - od_frac_k(j)
      ELSE
        od_surf_k(j) = od_surf_k(j) + od_frac_k(j)
        od_level_k(levsurf,j) = od_level_k(levsurf,j) - od_frac_k(j)
      ENDIF
      od_level_k(levsurf,j) = od_level_k(levsurf,j) + od_surf_k(j) * (1._jprb - aux%s(prof)%pfraction_surf)
      od_level_k(levsurf-1,j) = od_level_k(levsurf-1,j) + od_surf_k(j) * aux%s(prof)%pfraction_surf
      aux_k%s(j)%pfraction_surf = &
          aux_k%s(j)%pfraction_surf + od_surf_k(j) * (path2%od_level(levsurf-1,j) - path2%od_level(levsurf,j))
      od_surf_k(j) = 0._jprb
    ENDIF
  ENDDO

!-------------------------------------------
!K of assemble layer optical depths
!-------------------------------------------
  DO j = 1, nchanprof
    IF (solar(j)) THEN
      prof = chanprof(j)%prof
      tau_level_k(1,j) = 0._jprb
      DO lev = 2, nlevels
        lay = lev - 1

        pathsat_patheff_k(lay,j) = pathsat_patheff_k(lay,j) + &
            path1%tau_level(lev,j) * path2%od_level(lev,j) * tau_level_k(lev,j)

        od_level_k(lev,j) = od_level_k(lev,j) + &
            path1%tau_level(lev,j) * tau_level_k(lev,j) * &
            raytracing%pathsat(lay, prof) / raytracing%patheff(lay, prof)
      ENDDO
    ENDIF
  ENDDO

  DO j = 1, nchanprof
    IF (solar(j)) THEN
      DO lay = 1, nlayers
        prof = chanprof(j)%prof

        zpatheff = 1._jprb / raytracing%patheff(lay, prof)

        raytracing_k%pathsun(lay,j) = raytracing_k%pathsun(lay,j) + &
              pathsun_patheff_k(lay,j) * zpatheff

        raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) - &
              raytracing%pathsun(lay, prof) * pathsun_patheff_k(lay,j) * zpatheff ** 2_jpim

        raytracing_k%pathsat(lay,j) = raytracing_k%pathsat(lay,j) + &
              pathsat_patheff_k(lay,j) * zpatheff

        raytracing_k%patheff(lay,j) = raytracing_k%patheff(lay,j) - &
              raytracing%pathsat(lay, prof) * pathsat_patheff_k(lay,j) * zpatheff ** 2_jpim
      ENDDO
    ENDIF
  ENDDO

  DO j = 1, nchanprof
    IF (solar(j)) THEN
      chan = chanprof(j)%chan
      od_level_k(:,j) = od_level_k(:,j) + coef%ff_gam(chan) * tausun_level_k(:,j) * path2%tau_level(:,j)
      opdp_path_k%sun_level_path2(:,j) = opdp_path_k%sun_level_path2(:,j) + od_level_k(:,j)
    ENDIF
  ENDDO
  DO j = 1, nchanprof
    IF (solar(j)) THEN
      chan = chanprof(j)%chan
      DO lay = nlayers, 1, -1
        opdp_path_k%sun_level_path2(lay,j)   = &
          opdp_path_k%sun_level_path2(lay,j) + coef%ff_gam(chan) * od_singlelayer_k(lay,j)
        opdp_path_k%sun_level_path2(lay+1,j) = &
          opdp_path_k%sun_level_path2(lay+1,j) - coef%ff_gam(chan) * od_singlelayer_k(lay,j)
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_SOLAR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_solar_k
