! Description:
!> @file
!!   K of optical depth and transmittance calculations for emitted (thermal) radiation
!
!> @brief
!!   K of optical depth and transmittance calculations for emitted (thermal) radiation
!!
!!
!! @param[in]     opts                           RTTOV options structure
!! @param[in]     do_lambertian                  flag indicating whether Lambertian surface is active for each channel
!! @param[in]     nlayers                        number of layers in input profile
!! @param[in]     chanprof                       specifies channels and profiles to simulate
!! @param[in]     chanflag                       flags indicating which channels have emissive contribution
!! @param[in]     aux                            auxiliary profile variables
!! @param[in,out] aux_k                          auxiliary profile variable increments
!! @param[in]     coef                           optical depth coefficients structure
!! @param[in]     ircld                          computed cloud column data
!! @param[in]     geometry                       geometry structure
!! @param[in,out] opdp_path_k                    optical depth increments
!! @param[in]     od_level                       level-to-space optical depths
!! @param[in,out] transmission_k                 input gradient wrt transmittances (usually zero)
!! @param[in]     transmission_aux_path          auxiliary transmittances and optical depths for thermal radiation
!! @param[in,out] transmission_aux_path_k        auxiliary transmittance and optical depth increments for thermal radiation
!! @param[in,out] transmission_scatt_ir_k        cloud/aerosol layer optical depth increments
!! @param[in]     transmission_scatt_ir_dyn      cloud/aerosol level-to-space optical depths
!! @param[in,out] transmission_scatt_ir_dyn_k    cloud/aerosol level-to-space optical depth increments
!! @param[in]     tau_surf                       surface-to-space transmittances used by TL/AD/K
!! @param[in]     tau_level                      level-to-space transmittances used by TL/AD/K
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
SUBROUTINE rttov_transmit_k( &
              opts,                           &
              do_lambertian,                  &
              nlayers,                        &
              chanprof,                       &
              chanflag,                       &
              aux,                            &
              aux_k,                          &
              coef,                           &
              ircld,                          &
              geometry,                       &
              opdp_path_k,                    &
              od_level,                       &
              transmission_k,                 &
              transmission_aux_path,          &
              transmission_aux_path_k,        &
              transmission_scatt_ir_k,        &
              transmission_scatt_ir_dyn,      &
              transmission_scatt_ir_dyn_k,    &
              tau_surf,                       &
              tau_level)

  USE rttov_types, ONLY :  &
       rttov_options,               &
       rttov_chanprof,              &
       rttov_coef,                  &
       rttov_transmission,          &
       rttov_path_transmission,     &
       rttov_transmission_scatt_ir, &
       rttov_profile_aux,           &
       rttov_ircld,                 &
       rttov_geometry
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_math_mod, ONLY : exponential, exponential_ad
  USE rttov_const, ONLY : sec_theta_eff, ir_scatt_chou, ir_scatt_dom
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_lambertian(:)
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(INOUT) :: aux_k
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_geometry),              INTENT(in)    :: geometry(:)
  REAL(KIND=jprb),                   INTENT(INOUT) :: opdp_path_k(:,:)
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission_k
  TYPE(rttov_path_transmission),     INTENT(IN)    :: transmission_aux_path
  TYPE(rttov_path_transmission),     INTENT(INOUT) :: transmission_aux_path_k
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_dyn_k
  REAL(KIND=jprb),                   INTENT(IN)    :: od_level(nlayers+1, SIZE(chanprof))
  REAL(KIND=jprb),                   INTENT(IN)    :: tau_surf(SIZE(chanprof))
  REAL(KIND=jprb),                   INTENT(IN)    :: tau_level(nlayers+1, SIZE(chanprof))
!INTF_END

  REAL   (KIND=jprb) :: od_surf, od_surf_ac
  REAL   (KIND=jprb) :: od_surf_nadir, od_surf_nadir_ref, theta_eff, ref_power
  REAL   (KIND=jprb) :: od_surf_nadir_k, theta_eff_k, ref_power_k
  REAL   (KIND=jprb) :: od_level_k(nlayers+1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_k(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_k(SIZE(chanprof)), od_frac_ac_k
  REAL   (KIND=jprb) :: tau_surf_k(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_k(nlayers+1, SIZE(chanprof))
 
  INTEGER(KIND=jpim) :: lev, lay, chan, prof, j, levsurf, col
  INTEGER(KIND=jpim) :: nlevels, nchanprof

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  REAL   (KIND=jprb) :: ztemp_k(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(KIND=jplm) :: do_scatt, do_chou, do_dom

  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_K', 0_jpim, ZHOOK_HANDLE)
  nchanprof             = SIZE(chanprof)
  nlevels               = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_chou = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_chou)
  do_dom = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_dom)

!-----------------------------------------------------
!K of store transmittances for other columns
!-----------------------------------------------------
  od_level_k(:,:)       = 0._jprb
  tau_level_k(:,:)      = 0.0_jprb
  tau_surf_k(:)         = 0.0_jprb
  od_frac_k(:)          = 0.0_jprb
  od_surf_k(:)          = 0.0_jprb
  od_surf_ac_k          = 0.0_jprb
  od_frac_ac_k          = 0.0_jprb

  IF (do_chou) THEN
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof

        DO col = 0, ircld%ncolumn(prof)
          CALL exponential(-transmission_scatt_ir_dyn%opdpac(:,col,j), &
            ztemp(:,col))
        ENDDO
!start AD

        transmission_scatt_ir_k%opdpacl(:,:,j) = &
          transmission_scatt_ir_k%opdpacl(:,:,j) + &
          transmission_aux_path_k%od_singlelayer(:,:,j)

        DO col = 0, ircld%ncolumn(prof)
          tau_level_k(:,j) = tau_level_k(:,j) + &
            transmission_aux_path_k%tau_level(:,col,j) * ztemp(:,col)

          ztemp_k(:,col) = &!ztemp_k(:,col) + &
            tau_level(:,j) * transmission_aux_path_k%tau_level(:,col,j)
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          CALL exponential_ad(-ztemp(:,col), &
                             transmission_scatt_ir_dyn_k%opdpac(:,col,j), &
                             ztemp_k(:,col), acc = .TRUE._jplm)
        ENDDO
      ENDIF
    ENDDO
  ELSE
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        IF (.NOT. do_dom) THEN
          col = 0
          tau_level_k(:,j) = tau_level_k(:,j) + transmission_aux_path_k%tau_level(:,col,j)
          od_frac_k(j)  = od_frac_k(j) - transmission_aux_path_k%od_sfrac(col,j)
          tau_surf_k(j) = tau_surf_k(j) + transmission_aux_path_k%tau_surf(col,j)
        ENDIF

        tau_level_k(:,j) = tau_level_k(:,j) + transmission_k%tau_levels(:,j)
        transmission_k%tau_levels(:,j) = 0
        tau_surf_k(j) = tau_surf_k(j) + transmission_k%tau_total(j)
        transmission_k%tau_total(j) = 0
      ENDIF
    ENDDO
  ENDIF

  ! Not used anymore but no need to 0
  ! DO j = 1, nchanprof
  !   IF (chanflag(j)) THEN
  !     DO col = 0, ircld%ncolumn(prof)
  !       transmission_aux_path_k%tau_level(:, col,j)      = 0.0_jprb
  !       transmission_aux_path_k%od_singlelayer(:, col,j) = 0.0_jprb
  !     ENDDO
  !   ENDIF
  ! ENDDO

  IF (ANY(do_lambertian) .AND. .NOT. do_dom) THEN
    DO j = 1, nchanprof
      IF (chanflag(j) .AND. do_lambertian(j)) THEN
        prof = chanprof(j)%prof

        IF (opts%rt_all%lambertian_fixed_angle) THEN
          ref_power = sec_theta_eff * geometry(prof)%coszen
        ELSE
          od_surf = LOG(tau_surf(j))
          od_surf_nadir_ref = -od_surf * geometry(prof)%coszen
          od_surf_nadir = MAX(MIN(od_surf_nadir_ref, 4._jprb), 0.05_jprb)
          theta_eff = 1.029024_jprb - 0.367866_jprb * od_surf_nadir + 0.344010_jprb * od_surf_nadir**2 &
                      -0.219791_jprb * od_surf_nadir**3 + 0.078976_jprb * od_surf_nadir**4 &
                      -0.014515_jprb * od_surf_nadir**5 + 0.001061_jprb * od_surf_nadir**6
          ref_power = geometry(prof)%coszen / COS(theta_eff)

          ref_power_k = 0._jprb
          theta_eff_k = 0._jprb
          od_surf_nadir_k = 0._jprb
        ENDIF

        IF (do_chou) THEN
          DO col = 0, ircld%ncolumn(prof)
            od_surf_k(j) = od_surf_k(j) + &
              ref_power * transmission_aux_path%tau_surf_p(col,j) * transmission_aux_path_k%tau_surf_p(col,j)

            od_surf_ac_k(col,j) = od_surf_ac_k(col,j) - &
              ref_power * transmission_aux_path%tau_surf_p(col,j) * transmission_aux_path_k%tau_surf_p(col,j)

            ! must move ref_power to input because can't have it with intent inout
            CALL exponential_ad(transmission_aux_path%tau_level_p(:,col,j) * ref_power, &
              od_level_k(:,j), &
              transmission_aux_path_k%tau_level_p(:,col,j), acc=.TRUE._jplm)

            CALL exponential_ad(-transmission_aux_path%tau_level_p(:,col,j) * ref_power, &
              transmission_scatt_ir_dyn_k%opdpac(:,col,j), &
              transmission_aux_path_k%tau_level_p(:,col,j), acc=.TRUE._jplm)

            IF (.NOT. opts%rt_all%lambertian_fixed_angle) THEN
              levsurf = aux%s(prof)%nearestlev_surf
              od_surf_ac = transmission_scatt_ir_dyn%opdpac(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
                  (transmission_scatt_ir_dyn%opdpac(levsurf-1,col,j) - &
                   transmission_scatt_ir_dyn%opdpac(levsurf,col,j))

              ref_power_k = ref_power_k + (od_surf - od_surf_ac) * &
                transmission_aux_path%tau_surf_p(col,j) * transmission_aux_path_k%tau_surf_p(col,j)

              ref_power_k = ref_power_k + SUM((od_level(:,j) - transmission_scatt_ir_dyn%opdpac(:,col,j)) * &
                transmission_aux_path%tau_level_p(:,col,j) * transmission_aux_path_k%tau_level_p(:,col,j))
            ENDIF
          ENDDO
        ELSE
          od_surf_k(j) = od_surf_k(j) + &
            ref_power * transmission_aux_path%tau_surf_p(0,j) * transmission_aux_path_k%tau_surf_p(0,j)

          ! must move ref_power to input because can't have it with intent inout
          CALL exponential_ad(transmission_aux_path%tau_level_p(:,0,j) * ref_power, &
            od_level_k(:,j), &
            transmission_aux_path_k%tau_level_p(:,0,j), acc=.TRUE._jplm)

          IF (.NOT. opts%rt_all%lambertian_fixed_angle) THEN
            ref_power_k = ref_power_k + &
              od_surf * transmission_aux_path%tau_surf_p(0,j) * transmission_aux_path_k%tau_surf_p(0,j)

            ref_power_k = ref_power_k + SUM(transmission_aux_path%tau_level_p(:,0,j) * &
                                            od_level(:,j) * transmission_aux_path_k%tau_level_p(:,0,j))
          ENDIF
        ENDIF

        IF (.NOT. opts%rt_all%lambertian_fixed_angle) THEN
          IF (ref_power_k /= 0._jprb) THEN
            theta_eff_k = theta_eff_k + ref_power_k * TAN(theta_eff) * ref_power
            ref_power_k = 0._jprb

            od_surf_nadir_k = od_surf_nadir_k + theta_eff_k * &
                (-0.367866_jprb + 2 * 0.344010_jprb * od_surf_nadir &
                 -3 * 0.219791_jprb * od_surf_nadir**2 + 4 * 0.078976_jprb * od_surf_nadir**3 &
                 -5 * 0.014515_jprb * od_surf_nadir**4 + 6 * 0.001061_jprb * od_surf_nadir**5)
            theta_eff_k = 0._jprb
          ELSE
            od_surf_nadir_k = 0._jprb
          ENDIF

          IF (od_surf_nadir_ref < 0.05_jprb) THEN
            od_surf_nadir_k = 0._jprb
          ELSEIF (od_surf_nadir_ref > 4._jprb) THEN
            od_surf_nadir_k = 0._jprb
          ELSE
            od_surf_k(j) = od_surf_k(j) - od_surf_nadir_k * geometry(prof)%coszen
            od_surf_nadir_k = 0._jprb
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!K of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  IF (do_chou) THEN
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof    = chanprof(j)%prof
        chan    = chanprof(j)%chan
        levsurf = aux%s(prof)%nearestlev_surf
        DO col = ircld%ncolumn(prof), 0, -1
          transmission_aux_path_k%tau_surf_ac(col,j) =     &
              transmission_aux_path_k%tau_surf_ac(col,j) + &
              transmission_aux_path_k%tau_surf(col,j) * tau_surf(j)
          tau_surf_k(j)                          =      &
              tau_surf_k(j) + transmission_aux_path_k%tau_surf(col,j) * &
              transmission_aux_path%tau_surf_ac(col,j)
          od_frac_k(j) = od_frac_k(j) - transmission_aux_path_k%od_sfrac(col,j)
          od_frac_ac_k = od_frac_ac_k + transmission_aux_path_k%od_sfrac(col,j)
          od_surf_ac_k(col,j) = &
              od_surf_ac_k(col,j) - transmission_aux_path_k%tau_surf_ac(col,j) * &
              transmission_aux_path%tau_surf_ac(col,j)
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_surf_ac_k(col,j) = od_surf_ac_k(col,j) + od_frac_ac_k
            transmission_scatt_ir_dyn_k%opdpac(levsurf-1, col,j) =     &
                transmission_scatt_ir_dyn_k%opdpac(levsurf-1, col,j) - od_frac_ac_k
          ELSE
            od_surf_ac_k(col,j) = od_surf_ac_k(col,j) + od_frac_ac_k
            transmission_scatt_ir_dyn_k%opdpac(levsurf, col,j) =      &
                transmission_scatt_ir_dyn_k%opdpac(levsurf, col,j) - od_frac_ac_k
          ENDIF
          transmission_scatt_ir_dyn_k%opdpac(levsurf, col,j)     =      &
              transmission_scatt_ir_dyn_k%opdpac(levsurf, col,j) +      &
              od_surf_ac_k(col,j) * (1._jprb - aux%s(prof)%pfraction_surf)
          transmission_scatt_ir_dyn_k%opdpac(levsurf-1, col,j) =      &
              transmission_scatt_ir_dyn_k%opdpac(levsurf-1, col,j) + &
              od_surf_ac_k(col,j) * aux%s(prof)%pfraction_surf
          aux_k%s(j)%pfraction_surf = &
              aux_k%s(j)%pfraction_surf + od_surf_ac_k(col,j) * &
              (transmission_scatt_ir_dyn%opdpac(levsurf-1, col,j) - &
               transmission_scatt_ir_dyn%opdpac(levsurf, col,j))
!           od_surf_ac_k(col,j) = 0._jprb
          od_frac_ac_k = 0._jprb
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof         = chanprof(j)%prof
      levsurf      = aux%s(prof)%nearestlev_surf
      od_surf_k(j) = od_surf_k(j) + tau_surf_k(j) * tau_surf(j)
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_surf_k(j)            = od_surf_k(j) + od_frac_k(j)
        od_level_k(levsurf-1,j) = od_level_k(levsurf-1,j) - od_frac_k(j)
      ELSE
        od_surf_k(j)          = od_surf_k(j) + od_frac_k(j)
        od_level_k(levsurf,j) = od_level_k(levsurf,j) - od_frac_k(j)
      ENDIF
      od_level_k(levsurf,j)   = od_level_k(levsurf,j) + od_surf_k(j) * (1._jprb - aux%s(prof)%pfraction_surf)
      od_level_k(levsurf-1,j) = od_level_k(levsurf-1,j) + od_surf_k(j) * aux%s(prof)%pfraction_surf
      aux_k%s(j)%pfraction_surf  =      &
          aux_k%s(j)%pfraction_surf + od_surf_k(j) * (od_level(levsurf-1,j) - od_level(levsurf,j))
      od_surf_k(j) = 0._jprb
    ENDIF
  ENDDO

!-------------------------------------------
!K of assemble layer optical depths
!-------------------------------------------

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      DO lev = 1, nlevels
        od_level_k(lev,j)  = od_level_k(lev,j) + tau_level_k(lev,j) * tau_level(lev,j)
        od_level_k(lev,j)  = coef%ff_gam(chan) * od_level_k(lev,j)
        opdp_path_k(lev,j) = opdp_path_k(lev,j) + od_level_k(lev,j)
      ENDDO
    ENDIF
  ENDDO
! ad of level to space optical depths
! opdp_path_k(:,:) = opdp_path_k(:,:) + od_level_k(:,:)
! ad of single layer optical depths
  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      DO lay = nlayers, 1, -1
        opdp_path_k(lay,j)   = &
          opdp_path_k(lay,j) + coef%ff_gam(chan) * SUM(transmission_aux_path_k%od_singlelayer(:,lay,j))
        opdp_path_k(lay+1,j) = &
          opdp_path_k(lay+1,j) - coef%ff_gam(chan) * SUM(transmission_aux_path_k%od_singlelayer(:,lay,j))
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_k
