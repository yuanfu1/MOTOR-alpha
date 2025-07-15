! Description:
!> @file
!!   TL of optical depth and transmittance calculations for emitted (thermal) radiation
!
!> @brief
!!   TL of optical depth and transmittance calculations for emitted (thermal) radiation
!!
!!
!! @param[in]     opts                           RTTOV options structure
!! @param[in]     do_lambertian                  flag indicating whether Lambertian surface is active for each channel
!! @param[in]     nlayers                        number of layers in input profile
!! @param[in]     chanprof                       specifies channels and profiles to simulate
!! @param[in]     chanflag                       flags indicating which channels have emissive contribution
!! @param[in]     aux                            auxiliary profile variables
!! @param[in]     aux_tl                         auxiliary profile variable perturbations
!! @param[in]     coef                           optical depth coefficients structure
!! @param[in]     ircld                          computed cloud column data
!! @param[in]     geometry                       geometry structure
!! @param[in]     opdp_path_tl                   optical depth perturbations
!! @param[in]     od_level                       level-to-space optical depths
!! @param[in,out] transmission_tl                output transmittance perturbations
!! @param[in]     transmission_aux_path          auxiliary transmittances and optical depths for thermal radiation
!! @param[in,out] transmission_aux_path_tl       auxiliary transmittance and optical depth perturbations for thermal radiation
!! @param[in]     transmission_scatt_ir_tl       cloud/aerosol layer optical depth perturbations
!! @param[in]     transmission_scatt_ir_dyn      cloud/aerosol level-to-space optical depths
!! @param[in]     transmission_scatt_ir_dyn_tl   cloud/aerosol level-to-space optical depth perturbations
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
SUBROUTINE rttov_transmit_tl( &
              opts,                            &
              do_lambertian,                   &
              nlayers,                         &
              chanprof,                        &
              chanflag,                        &
              aux,                             &
              aux_tl,                          &
              coef,                            &
              ircld,                           &
              geometry,                        &
              opdp_path_tl,                    &
              od_level,                        &
              transmission_tl,                 &
              transmission_aux_path,           &
              transmission_aux_path_tl,        &
              transmission_scatt_ir_tl,        &
              transmission_scatt_ir_dyn,       &
              transmission_scatt_ir_dyn_tl,    &
              tau_surf,                        &
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
  USE rttov_math_mod, ONLY : exponential, exponential_tl
  USE rttov_const, ONLY : sec_theta_eff, ir_scatt_chou, ir_scatt_dom
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_lambertian(:)
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux_tl
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_geometry),              INTENT(IN)    :: geometry(:)
  REAL(KIND=jprb),                   INTENT(IN)    :: opdp_path_tl(:,:)
  REAL(KIND=jprb),                   INTENT(IN)    :: od_level(nlayers+1, SIZE(chanprof))
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission_tl
  TYPE(rttov_path_transmission),     INTENT(IN)    :: transmission_aux_path
  TYPE(rttov_path_transmission),     INTENT(INOUT) :: transmission_aux_path_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_tl
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn_tl
  REAL(KIND=jprb),                   INTENT(IN)    :: tau_surf(SIZE(chanprof))
  REAL(KIND=jprb),                   INTENT(IN)    :: tau_level(nlayers+1, SIZE(chanprof))
!INTF_END

  REAL   (KIND=jprb) :: od_surf, od_surf_ac
  REAL   (KIND=jprb) :: od_surf_nadir, theta_eff, ref_power
  REAL   (KIND=jprb) :: od_surf_nadir_tl, theta_eff_tl, ref_power_tl
  REAL   (KIND=jprb) :: od_level_tl(nlayers+1, SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: od_surf_ac_tl(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac_tl(SIZE(chanprof)), od_frac_ac_tl
  REAL   (KIND=jprb) :: tau_surf_tl(SIZE(chanprof))
  REAL   (KIND=jprb) :: tau_level_tl(nlayers+1, SIZE(chanprof))

  INTEGER(KIND=jpim) :: lay, chan, prof, col, j, levsurf
  INTEGER(KIND=jpim) :: nchanprof

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  REAL   (KIND=jprb) :: ztemp_tl(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(KIND=jplm) :: do_scatt, do_chou, do_dom

  REAL   (KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 0_jpim, ZHOOK_HANDLE)
  nchanprof = SIZE(chanprof)

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_chou = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_chou)
  do_dom = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_dom)

!--------------------------------------------------------------
!1. Assemble layer optical depths
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      DO lay = 1, nlayers
        transmission_aux_path_tl%od_singlelayer(:,lay,j) =  &
          -coef%ff_gam(chan) * (opdp_path_tl(lay+1,j) - opdp_path_tl(lay,j))
      ENDDO
      od_level_tl(:,j)  = coef%ff_gam(chan) * opdp_path_tl(:,j)
      tau_level_tl(:,j) = od_level_tl(:,j) * tau_level(:,j)
    ENDIF
  ENDDO


!----------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!----------------------------------------------------------------------------------------
  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof = chanprof(j)%prof
      levsurf = aux%s(prof)%nearestlev_surf
      od_surf_tl(j) =  &
          od_level_tl(levsurf,j) + aux_tl%s(prof)%pfraction_surf * (od_level(levsurf-1,j) - od_level(levsurf,j)) +  &
          aux%s(prof)%pfraction_surf * (od_level_tl(levsurf-1,j) - od_level_tl(levsurf,j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf-1,j)
      ELSE
        od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf,j)
      ENDIF
      tau_surf_tl(j) = od_surf_tl(j) * tau_surf(j)
    ENDIF
  ENDDO

  IF (do_chou) THEN
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        levsurf = aux%s(prof)%nearestlev_surf
        DO col = 0, ircld%ncolumn(prof)
          od_surf_ac_tl(col,j) = transmission_scatt_ir_dyn_tl%opdpac(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn_tl%opdpac(levsurf-1,col,j) -                                              &
               transmission_scatt_ir_dyn_tl%opdpac(levsurf,col,j)) + aux_tl%s(prof)%pfraction_surf *               &
              (transmission_scatt_ir_dyn%opdpac(levsurf-1,col,j) -                                                 &
               transmission_scatt_ir_dyn%opdpac(levsurf,col,j))

          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_frac_ac_tl = od_surf_ac_tl(col,j) - transmission_scatt_ir_dyn_tl%opdpac(levsurf-1,col,j)
          ELSE
            od_frac_ac_tl = od_surf_ac_tl(col,j) - transmission_scatt_ir_dyn_tl%opdpac(levsurf,col,j)
          ENDIF
          transmission_aux_path_tl%tau_surf_ac(col,j) = - od_surf_ac_tl(col,j) * &
              transmission_aux_path%tau_surf_ac(col,j)
          transmission_aux_path_tl%od_sfrac(col,j)    = - od_frac_tl(j) + od_frac_ac_tl

          transmission_aux_path_tl%tau_surf(col,j) = tau_surf_tl(j) * &
              transmission_aux_path%tau_surf_ac(col,j) + &
              tau_surf(j) * transmission_aux_path_tl%tau_surf_ac(col,j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  IF (ANY(do_lambertian) .AND. .NOT. do_dom) THEN
    DO j = 1, nchanprof
      IF (chanflag(j) .AND. do_lambertian(j)) THEN
        prof = chanprof(j)%prof

        IF (opts%rt_all%lambertian_fixed_angle) THEN
          ref_power = sec_theta_eff * geometry(prof)%coszen
          ref_power_tl = 0._jprb
        ELSE
          od_surf = LOG(tau_surf(j))
          od_surf_nadir = -od_surf * geometry(prof)%coszen
          IF (od_surf_nadir < 0.05_jprb) THEN
            od_surf_nadir = 0.05_jprb
            od_surf_nadir_tl = 0._jprb
          ELSEIF (od_surf_nadir > 4._jprb) THEN
            od_surf_nadir = 4._jprb
            od_surf_nadir_tl = 0._jprb
          ELSE
            od_surf_nadir_tl = -od_surf_tl(j) * geometry(prof)%coszen
          ENDIF

          theta_eff = 1.029024_jprb - 0.367866_jprb * od_surf_nadir + 0.344010_jprb * od_surf_nadir**2 &
                      -0.219791_jprb * od_surf_nadir**3 + 0.078976_jprb * od_surf_nadir**4 &
                      -0.014515_jprb * od_surf_nadir**5 + 0.001061_jprb * od_surf_nadir**6
          ref_power = geometry(prof)%coszen / COS(theta_eff)
          IF (od_surf_nadir_tl /= 0._jprb) THEN
            theta_eff_tl = od_surf_nadir_tl * (-0.367866_jprb + 2 * 0.344010_jprb * od_surf_nadir &
                           -3 * 0.219791_jprb * od_surf_nadir**2 + 4 * 0.078976_jprb * od_surf_nadir**3 &
                           -5 * 0.014515_jprb * od_surf_nadir**4 + 6 * 0.001061_jprb * od_surf_nadir**5)
            ref_power_tl = theta_eff_tl * TAN(theta_eff) * ref_power
          ELSE
            ref_power_tl = 0._jprb
          ENDIF
        ENDIF

        IF (do_chou) THEN
          DO col = 0, ircld%ncolumn(prof)
            IF (opts%rt_all%lambertian_fixed_angle .OR. ref_power_tl == 0._jprb) THEN
              CALL exponential_tl((od_level_tl(:,j) - transmission_scatt_ir_dyn_tl%opdpac(:,col,j)) * ref_power, &
                transmission_aux_path%tau_level_p(:,col,j), &
                transmission_aux_path_tl%tau_level_p(:,col,j))

              transmission_aux_path_tl%tau_surf_p(col,j) = (od_surf_tl(j) - od_surf_ac_tl(col,j)) * ref_power * &
                                                           transmission_aux_path%tau_surf_p(col,j)
            ELSE
              CALL exponential_tl((od_level_tl(:,j) - transmission_scatt_ir_dyn_tl%opdpac(:,col,j)) * ref_power + &
                (od_level(:,j) - transmission_scatt_ir_dyn%opdpac(:,col,j)) * ref_power_tl, &
                transmission_aux_path%tau_level_p(:,col,j), &
                transmission_aux_path_tl%tau_level_p(:,col,j))

              levsurf = aux%s(prof)%nearestlev_surf
              od_surf_ac = transmission_scatt_ir_dyn%opdpac(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
                  (transmission_scatt_ir_dyn%opdpac(levsurf-1,col,j) - &
                   transmission_scatt_ir_dyn%opdpac(levsurf,col,j))

              transmission_aux_path_tl%tau_surf_p(col,j) = ((od_surf_tl(j) - od_surf_ac_tl(col,j)) * ref_power + &
                                                            (od_surf - od_surf_ac) * ref_power_tl) * &
                                                           transmission_aux_path%tau_surf_p(col,j)
            ENDIF
          ENDDO
        ELSE
          IF (opts%rt_all%lambertian_fixed_angle .OR. ref_power_tl == 0._jprb) THEN
            CALL exponential_tl(od_level_tl(:,j) * ref_power, &
              transmission_aux_path%tau_level_p(:,0,j), &
              transmission_aux_path_tl%tau_level_p(:,0,j))

            transmission_aux_path_tl%tau_surf_p(0,j) = od_surf_tl(j) * ref_power * transmission_aux_path%tau_surf_p(0,j)
          ELSE
            CALL exponential_tl(od_level_tl(:,j) * ref_power + od_level(:,j) * ref_power_tl, &
              transmission_aux_path%tau_level_p(:,0,j), &
              transmission_aux_path_tl%tau_level_p(:,0,j))

            transmission_aux_path_tl%tau_surf_p(0,j) = (od_surf_tl(j) * ref_power + od_surf * ref_power_tl) * &
                                                       transmission_aux_path%tau_surf_p(0,j)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!3. Store transmittances for other columns
!---------------------------------------------------------------------------------------
  IF (do_chou) THEN
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof

        DO col = 0, ircld%ncolumn(prof)
          CALL exponential(-transmission_scatt_ir_dyn%opdpac(:,col,j), &
            ztemp(:,col))
          CALL exponential_tl(ztemp(:,col) , &
                             -transmission_scatt_ir_dyn_tl%opdpac(:,col,j), &
                             ztemp_tl(:,col))
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          transmission_aux_path_tl%tau_level(:,col,j) = &
            tau_level_tl(:,j) * ztemp(:,col) + &
            tau_level(:,j) * ztemp_tl(:,col)
        ENDDO

        transmission_aux_path_tl%od_singlelayer(:,:,j) = &
          transmission_aux_path_tl%od_singlelayer(:,:,j) + &
          transmission_scatt_ir_tl%opdpacl(:,:,j)

!        WHERE (transmission_aux_path%od_singlelayer(:,:,j) < small_val)
!          transmission_aux_path_tl%od_singlelayer(:,:,j) = 0._jprb
!        ENDWHERE
      ENDIF
    ENDDO

  ELSE ! clear-sky or DOM

    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        transmission_tl%tau_total(j)    = tau_surf_tl(j)
        transmission_tl%tau_levels(:,j) = tau_level_tl(:,j)

        IF (.NOT. do_dom) THEN
          col = 0
          transmission_aux_path_tl%tau_level(:,col,j) = tau_level_tl(:,j)
          transmission_aux_path_tl%tau_surf(col,j) = tau_surf_tl(j)
          transmission_aux_path_tl%od_sfrac(col,j) = - od_frac_tl(j)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_tl
