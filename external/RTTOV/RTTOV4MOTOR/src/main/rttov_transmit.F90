! Description:
!> @file
!!   Calculate various optical depth and transmittance quantities for emitted
!!   (thermal) radiation
!
!> @brief
!!   Calculate various optical depth and transmittance quantities for emitted
!!   (thermal) radiation
!!
!!
!! @param[in]     opts                        RTTOV options structure
!! @param[in]     do_lambertian               flag indicating whether Lambertian surface is active for each channel
!! @param[in]     nlayers                     number of layers in input profile
!! @param[in]     chanprof                    specifies channels and profiles to simulate
!! @param[in]     chanflag                    flags indicating which channels have emissive contribution
!! @param[in]     aux                         auxiliary profile variables
!! @param[in]     coef                        optical depth coefficients structure
!! @param[in]     ircld                       computed cloud column data
!! @param[in]     geometry                    geometry structure
!! @param[in]     opdp_path                   optical depths calculated by rttov_opdep
!! @param[in,out] od_level                    clear-sky level-to-space optical depths
!! @param[in,out] transmission                output transmission structure
!! @param[in,out] transmission_aux            top-level auxiliary transmission structure
!! @param[in,out] transmission_aux_path       auxiliary transmittances and optical depths for thermal radiation
!! @param[in]     transmission_scatt_ir       cloud/aerosol layer optical depths
!! @param[in]     transmission_scatt_ir_dyn   cloud/aerosol level-to-space optical depths
!! @param[in,out] tau_surf                    clear-sky surface-to-space transmittances used by TL/AD/K
!! @param[in,out] tau_level                   clear-sky level-to-space transmittances used by TL/AD/K
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
SUBROUTINE rttov_transmit( &
              opts,                         &
              do_lambertian,                &
              nlayers,                      &
              chanprof,                     &
              chanflag,                     &
              aux,                          &
              coef,                         &
              ircld,                        &
              geometry,                     &
              opdp_path,                    &
              od_level,                     &
              transmission,                 &
              transmission_aux,             &
              transmission_aux_path,        &
              transmission_scatt_ir,        &
              transmission_scatt_ir_dyn,    &
              tau_surf,                     &
              tau_level)

  USE rttov_types, ONLY :  &
       rttov_options,               &
       rttov_chanprof,              &
       rttov_coef,                  &
       rttov_transmission,          &
       rttov_transmission_aux,      &
       rttov_path_transmission,     &
       rttov_transmission_scatt_ir, &
       rttov_profile_aux,           &
       rttov_ircld,                 &
       rttov_geometry

  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_math_mod, ONLY : exponential, reciprocal
  USE rttov_const, ONLY : &
    max_optical_depth, min_od, min_tau, &
    sec_theta_eff, ir_scatt_chou, ir_scatt_dom

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),               INTENT(IN)    :: opts
  LOGICAL(KIND=jplm),                INTENT(IN)    :: do_lambertian(:)
  INTEGER(KIND=jpim),                INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof),              INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm),                INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_profile_aux),           INTENT(IN)    :: aux
  TYPE(rttov_coef),                  INTENT(IN)    :: coef
  TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
  TYPE(rttov_geometry),              INTENT(in)    :: geometry(:)
  REAL(KIND=jprb),                   INTENT(IN)    :: opdp_path(:,:)
  REAL(KIND=jprb),                   INTENT(INOUT) :: od_level(nlayers+1, SIZE(chanprof))
  TYPE(rttov_transmission),          INTENT(INOUT) :: transmission
  TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux
  TYPE(rttov_path_transmission),     INTENT(INOUT) :: transmission_aux_path
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir_dyn
  REAL(KIND=jprb),                   INTENT(INOUT) :: tau_surf(SIZE(chanprof))
  REAL(KIND=jprb),                   INTENT(INOUT) :: tau_level(nlayers+1, SIZE(chanprof))
!INTF_END

  REAL   (KIND=jprb) :: od_surf(SIZE(chanprof))                 ! sat to surface optical depth
  REAL   (KIND=jprb) :: od_surf_ac(0:SIZE(transmission_aux_path%tau_level(1,:,1)), SIZE(chanprof))
  REAL   (KIND=jprb) :: od_frac(SIZE(chanprof)), od_frac_ac
  ! pgf90 doesn't like opdpac(nlayers+1) with do_lambertian + IR cld/aer simulations:
  REAL   (KIND=jprb) :: opdpac(SIZE(transmission_aux_path%tau_level(:,0,1)))
  REAL   (KIND=jprb) :: opdpsurfac
  REAL   (KIND=jprb) :: small_val
  REAL   (KIND=jprb) :: od_surf_nadir, theta_eff, ref_power
  INTEGER(KIND=jpim) :: lev, lay, chan, j, prof, col, coli, levsurf
  INTEGER(KIND=jpim) :: nchanprof

! allocate automatic arrays using SIZE(transmission_aux_path%tau_level(1,:,1) to keep pgf90 compiler happy
  REAL   (KIND=jprb) :: ztemp(nlayers+1, 0:SIZE(transmission_aux_path%tau_level(1,:,1)))
  LOGICAL(KIND=jplm) :: do_scatt, do_chou, do_dom

  REAL   (KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 0_jpim, ZHOOK_HANDLE)
  small_val = (TINY(1._jprb)) ** (0.333333_jprb) ! XLF doesn't like 1/3
  nchanprof = SIZE(chanprof)

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  do_chou = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_chou)
  do_dom = do_scatt .AND. (opts%rt_ir%ir_scatt_model == ir_scatt_dom)

!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that value of opticaldepth is sensible
!cdir NODEP
!cdir COLLAPSE
  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      chan = chanprof(j)%chan
      DO lay = 1, nlayers
        transmission_aux_path%od_singlelayer(:,lay,j) = &
          -coef%ff_gam(chan) * (opdp_path(lay+1,j) - opdp_path(lay,j))
      ENDDO
      od_level(:,j) = MAX(coef%ff_gam(chan) * opdp_path(:,j), -max_optical_depth)
    ENDIF
  ENDDO

  IF (ALL(chanflag)) THEN
    CALL exponential(od_level, tau_level)
    transmission%tau_levels(:,1:nchanprof) = tau_level(:,:)
  ELSE
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        CALL exponential(od_level(:,j), tau_level(:,j))
        transmission%tau_levels(:,j) = tau_level(:,j)
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof = chanprof(j)%prof
! as defined in rttov_calc_nearest_lev
      levsurf = aux%s(prof)%nearestlev_surf
! layer above this
! NB all od_level -ve
! if surface below nlevels, pfraction -ve
      od_surf(j) = od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf-1,j) - od_level(levsurf,j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac(j) = od_surf(j) - od_level(levsurf-1,j)
      ELSE
        od_frac(j) = od_surf(j) - od_level(levsurf,j)
      ENDIF
      tau_surf(j) = EXP(od_surf(j))
    ENDIF
  ENDDO

  IF (do_chou) THEN
    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        levsurf = aux%s(prof)%nearestlev_surf
! all arrays here based on layers, not levels
        DO col = 0, ircld%ncolumn(prof)
          od_surf_ac(col,j) = transmission_scatt_ir_dyn%opdpac(levsurf,col,j) + aux%s(prof)%pfraction_surf * &
              (transmission_scatt_ir_dyn%opdpac(levsurf-1,col,j) - &
               transmission_scatt_ir_dyn%opdpac(levsurf,col,j))

          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            od_frac_ac = od_surf_ac(col,j) - transmission_scatt_ir_dyn%opdpac(levsurf-1,col,j)
          ELSE
            od_frac_ac = od_surf_ac(col,j) - transmission_scatt_ir_dyn%opdpac(levsurf,col,j)
          ENDIF
          transmission_aux_path%tau_surf_ac(col,j) = EXP(-od_surf_ac(col,j))
          transmission_aux_path%od_sfrac(col,j)    = -od_frac(j) + od_frac_ac

          transmission_aux_path%tau_surf(col,j) = tau_surf(j) * &
              transmission_aux_path%tau_surf_ac(col,j)
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
        ELSE
          od_surf_nadir = MAX(MIN(-od_surf(j) * geometry(prof)%coszen, 4._jprb), 0.05_jprb)

          ! The following is a polynomial approximation to the actual formula arccos(-od/ln(2 * E_3(od))) for
          ! the effective zenith angle [rad], with less than 0.003 error for od between 0.05 and 4.

          theta_eff = 1.029024_jprb - 0.367866_jprb * od_surf_nadir + 0.344010_jprb * od_surf_nadir**2 &
                      -0.219791_jprb * od_surf_nadir**3 + 0.078976_jprb * od_surf_nadir**4 &
                      -0.014515_jprb * od_surf_nadir**5 + 0.001061_jprb * od_surf_nadir**6
          ref_power = geometry(prof)%coszen / COS(theta_eff)
        ENDIF

        ! Recall that od_level and od_surf are negative, but the aer/cld opdeps are positive
        IF (do_scatt) THEN
          DO col = 0, ircld%ncolumn(prof)
            WHERE (transmission_scatt_ir_dyn%opdpac(:,col,j) < max_optical_depth)
              opdpac = transmission_scatt_ir_dyn%opdpac(:,col,j)
            ELSEWHERE
              opdpac = max_optical_depth
            ENDWHERE
            CALL exponential((od_level(:,j) - opdpac) * ref_power, &
                             transmission_aux_path%tau_level_p(:,col,j))
            opdpsurfac = MIN(od_surf_ac(col,j), max_optical_depth)
            transmission_aux_path%tau_surf_p(col,j) = EXP((od_surf(j) - opdpsurfac) * ref_power)
          ENDDO
        ELSE
          CALL exponential(od_level(:,j) * ref_power, transmission_aux_path%tau_level_p(:,0,j))
          transmission_aux_path%tau_surf_p(0,j) = EXP(od_surf(j) * ref_power)
        ENDIF
        DO col = 0, ircld%ncolumn(prof)
          CALL reciprocal(transmission_aux_path%tau_level_p(:,col,j), transmission_aux_path%tau_level_p_r(:,col,j))
          transmission_aux_path%tau_surf_p_r(col,j) = 1._jprb / transmission_aux_path%tau_surf_p(col,j)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

!---------------------------------------------------------------------------------------
!3. Store transmittances for other columns and single column for o/p
!---------------------------------------------------------------------------------------
  IF (do_chou) THEN
    WHERE (chanflag(:)) transmission%tau_total(1:nchanprof) = transmission_aux_path%tau_surf(0,:)

    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof

        DO col = 0, ircld%ncolumn(prof)
          CALL exponential(-transmission_scatt_ir_dyn%opdpac(:,col,j), &
            ztemp(:,col))
        ENDDO

        WHERE (ztemp(:,0:ircld%ncolumn(prof)) < EXP(-max_optical_depth))
          ztemp(:,0:ircld%ncolumn(prof)) = EXP(-max_optical_depth)
        ENDWHERE

        DO col = 0, ircld%ncolumn(prof)
          transmission_aux_path%tau_level(:,col,j) = tau_level(:,j) * ztemp(:,col)
        ENDDO
        CALL reciprocal(transmission_aux_path%tau_level(:,0:ircld%ncolumn(prof),j), &
                        transmission_aux_path%tau_level_r(:,0:ircld%ncolumn(prof),j))

        transmission_aux_path%od_singlelayer(:,:,j) = &
          transmission_aux_path%od_singlelayer(:,:,j) + transmission_scatt_ir%opdpacl(:,:,j)

        WHERE (transmission_aux_path%od_singlelayer(:,:,j) < small_val)
          transmission_aux_path%od_singlelayer(:,:,j) = small_val
        ENDWHERE
        CALL reciprocal(transmission_aux_path%od_singlelayer(:,:,j), &
                        transmission_aux_path%od_singlelayer_r(:,:,j))

        transmission_aux_path%tau_surf(0:ircld%ncolumn(prof),j) = &
          MAX(small_val, transmission_aux_path%tau_surf(0:ircld%ncolumn(prof),j))
        CALL reciprocal(transmission_aux_path%tau_surf(0:ircld%ncolumn(prof),j), &
                        transmission_aux_path%tau_surf_r(0:ircld%ncolumn(prof),j))

        transmission%tau_levels(:,j) = transmission_aux_path%tau_level(:,0,j)
      ENDIF
    ENDDO

  ELSE ! clear-sky or DOM
    WHERE (chanflag(:)) transmission%tau_total(1:nchanprof) = tau_surf(:)

    DO j = 1, nchanprof
      IF (chanflag(j)) THEN
        prof = chanprof(j)%prof
        col = 0
        IF (.NOT. do_dom) transmission_aux_path%od_sfrac(col,j) = -od_frac(j)

! DAR stop changes in emissivity_ad/k for v. small tau and stop nag from complaining for under/overflows
        IF (.NOT. do_dom) THEN
          transmission_aux_path%tau_surf(col,j) = MAX(small_val, tau_surf(j))
          transmission_aux_path%tau_surf_r(0,j) = &
                              1._jprb / transmission_aux_path%tau_surf(0,j)
        ENDIF

        IF (.NOT. do_dom) THEN
          transmission_aux_path%tau_level(:,col,j) = tau_level(:,j)
          CALL reciprocal(transmission_aux_path%tau_level(:,0,j), &
                          transmission_aux_path%tau_level_r(:,0,j))
        ENDIF

        transmission_aux_path%od_singlelayer(0,:,j) = &
            MAX(small_val, transmission_aux_path%od_singlelayer(0,:,j))
        IF (.NOT. do_dom) CALL reciprocal(transmission_aux_path%od_singlelayer(0,:,j), &
                                          transmission_aux_path%od_singlelayer_r(0,:,j))
      ENDIF
    ENDDO
  ENDIF

  IF (do_dom) RETURN

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof = chanprof(j)%prof
! stop floating overflow
      transmission_aux_path%od_sfrac_r(0:ircld%ncolumn(prof),j) = &
        1.0_jprb / (transmission_aux_path%od_sfrac(0:ircld%ncolumn(prof),j) + small_val)
    ENDIF
  ENDDO

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof = chanprof(j)%prof
      DO col = 0, ircld%ncolumn(prof)
        DO lay = 1, nlayers
          coli = ircld%icldarr(col,lay,prof)
          lev = lay + 1

! DAR: separated fac1 and fac2 becuase fac2 can be used separately from fac1 in Phil Watts calcs
          IF (transmission_aux_path%od_singlelayer(coli,lay,j) < min_od) THEN
            transmission_aux%fac1(lay,col,j) = 0.0_jprb
          ELSE
            IF (opts%rt_all%dtau_test .AND. &
                transmission_aux_path%tau_level(lay,col,j) - &
                transmission_aux_path%tau_level(lev,col,j) < min_od) THEN
              transmission_aux%fac1(lay,col,j) = 0.0_jprb
            ELSE
              transmission_aux%fac1(lay,col,j) = 1.0_jprb
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      WHERE (transmission_aux_path%tau_level(:,0:ircld%ncolumn(prof),j) < min_tau)
        transmission_aux_path%fac2(:,0:ircld%ncolumn(prof),j) = 0.0_jprb
      ELSEWHERE
        transmission_aux_path%fac2(:,0:ircld%ncolumn(prof),j) = 1.0_jprb
      ENDWHERE
    ENDIF
  ENDDO

  DO j = 1, nchanprof
    IF (chanflag(j)) THEN
      prof = chanprof(j)%prof
      WHERE (transmission_aux_path%tau_surf(0:ircld%ncolumn(prof),j) > min_tau)
        transmission_aux%surf_fac(0:ircld%ncolumn(prof),j) = 1.0_jprb
      ELSEWHERE
        transmission_aux%surf_fac(0:ircld%ncolumn(prof),j) = 0.0_jprb
      ENDWHERE
    ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit
