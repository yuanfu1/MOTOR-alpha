! Description:
!> @file
!!   AD of integration of the radiative transfer equation
!
!> @brief
!!   AD of integration of the radiative transfer equation
!!
!!
!! @param[in]     opts                        RTTOV options structure
!! @param[in]     maxncolumns                 largest number of cloud columns across all profiles
!! @param[in]     chanprof                    specifies channels and profiles to simulate
!! @param[in]     emissivity                  input/output surface emissivities
!! @param[in,out] emissivity_ad               input/output surface emissivity increments
!! @param[in]     reflectance                 input/output surface reflectances for direct solar beam
!! @param[in,out] reflectance_ad              input/output surface reflectance increments for direct solar beam
!! @param[in]     refl_norm                   surface relfectance normalisation factors
!! @param[in]     diffuse_refl                surface reflectance for downwelling radiation
!! @param[in,out] diffuse_refl_ad             surface reflectance increments for downwelling radiation
!! @param[in]     do_lambertian               flag indicating whether Lambertian surface is active for each channel
!! @param[in]     do_mfasis                   flag to indicate MFASIS simulation
!! @param[in]     thermal                     per-channel flag to indicate if emissive simulations are being performed
!! @param[in]     dothermal                   flag to indicate if any emissive simulations are being performed
!! @param[in]     solar                       per-channel flag to indicate if solar simulations are being performed
!! @param[in]     dosolar                     flag to indicate if any solar simulations are being performed
!! @param[in]     do_rayleigh_ss              flag to indicate if Rayleigh single-scattering should be included
!! @param[in]     solar_spectrum              TOA solar irradiance for each channel
!! @param[in]     transmission_aux            top-level auxiliary transmission structure
!! @param[in,out] transmission_aux_ad         auxiliary transmission increments
!! @param[in]     transmission_scatt_ir       visible/IR cloud/aerosol scattering parameters
!! @param[in,out] transmission_scatt_ir_ad    visible/IR cloud/aerosol scattering parameter increments
!! @param[in]     profiles                    input atmospheric profiles and surface variables
!! @param[in,out] profiles_ad                 input atmospheric profile increments
!! @param[in]     profiles_dry                profiles in internal units
!! @param[in,out] profiles_dry_ad             profile increments in internal units
!! @param[in]     aux_prof                    auxiliary profile variables
!! @param[in,out] aux_prof_ad                 auxiliary profile variable increments
!! @param[in]     coef                        optical depth coefficients structure
!! @param[in]     raytracing                  raytracing structure
!! @param[in,out] raytracing_ad               raytracing increments
!! @param[in]     ircld                       computed cloud column data
!! @param[in,out] ircld_ad                    cloud column increments
!! @param[in]     rad                         primary output radiance structure
!! @param[in]     auxrad                      Planck radiances
!! @param[in,out] auxrad_ad                   Planck radiance increments
!! @param[in]     auxrad_column               internal radiance structure
!! @param[in,out] auxrad_column_ad            internal radiance increments
!! @param[in,out] rad_ad                      output radiance increments
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
SUBROUTINE rttov_integrate_ad(opts, maxncolumns, chanprof, emissivity, emissivity_ad, &
                              reflectance, reflectance_ad, refl_norm, diffuse_refl, diffuse_refl_ad, &
                              do_lambertian, do_mfasis, thermal, dothermal, solar, dosolar, do_rayleigh_ss, solar_spectrum, &
                              transmission_aux, transmission_aux_ad, transmission_scatt_ir, &
                              transmission_scatt_ir_ad, profiles, profiles_ad, profiles_dry, &
                              profiles_dry_ad, aux_prof, aux_prof_ad, coef, raytracing, &
                              raytracing_ad, ircld, ircld_ad, rad, auxrad, auxrad_ad, auxrad_column, &
                              auxrad_column_ad, rad_ad)

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_options, rttov_profile, rttov_profile_aux, &
                          rttov_transmission_aux, rttov_transmission_scatt_ir, rttov_radiance, &
                          rttov_ircld, rttov_raytracing, rttov_radiance_aux, rttov_emissivity, rttov_reflectance
!INTF_OFF
  USE rttov_const, ONLY : realtol, sensor_id_po, min_od, min_tau, pi_r, z4pi_r, &
                          deg2rad, gravity, na, rgc, Mh2o, Mair, &
                          overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
                          vis_scatt_dom, vis_scatt_single, ir_scatt_dom, &
                          ray_ps, ray_ts, ray_scs_wlm, &
                          ray_scs_a1x, ray_scs_b1, ray_scs_c1, ray_scs_d1, &
                          ray_scs_a2x, ray_scs_b2, ray_scs_c2, ray_scs_d2
! 
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),                INTENT(IN)              :: opts
  INTEGER(jpim),                      INTENT(IN)              :: maxncolumns
  TYPE(rttov_chanprof),               INTENT(IN)              :: chanprof(:)
  TYPE(rttov_emissivity),             INTENT(IN),    OPTIONAL :: emissivity(:)
  TYPE(rttov_emissivity),             INTENT(INOUT), OPTIONAL :: emissivity_ad(:)
  TYPE(rttov_reflectance),            INTENT(IN),    OPTIONAL :: reflectance(:)
  TYPE(rttov_reflectance),            INTENT(INOUT), OPTIONAL :: reflectance_ad(:)
  REAL(jprb),                         INTENT(IN)              :: refl_norm(:)
  REAL(jprb),                         INTENT(IN)              :: diffuse_refl(:)
  REAL(jprb),                         INTENT(INOUT)           :: diffuse_refl_ad(:)
  LOGICAL(jplm),                      INTENT(IN)              :: do_lambertian(:)
  LOGICAL(jplm),                      INTENT(IN)              :: do_mfasis
  LOGICAL(jplm),                      INTENT(IN)              :: thermal(:)
  LOGICAL(jplm),                      INTENT(IN)              :: dothermal
  LOGICAL(jplm),                      INTENT(IN)              :: solar(:)
  LOGICAL(jplm),                      INTENT(IN)              :: dosolar
  LOGICAL(jplm),                      INTENT(IN)              :: do_rayleigh_ss
  REAL(jprb),                         INTENT(IN)              :: solar_spectrum(:)
  TYPE(rttov_transmission_aux),       INTENT(IN)              :: transmission_aux
  TYPE(rttov_transmission_aux),       INTENT(INOUT)           :: transmission_aux_ad
  TYPE(rttov_transmission_scatt_ir),  INTENT(IN)              :: transmission_scatt_ir
  TYPE(rttov_transmission_scatt_ir),  INTENT(INOUT)           :: transmission_scatt_ir_ad
  TYPE(rttov_profile),                INTENT(IN)              :: profiles(:)
  TYPE(rttov_profile),                INTENT(INOUT)           :: profiles_ad(:)
  TYPE(rttov_profile),                INTENT(IN)              :: profiles_dry(:)
  TYPE(rttov_profile),                INTENT(INOUT)           :: profiles_dry_ad(:)
  TYPE(rttov_profile_aux),            INTENT(IN)              :: aux_prof
  TYPE(rttov_profile_aux),            INTENT(INOUT)           :: aux_prof_ad
  TYPE(rttov_coef),                   INTENT(IN)              :: coef
  TYPE(rttov_raytracing),             INTENT(IN)              :: raytracing
  TYPE(rttov_raytracing),             INTENT(INOUT)           :: raytracing_ad
  TYPE(rttov_ircld),                  INTENT(IN)              :: ircld
  TYPE(rttov_ircld),                  INTENT(INOUT)           :: ircld_ad
  TYPE(rttov_radiance),               INTENT(IN)              :: rad
  TYPE(rttov_radiance_aux),           INTENT(IN)              :: auxrad
  TYPE(rttov_radiance_aux),           INTENT(INOUT)           :: auxrad_ad
  TYPE(rttov_radiance_aux),           INTENT(IN)              :: auxrad_column
  TYPE(rttov_radiance_aux),           INTENT(INOUT)           :: auxrad_column_ad
  TYPE(rttov_radiance),               INTENT(INOUT)           :: rad_ad
!INTF_END

#include "rttov_calcrad_ad.interface"

  INTEGER(jpim) :: i, lev, col, lay, nchanprof, nlayers, nlevels
  REAL(jprb)    :: refl, refl_norm_scat
  LOGICAL(jplm) :: keyradonly, do_scatt, dom_ir, dom_vis, do_single_scatt

  INTEGER(jpim) :: iv2lay(SIZE(chanprof)), iv2lev(SIZE(chanprof))
  INTEGER(jpim) :: iv3lay(SIZE(chanprof)), iv3lev(SIZE(chanprof))
  INTEGER(jpim) :: pol_id(SIZE(chanprof))

  REAL(jprb)    :: cfraction(SIZE(chanprof)), cfraction_ad(SIZE(chanprof)), pfraction(SIZE(chanprof))
  LOGICAL(jplm) :: sateqsun(profiles(1)%nlayers,SIZE(profiles(:)))

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header ------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',0_jpim,ZHOOK_HANDLE)

#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

  ! auxrad_column_ad is initialised in rttov_ad

  !-------------------------------------------------------------------------------
  ! Initialise useful variables
  !-------------------------------------------------------------------------------

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  dom_ir = do_scatt .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom
  dom_vis = do_scatt .AND. opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom
  do_single_scatt = do_scatt .AND. dosolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single
  keyradonly = opts%rt_ir%addaerosl .OR. opts%rt_ir%pc%addpc .OR. (dom_ir .AND. dom_vis)

  IF (dosolar .AND. .NOT. do_mfasis) &
    sateqsun(:,:) = (ABS(raytracing%pathsat(:,:) - raytracing%pathsun(:,:)) < realtol)

  DO i = 1, nchanprof
    cfraction(i) = aux_prof%s(prof)%cfraction
    pfraction(i) = aux_prof%s(prof)%pfraction_surf
  ENDDO

  DO i = 1, nchanprof
    ! case-1: surf lies above lev=nlevels
    iv3lev(i) = aux_prof%s(prof)%nearestlev_surf - 1   ! lowest lev above surf
    ! case-2: surf lies below lev=nlevels
    IF (pfraction(i) < 0._jprb) iv3lev(i) = iv3lev(i) + 1  ! iv3lev=iv2lev=lowest lev above surf

    iv2lev(i) = aux_prof%s(prof)%nearestlev_surf       ! highest lev below surf
    iv2lay(i) = iv2lev(i) - 1                          ! same layer as that numbered by  iv2 in RTTOV-9
    iv3lay(i) = iv3lev(i) - 1                          ! same layer as that numbered by  iv3 in RTTOV-9
  ENDDO

  IF (coef%id_sensor == sensor_id_po) THEN
    DO i = 1, nchanprof
      pol_id(i) = coef%fastem_polar(chan) + 1_jpim
    ENDDO
  ELSE
    pol_id(:) = 0_jpim
  ENDIF

  DO i = 1, nchanprof
    IF (do_lambertian(i) .AND. .NOT. dom_ir) THEN
      transmission_aux_ad%thermal_path1%tau_level_p(:,:,i)   = 0._jprb
      transmission_aux_ad%thermal_path1%tau_surf_p(:,i)      = 0._jprb
    ENDIF
  ENDDO

  !-------------------------------------------------------------------------------
  ! Calculate total radiance
  !-------------------------------------------------------------------------------
  cfraction_ad(:) = 0._jprb

  IF (do_scatt) THEN
    !---------------------------------------------------
    ! Calculate complex cloudy radiances
    !---------------------------------------------------
    DO i = 1, nchanprof
      IF (solar(i) .AND. do_mfasis) CYCLE
      IF (thermal(i) .AND. dom_ir .AND. .NOT. solar(i)) CYCLE

      rad_ad%cloudy(i) = rad_ad%total(i)

      auxrad_column_ad%cloudy(0,i) = auxrad_column_ad%cloudy(0,i) + &
                                     rad_ad%cloudy(i) * ircld%xcolclr(prof)
      ircld_ad%xcolclr(prof) = ircld_ad%xcolclr(prof) + rad_ad%cloudy(i) * auxrad_column%cloudy(0,i)

      DO col = ircld%ncolumn(prof), 1, -1 !reverse loop nec. due to data dependence
        auxrad_column_ad%cloudy(col,i) = auxrad_column_ad%cloudy(col,i) + &
          rad_ad%cloudy(i) * (ircld%xcol(col+1,prof) - ircld%xcol(col,prof))

        ircld_ad%xcol(col+1,prof) = ircld_ad%xcol(col+1,prof) + &
          rad_ad%cloudy(i) * auxrad_column%cloudy(col,i)

        ircld_ad%xcol(col,prof) = ircld_ad%xcol(col,prof) - &
          rad_ad%cloudy(i) * auxrad_column%cloudy(col,i)
      ENDDO
    ENDDO
  ELSE
    !---------------------------------------------------
    ! Calculate total radiance (clear case/simple cloud)
    !---------------------------------------------------
    IF (opts%rt_ir%pc%addpc) THEN
      rad_ad%clear(:)  = rad_ad%total(:)
    ELSE
      ! Interpolate to given cloud-top pressures
      DO i = 1, nchanprof
        IF (solar(i) .AND. do_mfasis) CYCLE

        rad_ad%clear(i)  = rad_ad%clear(i) + (1._jprb - cfraction(i)) * rad_ad%total(i)
        rad_ad%cloudy(i) = cfraction(i) * rad_ad%total(i)
        cfraction_ad(i)  = (rad%cloudy(i) - rad%clear(i)) * rad_ad%total(i)

        lay = aux_prof%s(prof)%nearestlev_ctp - 1

        rad_ad%overcast(lay,i) = & !rad_ad%overcast(lay,i) + & !first use of rad_ad%overcast
          (1._jprb - aux_prof%s(prof)%pfraction_ctp) * rad_ad%cloudy(i)

        rad_ad%overcast(lay-1,i) = & !rad_ad%overcast(lay-1,i) + & !first use of rad_ad%overcast
          aux_prof%s(prof)%pfraction_ctp * rad_ad%cloudy(i)

        aux_prof_ad%s(prof)%pfraction_ctp = aux_prof_ad%s(prof)%pfraction_ctp + &
          (rad%overcast(lay-1,i) - rad%overcast(lay,i)) * rad_ad%cloudy(i)
      ENDDO
    ENDIF
  ENDIF
  ! rad_ad%cloudy done
  WHERE (.NOT. (solar(1:nchanprof) .AND. do_mfasis)) &
    auxrad_column_ad%cloudy(0,1:nchanprof) = auxrad_column_ad%cloudy(0,1:nchanprof) + rad_ad%clear(1:nchanprof)

  !-------------------------------------------------------------------------------
  ! Calculate overcast radiances
  !-------------------------------------------------------------------------------
  IF (.NOT. keyradonly) THEN
    col = 0_jpim
    !cdir nodep
    DO i = 1, nchanprof
      IF (thermal(i) .AND. .NOT. dom_ir) THEN
        DO lay = 1, nlayers
          ! Make exception for iv2lay because up(iv2lay) was replaced by meanrad_up in direct and
          ! overcast(iv2lay) was replaced by surface contribution
          IF (lay == iv2lay(i)) THEN
            auxrad_column_ad%meanrad_up(col,i) = auxrad_column_ad%meanrad_up(col,i) + rad_ad%overcast(lay,i)
            auxrad_ad%surfair(i) = auxrad_ad%surfair(i) + rad_ad%overcast(lay,i) * &
              transmission_aux%thermal_path1%tau_surf(col,i)
            transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
              transmission_aux_ad%thermal_path1%tau_surf(col,i) + &
              rad_ad%overcast(lay,i) * auxrad%surfair(i)
          ELSE
            auxrad_column_ad%up(lay,col,i) = auxrad_column_ad%up(lay,col,i) + rad_ad%overcast(lay,i)
            auxrad_ad%air(lay+1,i) = &!auxrad_ad%air(lay+1,i) + &
              rad_ad%overcast(lay,i) * transmission_aux%thermal_path1%tau_level(lay+1,col,i)
            transmission_aux_ad%thermal_path1%tau_level(lay+1,col,i) = &
              transmission_aux_ad%thermal_path1%tau_level(lay+1,col,i) + &
              auxrad%air(lay+1,i) * rad_ad%overcast(lay,i)
          ENDIF
        ENDDO
      ELSE IF (solar(i) .AND. .NOT. do_mfasis .AND. .NOT. dom_vis) THEN
        IF (reflectance(i)%refl_cloud_top > 0) THEN
          refl = reflectance(i)%refl_cloud_top
        ELSE
          IF (coef%ff_cwn(chan) > overcast_albedo_wvn) THEN
            refl = overcast_albedo1 * pi_r
          ELSE
            refl = overcast_albedo2 * pi_r
          ENDIF
        ENDIF
        DO lay = 1, nlayers
          lev = lay + 1
          transmission_aux_ad%solar_path2%Tau_level(lev,col,i) = &
            transmission_aux_ad%solar_path2%Tau_level(lev,col,i) + &
            solar_spectrum(i) * refl / raytracing%pathsun(lay,prof) * &
            rad_ad%overcast(lay,i)
          raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - &
            solar_spectrum(i) * refl / raytracing%pathsun(lay,prof)**2_jpim * &
            transmission_aux%solar_path2%Tau_level(lev,col,i) * rad_ad%overcast(lay,i)

          transmission_aux_ad%solar_path1%Tau_level(lev,col,i) = &
            transmission_aux_ad%solar_path1%Tau_level(lev,col,i) + &
            refl / raytracing%pathsat(lay,prof) * &
            auxrad_column%down_solar(lay,col,i) * &
            transmission_aux%solar_path1%Tau_level(lev,col,i) * &
            2_jpim * rad_ad%overcast(lay,i)
          raytracing_ad%pathsat(lay,prof) = &
            raytracing_ad%pathsat(lay,prof) - &
            refl / raytracing%pathsat(lay,prof) * &
            transmission_aux%solar_path1%Tau_level(lev,col,i)**2_jpim * &
            auxrad_column%down_solar(lay,col,i) / raytracing%pathsat(lay,prof) * &
            rad_ad%overcast(lay,i)
          auxrad_column_ad%down_solar(lay,col,i) = &
            auxrad_column_ad%down_solar(lay,col,i) + &
            refl / raytracing%pathsat(lay,prof) * &
            transmission_aux%solar_path1%Tau_level(lev,col,i)**2_jpim * &
            rad_ad%overcast(lay,i)

          ! Make exception for iv2lay because up_solar(iv2lay) was replaced by meanrad_up_solar in direct
          IF (lay == iv2lay(i)) THEN
            auxrad_column_ad%meanrad_up_solar(col,i) = &
              auxrad_column_ad%meanrad_up_solar(col,i) + rad_ad%overcast(lay,i)
          ELSE
            auxrad_column_ad%up_solar(lay,col,i) = &
              auxrad_column_ad%up_solar(lay,col,i) + rad_ad%overcast(lay,i)
          ENDIF
        ENDDO
      ELSE
        rad_ad%overcast(:,i) = 0._jprb
      ENDIF
    ENDDO
    ! rad_ad%overcast done
    ! rad_ad%clear done

    DO i = 1, nchanprof
      aux_prof_ad%s(prof)%cfraction = aux_prof_ad%s(prof)%cfraction + cfraction_ad(i)
    ENDDO
  ENDIF

  rad_ad%clear(:) = 0._jprb

  !-------------------------------------------------------------------------------
  ! Solar surface contribution
  !-------------------------------------------------------------------------------
  IF (dosolar .AND. .NOT. (dom_vis .OR. do_mfasis)) &
    CALL solar_surface_contribution_ad(transmission_aux, transmission_aux_ad, ircld, &
                                       reflectance%refl_out, reflectance_ad%refl_out, &
                                       auxrad_column_ad)

  !-------------------------------------------------------------------------------
  ! Calculate surface emission contribution
  !-------------------------------------------------------------------------------
  IF (dothermal .AND. .NOT. dom_ir) THEN
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        DO col = 0, ircld%ncolumn(prof) !rev loop not necessary
          emissivity_ad(i)%emis_out = emissivity_ad(i)%emis_out + &
            auxrad_column_ad%cloudy(col,i) * auxrad%skin(i) * &
            transmission_aux%thermal_path1%tau_surf(col,i)

          transmission_aux_ad%thermal_path1%tau_surf(col,i) = transmission_aux_ad%thermal_path1%tau_surf(col,i) + &
            auxrad_column_ad%cloudy(col,i) * auxrad%skin(i) * emissivity(i)%emis_out

          auxrad_ad%skin(i) = auxrad_ad%skin(i) + &
            auxrad_column_ad%cloudy(col,i) * emissivity(i)%emis_out * &
            transmission_aux%thermal_path1%tau_surf(col,i)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !-------------------------------------------------------------------------------
  ! Add thermal and solar atmospheric contributions to the clear and cloudy columns
  !-------------------------------------------------------------------------------
  !cdir nodep
  DO i = 1, nchanprof
    IF (thermal(i) .AND. .NOT. dom_ir) THEN
      DO col = 0, ircld%ncolumn(prof)
        ! add upward and downward parts
        auxrad_column_ad%meanrad_up(col,i) = auxrad_column_ad%meanrad_up(col,i) + &
          auxrad_column_ad%cloudy(col,i)

        IF (do_lambertian(i)) THEN
          auxrad_column_ad%meanrad_down_p(col,i) = &!auxrad_column_ad%meanrad_down_p(col,i) + &
            diffuse_refl(i) * (1._jprb - emissivity(i)%specularity) * &
            transmission_aux%thermal_path1%tau_surf_p(col,i) * &
            transmission_aux%thermal_path1%tau_surf(col,i) * auxrad_column_ad%cloudy(col,i)

          auxrad_column_ad%meanrad_down(col,i) = &!auxrad_column_ad%meanrad_down(col,i) + &
            diffuse_refl(i) * emissivity(i)%specularity * &
            transmission_aux%thermal_path1%tau_surf(col,i)**2 * auxrad_column_ad%cloudy(col,i)

          diffuse_refl_ad(i) = diffuse_refl_ad(i) + &
            (emissivity(i)%specularity * &
             auxrad_column%meanrad_down(col,i) * &
             transmission_aux%thermal_path1%tau_surf(col,i) + &
             (1._jprb - emissivity(i)%specularity) * &
             auxrad_column%meanrad_down_p(col,i) * &
             transmission_aux%thermal_path1%tau_surf_p(col,i)) * &
            transmission_aux%thermal_path1%tau_surf(col,i) * auxrad_column_ad%cloudy(col,i)

          transmission_aux_ad%thermal_path1%tau_surf_p(col,i) = &
            transmission_aux_ad%thermal_path1%tau_surf_p(col,i) + &
            auxrad_column%meanrad_down_p(col,i) * &
            diffuse_refl(i) * (1._jprb - emissivity(i)%specularity) * &
            transmission_aux%thermal_path1%tau_surf(col,i) * auxrad_column_ad%cloudy(col,i)

          transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
            transmission_aux_ad%thermal_path1%tau_surf(col,i) + &
            auxrad_column%meanrad_down(col,i) * &
            diffuse_refl(i) * emissivity(i)%specularity * &
            transmission_aux%thermal_path1%tau_surf(col,i) * auxrad_column_ad%cloudy(col,i)

          transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
            transmission_aux_ad%thermal_path1%tau_surf(col,i) + &
            (emissivity(i)%specularity * &
             auxrad_column%meanrad_down(col,i) * &
             transmission_aux%thermal_path1%tau_surf(col,i) + &
             (1._jprb - emissivity(i)%specularity) * &
             auxrad_column%meanrad_down_p(col,i) * &
             transmission_aux%thermal_path1%tau_surf_p(col,i)) * &
            diffuse_refl(i) * auxrad_column_ad%cloudy(col,i)

          emissivity_ad(i)%specularity = emissivity_ad(i)%specularity + &
            (auxrad_column%meanrad_down(col,i) * transmission_aux%thermal_path1%tau_surf(col,i) - &
             auxrad_column%meanrad_down_p(col,i) * transmission_aux%thermal_path1%tau_surf_p(col,i)) * &
            transmission_aux%thermal_path1%tau_surf(col,i) * &
            diffuse_refl(i) * auxrad_column_ad%cloudy(col,i)
        ELSE
          auxrad_column_ad%meanrad_down(col,i) = & ! auxrad_column_ad%meanrad_down(col,i) + &
            auxrad_column_ad%cloudy(col,i) * transmission_aux%thermal_path1%tau_surf(col,i)**2 * &
            diffuse_refl(i)

          transmission_aux_ad%thermal_path1%tau_surf(col,i) = transmission_aux_ad%thermal_path1%tau_surf(col,i) + &
            auxrad_column_ad%cloudy(col,i) * &
            2._jprb * transmission_aux%thermal_path1%tau_surf(col,i) * diffuse_refl(i) * &
            auxrad_column%meanrad_down(col,i)

          diffuse_refl_ad(i) = diffuse_refl_ad(i) + &
            auxrad_column_ad%cloudy(col,i) * &
            transmission_aux%thermal_path1%tau_surf(col,i)**2 * &
            auxrad_column%meanrad_down(col,i)
        ENDIF
      ENDDO
    ENDIF

    IF (solar(i) .AND. .NOT. do_mfasis) THEN
      refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad) * pi_r
      DO col = 0, ircld%ncolumn(prof)

        auxrad_column_ad%meanrad_up_solar(col,i) = auxrad_column_ad%meanrad_up_solar(col,i) + &
          auxrad_column_ad%cloudy(col,i)

        diffuse_refl_ad(i) = diffuse_refl_ad(i) + &
          transmission_aux%solar_path1%Tau_surf(col,i)**2_jpim * refl_norm_scat * &
          auxrad_column%meanrad_down_solar(col,i) * auxrad_column_ad%cloudy(col,i)

        auxrad_column_ad%meanrad_down_solar(col,i) = &! meanrad_down_solar_ad(col,i) + &
          transmission_aux%solar_path1%Tau_surf(col,i)**2_jpim * &
          diffuse_refl(i) * refl_norm_scat * auxrad_column_ad%cloudy(col,i)

        transmission_aux_ad%solar_path1%Tau_surf(col,i) = transmission_aux_ad%solar_path1%Tau_surf(col,i) + &
          transmission_aux%solar_path1%Tau_surf(col,i) * &
          diffuse_refl(i) * refl_norm_scat * 2._jprb * auxrad_column_ad%cloudy(col,i) * &
          auxrad_column%meanrad_down_solar(col,i)
      ENDDO
    ENDIF
  ENDDO

  !-------------------------------------------------------------------------------
  ! Calculate clear-sky Rayleigh scattering contribution
  !-------------------------------------------------------------------------------
  IF (do_rayleigh_ss .AND. .NOT. do_mfasis) &
    CALL solar_rayleigh_ad(raytracing, raytracing_ad, ircld, profiles, profiles_ad, &
                           profiles_dry, profiles_dry_ad, transmission_aux, &
                           transmission_aux_ad, auxrad_column_ad)

  !-------------------------------------------------------------------------------
  ! Calculate near-surface layer contribution
  !-------------------------------------------------------------------------------
  IF (do_single_scatt) &
    CALL solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                       transmission_scatt_ir_ad, diffuse_refl, auxrad_column, auxrad_column_ad)

  IF (dothermal .AND. .NOT. dom_ir) &
    CALL calc_near_surf_contribution_ad(transmission_aux, transmission_aux_ad, auxrad, auxrad_ad, ircld, &
                                        auxrad_column_ad)

  !-------------------------------------------------------------------------------
  ! Calculate atmospheric contribution
  !-------------------------------------------------------------------------------
  IF (do_single_scatt) &
    CALL solar_scattering_air_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                 transmission_scatt_ir_ad, diffuse_refl, auxrad_column, auxrad_column_ad)

  IF (dothermal .AND. .NOT. dom_ir) &
    CALL calc_atmospheric_radiance_ad(transmission_aux, transmission_aux_ad, auxrad, &
                                      auxrad_ad, ircld, auxrad_column_ad)

  !-------------------------------------------------------------------------------
  ! Calculate layer radiances
  !-------------------------------------------------------------------------------
  IF (dothermal) CALL rttov_calcrad_ad(chanprof, profiles, profiles_ad, coef, thermal, auxrad, auxrad_ad)

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',1_jpim,ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE calc_atmospheric_radiance_ad(transmission_aux, transmission_aux_ad, &
                                          auxrad, auxrad_ad, ircld, auxrad_column_ad)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux), INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_radiance_aux),     INTENT(IN)    :: auxrad
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_ad
    TYPE(rttov_ircld),            INTENT(IN)    :: ircld
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column_ad

    INTEGER(jpim) :: i, col, coli, lay, lev, levm1
    REAL(jprb)    :: b1_ad, b2_ad, b3_ad, temp, ztemp
    REAL(jprb)    :: rad_air_avg(nlayers), rad_air_diff(nlayers)
    REAL(jprb)    :: dtau(nlayers)

#define tau_lev transmission_aux%thermal_path1%Tau_level(lay+1,col,i)
#define tau_levm1 transmission_aux%thermal_path1%Tau_level(lay,col,i)
#define tau_lev_r transmission_aux%thermal_path1%Tau_level_r(lay+1,col,i)
#define tau_levm1_r transmission_aux%thermal_path1%Tau_level_r(lay,col,i)
#define tausun_lev transmission_aux%solar_path2%Tau_level(lay+1,col,i)
#define tausun_levm1 transmission_aux%solar_path2%Tau_level(lay,col,i)
#define B1_2 rad_air_diff(lay) * dtau(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)

    DO i = 1, nchanprof
      IF (.NOT. thermal(i)) CYCLE

      IF (do_lambertian(i) .OR. .NOT. opts%rt_all%rad_down_lin_tau) THEN
        rad_air_avg = 0.5_jprb * (auxrad%air(1:nlayers,i) + auxrad%air(2:nlevels,i))
      ENDIF
      rad_air_diff = auxrad%air(2:nlevels,i) - auxrad%air(1:nlayers,i)

      DO col = 0, ircld%ncolumn(prof)

        dtau = transmission_aux%thermal_path1%Tau_level(1:nlayers,col,i) - &
               transmission_aux%thermal_path1%Tau_level(2:nlevels,col,i)

        IF (do_lambertian(i)) THEN
          WHERE (auxrad_column%down_p_ref(:,col,i) < 0._jprb)
            auxrad_column_ad%down_p(:,col,i) = 0._jprb
          ENDWHERE
          DO lay = nlayers, 2, -1
            IF (auxrad_column%down_p_ref(lay,col,i) < 0._jprb) auxrad_column_ad%down_p(lay,col,i) = 0._jprb
            auxrad_column_ad%down_p(lay-1,col,i) = auxrad_column_ad%down_p(lay-1,col,i) + auxrad_column_ad%down_p(lay,col,i)
          ENDDO
        ENDIF

        WHERE (auxrad_column%down_ref(:,col,i) < 0._jprb)
          auxrad_column_ad%down(:,col,i) = 0._jprb
        ENDWHERE
        DO lay = nlayers, 2, -1
          auxrad_column_ad%up(lay-1,col,i) = auxrad_column_ad%up(lay-1,col,i) + auxrad_column_ad%up(lay,col,i)
          auxrad_column_ad%down(lay-1,col,i) = auxrad_column_ad%down(lay-1,col,i) + auxrad_column_ad%down(lay,col,i)
        ENDDO

        DO lay = 1, nlayers
          levm1 = lay
          lev = lay + 1
          coli = ircld%icldarr(col,lay,prof)
          IF (transmission_aux%fac1(lay,col,i) > 0._jprb) THEN
            IF (do_lambertian(i)) THEN
              temp = transmission_aux%thermal_path1%fac2(lay+1,col,i) * &
                     auxrad_column_ad%down_p(lay,col,i)

              ztemp = 0.5_jprb * temp * &
                (transmission_aux%thermal_path1%tau_level_p_r(lev,col,i) - &
                 transmission_aux%thermal_path1%tau_level_p_r(lay,col,i))

              auxrad_ad%air(lay,i) = auxrad_ad%air(lay,i) + ztemp
              auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) + ztemp

              transmission_aux_ad%thermal_path1%tau_level_p(lev,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level_p(lev,col,i) - &
                transmission_aux%thermal_path1%tau_level_p_r(lev,col,i)**2 * temp * rad_air_avg(lay)
              transmission_aux_ad%thermal_path1%tau_level_p(lay,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level_p(lay,col,i) + &
                transmission_aux%thermal_path1%tau_level_p_r(lay,col,i)**2 * temp * rad_air_avg(lay)
            ENDIF

            IF (opts%rt_all%rad_down_lin_tau) THEN
              ! Upwelling and linear-in-tau specular downwelling
              temp = transmission_aux%thermal_path1%fac2(lay+1,col,i) * &
                     auxrad_column_ad%down(lay,col,i)
              B1_AD = temp  * (tau_levm1_r * tau_lev_r)
              B2_AD = temp * (tau_lev_r)**2
              B3_AD = -temp * (tau_levm1_r * tau_lev_r)

              B1_AD = B1_AD + auxrad_column_ad%up(lay,col,i)
              B2_AD = B2_AD - auxrad_column_ad%up(lay,col,i)
              B3_AD = B3_AD + auxrad_column_ad%up(lay,col,i)

              transmission_aux_ad%thermal_path1%tau_level(levm1,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level(levm1,col,i) + &
                B1_AD * auxrad%air(levm1,i) - &
                temp * &
                (auxrad%air(levm1,i) * dtau(lay) - &
                B1_2) * &
                (tau_lev_r * &
                tau_levm1_r**2) + &
                B3_AD * rad_air_diff(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)

              transmission_aux_ad%thermal_path1%tau_level(lev,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level(lev,col,i) - &
                temp * &
                (auxrad%air(levm1,i) * dtau(lay) - &
                B1_2) * &
                (tau_lev_r**2 * tau_levm1_r) - &
                B3_AD * rad_air_diff(lay) * &
                transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i) - &
                temp * 2._jprb * rad_air_diff(lay) * tau_lev * &
                (tau_lev_r**3) - &
                B1_AD * auxrad%air(levm1,i) + &
                B2_AD * rad_air_diff(lay)

              transmission_aux_ad%thermal_path1%od_singlelayer(coli,lay,i) = &
                transmission_aux_ad%thermal_path1%od_singlelayer(coli,lay,i) - &
                B3_AD * &
                B1_2 * &
                transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)

              auxrad_ad%air(levm1,i) = auxrad_ad%air(levm1,i) - &
                B3_AD * dtau(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i) + &
                B1_AD * dtau(lay) - &
                B2_AD * tau_lev

              auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) + &
                B3_AD * dtau(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i) + &
                B2_AD * tau_lev
            ELSE
              ! Upwelling and layer-average specular downwelling
              temp = transmission_aux%thermal_path1%fac2(lay+1,col,i) * &
                     auxrad_column_ad%down(lay,col,i)

              ztemp = 0.5_jprb * temp * &
                (transmission_aux%thermal_path1%tau_level_r(lev,col,i) - &
                 transmission_aux%thermal_path1%tau_level_r(lay,col,i))

              transmission_aux_ad%thermal_path1%tau_level(levm1,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level(levm1,col,i) + &
                auxrad_column_ad%up(lay,col,i) * (auxrad%air(levm1,i) + &
                rad_air_diff(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)) + &
                transmission_aux%thermal_path1%tau_level_r(lay,col,i)**2 * temp * rad_air_avg(lay)

              transmission_aux_ad%thermal_path1%tau_level(lev,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level(lev,col,i) - &
                auxrad_column_ad%up(lay,col,i) * (rad_air_diff(lay) * &
                 (transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i) + 1._jprb) + &
                 auxrad%air(levm1,i)) - &
                transmission_aux%thermal_path1%tau_level_r(lev,col,i)**2 * temp * rad_air_avg(lay)

              transmission_aux_ad%thermal_path1%od_singlelayer(coli,lay,i) = &
                transmission_aux_ad%thermal_path1%od_singlelayer(coli,lay,i) - &
                auxrad_column_ad%up(lay,col,i) * &
                B1_2 * &
                transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)

              auxrad_ad%air(levm1,i) = auxrad_ad%air(levm1,i) + &
                auxrad_column_ad%up(lay,col,i) * (dtau(lay) * &
                 (1._jprb - transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i)) + &
                 tau_lev) + ztemp

              auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) + &
                auxrad_column_ad%up(lay,col,i) * &
                (dtau(lay) * transmission_aux%thermal_path1%od_singlelayer_r(coli,lay,i) - &
                 tau_lev) + ztemp
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE calc_atmospheric_radiance_ad

  SUBROUTINE calc_near_surf_contribution_ad(transmission_aux, transmission_aux_ad, auxrad, auxrad_ad, ircld, &
                                            auxrad_column_ad)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux), INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_radiance_aux),     INTENT(IN)    :: auxrad
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_ad
    TYPE(rttov_ircld),            INTENT(IN)    :: ircld
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column_ad

    INTEGER(jpim) :: i, col, lay, lev
    REAL(jprb)    :: b1_3, b2_3, b3_3, b1_ad, b2_ad, b3_ad, ztemp
    REAL(jprb)    :: tau_surf_p_r_ad, tau_lev_p_r_ad
    REAL(jprb)    :: tau_surf_r_ad, tau_lev_r_ad
    REAL(jprb)    :: rad_air_avg, rad_air_diff

    DO i = 1, nchanprof
      IF (.NOT. thermal(i)) CYCLE
      lay = iv3lay(i)
      lev = iv3lev(i)

      rad_air_avg = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i))
      rad_air_diff = auxrad%surfair(i) - auxrad%air(lev,i)

      DO col = 0, ircld%ncolumn(prof) !rev loop not necessary
        ! assume there is no atmospheric source term for 3rd/4th stokes vector elements
        IF (pol_id(i) >= 6_jpim) auxrad_column_ad%meanrad_up(col,i) = 0._jprb

        auxrad_column_ad%up(lay,col,i) = auxrad_column_ad%up(lay,col,i) + &
                            auxrad_column_ad%meanrad_up(col,i)

        auxrad_column_ad%down(lay,col,i) = auxrad_column_ad%down(lay,col,i) + &
                              auxrad_column_ad%meanrad_down(col,i)

        IF (do_lambertian(i)) THEN
          auxrad_column_ad%down_p(lay,col,i) = auxrad_column_ad%down_p(lay,col,i) + &
                                auxrad_column_ad%meanrad_down_p(col,i)
        ENDIF

        IF (transmission_aux%thermal_path1%od_sfrac(col,i) < min_od .OR. &
            (opts%rt_all%dtau_test .AND. &
             (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
              transmission_aux%thermal_path1%tau_surf(col,i)) < min_od)) THEN
        ELSE
          IF (transmission_aux%thermal_path1%tau_surf(col,i) > min_tau) THEN
            IF (do_lambertian(i)) THEN
              tau_surf_p_r_ad = rad_air_avg * auxrad_column_ad%meanrad_down_p(col,i)
              tau_lev_p_r_ad = -tau_surf_p_r_ad

              ztemp = 0.5_jprb * auxrad_column_ad%meanrad_down_p(col,i) * &
                (transmission_aux%thermal_path1%tau_surf_p_r(col,i) - &
                 transmission_aux%thermal_path1%tau_level_p_r(lev,col,i))

              auxrad_ad%surfair(i) = auxrad_ad%surfair(i) + ztemp
              auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) + ztemp

              transmission_aux_ad%thermal_path1%tau_surf_p(col,i) = &
                transmission_aux_ad%thermal_path1%tau_surf_p(col,i) - &
                transmission_aux%thermal_path1%tau_surf_p_r(col,i)**2 * tau_surf_p_r_ad
              transmission_aux_ad%thermal_path1%tau_level_p(lev,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level_p(lev,col,i) - &
                transmission_aux%thermal_path1%tau_level_p_r(lev,col,i)**2 * tau_lev_p_r_ad

              B1_AD = 0._jprb
              B2_AD = 0._jprb
              B3_AD = 0._jprb
            ENDIF

            IF (opts%rt_all%rad_down_lin_tau) THEN
              !direct calc
              B1_3 = auxrad%air(lev,i) * &
                    (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
                    transmission_aux%thermal_path1%tau_surf(col,i))

              B2_3 = rad_air_diff * transmission_aux%thermal_path1%tau_surf(col,i)

              B3_3 = rad_air_diff * &
                     (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
                     transmission_aux%thermal_path1%tau_surf(col,i)) * &
                     (transmission_aux%thermal_path1%od_sfrac_r(col,i))
              !first use B1,b2,b3c
              B1_AD = auxrad_column_ad%meanrad_down(col,i) * &
                      (transmission_aux%thermal_path1%tau_level_r(lev,col,i) * &
                      transmission_aux%thermal_path1%tau_surf_r(col,i))
              B3_AD = -auxrad_column_ad%meanrad_down(col,i) * &
                      (transmission_aux%thermal_path1%tau_level_r(lev,col,i) * &
                      transmission_aux%thermal_path1%tau_surf_r(col,i))
              B2_AD = auxrad_column_ad%meanrad_down(col,i) * &
                      (transmission_aux%thermal_path1%tau_surf_r(col,i))**2

              transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
                    transmission_aux_ad%thermal_path1%tau_surf(col,i) - &
                    auxrad_column_ad%meanrad_down(col,i) * 2._jprb * B2_3 * &
                    transmission_aux%thermal_path1%tau_surf_r(col,i)**3_jpim

              transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
                    transmission_aux_ad%thermal_path1%tau_surf(col,i) - &
                    auxrad_column_ad%meanrad_down(col,i) * (B1_3 - B3_3) * &
                    transmission_aux%thermal_path1%tau_level_r(lev,col,i) * &
                    transmission_aux%thermal_path1%tau_surf_r(col,i)**2_jpim

              transmission_aux_ad%thermal_path1%tau_level(lev,col,i) = &
                    transmission_aux_ad%thermal_path1%tau_level(lev,col,i) - &
                    auxrad_column_ad%meanrad_down(col,i) * (B1_3 - B3_3) * &
                    transmission_aux%thermal_path1%tau_level_r(lev,col,i)**2_jpim * &
                    transmission_aux%thermal_path1%tau_surf_r(col,i)
            ELSE
              ! Downwelling

              tau_surf_r_ad = rad_air_avg * auxrad_column_ad%meanrad_down(col,i)
              tau_lev_r_ad = -tau_surf_r_ad

              ztemp = 0.5_jprb * auxrad_column_ad%meanrad_down(col,i) * &
                (transmission_aux%thermal_path1%tau_surf_r(col,i) - &
                 transmission_aux%thermal_path1%tau_level_r(lev,col,i))

              auxrad_ad%surfair(i) = auxrad_ad%surfair(i) + ztemp
              auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) + ztemp

              transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
                transmission_aux_ad%thermal_path1%tau_surf(col,i) - &
                transmission_aux%thermal_path1%tau_surf_r(col,i)**2 * tau_surf_r_ad
              transmission_aux_ad%thermal_path1%tau_level(lev,col,i) = &
                transmission_aux_ad%thermal_path1%tau_level(lev,col,i) - &
                transmission_aux%thermal_path1%tau_level_r(lev,col,i)**2 * tau_lev_r_ad

              B1_AD = 0._jprb
              B2_AD = 0._jprb
              B3_AD = 0._jprb
            ENDIF
          ELSE
            B1_AD = 0._jprb
            B2_AD = 0._jprb
            B3_AD = 0._jprb
          ENDIF

          B1_AD = B1_AD + auxrad_column_ad%meanrad_up(col,i)
          B2_AD = B2_AD - auxrad_column_ad%meanrad_up(col,i)
          B3_AD = B3_AD + auxrad_column_ad%meanrad_up(col,i)

          transmission_aux_ad%thermal_path1%od_sfrac(col,i) = &
                                                transmission_aux_ad%thermal_path1%od_sfrac(col,i) - &
                                                B3_AD * rad_air_diff * &
                                               (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
                                                transmission_aux%thermal_path1%tau_surf(col,i)) * &
                                               (transmission_aux%thermal_path1%od_sfrac_r(col,i)**2)

          transmission_aux_ad%thermal_path1%tau_surf(col,i) = &
                                                transmission_aux_ad%thermal_path1%tau_surf(col,i) - &
                                                B3_AD * rad_air_diff * &
                                                transmission_aux%thermal_path1%od_sfrac_r(col,i) + &
                                                B2_AD * rad_air_diff - &
                                                B1_AD * auxrad%air(lev,i)

          transmission_aux_ad%thermal_path1%tau_level(lev,col,i) = &
                                                transmission_aux_ad%thermal_path1%tau_level(lev,col,i) + &
                                                B3_AD * rad_air_diff * &
                                                transmission_aux%thermal_path1%od_sfrac_r(col,i) + &
                                                B1_AD * auxrad%air(lev,i)

          auxrad_ad%air(lev,i) = auxrad_ad%air(lev,i) - &
               B3_AD * (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
               transmission_aux%thermal_path1%tau_surf(col,i)) * &
               transmission_aux%thermal_path1%od_sfrac_r(col,i) - &
               B2_AD * transmission_aux%thermal_path1%tau_surf(col,i) + &
               B1_AD * (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
               transmission_aux%thermal_path1%tau_surf(col,i))

          auxrad_ad%surfair(i) = auxrad_ad%surfair(i) + &
               B3_AD * (transmission_aux%thermal_path1%tau_level(lev,col,i) - &
               transmission_aux%thermal_path1%tau_surf(col,i)) * &
               transmission_aux%thermal_path1%od_sfrac_r(col,i) + &
               B2_AD * transmission_aux%thermal_path1%tau_surf(col,i)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_near_surf_contribution_ad

  SUBROUTINE solar_scattering_air_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                     transmission_scatt_ir_ad, refl, auxrad_column, auxrad_column_ad)

    TYPE(rttov_transmission_aux),      INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
    TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
    TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_ad
    TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_ad
    REAL(jprb),                        INTENT(IN)    :: refl(:)
    TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad_column
    TYPE(rttov_radiance_aux),          INTENT(INOUT) :: auxrad_column_ad

    INTEGER(jpim) :: i, col, coli, lay, su
    REAL(jprb)    :: temp_ad(1:nlayers,0:maxncolumns), tmpval
    REAL(jprb)    :: fac1_2_ad(0:1,1:nlayers), fac3_2_ad(1:nlayers), fac4_2_ad(0:1,1:nlayers), &
                     fac5_2_ad(0:1,1:nlayers), fac6_2_ad(0:1,1:nlayers), fac7_2_ad(1:nlayers)

    su = 0
    IF (opts%rt_ir%addclouds) su = 1

    DO i = 1, nchanprof
      IF (.NOT. solar(i)) CYCLE

      Fac1_2_ad = 0._jprb
      Fac3_2_ad = 0._jprb
      Fac4_2_ad = 0._jprb
      Fac5_2_ad = 0._jprb
      Fac6_2_ad = 0._jprb
      Fac7_2_ad = 0._jprb

      !-------------------Downward single scattering of the solar beam------------------

      IF (refl(i) > 0._jprb) THEN

        temp_ad = 0._jprb
        DO col = 0, ircld%ncolumn(prof)
          DO lay = nlayers, 2, -1
            IF (auxrad_column%down_ref_solar(lay,col,i) < 0._jprb) THEN
              auxrad_column_ad%down_solar(lay,col,i) = 0._jprb
            ELSE
              temp_ad(lay,col) = temp_ad(lay,col) + auxrad_column_ad%down_solar(lay,col,i)
            ENDIF
            temp_ad(lay-1,col) = temp_ad(lay-1,col) + temp_ad(lay,col)
          ENDDO

          lay = 1_jpim
          IF (auxrad_column%down_ref_solar(lay,col,i) < 0._jprb) THEN
            auxrad_column_ad%down_solar(lay,col,i) = 0._jprb
          ELSE
            temp_ad(lay,col) = temp_ad(lay,col) + auxrad_column_ad%down_solar(lay,col,i)
          ENDIF
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          DO lay = nlayers, 1, -1
            coli = ircld%icldarr(col,lay,prof)

            tmpval = &
                transmission_aux%solar_path1%fac2(lay+1,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay+1,col,i)

            IF (.NOT. sateqsun(lay,prof)) THEN
              Fac6_2_ad(coli,lay) = Fac6_2_ad(coli,lay) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i)) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              transmission_scatt_ir_ad%ssa_solar(coli,lay,i) = &
                transmission_scatt_ir_ad%ssa_solar(coli,lay,i) + temp_ad(lay,col) * tmpval * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i)) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              Fac7_2_ad(lay) = Fac7_2_ad(lay) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i)) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              Fac5_2_ad(coli,lay) = Fac5_2_ad(coli,lay) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              Fac4_2_ad(coli,lay) = Fac4_2_ad(coli,lay) - temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              transmission_aux_ad%solar_path2%Tau_level(lay,col,i) = &
                transmission_aux_ad%solar_path2%Tau_level(lay,col,i) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i))

              transmission_aux_ad%solar_path1%Tau_level(lay,col,i) = &
                transmission_aux_ad%solar_path1%Tau_level(lay,col,i) - temp_ad(lay,col) * tmpval * &
                transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay+1,col,i) * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i)) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level(lay+1,col,i)

              transmission_aux_ad%solar_path1%Tau_level(lay+1,col,i) = &
                transmission_aux_ad%solar_path1%Tau_level(lay+1,col,i) - temp_ad(lay,col) * tmpval * &
                transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay+1,col,i) * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac7_2(lay,i) * &
                (auxrad_column%Fac5_2(coli,lay,i) - auxrad_column%Fac4_2(coli,lay,i)) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level(lay,col,i)

            ELSE
              Fac6_2_ad(coli,lay) = Fac6_2_ad(coli,lay) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              transmission_scatt_ir_ad%ssa_solar(coli,lay,i) = &
                transmission_scatt_ir_ad%ssa_solar(coli,lay,i) + temp_ad(lay,col) * tmpval * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              transmission_aux_ad%solar_path2%Od_singlelayer(coli,lay,i) = &
                transmission_aux_ad%solar_path2%Od_singlelayer(coli,lay,i) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              Fac4_2_ad(coli,lay) = Fac4_2_ad(coli,lay) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i)

              transmission_aux_ad%solar_path2%Tau_level(lay,col,i) = &
                transmission_aux_ad%solar_path2%Tau_level(lay,col,i) + temp_ad(lay,col) * tmpval * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i)

              transmission_aux_ad%solar_path1%Tau_level(lay,col,i) = &
                transmission_aux_ad%solar_path1%Tau_level(lay,col,i) - temp_ad(lay,col) * tmpval * &
                transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay+1,col,i) * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level(lay+1,col,i)

              transmission_aux_ad%solar_path1%Tau_level(lay+1,col,i) = &
                transmission_aux_ad%solar_path1%Tau_level(lay+1,col,i) - temp_ad(lay,col) * tmpval * &
                transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level_r(lay+1,col,i) * &
                transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                transmission_aux%solar_path2%Od_singlelayer(coli,lay,i) * &
                auxrad_column%Fac6_2(coli,lay,i) * &
                auxrad_column%Fac4_2(coli,lay,i) * &
                transmission_aux%solar_path2%Tau_level(lay,col,i) * &
                transmission_aux%solar_path1%Tau_level(lay,col,i)
            ENDIF
          ENDDO
        ENDDO

        DO lay = 1, nlayers
          IF (.NOT. sateqsun(lay,prof)) THEN
            raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + Fac7_2_ad(lay) * &
              raytracing%pathsun(lay,prof) / (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2_jpim
            raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - Fac7_2_ad(lay) * &
              raytracing%pathsat(lay,prof) / (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2_jpim
          ENDIF
        ENDDO
      ENDIF ! refl(i) > 0.

      !----------------Upward single scattering of the solar beam-----------------------

      DO col = 0, ircld%ncolumn(prof)
        DO lay = nlayers, 2, -1
          auxrad_column_ad%up_solar(lay-1,col,i) = auxrad_column_ad%up_solar(lay-1,col,i) + &
                                                   auxrad_column_ad%up_solar(lay,col,i)
        ENDDO
      ENDDO

      DO col = 0, ircld%ncolumn(prof)
        DO lay = nlayers, 1, -1
          coli = ircld%icldarr(col,lay,prof)

          Fac1_2_ad(coli,lay) = Fac1_2_ad(coli,lay) + auxrad_column_ad%up_solar(lay,col,i) * &
            (transmission_aux%solar_path2%Tau_level(lay,col,i) - &
             transmission_aux%solar_path2%Tau_level(lay+1,col,i)) * &
            transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%Fac3_2(lay,i)

          transmission_scatt_ir_ad%ssa_solar(coli,lay,i) = &
            transmission_scatt_ir_ad%ssa_solar(coli,lay,i) + auxrad_column_ad%up_solar(lay,col,i) * &
            (transmission_aux%solar_path2%Tau_level(lay,col,i) - &
             transmission_aux%solar_path2%Tau_level(lay+1,col,i)) * &
            auxrad_column%Fac1_2(coli,lay,i) * auxrad_column%Fac3_2(lay,i)

          Fac3_2_ad(lay) = Fac3_2_ad(lay) + auxrad_column_ad%up_solar(lay,col,i) * &
            (transmission_aux%solar_path2%Tau_level(lay,col,i) - &
             transmission_aux%solar_path2%Tau_level(lay+1,col,i)) * &
            transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%Fac1_2(coli,lay,i)

          transmission_aux_ad%solar_path2%Tau_level(lay,col,i) = &
            transmission_aux_ad%solar_path2%Tau_level(lay,col,i) + auxrad_column_ad%up_solar(lay,col,i) * &
            auxrad_column%Fac1_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%Fac3_2(lay,i)

          transmission_aux_ad%solar_path2%Tau_level(lay+1,col,i) = &
            transmission_aux_ad%solar_path2%Tau_level(lay+1,col,i) - auxrad_column_ad%up_solar(lay,col,i) * &
            auxrad_column%Fac1_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%Fac3_2(lay,i)
        ENDDO
      ENDDO

      raytracing_ad%pathsat(:,prof) = raytracing_ad%pathsat(:,prof) + &
        Fac3_2_ad(:) / raytracing%patheff(:,prof)
      raytracing_ad%patheff(:,prof) = raytracing_ad%patheff(:,prof) - &
        Fac3_2_ad(:) * raytracing%pathsat(:,prof) / raytracing%patheff(:,prof)**2_jpim

      transmission_aux_ad%solar_path1%Od_singlelayer(0:su,:,i) = &
        transmission_aux_ad%solar_path1%Od_singlelayer(0:su,:,i) - &
        Fac5_2_ad(0:su,:) * auxrad_column%Fac5_2(0:su,:,i)

      transmission_aux_ad%solar_path2%Od_singlelayer(0:su,:,i) = &
        transmission_aux_ad%solar_path2%Od_singlelayer(0:su,:,i) - &
        Fac4_2_ad(0:su,:) * auxrad_column%Fac4_2(0:su,:,i)

      transmission_scatt_ir_ad%phup(0:su,:,i) = &
        transmission_scatt_ir_ad%phup(0:su,:,i) + &
        solar_spectrum(i) * z4pi_r * Fac1_2_ad(0:su,:)

      transmission_scatt_ir_ad%phdo(0:su,:,i) = &
        transmission_scatt_ir_ad%phdo(0:su,:,i) + &
        solar_spectrum(i) * z4pi_r * Fac6_2_ad(0:su,:)
    ENDDO

  END SUBROUTINE solar_scattering_air_ad

  SUBROUTINE solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, ircld, raytracing, raytracing_ad, &
                                           transmission_scatt_ir_ad, refl, auxrad_column, auxrad_column_ad)

    TYPE(rttov_transmission_aux),      INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux),      INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_ircld),                 INTENT(IN)    :: ircld
    TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
    TYPE(rttov_raytracing),            INTENT(INOUT) :: raytracing_ad
    TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: transmission_scatt_ir_ad
    REAL(jprb),                        INTENT(IN)    :: refl(:)
    TYPE(rttov_radiance_aux),          INTENT(IN)    :: auxrad_column
    TYPE(rttov_radiance_aux),          INTENT(INOUT) :: auxrad_column_ad

    INTEGER(jpim) :: i, col, coli, lay, lev, lay1
    REAL(jprb)    :: fac_3_ad(7), temp(4,0:maxncolumns)

    DO i = 1, nchanprof
      IF (.NOT. solar(i)) CYCLE
      lay = iv3lay(i)
      lev = iv3lev(i)

      IF (pfraction(i) < 0._jprb) THEN
         lay1 = iv3lay(i)
      ELSE
         lay1 = iv3lay(i) + 1
      ENDIF

      DO col = 0, ircld%ncolumn(prof) !rev loop not necessary
        coli = ircld%icldarr(col,lay,prof)

        !-------------------Downward single scattering of the solar beam------------------

        IF (refl(i) > 0._jprb) THEN

          IF (auxrad_column%meanrad_down_solar(col,i) < 0._jprb) THEN
            auxrad_column_ad%meanrad_down_solar(col,i) = 0._jprb
          ENDIF

          auxrad_column_ad%down_solar(lay,col,i) = auxrad_column_ad%down_solar(lay,col,i) + &
                                                   auxrad_column_ad%meanrad_down_solar(col,i)

          fac_3_ad(:) = 0._jprb
          temp(1,col) = auxrad_column%fac6_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                        auxrad_column_ad%meanrad_down_solar(col,i)
          temp(2,col) = transmission_aux%solar_path1%tau_level_r(lev,col,i) * &
                        transmission_aux%solar_path1%tau_surf_r(col,i)
          ! By doing the multiplication in the following way we avoid overflows
          temp(3,col) = (transmission_aux%solar_path1%fac2(lev,col,i) * temp(2,col)) * temp(2,col)

          IF (.NOT. sateqsun(lay,prof)) THEN
            temp(4,col) = (auxrad_column%fac5_3(col,i) - auxrad_column%fac4_3(col,i)) * auxrad_column%fac7_2(lay1,i)

            fac_3_ad(4) = &!fac4_3_ad &
                        - temp(1,col) * temp(2,col) * auxrad_column%fac7_2(lay1,i) * &
                          transmission_aux%solar_path2%tau_level(lev,col,i)

            fac_3_ad(5) = &!fac5_3_ad + &
                        + temp(1,col) * temp(2,col) * auxrad_column%fac7_2(lay1,i) * &
                          transmission_aux%solar_path2%tau_level(lev,col,i)

            fac_3_ad(7) = &!fac7_3_ad + &
                        temp(1,col) * temp(2,col) * &
                       (auxrad_column%fac5_3(col,i) - auxrad_column%fac4_3(col,i)) * &
                        transmission_aux%solar_path2%tau_level(lev,col,i)

            raytracing_ad%pathsat(lay1,prof) = raytracing_ad%pathsat(lay1,prof) + &
                                               fac_3_ad(7) * raytracing%pathsun(lay1,prof) / &
                                              (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2

            raytracing_ad%pathsun(lay1,prof) = raytracing_ad%pathsun(lay1,prof) - &
                                               fac_3_ad(7) * raytracing%pathsat(lay1,prof) / &
                                              (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2
          ELSE
            temp(4,col) = auxrad_column%fac4_3(col,i) * transmission_aux%solar_path2%od_sfrac(col,i)

            transmission_aux_ad%solar_path2%od_sfrac(col,i) = transmission_aux_ad%solar_path2%od_sfrac(col,i) + &
                                                     temp(1,col) * temp(2,col) * &
                                                     auxrad_column%fac4_3(col,i) * &
                                                     transmission_aux%solar_path2%tau_level(lev,col,i)

            fac_3_ad(4) = &!fac4_3_ad + &
                        temp(1,col) * temp(2,col) * &
                        transmission_aux%solar_path2%tau_level(lev,col,i) * &
                        transmission_aux%solar_path2%od_sfrac(col,i)
          ENDIF

          IF (transmission_aux%solar_path1%tau_level(lev,col,i) > min_tau) THEN

            transmission_aux_ad%solar_path1%tau_surf(col,i) = transmission_aux_ad%solar_path1%tau_surf(col,i) - &
                                                  temp(1,col) * temp(3,col) * temp(4,col) * &
                                                  transmission_aux%solar_path2%tau_level(lev,col,i) * &
                                                  transmission_aux%solar_path1%tau_level(lev,col,i)

            transmission_aux_ad%solar_path1%tau_level(lev,col,i) = &
                                                       transmission_aux_ad%solar_path1%tau_level(lev,col,i) - &
                                                       temp(1,col) * temp(3,col) * temp(4,col) * &
                                                       transmission_aux%solar_path2%tau_level(lev,col,i) * &
                                                       transmission_aux%solar_path1%tau_surf(col,i)

            transmission_aux_ad%solar_path2%tau_level(lev,col,i) = &
                                                       transmission_aux_ad%solar_path2%tau_level(lev,col,i) + &
                                                       temp(1,col) * temp(2,col) * temp(4,col)

            fac_3_ad(2) = &!fac2_3_ad + &
                        auxrad_column_ad%meanrad_down_solar(col,i) * auxrad_column%fac6_2(coli,lay,i) * &
                        temp(2,col) * temp(4,col) * transmission_aux%solar_path2%tau_level(lev,col,i)

            fac_3_ad(6) = &!fac6_3_ad + &
                        auxrad_column_ad%meanrad_down_solar(col,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                        temp(2,col) * temp(4,col) * transmission_aux%solar_path2%tau_level(lev,col,i)
          ENDIF

        ELSE
          fac_3_ad(2) = 0._jprb
          fac_3_ad(4) = 0._jprb
          fac_3_ad(5) = 0._jprb
          fac_3_ad(6) = 0._jprb
        ENDIF ! refl(i) > 0.

        !--------------Upward single scattering of the solar beam---------------------------------

        auxrad_column_ad%up_solar(lay,col,i) = auxrad_column_ad%up_solar(lay,col,i) + &
                                               auxrad_column_ad%meanrad_up_solar(col,i)

        transmission_aux_ad%solar_path2%tau_level(lev,col,i) = &
                   transmission_aux_ad%solar_path2%tau_level(lev,col,i) + &
                   auxrad_column_ad%meanrad_up_solar(col,i) * &
                   auxrad_column%fac1_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%fac3_2(lay,i)

        transmission_aux_ad%solar_path2%tau_surf(col,i) = &
                   transmission_aux_ad%solar_path2%tau_surf(col,i) - &
                   auxrad_column_ad%meanrad_up_solar(col,i) * &
                   auxrad_column%fac1_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%fac3_2(lay,i)

        fac_3_ad(3) = &!fac3_3_ad + &
                    auxrad_column_ad%meanrad_up_solar(col,i) * &
                    auxrad_column%fac1_2(coli,lay,i) * transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                   (transmission_aux%solar_path2%tau_level(lev,col,i) - &
                    transmission_aux%solar_path2%tau_surf(col,i))

        fac_3_ad(2) = fac_3_ad(2) + &
                    auxrad_column_ad%meanrad_up_solar(col,i) * &
                    auxrad_column%fac1_2(coli,lay,i) * auxrad_column%fac3_2(lay,i) * &
                   (transmission_aux%solar_path2%tau_level(lev,col,i) - &
                    transmission_aux%solar_path2%tau_surf(col,i))

        fac_3_ad(1) = &!fac1_3_ad + &
                    auxrad_column_ad%meanrad_up_solar(col,i) * &
                    transmission_scatt_ir%ssa_solar(coli,lay,i) * auxrad_column%fac3_2(lay,i) * &
                   (transmission_aux%solar_path2%tau_level(lev,col,i) - &
                    transmission_aux%solar_path2%tau_surf(col,i))

        transmission_scatt_ir_ad%phdo(coli,lay,i) = &
             transmission_scatt_ir_ad%phdo(coli,lay,i) + &
             fac_3_ad(6) * solar_spectrum(i) * z4pi_r

        transmission_aux_ad%solar_path1%od_sfrac(col,i) = transmission_aux_ad%solar_path1%od_sfrac(col,i) - &
             fac_3_ad(5) * auxrad_column%fac5_3(col,i)

        transmission_aux_ad%solar_path2%od_sfrac(col,i) = transmission_aux_ad%solar_path2%od_sfrac(col,i) - &
             fac_3_ad(4) * auxrad_column%fac4_3(col,i)

        raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
             fac_3_ad(3) / raytracing%patheff(lay,prof)

        raytracing_ad%patheff(lay,prof) = raytracing_ad%patheff(lay,prof) - &
             fac_3_ad(3) * raytracing%pathsat(lay,prof) / raytracing%patheff(lay,prof) ** 2

        transmission_scatt_ir_ad%ssa_solar(coli,lay,i) = transmission_scatt_ir_ad%ssa_solar(coli,lay,i) + &
             fac_3_ad(2)

        transmission_scatt_ir_ad%phup(coli,lay,i) = &
             transmission_scatt_ir_ad%phup(coli,lay,i) + &
             fac_3_ad(1) * solar_spectrum(i) * z4pi_r
      ENDDO
    ENDDO
  END SUBROUTINE solar_scattering_near_surf_ad

  SUBROUTINE solar_rayleigh_ad(raytracing, raytracing_ad, ircld, profiles, profiles_ad, profiles_dry, &
                               profiles_dry_ad, transmission_aux, transmission_aux_ad, auxrad_column_ad)

    TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
    TYPE(rttov_raytracing),       INTENT(INOUT) :: raytracing_ad
    TYPE(rttov_ircld),            INTENT(IN)    :: ircld
    TYPE(rttov_profile),          INTENT(IN)    :: profiles(:)
    TYPE(rttov_profile),          INTENT(INOUT) :: profiles_ad(:)
    TYPE(rttov_profile),          INTENT(IN)    :: profiles_dry(:)
    TYPE(rttov_profile),          INTENT(INOUT) :: profiles_dry_ad(:)
    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux), INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column_ad

    INTEGER(jpim) :: i, col, lay, lev
    REAL(jprb) :: wlm, ss_param
    REAL(jprb) :: v_h2o(nlayers), v_h2o_dry(nlayers), m(nlayers)
    REAL(jprb) :: v_h2o_surf, v_h2o_dry_surf, m_surf
    REAL(jprb) :: cossat, cossol, cosazi, cosscata_term1, cosscata_term2, cosscata
    REAL(jprb) :: ray_phase, solar_src, solar_src_updn
    REAL(jprb) :: v_h2o_ad(nlayers), v_h2o_dry_ad(nlayers), m_ad(nlayers)
    REAL(jprb) :: v_h2o_surf_ad, v_h2o_dry_surf_ad, m_surf_ad, plev_ad, dp_ad
    REAL(jprb) :: cosscata_term1_ad, cosscata_term2_ad, cosscata_ad
    REAL(jprb) :: ray_phase_ad, solar_src_ad, solar_src_updn_ad
    REAL(jprb) :: rayrad_up_ad(0:nlayers,0:maxncolumns), rayrad_dn_ad(0:nlayers,0:maxncolumns)

    DO i = 1, nchanprof
      wlm = 10000._jprb / coef%ff_cwn(chan)    ! Wavelength in microns
      IF (.NOT. solar(i) .OR. wlm > opts%rt_ir%rayleigh_max_wavelength) CYCLE

      ! Calculate layer-independent scattering parameter

      IF (coef%fmv_model_ver <= 9) THEN
        IF (wlm < ray_scs_wlm) THEN
          ss_param = ray_scs_a1x * wlm ** (ray_scs_b1 + ray_scs_c1 * wlm + ray_scs_d1/wlm)
        ELSE
          ss_param = ray_scs_a2x * wlm ** (ray_scs_b2 + ray_scs_c2 * wlm + ray_scs_d2/wlm)
        ENDIF
        ss_param = ss_param * 0.01_jprb ** 2_jpim * 1.E3_jprb * 100._jprb * na * z4pi_r / gravity

        v_h2o_dry = 0.5_jprb * (profiles_dry(prof)%q(1:nlevels-1) + profiles_dry(prof)%q(2:nlevels)) * 1.E-6_jprb
        v_h2o = v_h2o_dry / (1._jprb + v_h2o_dry)
        m = ((1._jprb  - v_h2o) * Mair + v_h2o * Mh2o)
      ELSE
        ss_param = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc * z4pi_r / gravity
        m = mair
        m_surf = mair
      ENDIF

      rayrad_up_ad(:,:) = 0._jprb
      rayrad_dn_ad(:,:) = 0._jprb

      ! Add the contribution from the part-layer above the surface

      IF (profiles(prof)%s2m%p >= opts%rt_ir%rayleigh_min_pressure) THEN

        IF (coef%fmv_model_ver <= 9) THEN
          IF (opts%rt_all%use_q2m) THEN
            v_h2o_dry_surf = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-6_jprb
          ELSE
            v_h2o_dry_surf = profiles_dry(prof)%q(iv3lev(i)) * 1.E-6_jprb
          ENDIF
          v_h2o_surf = v_h2o_dry_surf / (1._jprb + v_h2o_dry_surf)
          m_surf = ((1._jprb  - v_h2o_surf) * Mair + v_h2o_surf * Mh2o)
        ENDIF

        solar_src = solar_spectrum(i) * &
          ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * &
          raytracing%pathsat(iv2lay(i), prof) * ss_param / m_surf

        cossat = 1._jprb - raytracing%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof)
        cossol = 1._jprb - raytracing%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof)
        cosscata_term1 = SQRT(cossat * cossol)
        cosazi = COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
        cosscata_term2 = raytracing%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * cosazi

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1._jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        cosscata_term1_ad = 0._jprb
        cosscata_term2_ad = 0._jprb
        solar_src_updn_ad = 0._jprb

        DO col = 0, ircld%ncolumn(prof)
          IF (transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) > min_tau) THEN
            solar_src_updn_ad = solar_src_updn_ad + &
              auxrad_column_ad%meanrad_down_solar(col,i) * &
              transmission_aux%solar_path2%Tau_level(iv3lev(i),col,i) / &
              transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) ** 3_jpim

            transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),col,i) = &
              transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),col,i) + &
              auxrad_column_ad%meanrad_down_solar(col,i) * &
              solar_src_updn / transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) ** 3_jpim

            transmission_aux_ad%solar_path1%Tau_level(iv3lev(i),col,i) = &
              transmission_aux_ad%solar_path1%Tau_level(iv3lev(i),col,i) - &
              auxrad_column_ad%meanrad_down_solar(col,i) * &
              solar_src_updn * 3_jpim * transmission_aux%solar_path2%Tau_level(iv3lev(i),col,i) / &
              transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) ** 4_jpim
          ENDIF

          solar_src_updn_ad = solar_src_updn_ad + auxrad_column_ad%meanrad_up_solar(col,i) * &
            transmission_aux%solar_path2%Tau_level(iv3lev(i),col,i)

          transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),col,i) = &
            transmission_aux_ad%solar_path2%Tau_level(iv3lev(i),col,i) + &
            solar_src_updn * auxrad_column_ad%meanrad_up_solar(col,i)
        ENDDO

        ray_phase_ad = solar_src * solar_src_updn_ad  ! First use of ray_phase_ad, no accumulation
        solar_src_ad = solar_src_updn_ad * ray_phase  ! First use of solar_src_ad
        cosscata_ad = 2._jprb * 0.75_jprb * ray_phase_ad * cosscata ! First use of cossacata_ad
        cosscata_term1_ad = cosscata_term1_ad - cosscata_ad
        cosscata_term2_ad = cosscata_term2_ad - cosscata_ad

        raytracing_ad%zasat(iv2lay(i), prof) = raytracing_ad%zasat(iv2lay(i), prof) + &
          cosscata_term2_ad * raytracing%zasun(iv2lay(i), prof) * cosazi
        raytracing_ad%zasun(iv2lay(i), prof) = raytracing_ad%zasun(iv2lay(i), prof) + &
          cosscata_term2_ad * raytracing%zasat(iv2lay(i), prof) * cosazi

        raytracing_ad%zasat(iv2lay(i), prof) = raytracing_ad%zasat(iv2lay(i), prof) - &
          cosscata_term1_ad * raytracing%zasat(iv2lay(i), prof) * cossol / cosscata_term1
        raytracing_ad%zasun(iv2lay(i), prof) = raytracing_ad%zasun(iv2lay(i), prof) - &
          cosscata_term1_ad * raytracing%zasun(iv2lay(i), prof) * cossat / cosscata_term1

        raytracing_ad%pathsat(iv2lay(i), prof) = raytracing_ad%pathsat(iv2lay(i), prof) + &
            solar_src_ad * solar_spectrum(i) * (ss_param / m_surf) * &
            ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i)))

        dp_ad = solar_src_ad * solar_spectrum(i) * (ss_param / m_surf) * &
                raytracing%pathsat(iv2lay(i), prof)

        IF (profiles(prof)%s2m%p >= profiles(prof)%p(iv3lev(i))) THEN
          profiles_ad(prof)%s2m%p = profiles_ad(prof)%s2m%p + dp_ad
          plev_ad = -dp_ad
        ELSE
          profiles_ad(prof)%s2m%p = profiles_ad(prof)%s2m%p - dp_ad
          plev_ad = dp_ad
        ENDIF

        IF (opts%interpolation%lgradp) THEN
          profiles_ad(prof)%p(iv3lev(i)) = profiles_ad(prof)%p(iv3lev(i)) + plev_ad
        ENDIF

        IF (coef%fmv_model_ver <= 9) THEN
          m_surf_ad = -solar_src_ad * solar_src / m_surf
          v_h2o_surf_ad = m_surf_ad * (Mh2o - Mair)
          v_h2o_dry_surf_ad = v_h2o_surf_ad * (1._jprb - v_h2o_surf) / (1._jprb + v_h2o_dry_surf)

          IF (opts%rt_all%use_q2m) THEN
            profiles_dry_ad(prof)%q(iv3lev(i)) = profiles_dry_ad(prof)%q(iv3lev(i)) + &
              0.5_jprb * 1.E-6_jprb * v_h2o_dry_surf_ad
            profiles_dry_ad(prof)%s2m%q = profiles_dry_ad(prof)%s2m%q + &
              0.5_jprb * 1.E-6_jprb * v_h2o_dry_surf_ad
          ELSE
            profiles_dry_ad(prof)%q(iv3lev(i)) = profiles_dry_ad(prof)%q(iv3lev(i)) + &
              v_h2o_dry_surf_ad * 1.E-6_jprb
          ENDIF
        ENDIF
      ENDIF

      ! Atmospheric contribution
      DO col = 0, ircld%ncolumn(prof)
        rayrad_dn_ad(iv3lay(i),col) = rayrad_dn_ad(iv3lay(i),col) + auxrad_column_ad%meanrad_down_solar(col,i)
      
        rayrad_dn_ad(1:nlayers,col) = rayrad_dn_ad(1:nlayers,col) + auxrad_column_ad%down_solar(:,col,i)
      
        rayrad_up_ad(iv3lay(i),col) = rayrad_up_ad(iv3lay(i),col) + auxrad_column_ad%meanrad_up_solar(col,i)
      
        rayrad_up_ad(1:nlayers,col) = rayrad_up_ad(1:nlayers,col) + auxrad_column_ad%up_solar(:,col,i)
      ENDDO

      DO lev = nlevels, 2, -1
        ! Skip calculation for whole layers above min pressure
        IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) CYCLE

        lay = lev - 1

        solar_src = solar_spectrum(i) * &
              (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * &
              raytracing%pathsat(lay, prof) * ss_param / m(lev-1)

        cossat = 1._jprb - raytracing%zasat(lay, prof) * raytracing%zasat(lay, prof)
        cossol = 1._jprb - raytracing%zasun(lay, prof) * raytracing%zasun(lay, prof)
        cosscata_term1 = SQRT(cossat * cossol)
        cosazi = COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)
        cosscata_term2 = raytracing%zasat(lay, prof) * raytracing%zasun(lay, prof) * cosazi

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1._jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        cosscata_term1_ad = 0._jprb
        cosscata_term2_ad = 0._jprb
        solar_src_updn_ad = 0._jprb

        DO col = 0, ircld%ncolumn(prof)
          IF (transmission_aux%solar_path1%Tau_level(lev-1,col,i) > min_tau) THEN
            solar_src_updn_ad = solar_src_updn_ad + rayrad_dn_ad(lay,col) * &
              transmission_aux%solar_path2%Tau_level(lev-1,col,i) / &
              transmission_aux%solar_path1%Tau_level(lev-1,col,i) ** 3_jpim

            transmission_aux_ad%solar_path2%Tau_level(lev-1,col,i) = &
              transmission_aux_ad%solar_path2%Tau_level(lev-1,col,i) + &
              rayrad_dn_ad(lay,col) * solar_src_updn / &
              transmission_aux%solar_path1%Tau_level(lev-1,col,i) ** 3_jpim

            transmission_aux_ad%solar_path1%Tau_level(lev-1,col,i) = &
              transmission_aux_ad%solar_path1%Tau_level(lev-1,col,i) - &
              rayrad_dn_ad(lay,col) * solar_src_updn * &
              3_jpim * transmission_aux%solar_path2%Tau_level(lev-1,col,i) / &
              transmission_aux%solar_path1%Tau_level(lev-1,col,i) ** 4_jpim
          ENDIF

          rayrad_dn_ad(lay-1,col) = rayrad_dn_ad(lay-1,col) + rayrad_dn_ad(lay,col)

          solar_src_updn_ad = solar_src_updn_ad + rayrad_up_ad(lay,col) * &
            transmission_aux%solar_path2%Tau_level(lev-1,col,i)

          transmission_aux_ad%solar_path2%Tau_level(lev-1,col,i) = &
            transmission_aux_ad%solar_path2%Tau_level(lev-1,col,i) + &
            solar_src_updn * rayrad_up_ad(lay,col)

          rayrad_up_ad(lay-1,col) = rayrad_up_ad(lay-1,col) + rayrad_up_ad(lay,col)
        ENDDO

        ray_phase_ad = solar_src * solar_src_updn_ad  ! First use of ray_phase_ad, no accumulation
        solar_src_ad = solar_src_updn_ad * ray_phase  ! First use of solar_src_ad
        cosscata_ad = 2._jprb * 0.75_jprb * ray_phase_ad * cosscata ! First use of cossacata_ad
        cosscata_term1_ad = cosscata_term1_ad - cosscata_ad
        cosscata_term2_ad = cosscata_term2_ad - cosscata_ad

        raytracing_ad%zasat(lay, prof) = raytracing_ad%zasat(lay, prof) + &
          cosscata_term2_ad * raytracing%zasun(lay, prof) * cosazi
        raytracing_ad%zasun(lay, prof) = raytracing_ad%zasun(lay, prof) + &
          cosscata_term2_ad * raytracing%zasat(lay, prof) * cosazi

        raytracing_ad%zasat(lay, prof) = raytracing_ad%zasat(lay, prof) - &
          cosscata_term1_ad * raytracing%zasat(lay, prof) * cossol / cosscata_term1
        raytracing_ad%zasun(lay, prof) = raytracing_ad%zasun(lay, prof) - &
          cosscata_term1_ad * raytracing%zasun(lay, prof) * cossat / cosscata_term1

        IF (opts%interpolation%lgradp) THEN
          profiles_ad(prof)%p(lev) = profiles_ad(prof)%p(lev) + &
              solar_src_ad * solar_spectrum(i) * (ss_param / m(lev-1)) * &
              raytracing%pathsat(lay, prof)
          profiles_ad(prof)%p(lev-1) = profiles_ad(prof)%p(lev-1) - &
              solar_src_ad * solar_spectrum(i) * (ss_param / m(lev-1)) * &
              raytracing%pathsat(lay, prof)
        ENDIF

        raytracing_ad%pathsat(lay, prof) = raytracing_ad%pathsat(lay, prof) + &
            solar_src_ad * solar_spectrum(i) * (ss_param / m(lev-1)) * &
            (profiles(prof)%p(lev) - profiles(prof)%p(lev-1))

        m_ad(lev-1) = -solar_src_ad * solar_src / m(lev-1)
      ENDDO

      IF (coef%fmv_model_ver <= 9) THEN
        v_h2o_ad = m_ad * (Mh2o - Mair)
        v_h2o_dry_ad = v_h2o_ad * (1._jprb - v_h2o) / (1._jprb + v_h2o_dry)

        profiles_dry_ad(prof)%q(1:nlevels-1) = profiles_dry_ad(prof)%q(1:nlevels-1) + &
          0.5_jprb * v_h2o_dry_ad * 1.E-6_jprb
        profiles_dry_ad(prof)%q(2:nlevels) = profiles_dry_ad(prof)%q(2:nlevels) + &
          0.5_jprb * v_h2o_dry_ad * 1.E-6_jprb
      ENDIF
    ENDDO
  END SUBROUTINE solar_rayleigh_ad

  SUBROUTINE solar_surface_contribution_ad(transmission_aux, transmission_aux_ad, ircld, &
                                           reflectance, reflectance_ad, auxrad_column_ad)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_transmission_aux), INTENT(INOUT) :: transmission_aux_ad
    TYPE(rttov_ircld),            INTENT(IN)    :: ircld
    REAL(jprb),                   INTENT(IN)    :: reflectance(:)
    REAL(jprb),                   INTENT(INOUT) :: reflectance_ad(:)
    TYPE(rttov_radiance_aux),     INTENT(IN)    :: auxrad_column_ad

    INTEGER(jpim) :: i, col

    DO i = 1, nchanprof
      IF (solar(i)) THEN
        DO col = 0, ircld%ncolumn(prof)
          reflectance_ad(i) = reflectance_ad(i) + solar_spectrum(i) * refl_norm(i) * &
            auxrad_column_ad%cloudy(col,i) * transmission_aux%solar_path2%Tau_surf(col,i)

          transmission_aux_ad%solar_path2%Tau_surf(col,i) = transmission_aux_ad%solar_path2%Tau_surf(col,i) + &
            solar_spectrum(i) * refl_norm(i) * auxrad_column_ad%cloudy(col,i) * reflectance(i)
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE solar_surface_contribution_ad

END SUBROUTINE rttov_integrate_ad
