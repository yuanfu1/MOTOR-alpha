! Description:
!> @file
!!   Integrate the radiative transfer equation
!
!> @brief
!!   Integrate the radiative transfer equation
!!
!! @details
!!   The layer emission terms are calculated using the linear-in-tau
!!   approximation, although there is an option to use the simple layer-average
!!   source term for the downwelling emission which is slightly more efficient.
!!
!!   For v9 predictor coefficients, Rayleigh extinction is included in the
!!   gas optical depths. For v13 predictor coefficients, the Rayleigh
!!   extinction is computed separately within RTTOV. In either case a Rayleigh
!!   single-scattering contribution is computed here if requested. The Rayleigh
!!   single-scattering is calculated even when the DOM solver is used, unless
!!   the user sets the dom_rayleigh option. Rayleigh scattering is computed
!!   using the parameterisation in Bucholtz (1995). For v9 predictors, the
!!   parameterised scattering cross-section is used, and the calculation takes
!!   account of variable water vapour. In practice the impact of this is
!!   negligible. For v13 predictors, the parameterised extinction coefficient
!!   is used (as in rttov_rayleigh_extinction.F90), and water vapour is ignored.
!!
!!   Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
!!   systems. ECMWF Research Dept. Tech. Memo. 176
!!
!!   Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
!!   Transfer Model for Assimilation of Satellite Radiance Observations.
!!   QJRMS, 125, 1407-1425.
!!
!!   Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
!!   transfer model for the infrared atmospheric sounding interferometer.
!!   ECMWF Research Dept. Tech. Memo. 425
!!
!!   Matricardi, M. 2005 The inclusion of aerosols and clouds in RTIASI,the
!!   ECMWF radiative transfer model for the infrared atmospheric sounding
!!   interferometer. ECMWF Research Dept. Tech. Memo. 474
!!
!!   Bucholtz, A., 1995: Rayleigh-scattering calculations for the terrestrial
!!   atmosphere. Applied Optics, 34, 15, 2765-2773.
!!
!!
!! @param[in]     addcosmic                   flag for inclusion of CMBR term in MW simulations
!! @param[in]     opts                        RTTOV options structure
!! @param[in]     maxncolumns                 largest number of cloud columns across all profiles
!! @param[in]     chanprof                    specifies channels and profiles to simulate
!! @param[in]     emissivity                  input/output surface emissivities
!! @param[in]     reflectance                 input/output surface reflectances for direct solar beam
!! @param[in]     refl_norm                   surface relfectance normalisation factors
!! @param[in]     diffuse_refl                surface reflectance for downwelling radiation
!! @param[in]     do_lambertian               flag indicating whether Lambertian surface is active for each channel
!! @param[in]     do_mfasis                   flag to indicate MFASIS simulation
!! @param[in]     thermal                     per-channel flag to indicate if emissive simulations are being performed
!! @param[in]     dothermal                   flag to indicate if any emissive simulations are being performed
!! @param[in]     solar                       per-channel flag to indicate if solar simulations are being performed
!! @param[in]     dosolar                     flag to indicate if any solar simulations are being performed
!! @param[in]     do_rayleigh_ss              flag to indicate if Rayleigh single-scattering should be included
!! @param[in]     solar_spectrum              TOA solar irradiance for each channel
!! @param[in]     transmission_aux            top-level auxiliary transmission structure
!! @param[in]     transmission_scatt_ir       visible/IR cloud/aerosol scattering parameters
!! @param[in]     profiles                    input atmospheric profiles and surface variables
!! @param[in]     profiles_dry                profiles in internal units
!! @param[in]     aux_prof                    auxiliary profile variables
!! @param[in]     coef                        optical depth coefficients structure
!! @param[in]     raytracing                  raytracing structure
!! @param[in]     ircld                       computed cloud column data
!! @param[in,out] rad                         primary output radiance structure
!! @param[in,out] rad2                        secondary output radiance structure
!! @param[in,out] auxrad                      Planck radiances
!! @param[in,out] auxrad_column               internal radiance structure
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
SUBROUTINE rttov_integrate(addcosmic, opts, maxncolumns, chanprof, &
                           emissivity, reflectance, refl_norm, diffuse_refl, do_lambertian, &
                           do_mfasis, thermal, dothermal, solar, dosolar, do_rayleigh_ss, solar_spectrum, &
                           transmission_aux, transmission_scatt_ir, &
                           profiles, profiles_dry, aux_prof, coef, raytracing, ircld, &
                           rad, rad2, auxrad, auxrad_column)

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, rttov_options, rttov_profile, rttov_profile_aux, &
                          rttov_transmission_aux, rttov_transmission_scatt_ir, rttov_radiance, rttov_radiance2, &
                          rttov_ircld, rttov_raytracing, rttov_radiance_aux, rttov_emissivity, rttov_reflectance
!INTF_OFF
  USE rttov_const, ONLY : realtol, sensor_id_po, min_od, min_tau, pi_r, z4pi_r, &
                          deg2rad, gravity, na, rgc, Mh2o, Mair, &
                          overcast_albedo_wvn, overcast_albedo1, overcast_albedo2, &
                          vis_scatt_dom, vis_scatt_single, ir_scatt_dom, &
                          ray_ps, ray_ts, ray_scs_wlm, &
                          ray_scs_a1x, ray_scs_b1, ray_scs_c1, ray_scs_d1, &
                          ray_scs_a2x, ray_scs_b2, ray_scs_c2, ray_scs_d2

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  LOGICAL(jplm),                      INTENT(IN)              :: addcosmic
  TYPE(rttov_options),                INTENT(IN)              :: opts
  INTEGER(jpim),                      INTENT(IN)              :: maxncolumns
  TYPE(rttov_chanprof),               INTENT(IN)              :: chanprof(:)
  TYPE(rttov_emissivity),             INTENT(IN),    OPTIONAL :: emissivity(:)
  TYPE(rttov_reflectance),            INTENT(IN),    OPTIONAL :: reflectance(:)
  REAL(jprb),                         INTENT(IN)              :: refl_norm(:)
  REAL(jprb),                         INTENT(IN)              :: diffuse_refl(:)
  LOGICAL(jplm),                      INTENT(IN)              :: do_lambertian(:)
  LOGICAL(jplm),                      INTENT(IN)              :: do_mfasis
  LOGICAL(jplm),                      INTENT(IN)              :: thermal(:)
  LOGICAL(jplm),                      INTENT(IN)              :: dothermal
  LOGICAL(jplm),                      INTENT(IN)              :: solar(:)
  LOGICAL(jplm),                      INTENT(IN)              :: dosolar
  LOGICAL(jplm),                      INTENT(IN)              :: do_rayleigh_ss
  REAL(jprb),                         INTENT(IN)              :: solar_spectrum(:)
  TYPE(rttov_transmission_aux),       INTENT(IN)              :: transmission_aux
  TYPE(rttov_transmission_scatt_ir),  INTENT(IN)              :: transmission_scatt_ir
  TYPE(rttov_profile),                INTENT(IN)              :: profiles(:)
  TYPE(rttov_profile),                INTENT(IN)              :: profiles_dry(:)
  TYPE(rttov_profile_aux),            INTENT(IN)              :: aux_prof
  TYPE(rttov_coef),                   INTENT(IN)              :: coef
  TYPE(rttov_raytracing),             INTENT(IN)              :: raytracing
  TYPE(rttov_ircld),                  INTENT(IN)              :: ircld
  TYPE(rttov_radiance),               INTENT(INOUT)           :: rad
  TYPE(rttov_radiance2),              INTENT(INOUT), OPTIONAL :: rad2
  TYPE(rttov_radiance_aux),           INTENT(INOUT)           :: auxrad
  TYPE(rttov_radiance_aux),           INTENT(INOUT)           :: auxrad_column
!INTF_END

#include "rttov_calcrad.interface"

  INTEGER(jpim) :: i, lev, col, lay, nchanprof, nlayers, nlevels
  REAL(jprb)    :: refl, refl_norm_scat
  LOGICAL(jplm) :: keyradonly, do_scatt, dom_ir, dom_vis, do_single_scatt

  INTEGER(jpim) :: iv2lev(SIZE(chanprof)), iv2lay(SIZE(chanprof))
  INTEGER(jpim) :: iv3lev(SIZE(chanprof)), iv3lay(SIZE(chanprof))
  INTEGER(jpim) :: pol_id(SIZE(chanprof))

  REAL(jprb)    :: cfraction(SIZE(chanprof)), pfraction(SIZE(chanprof))
  LOGICAL(jplm) :: sateqsun(profiles(1)%nlayers,SIZE(profiles(:)))

  REAL(jprb) :: ZHOOK_HANDLE
!- End of header ------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',0_jpim,ZHOOK_HANDLE)

! Define macros for commonly used variables
#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

#define tau_surf_p_r transmission_aux%thermal_path1%Tau_surf_p_r(col,i)
#define tau_surf_p transmission_aux%thermal_path1%Tau_surf_p(col,i)
#define tau_surf_r transmission_aux%thermal_path1%Tau_surf_r(col,i)
#define tau_surf transmission_aux%thermal_path1%Tau_surf(col,i)
#define tau_layer_p_r transmission_aux%thermal_path1%Tau_level_p_r(lay,col,i)
#define tau_layer_p transmission_aux%thermal_path1%Tau_level_p(lay,col,i)
#define tau_layer_r transmission_aux%thermal_path1%Tau_level_r(lay,col,i)
#define tau_layer transmission_aux%thermal_path1%Tau_level(lay,col,i)
#define tau_level_p_r transmission_aux%thermal_path1%Tau_level_p_r(lay+1,col,i)
#define tau_level_p transmission_aux%thermal_path1%Tau_level_p(lay+1,col,i)
#define tau_level_r transmission_aux%thermal_path1%Tau_level_r(lay+1,col,i)
#define tau_level transmission_aux%thermal_path1%Tau_level(lay+1,col,i)
#define od_singlelayer_r transmission_aux%thermal_path1%Od_singlelayer_r(coli,lay,i)
#define od_singlelayer transmission_aux%thermal_path1%Od_singlelayer(coli,lay,i)
#define fac1 transmission_aux%Fac1(lay,col,i)
#define fac2_thermal_path1 transmission_aux%thermal_path1%fac2(lay+1,col,i)
#define fac2_solar_path1 transmission_aux%solar_path1%fac2(lay+1,col,i)
#define surf_fac transmission_aux%Surf_fac(col,i)

  !-------------------------------------------------------------------------------
  ! Initialise useful variables
  !-------------------------------------------------------------------------------

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  nlevels = nlayers + 1

  do_scatt = opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl
  dom_ir = do_scatt .AND. opts%rt_ir%ir_scatt_model == ir_scatt_dom
  dom_vis = do_scatt .AND. opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom ! DOM vis *selected*
  do_single_scatt = do_scatt .AND. dosolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_single  ! Single-scatt *actually used*
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

  auxrad_column%cloudy = 0._jprb
  IF (dosolar .AND. .NOT. do_mfasis) THEN
    auxrad_column%up_solar = 0._jprb
    auxrad_column%meanrad_up_solar = 0._jprb
    auxrad_column%down_solar = 0._jprb
    auxrad_column%meanrad_down_solar = 0._jprb
  ENDIF

  !-------------------------------------------------------------------------------
  ! Calculate layer radiances
  !-------------------------------------------------------------------------------
  IF (dothermal) CALL rttov_calcrad(addcosmic, chanprof, profiles, coef, thermal, auxrad)

  !-------------------------------------------------------------------------------
  ! Calculate atmospheric contribution from layers
  !-------------------------------------------------------------------------------
  IF (dothermal .AND. .NOT. dom_ir) &
    CALL calc_atmospheric_radiance(transmission_aux, auxrad, auxrad_column)

  ! Scattering of the solar beam
  IF (do_single_scatt) &
    CALL  solar_scattering_air(transmission_aux, diffuse_refl, raytracing, &
                               transmission_scatt_ir, auxrad_column)

  !-------------------------------------------------------------------------------
  ! Calculate near-surface layer contribution
  !-------------------------------------------------------------------------------
  IF (dothermal .AND. .NOT. dom_ir) &
    CALL calc_near_surf_contribution(transmission_aux, auxrad, auxrad_column)

  ! Scattering of the solar beam
  IF (do_single_scatt) &
    CALL solar_scattering_near_surf(transmission_aux, diffuse_refl, auxrad_column)

  !-------------------------------------------------------------------------------
  ! Calculate clear-sky Rayleigh scattering contribution
  !-------------------------------------------------------------------------------
  IF (do_rayleigh_ss .AND. .NOT. do_mfasis) &
    CALL solar_rayleigh(raytracing, profiles, profiles_dry, transmission_aux, auxrad_column)

  !-------------------------------------------------------------------------------
  ! Cosmic temperature correction - direct model only since tcosmic is fixed
  !-------------------------------------------------------------------------------
  ! Add Planck source corresponding to tcosmic = 2.7k for microwave sensors only
  IF (addcosmic) THEN
    col = 0
    DO i = 1, nchanprof
      auxrad_column%meanrad_down(col,i) = auxrad_column%meanrad_down(col,i) + auxrad%cosmic(i)
      IF (do_lambertian(i)) THEN
        auxrad_column%meanrad_down_p(col,i) = auxrad_column%meanrad_down_p(col,i) + auxrad%cosmic(i)
      ENDIF
    ENDDO
  ENDIF

  !-------------------------------------------------------------------------------
  ! Add thermal and solar atmospheric contributions to the clear and cloudy columns
  !-------------------------------------------------------------------------------
  DO i = 1, nchanprof
    IF (thermal(i) .AND. .NOT. dom_ir) THEN
      DO col = 0, ircld%ncolumn(prof)
        IF (do_lambertian(i)) THEN
          auxrad_column%cloudy(col,i) = &
              auxrad_column%cloudy(col,i) + auxrad_column%meanrad_up(col,i) + &
              (emissivity(i)%specularity * &
               auxrad_column%meanrad_down(col,i) * tau_surf + &
               (1._jprb - emissivity(i)%specularity) * &
               auxrad_column%meanrad_down_p(col,i) * tau_surf_p) * &
              diffuse_refl(i) * tau_surf
        ELSE
          auxrad_column%cloudy(col,i) = &
              auxrad_column%cloudy(col,i) + auxrad_column%meanrad_up(col,i) + &
              auxrad_column%meanrad_down(col,i) * diffuse_refl(i) * &
              tau_surf**2_jpim
        ENDIF
      ENDDO

      ! Replace the upward radiances from the level at the bottom of the layer
      ! containing the surface with the calculated surface->ToA radiance (only for column 0)
      ! (for overcast and secondary radiances)
      auxrad_column%up(iv2lay(i),0,i) = auxrad_column%meanrad_up(0,i)
    ENDIF

    IF (solar(i) .AND. .NOT. do_mfasis) THEN
      ! Downward-scattered component: this radiation is travelling along the
      ! satellite line-of-sight. diffuse_refl must be divided by pi to give
      ! a BRDF and multiplied by cos(sat_zen_angle)
      refl_norm_scat = COS(profiles(prof)%zenangle * deg2rad) * pi_r

      DO col = 0, ircld%ncolumn(prof)
        auxrad_column%cloudy(col,i) = auxrad_column%cloudy(col,i) + &
                                      auxrad_column%meanrad_up_solar(col,i) + &
                                      auxrad_column%meanrad_down_solar(col,i) * &
                                      diffuse_refl(i) * refl_norm_scat * &
                                      transmission_aux%solar_path1%Tau_surf(col,i)**2_jpim
      ENDDO

      ! Replace the upward radiances from the level at the bottom of the layer
      ! containing the surface with the calculated surface->ToA radiance (only for column 0)
      ! (for overcast and secondary radiances)
      auxrad_column%up_solar(iv2lay(i),0,i) = auxrad_column%meanrad_up_solar(0,i)
    ENDIF
  ENDDO

  !-------------------------------------------------------------------------------
  ! Calculate surface emission contribution
  !-------------------------------------------------------------------------------
  IF (dothermal .AND. .NOT. dom_ir) THEN
  !cdir nodep
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
  !cdir nodep
        DO col = 0, ircld%ncolumn(prof)
          auxrad_column%cloudy(col,i) = auxrad_column%cloudy(col,i) + &
                                        auxrad%skin(i) * emissivity(i)%emis_out * tau_surf
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  !-------------------------------------------------------------------------------
  ! Solar surface contribution
  !-------------------------------------------------------------------------------
  IF (dosolar .AND. .NOT. (dom_vis .OR. do_mfasis)) &
    CALL solar_surface_contribution(transmission_aux, reflectance%refl_out, auxrad_column)

  !-------------------------------------------------------------------------------
  ! Calculate secondary radiances
  !-------------------------------------------------------------------------------
  ! These direct-model-only outputs include thermal contributions (no solar) and are calculated
  ! only for non-aerosol, non-DOM and non-PC simulations if the rad2 parameter is present

  IF (PRESENT(rad2) .AND. dothermal .AND. .NOT. keyradonly) then
    !cdir nodep
    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        col = 0

        ! Clear-sky upwelling atmospheric emission at TOA
        ! Note the surface layer value in auxrad_column%up has been modified above
        rad2%up(:,i) = auxrad_column%up(:,col,i)

        ! Clear-sky upwelling radiance at TOA (without surface reflected term)
        rad2%upclear(i) = auxrad_column%meanrad_up(col,i) + &
                          auxrad%skin(i) * emissivity(i)%emis_out * tau_surf

        ! Clear-sky downwelling radiance at surface (before reflection)
        IF (do_lambertian(i)) THEN  ! Lambertian/specular mixed
          rad2%dnclear(i) = emissivity(i)%specularity * &
                            auxrad_column%meanrad_down(col,i) * tau_surf + &
                            (1._jprb - emissivity(i)%specularity) * &
                            auxrad_column%meanrad_down_p(col,i) * tau_surf_p
        ELSE                        ! Specular
          rad2%dnclear(i) = auxrad_column%meanrad_down(col,i) * tau_surf
        ENDIF

        ! Reflected clear-sky downwelling radiance at TOA
        rad2%refldnclear(i) = rad2%dnclear(i) * diffuse_refl(i) * tau_surf

        ! Planck emission
        rad2%surf(2:nlayers,i) = auxrad%air(3:nlevels,i)
        rad2%surf(iv2lay(i),i) = auxrad%skin(i)

        ! Clear-sky downwelling radiance from TOA down to bottom of each layer
        ! at the level bounding the bottom of the layer. Note that tcosmic is
        ! included here (it is set to zero if not applicable).
        DO lay = 2, nlayers
          IF (tau_level > min_tau) THEN
            IF (do_lambertian(i)) THEN
              rad2%down(lay,i) = emissivity(i)%specularity * &
                                 (auxrad_column%down(lay,col,i) + auxrad%cosmic(i)) * &
                                 tau_level + &
                                 (1._jprb - emissivity(i)%specularity) * &
                                 (auxrad_column%down_p(lay,col,i) + auxrad%cosmic(i)) * &
                                 tau_level_p
            ELSE
              rad2%down(lay,i) = (auxrad_column%down(lay,col,i) + auxrad%cosmic(i)) * &
                                 tau_level
            ENDIF
          ELSE
            rad2%down(lay,i) = rad2%down(lay-1,i)
          ENDIF
        ENDDO
        IF (tau_surf > min_tau) THEN
          rad2%down(iv2lay(i),i) = rad2%dnclear(i)
        ELSE
          rad2%down(iv2lay(i),i) = rad2%down(iv3lay(i),i)
        ENDIF

      ENDIF
    ENDDO
  ENDIF

  !-------------------------------------------------------------------------------
  ! Calculate overcast radiances
  !-------------------------------------------------------------------------------
  ! Overcast radiances only calculated for non-aerosol, non-DOM, non-PC simulations

  IF (.NOT. keyradonly) THEN
    ! rad%overcast includes solar contribution ONLY for pure-solar channels because assumed
    !   cloud emissivity is 1.0 (i.e. no reflection) for all thermal channels
    col = 0
    DO i = 1, nchanprof
      IF (thermal(i) .AND. .NOT. dom_ir) THEN
        DO lay = 1, nlayers
          lev = lay + 1
          ! Overcast radiances at given cloud top
          rad%overcast(lay,i) = auxrad_column%up(lay,col,i) + auxrad%air(lev,i) * tau_level
        ENDDO
      ELSEIF (solar(i) .AND. .NOT. do_mfasis .AND. .NOT. dom_vis) THEN
        ! Very crude model: assumes clouds are Lambertian reflectors with fixed albedo
        ! Use input cloud top BRDF if user has supplied it, otherwise use default BRDF
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
          ! Overcast radiances at given cloud top
          rad%overcast(lay,i) = solar_spectrum(i) * refl / raytracing%pathsun(lay,prof) * &
                                transmission_aux%solar_path2%Tau_level(lev,col,i) + &
                                auxrad_column%up_solar(lay,col,i) + &
                                auxrad_column%down_solar(lay,col,i) * refl / raytracing%pathsat(lay,prof) * &
                                transmission_aux%solar_path1%Tau_level(lev,col,i) ** 2_jpim
        ENDDO
      ENDIF
    ENDDO

    IF (dothermal) THEN
      ! Add surface component to thermal overcast radiances in near-surface layer
      DO i = 1, nchanprof
        IF (thermal(i) .AND. .NOT. dom_ir) THEN
          lay = iv2lay(i)
          rad%overcast(lay,i) = auxrad_column%up(lay,col,i) + tau_surf * auxrad%surfair(i)
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  !-------------------------------------------------------------------------------
  ! Calculate total radiance
  !-------------------------------------------------------------------------------
  WHERE (.NOT. (solar(1:nchanprof) .AND. do_mfasis)) &
    rad%clear(1:nchanprof) = auxrad_column%cloudy(0,1:nchanprof)

  ! The simple cloudy scheme is not applied to aerosol-affected radiances
  IF (do_scatt) THEN
    !---------------------------------------------------
    ! Calculate complex cloudy radiances
    !---------------------------------------------------
    DO i = 1, nchanprof
      ! Skip MFASIS channels
      IF (solar(i) .AND. do_mfasis) CYCLE
      ! For thermal channels with DOM the whole radiance is calculated in rttov_dom.
      ! For solar-affected channels the Rayleigh contribution and/or the solar
      ! single-scattering are calculated here and need to be accumulated.
      IF (thermal(i) .AND. dom_ir .AND. .NOT. solar(i)) CYCLE
      DO col = 1, ircld%ncolumn(prof)
        rad%cloudy(i) = rad%cloudy(i) + &
                        auxrad_column%cloudy(col,i) * &
                        (ircld%xcol(col+1,prof) - ircld%xcol(col,prof))
      ENDDO

      rad%cloudy(i) = rad%cloudy(i) + auxrad_column%cloudy(0,i) * ircld%xcolclr(prof)

      rad%total(i) = rad%cloudy(i)
    ENDDO

  ELSE
    !---------------------------------------------------
    ! Calculate total radiance (clear case/simple cloud)
    !---------------------------------------------------
    IF (opts%rt_ir%pc%addpc) THEN
      rad%total(1:nchanprof) = rad%clear(1:nchanprof)
    ELSE
      ! Interpolate to given cloud-top pressures
      DO i = 1, nchanprof
        ! Skip MFASIS channels
        IF (solar(i) .AND. do_mfasis) CYCLE

        lay = aux_prof%s(prof)%nearestlev_ctp - 1
        rad%cloudy(i) = rad%overcast(lay,i) * &
                        (1._jprb - aux_prof%s(prof)%pfraction_ctp) + &
                        rad%overcast(lay-1,i) * aux_prof%s(prof)%pfraction_ctp

        rad%total(i) = rad%clear(i) + cfraction(i) * (rad%cloudy(i) - rad%clear(i))
      ENDDO
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',1_jpim,ZHOOK_HANDLE)

CONTAINS

! DAR: This subroutine calculates the individual layer clear-sky radiances and then does a cumulative sum to determine the
!      radiance observed from the top of atmosphere to a particular layer.
! DAR: As expected, this routine consumes most of the time (loops over channels, columns and levels)
!      and has been most heavily optimised (for IBM only so far, Intel shows neutral impact - will look into this)
!      I have taken out the code that switches the order of the loops on the NEC and will let MF test this impact
!      See subroutine comments for more details of individual changes.

  SUBROUTINE calc_atmospheric_radiance(transmission_aux, auxrad, auxrad_column)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_radiance_aux),     INTENT(IN)    :: auxrad
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column

    INTEGER(jpim) :: i, col, coli, lay
    REAL(jprb)    :: up_laym1, down_laym1, down_p_laym1
    REAL(jprb)    :: rad_air_avg(nlayers), rad_air_diff(nlayers)
    REAL(jprb)    :: dtau(nlayers)

! DAR: fac and fac2 contain real 1 or 0 depending on whether a calculation should be performed or not.
!      IBM performance was suffering as a result of doing lots of branching (mispredicts?) so these arrays are populated in advance
!      and the calculation is performed regardless.

! DAR: I've removed a lot of the 'temporary' variables that were eventually summed together because they were taking up a
!      signficant amount of time on the IBM. Now the big sum is done with all the variables in full because the 16 prefetch streams
!      can handle this (POWER6) - maybe will have to split this up for POWER7 (only 12). Will do more testing on Intel with vtune
!      to check performance impact.

! DAR: I've also removed od_singlelayer_r and added it to transmission_aux so it can be reused the TL/AD/K code when ready
!      and doesn't have to be recalculated - This should be faster...

! DAR: cumulative sum is done at same time as it's significantly quicker than doing it later.
!      This code is still very hard to read though.
    DO i = 1, nchanprof
      IF (thermal(i)) THEN

        IF (do_lambertian(i) .OR. .NOT. opts%rt_all%rad_down_lin_tau) THEN
          rad_air_avg = 0.5_jprb * (auxrad%air(1:nlayers,i) + auxrad%air(2:nlevels,i))
        ENDIF
        rad_air_diff = auxrad%air(2:nlevels,i) - auxrad%air(1:nlayers,i)

        DO col = 0, ircld%ncolumn(prof)
          up_laym1 = 0._jprb
          down_laym1 = 0._jprb
          down_p_laym1 = 0._jprb

          dtau = transmission_aux%thermal_path1%Tau_level(1:nlayers,col,i) - &
                 transmission_aux%thermal_path1%Tau_level(2:nlevels,col,i)

          DO lay = 1, nlayers
            coli = ircld%icldarr(col,lay,prof)
            auxrad_column%up(lay,col,i) = up_laym1 + &
              fac1 * &
              (dtau(lay) * &
              (auxrad%air(lay,i) + &
              rad_air_diff(lay) * &
              od_singlelayer_r) - &
              rad_air_diff(lay) * &
              tau_level)
            IF (do_lambertian(i)) THEN
              ! Lambertian reflected downwelling layer-average
              auxrad_column%down_p(lay,col,i) = down_p_laym1 + &
                fac1 * &
                fac2_thermal_path1 * rad_air_avg(lay) * &
                ((tau_level_p_r) - (tau_layer_p_r))
            ENDIF
            IF (opts%rt_all%rad_down_lin_tau) THEN
              ! Specular reflected downwelling linear-in-tau
              auxrad_column%down(lay,col,i) = down_laym1 + &
                fac1 * &
                fac2_thermal_path1 * &
              ((dtau(lay) * &
                (auxrad%air(lay,i) - &
                rad_air_diff(lay) * &
                od_singlelayer_r) * &
                (tau_level_r * &
                tau_layer_r)) + &
                rad_air_diff(lay) * &
                tau_level_r)
            ELSE
              ! Specular reflected downwelling layer-average
              auxrad_column%down(lay,col,i) = down_laym1 + &
                fac1 * &
                fac2_thermal_path1 * rad_air_avg(lay) * &
                ((tau_level_r) - (tau_layer_r))
            ENDIF

            up_laym1   = auxrad_column%up(lay,col,i)
            down_laym1 = auxrad_column%down(lay,col,i)
            IF (do_lambertian(i)) down_p_laym1 = auxrad_column%down_p(lay,col,i)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    DO i = 1, nchanprof
      IF (thermal(i)) THEN
        DO col = 0, ircld%ncolumn(prof)
          auxrad_column%down_ref(:,col,i) = auxrad_column%down(:,col,i)
          auxrad_column%down(:,col,i) = MAX(auxrad_column%down_ref(:,col,i), 0._jprb)
        ENDDO
        IF (do_lambertian(i)) THEN
          DO col = 0, ircld%ncolumn(prof)
            auxrad_column%down_p_ref(:,col,i) = auxrad_column%down_p(:,col,i)
            auxrad_column%down_p(:,col,i) = MAX(auxrad_column%down_p_ref(:,col,i), 0._jprb)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  END SUBROUTINE calc_atmospheric_radiance

  ! DAR: This is the next biggest consumer of CPU time and should be looked at next. The trouble seems to be that you have to change a
  !      lot of non-consecutive data.

  SUBROUTINE calc_near_surf_contribution(transmission_aux, auxrad, auxrad_column)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_radiance_aux),     INTENT(IN)    :: auxrad
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column

    INTEGER(jpim) :: i, col, lay, lev
    REAL(jprb)    :: rad_air_avg, rad_air_diff

#define B1_3 auxrad%air(lev,i) * (tau_level  - tau_surf)
#define B2_3 rad_air_diff * tau_surf
#define B3_3 rad_air_diff * (tau_level - tau_surf) * (transmission_aux%thermal_path1%od_sfrac_r(col,i))

    DO i = 1, nchanprof
      IF (.NOT. thermal(i)) CYCLE
      lay = iv3lay(i)
      lev = iv3lev(i)

      rad_air_avg = 0.5_jprb * (auxrad%surfair(i) + auxrad%air(lev,i))
      rad_air_diff = auxrad%surfair(i) - auxrad%air(lev,i)

      DO col = 0, ircld%ncolumn(prof)

        IF (transmission_aux%thermal_path1%od_sfrac(col,i) < min_od .OR. &
            (opts%rt_all%dtau_test .AND. &
             (tau_level - tau_surf) < min_od)) THEN
          ! small optical depth or optical depth change set radiance to zero
          auxrad_column%meanrad_up(col,i) = 0._jprb
          auxrad_column%meanrad_down(col,i) = 0._jprb
          IF (do_lambertian(i)) auxrad_column%meanrad_down_p(col,i) = 0._jprb
        ELSE
          ! Upwelling linear-in-tau
          auxrad_column%meanrad_up(col,i) = &
            B1_3 - &
            B2_3 + &
            B3_3
          IF (do_lambertian(i)) THEN 
            ! Lambertian reflected downwelling layer-average
            auxrad_column%meanrad_down_p(col,i) = &
              surf_fac * rad_air_avg * &
              ((tau_surf_p_r) - (tau_level_p_r))
          ENDIF
          IF (opts%rt_all%rad_down_lin_tau) THEN
            ! Specular reflected downwelling linear-in-tau
             auxrad_column%meanrad_down(col,i) = &
                surf_fac * &
                (B1_3 - &
                B3_3) * &
                tau_level_r * tau_surf_r + &
                B2_3 * tau_surf_r**2
          ELSE
            ! Specular reflected downwelling layer-average
            auxrad_column%meanrad_down(col,i) = &
              surf_fac * rad_air_avg * &
              ((tau_surf_r) - (tau_level_r))
          ENDIF
        ENDIF

        IF (pol_id(i) >= 6_jpim) THEN
          auxrad_column%meanrad_up(col,i) = 0._jprb
        ELSE
          auxrad_column%meanrad_up(col,i) = auxrad_column%meanrad_up(col,i) + auxrad_column%up(lay,col,i)
        ENDIF

        auxrad_column%meanrad_down(col,i) = auxrad_column%meanrad_down(col,i) + auxrad_column%down(lay,col,i)
        auxrad_column%meanrad_down(col,i) = MAX(auxrad_column%meanrad_down(col,i), 0._jprb)
        IF (do_lambertian(i)) THEN
          auxrad_column%meanrad_down_p(col,i) = auxrad_column%meanrad_down_p(col,i) + auxrad_column%down_p(lay,col,i)
          auxrad_column%meanrad_down_p(col,i) = MAX(auxrad_column%meanrad_down_p(col,i), 0._jprb)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE calc_near_surf_contribution

  SUBROUTINE solar_scattering_air(transmission_aux, refl, raytracing, transmission_scatt_ir, auxrad_column)

    TYPE(rttov_transmission_aux),      INTENT(IN)    :: transmission_aux
    REAL(jprb),                        INTENT(IN)    :: refl(:)
    TYPE(rttov_raytracing),            INTENT(IN)    :: raytracing
    TYPE(rttov_transmission_scatt_ir), INTENT(IN)    :: transmission_scatt_ir
    TYPE(rttov_radiance_aux),          INTENT(INOUT) :: auxrad_column

    REAL(jprb)    :: temp(nlayers,0:maxncolumns)
    INTEGER(jpim) :: i, col, coli, lay, su

    ! The solar_path1 and solar_path2 quantities are completely consistent (i.e.
    ! they are based on the same optical depth regression, namely the solar one).
    ! Note that solar_path2%tau_level is on the sun-surface-satellite path
    ! but solar_path2%od_single_layer is on the sun-surface path.

#define fac1_2 auxrad_column%Fac1_2(coli,lay,i)
#define fac3_2 auxrad_column%Fac3_2(lay,i)
#define fac4_2 auxrad_column%Fac4_2(coli,lay,i)
#define fac5_2 auxrad_column%Fac5_2(coli,lay,i)
#define fac6_2 auxrad_column%Fac6_2(coli,lay,i)
#define fac7_2 auxrad_column%Fac7_2(lay,i)
#define dfac54_2 (fac5_2 - fac4_2)
#define tausun_layer transmission_aux%solar_path2%Tau_level(lay,col,i)

    su = 0
    IF (opts%rt_ir%addclouds) su = 1

    DO i = 1, nchanprof
      IF (.NOT. solar(i)) CYCLE

      auxrad_column%Fac6_2(0:su,:,i) = solar_spectrum(i) * z4pi_r * transmission_scatt_ir%phdo(0:su,:,i)
      auxrad_column%Fac1_2(0:su,:,i) = solar_spectrum(i) * z4pi_r * transmission_scatt_ir%phup(0:su,:,i)

      auxrad_column%Fac3_2(:,i) = raytracing%pathsat(:,prof) / raytracing%patheff(:,prof)

      auxrad_column%Fac4_2(0:su,:,i) = EXP(-transmission_aux%solar_path2%Od_singlelayer(0:su,:,i))
      auxrad_column%Fac5_2(0:su,:,i) = EXP(-transmission_aux%solar_path1%Od_singlelayer(0:su,:,i))

      DO lay = 1, nlayers
        IF (.NOT. sateqsun(lay,prof)) THEN
          auxrad_column%Fac7_2(lay,i) = raytracing%pathsat(lay,prof) / &
                                        (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))
        ENDIF
      ENDDO

      !----------------Upward single scattering of the solar beam-----------------------

      DO col = 0, ircld%ncolumn(prof)
        DO lay = 1, nlayers
          coli = ircld%icldarr(col,lay,prof)
          auxrad_column%up_solar(lay,col,i) = &
                        (auxrad_column%Fac1_2(coli,lay,i) * &
                        transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                        auxrad_column%Fac3_2(lay,i) * &
                        (transmission_aux%solar_path2%Tau_level(lay,col,i) - &
                        transmission_aux%solar_path2%Tau_level(lay+1,col,i)))
        ENDDO
      ENDDO

      DO col = 0, ircld%ncolumn(prof)
        DO lay = 2, nlayers
          auxrad_column%up_solar(lay,col,i) = auxrad_column%up_solar(lay,col,i) + auxrad_column%up_solar(lay-1,col,i)
        ENDDO
      ENDDO

      !-------------------Downward single scattering of the solar beam------------------

      IF (refl(i) > 0._jprb) THEN
        DO col = 0, ircld%ncolumn(prof)
          DO lay = 1, nlayers
            coli = ircld%icldarr(col,lay,prof)

            temp(lay,col) = fac2_solar_path1 * &
                            fac6_2 * &
                            transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                            tausun_layer * &
                            transmission_aux%solar_path1%Tau_level_r(lay,col,i) * &
                            transmission_aux%solar_path1%Tau_level_r(lay+1,col,i)

            IF (.NOT. sateqsun(lay,prof)) THEN
              temp(lay,col) = temp(lay,col) * &
                              fac7_2 * &
                              dfac54_2
            ELSE
              temp(lay,col) = temp(lay,col) * &
                              fac4_2 * &
                              transmission_aux%solar_path2%Od_singlelayer(coli,lay,i)
            ENDIF
          ENDDO
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          lay = 1_jpim
          auxrad_column%down_ref_solar(lay,col,i) = temp(lay,col)
          auxrad_column%down_solar(lay,col,i) = MAX(auxrad_column%down_ref_solar(lay,col,i), 0._jprb)
          DO lay = 2, nlayers
            temp(lay,col) = temp(lay-1,col) + temp(lay,col)
            auxrad_column%down_ref_solar(lay,col,i) = temp(lay,col)
            auxrad_column%down_solar(lay,col,i) = MAX(auxrad_column%down_ref_solar(lay,col,i), 0._jprb)
          ENDDO
        ENDDO
      ENDIF ! refl(i) > 0.
    ENDDO
  END SUBROUTINE solar_scattering_air

  SUBROUTINE solar_scattering_near_surf(transmission_aux, refl, auxrad_column)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    REAL(jprb),                   INTENT(IN)    :: refl(:)
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column

    INTEGER(jpim) :: i, col, coli, lay, lev, lay1, ncolumns
    REAL(jprb) :: temp(0:maxncolumns)

#define fac4_3 auxrad_column%Fac4_3(col,i)
#define fac5_3 auxrad_column%Fac5_3(col,i)
#define dfac54_3 (fac5_3 - fac4_3)
#define dtausun_surf (transmission_aux%solar_path2%Tau_level(lev,col,i) - transmission_aux%solar_path2%Tau_surf(col,i))
#define tausun_level transmission_aux%solar_path2%Tau_level(lev,col,i)

    DO i = 1, nchanprof
      IF (.NOT. solar(i)) CYCLE

      ! lay is the layer above the one containing the surface
      ! lev is the nearest layer above the surface
      lay = iv3lay(i)
      lev = iv3lev(i)
      ncolumns = ircld%ncolumn(prof)

      ! lay1 is the layer containing the surface or the bottom layer
      !   if the surface lies below the bottom of the profile
      IF (pfraction(i) < 0._jprb) THEN
        lay1 = lay
      ELSE
        lay1 = lay + 1
      ENDIF

      auxrad_column%Fac4_3(0:ncolumns,i) = EXP(-transmission_aux%solar_path2%od_sfrac(0:ncolumns,i))
      auxrad_column%Fac5_3(0:ncolumns,i) = EXP(-transmission_aux%solar_path1%od_sfrac(0:ncolumns,i))

      !--------------Upward single scattering of the solar beam-------------------------

      DO col = 0, ircld%ncolumn(prof)
        coli = ircld%icldarr(col,lay,prof)

        auxrad_column%meanrad_up_solar(col,i) = &
          fac1_2 * &
          transmission_scatt_ir%ssa_solar(coli,lay,i) * &
          fac3_2 * &
          dtausun_surf
      ENDDO

      DO col = 0, ircld%ncolumn(prof)
        auxrad_column%meanrad_up_solar(col,i) = auxrad_column%meanrad_up_solar(col,i) + &
                                                auxrad_column%up_solar(lay,col,i)
      ENDDO

      !--------------Downward single scattering of the solar beam-----------------------

      IF (refl(i) > 0._jprb) THEN
        DO col = 0, ircld%ncolumn(prof)
          temp(col) = fac2_solar_path1 * &
                      fac6_2 * &
                      transmission_scatt_ir%ssa_solar(coli,lay,i) * &
                      tausun_level * &
                      (transmission_aux%solar_path1%Tau_level_r(lay+1,col,i) * &
                      transmission_aux%solar_path1%Tau_surf_r(col,i))
        ENDDO

        IF (.NOT. sateqsun(lay1,prof)) THEN
          DO col = 0, ircld%ncolumn(prof)
            auxrad_column%meanrad_down_solar(col,i) = &
                    temp(col) * &
                    auxrad_column%Fac7_2(lay1,i) * &
                    dfac54_3
          ENDDO
        ELSE
          DO col = 0, ircld%ncolumn(prof)
            auxrad_column%meanrad_down_solar(col,i) = &
                                temp(col) * fac4_3 * &
                                transmission_aux%solar_path2%od_sfrac(col,i)
          ENDDO
        ENDIF

        DO col = 0, ircld%ncolumn(prof)
          auxrad_column%meanrad_down_solar(col,i) = auxrad_column%meanrad_down_solar(col,i) + &
                                                    auxrad_column%down_solar(lay,col,i)

          auxrad_column%meanrad_down_solar(col,i) = MAX(auxrad_column%meanrad_down_solar(col,i), 0._jprb)
        ENDDO
      ENDIF ! refl(i) > 0.
    ENDDO
  END SUBROUTINE solar_scattering_near_surf

  SUBROUTINE solar_rayleigh(raytracing, profiles, profiles_dry, transmission_aux, auxrad_column)

    TYPE(rttov_raytracing),       INTENT(IN)    :: raytracing
    TYPE(rttov_profile),          INTENT(IN)    :: profiles(:)
    TYPE(rttov_profile),          INTENT(IN)    :: profiles_dry(:)
    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column

    INTEGER(jpim) :: i, col, lay, lev
    REAL(jprb) :: wlm, ss_param, v_h2o(nlayers), m(nlayers), v_h2o_surf, m_surf
    REAL(jprb) :: cosscata_term1, cosscata_term2, cosscata
    REAL(jprb) :: ray_phase, solar_src, solar_src_updn
    REAL(jprb) :: rayrad_up(0:nlayers,0:maxncolumns), rayrad_dn(0:nlayers,0:maxncolumns)

    ! The phase function could be made to account for depolarisation.
    ! The solar geometry is approximated using the sun-surface path values,
    ! and similarly for the sun-level-satellite transmittances.

    ! Rayleigh scattering only included for channels less than rayleigh_max_wavelength option.

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

        ! Layer H2O by volume as fraction:
        v_h2o = 0.5_jprb * (profiles_dry(prof)%q(1:nlevels-1) + profiles_dry(prof)%q(2:nlevels)) * 1.E-6_jprb
        ! Convert ppmv dry to ppmv wet
        v_h2o = v_h2o / (1._jprb + v_h2o)
      
        ! Average molar weight of wet air for the layer (g)
        m = ((1._jprb - v_h2o) * Mair + v_h2o * Mh2o)
      ELSE
        ss_param = coef%ss_rayleigh_ext(chan) * (ray_ts / ray_ps) * rgc * z4pi_r / gravity
        m = mair
        m_surf = mair
      ENDIF

      rayrad_up(0,:) = 0._jprb
      rayrad_dn(0,:) = 0._jprb

      ! Sum contributions from atmospheric layers
      DO lev = 2, nlevels
        lay = lev - 1

        ! Skip calculation for whole layers above min pressure
        IF (profiles(prof)%p(lev) < opts%rt_ir%rayleigh_min_pressure) THEN
          rayrad_up(lay,:) = 0._jprb
          rayrad_dn(lay,:) = 0._jprb
          CYCLE
        ENDIF

        solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
              (profiles(prof)%p(lev) - profiles(prof)%p(lev-1)) * &
              raytracing%pathsat(lay, prof) * ss_param / m(lev-1)

        ! cosine of scattering angle - raytracing%zasat/zasun contain the sine of the angles
        cosscata_term1 = SQRT((1._jprb - raytracing%zasat(lay, prof) * raytracing%zasat(lay, prof)) * &
                         (1._jprb - raytracing%zasun(lay, prof) * raytracing%zasun(lay, prof)))
        cosscata_term2 = raytracing%zasat(lay, prof) * raytracing%zasun(lay, prof) * &
                         COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1._jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        ! Phase function symmetry means upwelling and downwelling radiances are the same
        DO col = 0, ircld%ncolumn(prof)
          rayrad_up(lay,col) = rayrad_up(lay-1,col) + solar_src_updn * &
                               transmission_aux%solar_path2%Tau_level(lev-1,col,i)
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          IF (transmission_aux%solar_path1%Tau_level(lev-1,col,i) > min_tau) THEN
            rayrad_dn(lay,col) = rayrad_dn(lay-1,col) + solar_src_updn * &
                                 transmission_aux%solar_path2%Tau_level(lev-1,col,i) / &
                                 transmission_aux%solar_path1%Tau_level(lev-1,col,i) ** 3_jpim
          ELSE
            rayrad_dn(lay,col) = rayrad_dn(lay-1,col)
          ENDIF
        ENDDO
      ENDDO

      ! Add Rayleigh contributions to radiance totals
      DO col = 0, ircld%ncolumn(prof)
        auxrad_column%up_solar(:,col,i) = auxrad_column%up_solar(:,col,i) + &
                                          rayrad_up(1:nlayers,col)

        auxrad_column%meanrad_up_solar(col,i) = auxrad_column%meanrad_up_solar(col,i) + &
                                                rayrad_up(iv3lay(i),col)

        auxrad_column%down_solar(:,col,i) = auxrad_column%down_solar(:,col,i) + &
                                            rayrad_dn(1:nlayers,col)

        auxrad_column%meanrad_down_solar(col,i) = auxrad_column%meanrad_down_solar(col,i) + &
                                                  rayrad_dn(iv3lay(i),col)
      ENDDO

      ! Calculate the contribution from the part-layer above the surface

      ! Skip calculation for whole layers above min pressure
      IF (profiles(prof)%s2m%p >= opts%rt_ir%rayleigh_min_pressure) THEN

        IF (coef%fmv_model_ver <= 9) THEN
          ! Layer H2O by volume as fraction:
          IF (opts%rt_all%use_q2m) THEN
            v_h2o_surf = 0.5_jprb * (profiles_dry(prof)%q(iv3lev(i)) + profiles_dry(prof)%s2m%q) * 1.E-6_jprb
          ELSE
            v_h2o_surf = profiles_dry(prof)%q(iv3lev(i)) * 1.E-6_jprb
          ENDIF
          ! Convert ppmv dry to ppmv wet
          v_h2o_surf = v_h2o_surf / (1._jprb + v_h2o_surf)
        
          ! Average molar weight of wet air for the layer (g):
          m_surf = ((1._jprb - v_h2o_surf) * Mair + v_h2o_surf * Mh2o)
        ENDIF

        solar_src = solar_spectrum(i) * & ! mW m^-2 (cm^-1)^-1
              ABS(profiles(prof)%s2m%p - profiles(prof)%p(iv3lev(i))) * &
              raytracing%pathsat(iv2lay(i), prof) * ss_param / m_surf

        ! cosine of scattering angle
        cosscata_term1 = SQRT((1._jprb - raytracing%zasat(iv2lay(i), prof) * raytracing%zasat(iv2lay(i), prof)) * &
                          (1._jprb - raytracing%zasun(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof)))
        cosscata_term2 = raytracing%zasat(iv2lay(i), prof) * raytracing%zasun(iv2lay(i), prof) * &
                          COS((profiles(prof)%azangle - profiles(prof)%sunazangle)*deg2rad)

        cosscata = - cosscata_term1 - cosscata_term2
        ray_phase = 0.75_jprb * (1._jprb + cosscata * cosscata)
        solar_src_updn = solar_src * ray_phase

        ! Add near-surface contributions to the total radiances
        DO col = 0, ircld%ncolumn(prof)
          auxrad_column%meanrad_up_solar(col,i) = auxrad_column%meanrad_up_solar(col,i) + solar_src_updn * &
                                                  transmission_aux%solar_path2%Tau_level(iv3lev(i),col,i)
        ENDDO

        DO col = 0, ircld%ncolumn(prof)
          IF (transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) > min_tau) THEN
            auxrad_column%meanrad_down_solar(col,i) = auxrad_column%meanrad_down_solar(col,i) + solar_src_updn * &
                                                      transmission_aux%solar_path2%Tau_level(iv3lev(i),col,i) / &
                                                      transmission_aux%solar_path1%Tau_level(iv3lev(i),col,i) ** 3_jpim
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE solar_rayleigh

  SUBROUTINE solar_surface_contribution(transmission_aux, reflectance, auxrad_column)

    TYPE(rttov_transmission_aux), INTENT(IN)    :: transmission_aux
    REAL(jprb),                   INTENT(IN)    :: reflectance(:)
    TYPE(rttov_radiance_aux),     INTENT(INOUT) :: auxrad_column

    INTEGER(jpim) :: i, ncolumns

    DO i = 1, nchanprof
      IF (solar(i)) THEN
        ncolumns = ircld%ncolumn(prof)
        auxrad_column%cloudy(0:ncolumns,i) = auxrad_column%cloudy(0:ncolumns,i) + &
               solar_spectrum(i) * reflectance(i) * refl_norm(i) * &
               transmission_aux%solar_path2%Tau_surf(0:ncolumns,i)
      ENDIF
    ENDDO
  END SUBROUTINE solar_surface_contribution

END SUBROUTINE rttov_integrate
