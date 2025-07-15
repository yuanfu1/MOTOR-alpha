! Description:
!> @file
!!   Allocate and initialise/deallocate "static" (i.e. direct-model only)
!!   internal data structure.
!
!> @brief
!!   Allocate and initialise/deallocate "static" (i.e. direct-model only)
!!   internal data structure.
!!
!! @details
!!   The static trajectory structure contains internal direct model variables
!!   which are re-used by the TL, AD and K models for efficiency. This
!!   subroutine is used to allocate/deallocate these variables.
!!
!!   When allocating this subroutine also initialises various flags which
!!   are subsequently used to determine which other structures are required:
!!   - thermal and solar flags for each channel indicating whether emitted
!!     and/or solar radiation should be included
!!   - overall dothermal and dosolar flags which indicate if any emission
!!     or solar calculations are required at all
!!   - flags indicating whether MFASIS or DOM have been selected as
!!     scattering solvers
!!
!!   Some other variables/flags are also initialised here:
!!   - the plane_parallel flag which can be enabled manually in the options
!!     structure or in the lbl_check argument (intended for developers only)
!!     but is enforced when the DOM solver is being used
!!   - ensures the number of DOM streams is valid
!!   - do_lambertian array holding flags to indicate if the Lambertian option
!!     is valid for each channel (cannot be used with PC-RTTOV or with any sea
!!     surface emissivity model)
!!
!!   It is possible to gain access to some of the internal RTTOV state by
!!   passing the (unallocated) traj_sta argument into rttov_direct: on exit
!!   the structure is allocated and populated. This could potentially be
!!   useful for debugging purposes.
!!
!! @param[out]    err            status on exit
!! @param[in,out] traj_sta       RTTOV static trajectory structure
!! @param[in]     opts           options to configure the simulations
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     npcscores      number of PC scores requested, optional
!! @param[in]     channels_rec   list of channels for which to calculate reconstructed radiances, optional
!! @param[in]     calcemis       flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] lbl_check      used for coef verification, optional, not intended for general use, optional
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_alloc_traj_sta(err, traj_sta, opts, coefs, chanprof, profiles, asw, &
                                npcscores, channels_rec, calcemis, lbl_check)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : &
      rttov_traj_sta, &
      rttov_options,  &
      rttov_coefs,    &
      rttov_chanprof, &
      rttov_profile,  &
      rttov_lbl_check
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : &
      max_sol_zen,      &
      sensor_id_ir,     &
      sensor_id_hi,     &
      vis_scatt_mfasis, &
      vis_scatt_dom,    &
      ir_scatt_dom,     &
      dom_min_nstr,     &
      surftype_sea
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),         INTENT(OUT)          :: err
  TYPE(rttov_traj_sta),  INTENT(INOUT)        :: traj_sta
  TYPE(rttov_options),   INTENT(IN)           :: opts
  TYPE(rttov_coefs),     INTENT(IN)           :: coefs
  TYPE(rttov_chanprof),  INTENT(IN)           :: chanprof(:)
  TYPE(rttov_profile),   INTENT(IN)           :: profiles(:)
  INTEGER(jpim),         INTENT(IN)           :: asw
  INTEGER(jpim),         INTENT(IN), OPTIONAL :: npcscores
  INTEGER(jpim),         INTENT(IN), OPTIONAL :: channels_rec(:)
  LOGICAL(jplm),         INTENT(IN), OPTIONAL :: calcemis(SIZE(chanprof))
  TYPE(rttov_lbl_check), INTENT(IN), OPTIONAL :: lbl_check
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_pc_dimensions.interface"
#include "rttov_calc_solar_spec_esd.interface"

  INTEGER(jpim)       :: nlevels, nlayers, nchanprof, nprofiles
  INTEGER(jpim)       :: nlevels_coef, nlayers_coef, j, prof, chan
  LOGICAL(jplm)       :: do_dom, do_rayleigh
  TYPE(rttov_options) :: opts_coef

  TRY

  nlevels = profiles(1)%nlevels
  nlayers = nlevels - 1
  nchanprof = SIZE(chanprof)
  nprofiles = SIZE(profiles)

  nlevels_coef = coefs%coef%nlevels
  nlayers_coef = coefs%coef%nlayers

  opts_coef = opts
  opts_coef%rt_ir%addclouds = .FALSE.
  opts_coef%rt_ir%addaerosl = .FALSE.

  ! ---------------------------------------------------------------------------
  ! Allocation/initialisation
  ! ---------------------------------------------------------------------------
  IF (asw == 1_jpim) THEN

    ALLOCATE(traj_sta%thermal(nchanprof),        &
             traj_sta%solar(nchanprof),          &
             traj_sta%do_lambertian(nchanprof),  &
             traj_sta%angles(nprofiles),         &
             ! These are passed into rttov_integrate so are required even when addsolar is false:
             traj_sta%solar_spec_esd(nchanprof), &
             traj_sta%refl_norm(nchanprof),      &
             STAT = err)
    THROWM(err.NE.0,"Allocation of traj_sta failed")

    ! ---------------------------------------------------------------------------
    ! Set flags indicating if MFASIS or DOM have been selected
    ! ---------------------------------------------------------------------------
    ! NB This tests whether MFASIS/DOM simulations MIGHT be carried out (depending on channel
    ! selection). It does NOT check whether MFASIS/DOM simulations are ACTUALLY being carried
    ! out using dothermal/dosolar: if RTTOV is run through the parallel interface the
    ! latter test can cause some profiles to be run with plane_parallel and others not for
    ! MFASIS/DOM depending on how profiles/channels are grouped.
    traj_sta%do_mfasis = (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
                         opts%rt_ir%vis_scatt_model == vis_scatt_mfasis   .AND. &
                         opts%rt_ir%addsolar
    ! The DOM flag is only used locally, it is not stored
    do_dom = (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) .AND. &
             ((opts%rt_ir%vis_scatt_model == vis_scatt_dom .AND. opts%rt_ir%addsolar) .OR. &
               opts%rt_ir%ir_scatt_model == ir_scatt_dom)

    ! ---------------------------------------------------------------------------
    ! Set thermal flag for emissive calculations
    ! ---------------------------------------------------------------------------
    traj_sta%thermal(:) = .TRUE._jplm

    ! Turn off thermal calculations for any pure-solar or invalid channels
    ! ss_val_chn = 0 => thermal-only
    ! ss_val_chn = 1 => thermal + solar
    ! ss_val_chn = 2 => solar-only
    DO j = 1, nchanprof
      chan = chanprof(j)%chan
      traj_sta%thermal(j) = .NOT. (coefs%coef%ss_val_chn(chan) == 2_jpim) .AND. &
                            .NOT. (coefs%coef%ff_val_chn(chan) == 0_jpim)
    ENDDO

    traj_sta%dothermal = ANY(traj_sta%thermal(:))

    ! ---------------------------------------------------------------------------
    ! Set solar flag to include solar radiation
    ! ---------------------------------------------------------------------------
    traj_sta%solar(:) = .FALSE._jplm
    traj_sta%dosolar  = .FALSE._jplm
    IF (opts%rt_ir%addsolar) THEN ! If addsolar is FALSE no solar calculations will be done at all

      IF ((coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) .AND. &
           coefs%coef%fmv_model_ver >= 9) THEN

        ! Turn on solar calculations for any solar-affected, valid channels
        DO j = 1, nchanprof
          prof = chanprof(j)%prof
          chan = chanprof(j)%chan

          IF (profiles(prof)%sunzenangle >= 0.0_jprb .AND. &
              profiles(prof)%sunzenangle < max_sol_zen) THEN

            ! Ensure channels are solar-enabled and not invalid
            traj_sta%solar(j) = .NOT. (coefs%coef%ss_val_chn(chan) == 0_jpim) .AND. &
                                .NOT. (coefs%coef%ff_val_chn(chan) == 0_jpim)

            ! For MFASIS enable simulations only for channels in the LUT (these might
            ! only be a subset of the solar channels in the rtcoef file)
            IF (opts%rt_ir%vis_scatt_model == vis_scatt_mfasis) THEN
              IF (opts%rt_ir%addclouds) THEN
                traj_sta%solar(j) = traj_sta%solar(j) .AND. &
                                    coefs%coef_mfasis_cld%channel_lut_index(chan) > 0
              ENDIF
              IF (opts%rt_ir%addaerosl) THEN
                traj_sta%solar(j) = traj_sta%solar(j) .AND. &
                                    coefs%coef_mfasis_aer%channel_lut_index(chan) > 0
              ENDIF
            ENDIF

          ENDIF
        ENDDO

        traj_sta%dosolar = ANY(traj_sta%solar(:))

      ENDIF
    ENDIF

    ! ---------------------------------------------------------------------------
    ! Populate ToA solar irradiance array
    ! ---------------------------------------------------------------------------
    IF (traj_sta%dosolar) CALL rttov_calc_solar_spec_esd(coefs%coef, chanprof, profiles, traj_sta%solar_spec_esd)

    ! ---------------------------------------------------------------------------
    ! Rayleigh scattering flags
    ! ---------------------------------------------------------------------------
    ! Rayleigh flag true if solar active and minimum channel wavelength is below
    ! the user-specified Rayleigh max wavelength (and not MFASIS)
    do_rayleigh = .NOT. traj_sta%do_mfasis .AND. traj_sta%dosolar .AND. &
      10000._jprb / MAXVAL(coefs%coef%ff_cwn) < opts%rt_ir%rayleigh_max_wavelength

    ! The Rayleigh extinction parameterisation is turned on for v13 predictors
    traj_sta%do_rayleigh_param = do_rayleigh .AND. &
                                 coefs%coef%fmv_model_ver >= 13

    ! DOM Rayleigh is activated if solar DOM calculations are active and user
    ! selects the rayleigh_dom option with v13 predictor rtcoef file
    traj_sta%do_rayleigh_dom = traj_sta%do_rayleigh_param .AND. &
                               opts%rt_ir%dom_rayleigh .AND. do_dom .AND. &
                               opts%rt_ir%vis_scatt_model == vis_scatt_dom

    ! Rayleigh single-scattering is activated the single-scattering option is
    ! true and DOM Rayleigh is not active
    traj_sta%do_rayleigh_ss = do_rayleigh .AND. &
                              opts%rt_ir%rayleigh_single_scatt .AND. &
                              .NOT. traj_sta%do_rayleigh_dom

    ! ---------------------------------------------------------------------------
    ! Set RTTOV gas optical depth calculation switches
    ! ---------------------------------------------------------------------------
    traj_sta%do_opdep_calc = traj_sta%dothermal .OR. .NOT. traj_sta%do_mfasis
    traj_sta%do_solar_opdep_calc = traj_sta%dosolar .AND. .NOT. traj_sta%do_mfasis

    ! If the developer option is false then turn off the opdep calculation, but
    ! don't modify the do_opdep_calc switch if the developer option is true
    IF (.NOT. opts%dev%do_opdep_calc) traj_sta%do_opdep_calc = .FALSE.

    ! Note that do_opdep_calc enables/disables the whole optical depth
    ! parameterisation, but all other calculations on user levels (including
    ! transmittances) are carried out.

    ! The do_solar_opdep_calc flag is specifically for disabling the solar
    ! optical depth/transmittance/ToA refl calculations for solar channels if
    ! MFASIS is being run. Therefore we do not set do_solar_opdep_calc to false
    ! if the opts%dev%do_opdep_calc has been set to false.

    ! ---------------------------------------------------------------------------
    ! Set plane-parallel geometry switch
    ! ---------------------------------------------------------------------------
    IF (PRESENT(lbl_check)) THEN
      ! For lbl checks we may want to test in a plane geometry
      traj_sta%plane_parallel = lbl_check%plane_geometry
    ELSE
      traj_sta%plane_parallel = do_dom .OR. opts%rt_all%plane_parallel
    ENDIF

    ! ---------------------------------------------------------------------------
    ! For DOM ensure the number of streams is even and within allowed bounds
    ! ---------------------------------------------------------------------------
    IF (do_dom) THEN
      traj_sta%dom_nstreams = MAX(dom_min_nstr, opts%rt_ir%dom_nstreams)
      traj_sta%dom_nstreams = traj_sta%dom_nstreams + MOD(traj_sta%dom_nstreams, 2)
    ELSE
      traj_sta%dom_nstreams = 0
    ENDIF

    ! ---------------------------------------------------------------------------
    ! Determine whether Lambertian option should be switched on for each channel:
    ! it is not used with PC-RTTOV or with RTTOV's sea-surface emissivity models
    ! ---------------------------------------------------------------------------
    IF (PRESENT(calcemis) .AND. opts%rt_all%do_lambertian) THEN
      DO j = 1, nchanprof
        traj_sta%do_lambertian(j) = (profiles(chanprof(j)%prof)%skin%surftype /= surftype_sea) .OR. .NOT. calcemis(j)
      ENDDO
    ELSE
      traj_sta%do_lambertian(:) = opts%rt_all%do_lambertian .AND. .NOT. opts%rt_ir%pc%addpc
    ENDIF

    ! ---------------------------------------------------------------------------
    ! Allocate the remaining structures
    ! ---------------------------------------------------------------------------
    IF (traj_sta%dothermal) THEN
      ALLOCATE(traj_sta%thermal_path1, STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for thermal path1 failed")
      ALLOCATE(traj_sta%thermal_path1%tau_surf(nchanprof),           &
               traj_sta%thermal_path1%od_level(nlevels, nchanprof),  &
               traj_sta%thermal_path1%tau_level(nlevels, nchanprof), &
               STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for thermal path1 failed")

      IF (traj_sta%do_opdep_calc) THEN 
        IF (coefs%coef%fmv_model_ver > 9) THEN
          CALL alloc_opdp_ref(err, traj_sta%thermal_path1%opdp_ref)
        ELSE
          ALLOCATE(traj_sta%thermal_path1%opdp_ref_coef(nlayers_coef, nchanprof), STAT = err)
        ENDIF
        THROWM(err.NE.0,"Allocation of traj_sta for thermal path1 failed")
      ENDIF
    ENDIF

    IF (traj_sta%dosolar .AND. .NOT. traj_sta%do_mfasis) THEN
      ALLOCATE(traj_sta%solar_path2, STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path2 failed")
      ALLOCATE(traj_sta%solar_path2%tau_level(nlevels, nchanprof),      &
               traj_sta%solar_path2%tau_surf(nchanprof),                &
               traj_sta%solar_path2%od_level(nlevels, nchanprof),       &
               traj_sta%solar_path2%od_singlelayer(nlayers, nchanprof), &
               traj_sta%solar_path2%od_frac(nchanprof),                 &
               STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path2 failed")

      IF (traj_sta%do_opdep_calc) THEN 
        IF (coefs%coef%fmv_model_ver > 9) THEN
          CALL alloc_opdp_ref(err, traj_sta%solar_path2%opdp_ref)
        ELSE
          ALLOCATE(traj_sta%solar_path2%opdp_ref_coef(nlayers_coef, nchanprof), STAT = err)
        ENDIF
        THROWM(err.NE.0,"Allocation of traj_sta for solar path1 failed")
      ENDIF

      ALLOCATE(traj_sta%solar_path1, STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path1 failed")
      ALLOCATE(traj_sta%solar_path1%tau_level(nlevels, nchanprof), &
               traj_sta%solar_path1%tau_surf(nchanprof),           &
               STAT = err)
      THROWM(err.NE.0,"Allocation of traj_sta for solar path1 failed")
    ENDIF

    IF (traj_sta%do_opdep_calc) THEN
      IF (opts%config%apply_reg_limits .OR. opts%interpolation%addinterp) THEN
        ALLOCATE(traj_sta%profiles_coef_ref(nprofiles), STAT = err)
        THROWM(err.NE.0,"Allocation of traj_sta%profiles_coef_ref failed")

        CALL rttov_alloc_prof(err, nprofiles, traj_sta%profiles_coef_ref, nlevels_coef, opts_coef, asw)
        THROW(err.NE.0)
      ENDIF
    ENDIF

    IF (opts%rt_ir%pc%addpc) THEN
      CALL rttov_alloc_pc_dimensions( &
              err,                   &
              opts,                  &
              npcscores,             &
              nprofiles,             &
              traj_sta%chanprof_in,  &
              traj_sta%chanprof_pc,  &
              asw,                   &
              channels_rec = channels_rec)
      THROW(err.NE.0)

      IF (opts%rt_ir%addaerosl) THEN
        ALLOCATE(traj_sta%pc_aer_ref(SIZE(profiles(1)%aerosols, DIM=1),nlayers,nprofiles), &
                 traj_sta%pc_aer_min(SIZE(profiles(1)%aerosols, DIM=1),nlayers,nprofiles), &
                 traj_sta%pc_aer_max(SIZE(profiles(1)%aerosols, DIM=1),nlayers,nprofiles), STAT = err)
        THROWM(err.NE.0,"Allocation of traj_sta%pc_aer_ref/min/max failed")
      ENDIF
    ENDIF

  ENDIF


  ! ---------------------------------------------------------------------------
  ! Deallocation
  ! ---------------------------------------------------------------------------
  IF (asw == 0_jpim) THEN

    IF (opts%rt_ir%pc%addpc) THEN
      CALL rttov_alloc_pc_dimensions( &
              err,                   &
              opts,                  &
              npcscores,             &
              nprofiles,             &
              traj_sta%chanprof_in,  &
              traj_sta%chanprof_pc,  &
              asw,                   &
              channels_rec = channels_rec)
      THROW(err.NE.0)

      IF (opts%rt_ir%addaerosl) THEN
        DEALLOCATE(traj_sta%pc_aer_ref, &
                   traj_sta%pc_aer_min, &
                   traj_sta%pc_aer_max, STAT = err)
        THROWM(err.NE.0,"Deallocation of traj_sta%profiles_aer_ref/min/max failed")
      ENDIF
    ENDIF

    IF (traj_sta%do_opdep_calc) THEN
      IF (opts%config%apply_reg_limits .OR. opts%interpolation%addinterp) THEN
        CALL rttov_alloc_prof (err, nprofiles, traj_sta%profiles_coef_ref, nlevels_coef, opts_coef, asw)
        THROW(err.NE.0)

        DEALLOCATE(traj_sta%profiles_coef_ref, STAT = err)
        THROWM(err.NE.0,"Deallocation of traj_sta%profiles_coef_ref failed")
      ENDIF
    ENDIF

    IF (traj_sta%dothermal) THEN
      IF (traj_sta%do_opdep_calc) THEN 
        IF (coefs%coef%fmv_model_ver > 9) THEN
          CALL alloc_opdp_ref(err, traj_sta%thermal_path1%opdp_ref)
        ELSE
          DEALLOCATE(traj_sta%thermal_path1%opdp_ref_coef, STAT = err)
        ENDIF
        THROWM(err.NE.0,"Deallocation of traj_sta for thermal path1 failed")
      ENDIF

      DEALLOCATE(traj_sta%thermal_path1%tau_surf,      &
                 traj_sta%thermal_path1%od_level,      &
                 traj_sta%thermal_path1%tau_level,     &
                 STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for thermal path1 failed")
      DEALLOCATE(traj_sta%thermal_path1, STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for thermal path1 failed")
    ENDIF

    IF (traj_sta%dosolar .AND. .NOT. traj_sta%do_mfasis) THEN
      IF (traj_sta%do_opdep_calc) THEN 
        IF (coefs%coef%fmv_model_ver > 9) THEN
          CALL alloc_opdp_ref(err, traj_sta%solar_path2%opdp_ref)
        ELSE
          DEALLOCATE(traj_sta%solar_path2%opdp_ref_coef, STAT = err)
        ENDIF
        THROWM(err.NE.0,"Deallocation of traj_sta for solar path1 failed")
      ENDIF

      DEALLOCATE(traj_sta%solar_path2%tau_level,      &
                 traj_sta%solar_path2%tau_surf,       &
                 traj_sta%solar_path2%od_level,       &
                 traj_sta%solar_path2%od_singlelayer, &
                 traj_sta%solar_path2%od_frac,        &
                 STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path2 failed")
      DEALLOCATE(traj_sta%solar_path2, STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path2 failed")

      DEALLOCATE(traj_sta%solar_path1%tau_level,    &
                 traj_sta%solar_path1%tau_surf,     &
                 STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path1 failed")
      DEALLOCATE(traj_sta%solar_path1, STAT = err)
      THROWM(err.NE.0,"Deallocation of traj_sta for solar path1 failed")
    ENDIF

    DEALLOCATE(traj_sta%do_lambertian,  &
               traj_sta%thermal,        &
               traj_sta%solar,          &
               traj_sta%angles,         &
               traj_sta%solar_spec_esd, &
               traj_sta%refl_norm,      &
               STAT = err)
    THROWM(err.NE.0,"Deallocation of traj_sta failed")

  ENDIF

  CATCH
  
CONTAINS

  SUBROUTINE alloc_opdp_ref(err, od_ref)
    USE rttov_types, ONLY : rttov_opdp_ref_coef
    INTEGER(jpim),             INTENT(OUT)   :: err
    TYPE(rttov_opdp_ref_coef), INTENT(INOUT) :: od_ref

    TRY
    IF (asw == 1) THEN

      ALLOCATE(od_ref%od_tot_ref(nlayers_coef, nchanprof), &
               od_ref%od_mg_ref(nlayers_coef, nchanprof),  &
               od_ref%od_wv_ref(nlayers_coef, nchanprof),  &
               od_ref%od_wvcont_ref(nlayers_coef, nchanprof), STAT = err)
      THROW(err.NE.0)
      IF (coefs%coef%nozone > 0) THEN
        ALLOCATE(od_ref%od_o3_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF
      IF (coefs%coef%nco2 > 0) THEN
        ALLOCATE(od_ref%od_co2_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF
      IF (coefs%coef%nn2o > 0) THEN
        ALLOCATE(od_ref%od_n2o_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF
      IF (coefs%coef%nco > 0) THEN
        ALLOCATE(od_ref%od_co_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF
      IF (coefs%coef%nch4 > 0) THEN
        ALLOCATE(od_ref%od_ch4_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF
      IF (coefs%coef%nso2 > 0) THEN
        ALLOCATE(od_ref%od_so2_ref(nlayers_coef, nchanprof), STAT = err)
        THROW(err.NE.0)
      ENDIF

    ELSEIF (asw == 0) THEN

      IF (ASSOCIATED(od_ref%od_tot_ref))    DEALLOCATE(od_ref%od_tot_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_mg_ref))     DEALLOCATE(od_ref%od_mg_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_wv_ref))     DEALLOCATE(od_ref%od_wv_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_wvcont_ref)) DEALLOCATE(od_ref%od_wvcont_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_o3_ref))     DEALLOCATE(od_ref%od_o3_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_co2_ref))    DEALLOCATE(od_ref%od_co2_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_n2o_ref))    DEALLOCATE(od_ref%od_n2o_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_co_ref))     DEALLOCATE(od_ref%od_co_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_ch4_ref))    DEALLOCATE(od_ref%od_ch4_ref, STAT = err)
      THROW(err.NE.0)
      IF (ASSOCIATED(od_ref%od_so2_ref))    DEALLOCATE(od_ref%od_so2_ref, STAT = err)
      THROW(err.NE.0)

      NULLIFY(od_ref%od_tot_ref, od_ref%od_mg_ref, od_ref%od_wv_ref, &
              od_ref%od_wvcont_ref, od_ref%od_o3_ref, od_ref%od_co2_ref, &
              od_ref%od_n2o_ref, od_ref%od_co_ref, od_ref%od_ch4_ref, &
              od_ref%od_so2_ref)

    ENDIF
    CATCH
  END SUBROUTINE alloc_opdp_ref
END SUBROUTINE
