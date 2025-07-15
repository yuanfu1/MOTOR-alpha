! Description:
!> @file
!!   Throw a fatal error for illegal or strictly incompatible options.
!
!> @brief
!!   Throw a fatal error for illegal or strictly incompatible options.
!!
!> @details
!!   This is called within RTTOV and gives fatal errors for strictly
!!   illegal inputs:
!!   - incompatible options/coefs
!!   - missing optional arguments which are mandatory given selected options
!!   - illegal values in options
!!
!!   It does not give warnings for dubious, but harmless settings of options
!!   (that is done by rttov_user_options_checkinput).
!!
!!   It checks for profile data which are incompatible with options or coefs
!!   but it does not otherwise check for illegal values: this is done by
!!   rttov_check_profiles and rttov_check_reg_limits or by
!!   rttov_user_profile_checkinput.
!!
!! @param[out]    err                  status on exit
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     coefs                coefficients structure for instrument to simulate
!! @param[in]     chanprof             specifies channels and profiles to simulate
!! @param[in]     profiles             input atmospheric and surface data
!! @param[in]     aer_opt_param        input aerosol optical parameters, optional
!! @param[in]     cld_opt_param        input cloud optical parameters, optional
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
SUBROUTINE rttov_check_options( &
       err,           &
       opts,          &
       coefs,         &
       chanprof,      &
       profiles,      &
       aer_opt_param, &
       cld_opt_param)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : &
      rttov_options,  &
      rttov_coefs,    &
      rttov_chanprof, &
      rttov_profile,  &
      rttov_opt_param

  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK

  USE parkind1, ONLY : jprb

  USE rttov_const, ONLY : &
      sensor_id_ir,             &
      sensor_id_hi,             &
      sensor_id_mw,             &
      sensor_id_po,             &
      ninterp_modes,            &
      max_solar_sea_brdf_model, &
      max_ir_sea_emis_model,    &
      max_fastem_version,       &
      max_mw_clw_scheme,        &
      max_ir_scatt_model,       &
      max_vis_scatt_model,      &
      vis_scatt_mfasis,         &
      vis_scatt_dom,            &
      ir_scatt_dom,             &
      ir_scatt_chou,            &
      dom_min_nstr,             &
      clw_scheme_opac,          &
      clw_scheme_deff,          &
      ice_scheme_baum,          &
      mfasis_cld,               &
      mfasis_aer,               &
      aer_id_opac,              &
      ncloud_overlap,           &
      cloud_overlap_simple,     &
      ray_max_wlm
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),         INTENT(OUT)          :: err
  TYPE(rttov_options),   INTENT(IN)           :: opts
  TYPE(rttov_coefs),     INTENT(IN)           :: coefs
  TYPE(rttov_chanprof),  INTENT(IN)           :: chanprof(:)
  TYPE(rttov_profile),   INTENT(IN)           :: profiles(:)
  TYPE(rttov_opt_param), INTENT(IN), OPTIONAL :: aer_opt_param
  TYPE(rttov_opt_param), INTENT(IN), OPTIONAL :: cld_opt_param
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim)      :: i, chan, prof, pol_id
  INTEGER(jpim)      :: chan_m1, chan_m2, chan_m3
  CHARACTER(LEN=256) :: msg
  REAL(jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_OPTIONS',0_jpim,ZHOOK_HANDLE)

  ! Ensure units are the same for all profiles
  IF (ANY(profiles(:)%gas_units /= profiles(1)%gas_units)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'All profiles(:)%gas_units must be the same')
  ENDIF
  IF ((opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) .OR. &
      (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param)) THEN
    IF (ANY(profiles(:)%mmr_cldaer .NEQV. profiles(1)%mmr_cldaer)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'All profiles(:)%mmr_cldaer must be the same')
    ENDIF
  ENDIF

  ! Check for illegal chanprof values
  IF (MAXVAL(chanprof(:)%chan) > coefs%coef%fmv_chn) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Channel index in chanprof out of range for coefficients')
  ENDIF
  IF (MAXVAL(chanprof(:)%prof) > SIZE(profiles)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Profile index in chanprof out of range for input profiles(:) array')
  ENDIF
  IF (ANY(chanprof(:)%chan < 1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Channel index in chanprof must be greater than 0')
  ENDIF
  IF (ANY(chanprof(:)%prof < 1)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'Profile index in chanprof must be greater than 0')
  ENDIF

  ! Check number of levels
  IF (.NOT. opts%interpolation%addinterp) THEN
    IF (profiles(1)%nlevels .NE. coefs%coef%nlevels) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, 'The interpolator must be used if the input levels differ to the coefficient levels')
    ENDIF
  ENDIF

  ! Check for invalid polarimetric channel selection
  IF (coefs%coef%id_sensor == sensor_id_po) THEN
    DO i = 1, SIZE(chanprof)
      prof = chanprof(i)%prof
      chan = chanprof(i)%chan
      pol_id = coefs%coef%fastem_polar(chan) + 1
      IF (pol_id == 6) THEN
        IF (i < 3) THEN
          err = errorstatus_fatal
        ELSE
          chan_m2 = chanprof(i-2)%chan
          chan_m1 = chanprof(i-1)%chan
          IF (chanprof(i-2)%prof /= prof .OR. &
              chanprof(i-1)%prof /= prof .OR. &
              chan_m2 /= chan-2 .OR. &
              chan_m1 /= chan-1 .OR. &
              coefs%coef%fastem_polar(chan_m2) + 1 /= 4 .OR. &
              coefs%coef%fastem_polar(chan_m1) + 1 /= 5) THEN
            err = errorstatus_fatal
          ENDIF
        ENDIF
      ELSEIF (pol_id == 7) THEN
        IF (i < 4) THEN
          err = errorstatus_fatal
        ELSE
          chan_m3 = chanprof(i-3)%chan
          chan_m2 = chanprof(i-2)%chan
          IF (chanprof(i-2)%prof /= prof .OR. &
              chanprof(i-1)%prof /= prof .OR. &
              chan_m3 /= chan-3 .OR. &
              chan_m2 /= chan-2 .OR. &
              coefs%coef%fastem_polar(chan_m3) + 1 /= 4 .OR. &
              coefs%coef%fastem_polar(chan_m2) + 1 /= 5) THEN
            err = errorstatus_fatal
          ENDIF
        ENDIF
      ENDIF
      THROWM(err.NE.0, 'Error in polarimetric channel selection')
    ENDDO
  ENDIF

  ! Ensure opt_param arguments are present when required
  IF (opts%rt_ir%user_aer_opt_param .AND. .NOT. PRESENT(aer_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'If opts%rt_ir%user_aer_opt_param is TRUE then aer_opt_param is mandatory')
  ENDIF
  IF (opts%rt_ir%user_cld_opt_param .AND. .NOT. PRESENT(cld_opt_param)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'If opts%rt_ir%user_cld_opt_param is TRUE then cld_opt_param is mandatory')
  ENDIF

  ! Check options for valid values and also check that required arrays have been allocated
  ! where required (e.g. cfrac for IR scattering)
  IF (opts%interpolation%addinterp) THEN
    IF (opts%interpolation%interp_mode < 1 .OR. &
        opts%interpolation%interp_mode > ninterp_modes) THEN
      err = errorstatus_fatal
      msg = 'Error in specified interpolation mode'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%ozone_data .AND. coefs%coef%nozone > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%o3)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for ozone, but ozone_data flag is true; '// &
            'opts%rt_all%ozone_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%co2_data .AND. coefs%coef%nco2 > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%co2)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for co2, but co2_data flag is true; '// &
            'opts%rt_all%co2_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%co_data .AND. coefs%coef%nco > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%co)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for co, but co_data flag is true; '// &
            'opts%rt_all%co_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%n2o_data .AND. coefs%coef%nn2o > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%n2o)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for n2o, but n2o_data flag is true; '// &
            'opts%rt_all%n2o_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%ch4_data .AND. coefs%coef%nch4 > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%ch4)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for ch4, but ch4_data flag is true; '// &
            'opts%rt_all%ch4_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%so2_data .AND. coefs%coef%nso2 > 0) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%so2)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for so2, but so2_data flag is true; '// &
            'opts%rt_all%so2_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_mw%clw_data) THEN
    IF (.NOT. ASSOCIATED(profiles(1)%clw)) THEN
      err = errorstatus_fatal
      msg = 'Profiles structure not allocated for clw, but clw_data flag is true; '// &
            'opts%rt_mw%clw_data must be true when calling rttov_alloc_prof'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF


  ! IR instruments
  IF (coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) THEN

    IF (opts%rt_ir%addsolar) THEN
      IF (opts%rt_ir%solar_sea_brdf_model < 1 .OR. &
          opts%rt_ir%solar_sea_brdf_model > max_solar_sea_brdf_model) THEN
        err = errorstatus_fatal
        msg = 'Error: invalid specified solar sea BRDF model'
        THROWM(err.NE.0, msg)
      ENDIF

      ! Rayleigh scattering options
      IF (opts%rt_ir%rayleigh_max_wavelength > ray_max_wlm) THEN
        err = errorstatus_fatal
        msg = 'Error: invalid rayleigh_max_wavelength (too large)'
        THROWM(err.NE.0, msg)
      ENDIF
    ENDIF

    ! For PC-RTTOV the sea surface emissivity model is automatically chosen
    IF (.NOT. opts%rt_ir%pc%addpc) THEN

      ! Sea emissivity model
      IF (opts%rt_ir%ir_sea_emis_model < 1 .OR. &
          opts%rt_ir%ir_sea_emis_model > max_ir_sea_emis_model) THEN
        err = errorstatus_fatal
        msg = 'Error: invalid specified IR sea emissivity model'
        THROWM(err.NE.0, msg)
      ENDIF

      IF (opts%rt_ir%ir_sea_emis_model == 1 .AND. coefs%coef%ssirem_ver == 0) THEN
        err = errorstatus_fatal
        msg = 'Coefficient file does not support ISEM (ir_sea_emis_model 1)'
        THROWM(err.NE.0, msg)
      ENDIF

      IF (opts%rt_ir%ir_sea_emis_model == 2 .AND. coefs%coef%iremis_version == 0) THEN
        err = errorstatus_fatal
        msg = 'Coefficient file does not support ir_sea_emis_model 2'
        THROWM(err.NE.0, msg)
      ENDIF

    ENDIF

    ! Check cloud/aerosol setup
    IF (opts%rt_ir%addclouds) THEN
      IF (opts%rt_ir%cldcol_threshold >= 1._jprb) THEN
        err = errorstatus_fatal
        msg = 'Error in specified cldcol_threshold, must be less than 1.'
        THROWM(err.NE.0, msg)
      ENDIF

      IF (opts%rt_ir%cloud_overlap < 1 .OR. opts%rt_ir%cloud_overlap > ncloud_overlap) THEN
        err = errorstatus_fatal
        msg = 'Error: invalid specified cloud_overlap scheme.'
        THROWM(err.NE.0, msg)
      ENDIF

      IF (.NOT. ASSOCIATED(profiles(1)%cfrac)) THEN
        err = errorstatus_fatal
        msg = 'Profiles structure not allocated for clouds; '// &
              'opts%rt_ir%addclouds must be true when calling rttov_alloc_prof'
        THROWM(err.NE.0, msg)
      ENDIF

      IF (opts%rt_ir%user_cld_opt_param .AND. opts%rt_ir%cloud_overlap == cloud_overlap_simple) THEN
        err = errorstatus_fatal
        msg = 'Simple cloud overlap scheme cannot be used with explicit optical properties'
        THROWM(err.NE.0, msg)
      ENDIF
    ENDIF

    IF (opts%rt_ir%addaerosl) THEN
      IF (.NOT. opts%rt_ir%user_aer_opt_param) THEN
        IF (.NOT. ASSOCIATED(profiles(1)%aerosols)) THEN
          err = errorstatus_fatal
          msg = 'Profiles structure not allocated for aerosols; '// &
                'opts%rt_ir%addaerosl must be true when calling rttov_alloc_prof'
          THROWM(err.NE.0, msg)
        ENDIF
      ENDIF
    ENDIF

    ! Check scattering model
    IF (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) THEN
      IF (opts%rt_ir%ir_scatt_model < 1 .OR. opts%rt_ir%ir_scatt_model > max_ir_scatt_model) THEN
        err = errorstatus_fatal
        msg = 'Error in specified IR scattering model'
        THROWM(err.NE.0, msg)
      ENDIF
      IF (opts%rt_ir%addsolar) THEN
        IF (opts%rt_ir%vis_scatt_model < 1 .OR. opts%rt_ir%vis_scatt_model > max_vis_scatt_model) THEN
          err = errorstatus_fatal
          msg = 'Error in specified VIS/NIR scattering model'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%vis_scatt_model == vis_scatt_dom .AND. opts%rt_ir%dom_rayleigh) THEN
          IF (coefs%coef%fmv_model_ver < 13) THEN
            err = errorstatus_fatal
            msg = 'Error: Rayleigh DOM calculation can only be used with v13 predictors'
            THROWM(err.NE.0, msg)
          ENDIF
        ENDIF
      ENDIF

      ! Checks for MFASIS
      IF (opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_mfasis) THEN
        IF (opts%rt_ir%addclouds .AND. opts%rt_ir%addaerosl) THEN
          err = errorstatus_fatal
          msg = 'MFASIS cannot be used with both cloud AND aerosol'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addclouds .AND. coefs%coef_mfasis_cld%file_type /= mfasis_cld) THEN
          err = errorstatus_fatal
          msg = 'MFASIS cloud LUT not loaded'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addaerosl .AND. coefs%coef_mfasis_aer%file_type /= mfasis_aer) THEN
          err = errorstatus_fatal
          msg = 'MFASIS aerosol LUT not loaded'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addclouds .AND. opts%rt_ir%user_cld_opt_param) THEN
          err = errorstatus_fatal
          msg = 'MFASIS is not compatible with explicit cloud optical property input'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addaerosl .AND. opts%rt_ir%user_aer_opt_param) THEN
          err = errorstatus_fatal
          msg = 'MFASIS is not compatible with explicit aerosol optical property input'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (ALL(coefs%coef%ss_val_chn == 0)) THEN
          err = errorstatus_fatal
          msg = 'MFASIS requires a solar-enabled optical depth coefficient file'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addclouds) THEN
          ! Check consistency of cloud liquid and ice water schemes in LUT and profiles
          IF (ANY(profiles(:)%clw_scheme /= coefs%coef_mfasis_cld%clw_scheme)) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0, 'Different clw_scheme in MFASIS LUT and profiles')
          ENDIF
          IF (ANY(profiles(:)%ice_scheme /= coefs%coef_mfasis_cld%ice_scheme)) THEN
            err = errorstatus_fatal
            THROWM(err.NE.0, 'Different ice_scheme in MFASIS LUT and profiles')
          ENDIF
        ENDIF
      ENDIF

      ! Check scattering coefficient file supports selected optical properties
      IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
        IF (ANY(profiles(:)%clw_scheme == clw_scheme_opac) .AND. &
            coefs%coef_scatt%optp_wcl_opac%nchan == 0) THEN
          err = errorstatus_fatal
          msg = 'OPAC clw_scheme optical properties are not present in sccld file'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (ANY(profiles(:)%clw_scheme == clw_scheme_deff) .AND. &
            coefs%coef_scatt%optp_wcl_deff%nchan == 0) THEN
          err = errorstatus_fatal
          msg = 'Deff clw_scheme optical properties are not present in sccld file'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (ANY(profiles(:)%ice_scheme == ice_scheme_baum) .AND. &
            coefs%coef_scatt%optp_icl_baum%nchan == 0) THEN
          err = errorstatus_fatal
          msg = 'Baum/SSEC ice_scheme optical properties are not present in sccld file'
          THROWM(err.NE.0, msg)
        ENDIF
      ENDIF

      ! Check DOM streams vs Legendre coefficients
      IF ((opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                                     opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        IF (opts%rt_ir%dom_nstreams < dom_min_nstr) THEN
          err = errorstatus_fatal
          msg = 'DOM nstreams value too small'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addaerosl) THEN
          IF (opts%rt_ir%user_aer_opt_param) THEN
            IF (opts%rt_ir%dom_nstreams > aer_opt_param%nmom) THEN
              err = errorstatus_fatal
              msg = 'DOM nstreams exceeds the number of Legendre coefficients in the input aerosol optical parameters'
              THROWM(err.NE.0, msg)
            ENDIF
          ELSE
            IF (opts%rt_ir%dom_nstreams > coefs%coef_scatt%optp_aer%maxnmom) THEN
              err = errorstatus_fatal
              msg = 'DOM nstreams exceeds the number of Legendre coefficients in the aerosol optical parameter file'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF
        ENDIF
        IF (opts%rt_ir%addclouds) THEN
          IF (opts%rt_ir%user_cld_opt_param) THEN
            IF (opts%rt_ir%dom_nstreams > cld_opt_param%nmom) THEN
              err = errorstatus_fatal
              msg = 'DOM nstreams exceeds the number of Legendre coefficients in the input cloud optical parameters'
              THROWM(err.NE.0, msg)
            ENDIF
          ELSE
            IF (ANY(profiles(:)%clw_scheme == clw_scheme_opac)) THEN
              IF (opts%rt_ir%dom_nstreams > coefs%coef_scatt%optp_wcl_opac%maxnmom) THEN
                err = errorstatus_fatal
                msg = 'DOM nstreams exceeds the number of Legendre coefficients '// &
                      'for OPAC CLW data in the cloud optical parameter file'
                THROWM(err.NE.0, msg)
              ENDIF
            ENDIF
            IF (ANY(profiles(:)%clw_scheme == clw_scheme_deff)) THEN
              IF (opts%rt_ir%dom_nstreams > coefs%coef_scatt%optp_wcl_deff%maxnmom) THEN
                err = errorstatus_fatal
                msg = 'DOM nstreams exceeds the number of Legendre coefficients '// &
                      'for Deff clw_scheme data in the cloud optical parameter file'
                THROWM(err.NE.0, msg)
              ENDIF
            ENDIF
            IF (ANY(profiles(:)%ice_scheme == ice_scheme_baum)) THEN
              IF (opts%rt_ir%dom_nstreams > coefs%coef_scatt%optp_icl_baum%maxnmom) THEN
                err = errorstatus_fatal
                msg = 'DOM nstreams exceeds the number of Legendre coefficients '// &
                      'for Baum/SSEC ice cloud data in the cloud optical parameter file'
                THROWM(err.NE.0, msg)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF ! cloud or aerosol

    ! Check phase functions exist for solar simulations and opt_param init routine called
    IF (opts%rt_ir%addsolar .AND. ANY(coefs%coef%ss_val_chn > 0)) THEN
      IF (opts%rt_ir%addaerosl) THEN
        IF (opts%rt_ir%user_aer_opt_param) THEN
          IF (SIZE(aer_opt_param%phangle) < 2) THEN
            err = errorstatus_fatal
            msg = 'Input aerosol optical parameters do not contain phase functions'
            THROWM(err.NE.0, msg)
          ENDIF
          IF (.NOT. ASSOCIATED(aer_opt_param%phasefn_int%cosphangle)) THEN
            err = errorstatus_fatal
            msg = 'For solar scattering simulations with explicit optical properties you must call '// &
                  'rttov_init_opt_param for the aerosol optical parameter structure'
            THROWM(err.NE.0, msg)
          ENDIF
        ELSE
          IF (.NOT. ASSOCIATED(coefs%coef_scatt%optp_aer%chan_pha_index)) THEN
            err = errorstatus_fatal
            msg = 'Aerosol optical parameter file does not support solar simulations'
            THROWM(err.NE.0, msg)
          ENDIF
          DO i = 1, SIZE(chanprof)
            chan = chanprof(i)%chan
            IF (coefs%coef%ss_val_chn(chan) > 0 .AND. &
              coefs%coef_scatt%optp_aer%chan_pha_index(chan) == 0) THEN
              err = errorstatus_fatal
              msg = 'Aerosol phase function not present for selected solar-affected channel(s)'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        IF (opts%rt_ir%user_cld_opt_param) THEN
          IF (SIZE(cld_opt_param%phangle) < 2) THEN
            err = errorstatus_fatal
            msg = 'Input cloud optical parameters do not contain phase functions'
            THROWM(err.NE.0, msg)
          ENDIF
          IF (.NOT. ASSOCIATED(cld_opt_param%phasefn_int%cosphangle)) THEN
            err = errorstatus_fatal
            msg = 'For solar scattering simulations with explicit optical properties you must call '// &
                  'rttov_init_opt_param for the cloud optical parameter structure'
            THROWM(err.NE.0, msg)
          ENDIF
        ELSE
          IF (ANY(profiles(:)%clw_scheme == clw_scheme_opac)) THEN
            IF (.NOT. ASSOCIATED(coefs%coef_scatt%optp_wcl_opac%chan_pha_index)) THEN
              err = errorstatus_fatal
              msg = 'Cloud optical parameter file does not support solar simulations for OPAC CLW data'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF
          IF (ANY(profiles(:)%clw_scheme == clw_scheme_deff)) THEN
            IF (.NOT. ASSOCIATED(coefs%coef_scatt%optp_wcl_deff%chan_pha_index)) THEN
              err = errorstatus_fatal
              msg = 'Cloud optical parameter file does not support solar simulations for Deff clw_scheme data'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF
          IF (ANY(profiles(:)%ice_scheme == ice_scheme_baum)) THEN
            IF (.NOT. ASSOCIATED(coefs%coef_scatt%optp_icl_baum%chan_pha_index)) THEN
              err = errorstatus_fatal
              msg = 'Cloud optical parameter file does not support solar simulations for Baum/SSEC ice cloud data'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF

          DO i = 1, SIZE(chanprof)
            chan = chanprof(i)%chan
            IF (coefs%coef%ss_val_chn(chan) > 0) THEN
              IF (ANY(profiles(:)%clw_scheme == clw_scheme_opac)) THEN
                IF (coefs%coef_scatt%optp_wcl_opac%chan_pha_index(chan) == 0) THEN
                  err = errorstatus_fatal
                  msg = 'OPAC water cloud phase function not present for selected solar-affected channel(s)'
                  THROWM(err.NE.0, msg)
                ENDIF
              ENDIF
              IF (ANY(profiles(:)%clw_scheme == clw_scheme_deff)) THEN
                IF (coefs%coef_scatt%optp_wcl_deff%chan_pha_index(chan) == 0) THEN
                  err = errorstatus_fatal
                  msg = 'Deff clw_scheme phase function not present for selected solar-affected channel(s)'
                  THROWM(err.NE.0, msg)
                ENDIF
              ENDIF
              IF (ANY(profiles(:)%ice_scheme == ice_scheme_baum)) THEN
                IF (coefs%coef_scatt%optp_icl_baum%chan_pha_index(chan) == 0) THEN
                  err = errorstatus_fatal
                  msg = 'Ice cloud phase function not present for selected solar-affected channel(s)'
                  THROWM(err.NE.0, msg)
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDIF ! solar

    ! NLTE
    IF (opts%rt_ir%do_nlte_correction .AND. coefs%coef%nltecoef) THEN
      ! Check for v11 NLTE coefs instead of v12
      IF (coefs%coef%nlte_coef%nsat > 1) THEN
        err = errorstatus_fatal
        msg = 'Coefficient file appears to contain v11 NLTE coefficients which are incompatible with v12'
        THROWM(err.NE.0, msg)
      ENDIF
    ENDIF

  ENDIF ! IR instrument


  ! PMC-shift coefficients
  IF (coefs%coef%pmc_shift .AND. .NOT. opts%rt_all%co2_data) THEN
    err = errorstatus_fatal
    msg = 'PMC shift coefficients require input CO2 profile (co2_data true)'
    THROWM(err.NE.0, msg)
  ENDIF


  ! MW instruments
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
    IF (opts%rt_ir%addclouds) THEN
      err = errorstatus_fatal
      msg = 'Error: addclouds not applicable to MW instruments, use rttov_scatt interface instead'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl) THEN
      err = errorstatus_fatal
      msg = 'Error: aerosol calculations not applicable to MW instruments'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_mw%fastem_version < 0 .OR. &
        opts%rt_mw%fastem_version > max_fastem_version) THEN
      err = errorstatus_fatal
      msg = 'Error in specified FASTEM version'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_mw%clw_data) THEN
      IF (opts%rt_mw%clw_scheme < 0 .OR. &
          opts%rt_mw%clw_scheme > max_mw_clw_scheme) THEN
        err = errorstatus_fatal
        msg = 'Error in specified MW CLW scheme'
        THROWM(err.NE.0, msg)
      ENDIF
    ENDIF
  ENDIF


  ! PC-RTTOV
  IF (opts%rt_ir%pc%addpc) THEN
    IF (coefs%coef%id_comp_pc == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: optical depth coefficient file is not compatible with PCs'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%do_nlte_correction .AND. coefs%coef_pccomp%fmv_pc_nlte == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: do_nlte_correction is TRUE, but PC coef file is not compatible with NLTE simulations'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addclouds .AND. coefs%coef_pccomp%fmv_pc_cld == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: addclouds is TRUE, but PC coef file is not compatible with cloudy simulations'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl .AND. coefs%coef_pccomp%fmv_pc_aer == 0_jpim) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: addaerosl is TRUE, but PC coef file is not compatible with aerosol simulations'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl .AND. coefs%coef_scatt%optp_aer%id /= aer_id_opac) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: must be used with an RTTOV OPAC aerosol optical property file'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%pc%npcscores < 1 .OR. opts%rt_ir%pc%npcscores > coefs%coef_pccomp%fmv_pc_mnum) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: opts%rt_ir%pc%npcscores must lie between 1 and the maximum for the PC coef file'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%pc%ipcbnd < 1 .OR. opts%rt_ir%pc%ipcbnd > coefs%coef_pccomp%fmv_pc_bands) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: invalid spectral band (opts%rt_ir%pc%ipcbnd)'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%pc%ipcreg < 1 .OR. opts%rt_ir%pc%ipcreg > coefs%coef_pccomp%fmv_pc_sets(opts%rt_ir%pc%ipcbnd)) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: invalid predictor set (opts%rt_ir%pc%ipcreg) for selected spectral band (opts%rt_ir%pc%ipcbnd)'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl .AND. opts%rt_ir%user_aer_opt_param) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with user input aerosol optical parameters'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addclouds .AND. opts%rt_ir%user_cld_opt_param) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with user input cloud optical parameters'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
      IF (opts%rt_ir%ir_scatt_model /= ir_scatt_chou) THEN
        err = errorstatus_fatal
        msg = 'PC-RTTOV: Chou-scaling must be selected for the ir_scatt_model option'
        THROWM(err.NE.0, msg)
      ENDIF
    ENDIF

    IF (opts%rt_ir%addsolar) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with addsolar TRUE'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_all%do_lambertian) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with do_lambertian TRUE'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_all%plane_parallel) THEN
      err = errorstatus_fatal
      msg = 'PC-RTTOV: cannot be run with plane_parallel TRUE'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_OPTIONS',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_CHECK_OPTIONS',1_jpim,ZHOOK_HANDLE)
END SUBROUTINE rttov_check_options
