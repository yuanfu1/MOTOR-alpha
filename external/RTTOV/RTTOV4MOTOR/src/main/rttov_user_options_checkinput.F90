! Description:
!> @file
!!   Check consistency of selected options and coefficient file.
!
!> @brief
!!   Check consistency of selected options and coefficient file,
!!   useful for "debugging" simulations.
!!
!!   This subroutine reports illegal option settings, and also dubious
!!   options settings which will not result in bad simulations, but may
!!   indicate that the simulation has not been configured as intended.
!!
!! @param[out]    err                  status on exit
!! @param[in]     opts                 options to configure the simulations
!! @param[in]     coefs                coefficients structure for instrument to simulate
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
SUBROUTINE rttov_user_options_checkinput(err, opts, coefs)

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coefs, rttov_options
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
      mfasis_cld,               &
      mfasis_aer,               &
      aer_id_opac,              &
      ncloud_overlap,           &
      cloud_overlap_simple,     &
      ray_max_wlm
!INTF_ON

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT) :: err
  TYPE(rttov_options), INTENT(IN)  :: opts
  TYPE(rttov_coefs),   INTENT(IN)  :: coefs
!INTF_END

#include "rttov_errorreport.interface"

  CHARACTER(LEN=256) :: msg
  REAL(KIND=jprb) :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  TRY

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',0_jpim,ZHOOK_HANDLE)

  IF (opts%interpolation%addinterp) THEN
    IF (opts%interpolation%interp_mode < 1 .OR. &
        opts%interpolation%interp_mode > ninterp_modes) THEN
      err = errorstatus_fatal
      msg = 'Error in specified interpolation mode'
      THROWM(err.NE.0, msg)
    ENDIF
  ELSE
    IF (opts%interpolation%lgradp) THEN
      err = errorstatus_fatal
      msg = 'Interpolation should be enabled if LGRADP is set to TRUE'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF

  IF (opts%rt_all%ozone_data .AND. coefs%coef%nozone == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable O3'
    THROWM(err.NE.0, msg)
  ENDIF

  IF (opts%rt_all%co2_data .AND. coefs%coef%nco2 == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CO2'
    THROWM(err.NE.0, msg)
  ENDIF

  IF (opts%rt_all%co_data .AND. coefs%coef%nco == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CO'
    THROWM(err.NE.0, msg)
  ENDIF

  IF (opts%rt_all%n2o_data .AND. coefs%coef%nn2o == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable N2O'
    THROWM(err.NE.0, msg)
  ENDIF

  IF (opts%rt_all%ch4_data .AND. coefs%coef%nch4 == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable CH4'
    THROWM(err.NE.0, msg)
  ENDIF

  IF (opts%rt_all%so2_data .AND. coefs%coef%nso2 == 0_jpim) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not allow variable SO2'
    THROWM(err.NE.0, msg)
  ENDIF


  IF (opts%rt_ir%addsolar .AND. coefs%coef%fmv_model_ver < 9) THEN
    err = errorstatus_fatal
    msg = 'Coefficient file does not support solar calculations'
    THROWM(err.NE.0, msg)
  ENDIF


  ! All IR instruments
  IF (coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) THEN
    IF (opts%rt_mw%clw_data) THEN
      err = errorstatus_fatal
      msg = 'Cloud liquid water only applicable to MW instruments'
      THROWM(err.NE.0, msg)
    ENDIF

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

      IF (opts%rt_ir%user_cld_opt_param .AND. opts%rt_ir%cloud_overlap == cloud_overlap_simple) THEN
        err = errorstatus_fatal
        msg = 'Simple cloud overlap scheme cannot be used with explicit optical properties'
        THROWM(err.NE.0, msg)
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
      ELSEIF (opts%rt_ir%vis_scatt_model == vis_scatt_mfasis) THEN
        err = errorstatus_fatal
        msg = 'MFASIS requires addsolar to be true'
        THROWM(err.NE.0, msg)
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
      ENDIF

      ! Check DOM streams vs Legendre coefficients
      IF ((opts%rt_ir%addsolar .AND. opts%rt_ir%vis_scatt_model == vis_scatt_dom) .OR. &
                                     opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
        IF (opts%rt_ir%dom_nstreams < dom_min_nstr) THEN
          err = errorstatus_fatal
          msg = 'DOM nstreams value too small'
          THROWM(err.NE.0, msg)
        ENDIF
        IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
          IF (opts%rt_ir%dom_nstreams > coefs%coef_scatt%optp_aer%maxnmom) THEN
            err = errorstatus_fatal
            msg = 'DOM nstreams exceeds the number of Legendre coefficients in the aerosol optical parameter file'
            THROWM(err.NE.0, msg)
          ENDIF
        ENDIF
        IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
          ! We don't know which optical properties the user will select in profiles
          ! and the sccld file may not contain properties for all cloud types so we
          ! can't be more specific than checking the largest maxnmom among cloud types
          IF (opts%rt_ir%dom_nstreams > MAX(coefs%coef_scatt%optp_wcl_opac%maxnmom, &
                                            coefs%coef_scatt%optp_wcl_deff%maxnmom, &
                                            coefs%coef_scatt%optp_icl_baum%maxnmom)) THEN
            err = errorstatus_fatal
            msg = 'DOM nstreams exceeds the number of Legendre coefficients in the cloud optical parameter file'
            THROWM(err.NE.0, msg)
          ENDIF
        ENDIF
      ENDIF

      IF (opts%rt_ir%addsolar .AND. ANY(coefs%coef%ss_val_chn > 0)) THEN
        IF (opts%rt_ir%addaerosl) THEN
          IF (.NOT. opts%rt_ir%user_aer_opt_param) THEN
            IF (.NOT. ASSOCIATED(coefs%coef_scatt%optp_aer%chan_pha_index)) THEN
              err = errorstatus_fatal
              msg = 'Aerosol optical parameter file does not support solar simulations'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF
        ENDIF
        IF (opts%rt_ir%addclouds) THEN
          IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN
            ! We don't know which optical properties the user will select in profiles
            ! and the sccld file may not contain properties for all cloud types so we
            ! can't be more specific than checking for any phase functions in the file
            IF (.NOT. (ASSOCIATED(coefs%coef_scatt%optp_wcl_opac%chan_pha_index) .OR. &
                       ASSOCIATED(coefs%coef_scatt%optp_wcl_deff%chan_pha_index) .OR. &
                       ASSOCIATED(coefs%coef_scatt%optp_icl_baum%chan_pha_index))) THEN
              err = errorstatus_fatal
              msg = 'Cloud optical parameter file does not support solar simulations'
              THROWM(err.NE.0, msg)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF ! clouds or aerosol
  ENDIF ! IR instrument


  ! NLTE
  IF (opts%rt_ir%do_nlte_correction) THEN
    IF (.NOT. coefs%coef%nltecoef) THEN
      err = errorstatus_fatal
      msg = 'Coefficient file does not support NLTE correction or channel selection does not include any NLTE-affected channels'
      THROWM(err.NE.0, msg)
    ENDIF

    ! Check for v11 NLTE coefs instead of v12
    IF (coefs%coef%nlte_coef%nsat > 1) THEN
      err = errorstatus_fatal
      msg = 'Coefficient file appears to contain v11 NLTE coefficients which are incompatible with v12'
      THROWM(err.NE.0, msg)
    ENDIF
  ENDIF


  ! PMC-shift coefficients
  IF (coefs%coef%pmc_shift .AND. .NOT. opts%rt_all%co2_data) THEN
    err = errorstatus_fatal
    msg = 'PMC shift coefficients require input CO2 profile (co2_data true)'
    THROWM(err.NE.0, msg)
  ENDIF


  ! MW instruments
  IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
    IF (opts%rt_ir%addsolar) THEN
      err = errorstatus_fatal
      msg = 'Solar calculations not applicable to MW instruments'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addclouds) THEN
      err = errorstatus_fatal
      msg = 'addclouds not applicable to MW instruments: use rttov_scatt interface instead'
      THROWM(err.NE.0, msg)
    ENDIF

    IF (opts%rt_ir%addaerosl) THEN
      err = errorstatus_fatal
      msg = 'Aerosol calculations not applicable to MW instruments'
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

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',1_jpim,ZHOOK_HANDLE)

  CATCH

  IF (LHOOK) CALL DR_HOOK('RTTOV_USER_OPTIONS_CHECKINPUT',1_jpim,ZHOOK_HANDLE)
END SUBROUTINE rttov_user_options_checkinput
