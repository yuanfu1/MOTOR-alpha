! Description:
!> @file
!!   Runs RTTOV direct model
!
!> @brief
!!   Runs RTTOV direct model
!!
!! @details
!!   Computes multi-channel level to space transmittances,
!!   top of atmosphere and level to space radiances, brightness
!!   temperatures and reflectances, and optionally surface
!!   emissivities and BRDFs for multiple profiles in a single call,
!!   for nadir-viewing satellite-based visible, infrared or microwave
!!   sensors. Requires a coefficient file for the sensor for which
!!   simulated radiances are requested.
!!
!!   The methodology is described in the following:
!!
!!   Eyre J.R. and H.M. Woolf  1988 Transmittance of atmospheric gases
!!   in the microwave region: a fast model. Applied Optics 27  3244-3249
!!
!!   Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
!!   systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
!!   librarian at ECMWF).
!!
!!   Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
!!   Transfer Model for Assimilation of Satellite Radiance Observations.
!!   QJRMS, 125, 1407-1425.
!!
!!   Matricardi, M., F. Chevallier and S. Tjemkes 2001 An improved general
!!   fast radiative transfer model for the assimilation of radiance
!!   observations. ECMWF Research Dept. Tech. Memo. 345
!!   (available from the librarian at ECMWF).
!!
!!   Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
!!   transfer model for the infrared atmospheric sounding interferometer.
!!   ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
!!
!!   Rochon Y.J., L. Garand, D.S. Turner and S. Polavarapu Jacobian mapping
!!   between  vertical coordinate systems in data assimilation (submitted QJRMS
!!   June 2006)
!!
!!   Matricardi, M. 2009: An Observation operator for the assimilation of
!!   principal component scores into a NWP system. Available from EUMETSAT
!!
!! @param[out]    errorstatus    status on exit
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     opts           options to configure the simulations
!! @param[in]     profiles       input atmospheric profiles and surface variables
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in,out] transmission   output transmittances
!! @param[in,out] radiance       output radiances and corresponding BTs and BRFs
!! @param[in,out] radiance2      secondary output radiances, optional
!! @param[in]     calcemis       flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity     input/output surface emissivities, optional
!! @param[in]     calcrefl       flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance    input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in]     aer_opt_param  input aerosol optical parameters, optional
!! @param[in]     cld_opt_param  input cloud optical parameters, optional
!! @param[in,out] traj           RTTOV internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp         output PC scores and radiances from PC-RTTOV, optional
!! @param[in]     channels_rec   list of channels for which to calculate reconstructed radiances, optional
!! @param[in,out] traj_sta       RTTOV internal state, optional, not intended for general use
!! @param[in,out] traj_dyn       RTTOV internal state, optional, not intended for general use
!! @param[in,out] lbl_check      used for coef verification, optional, not intended for general use
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
SUBROUTINE rttov_direct( &
              errorstatus,    &
              chanprof,       &
              opts,           &
              profiles,       &
              coefs,          &
              transmission,   &
              radiance,       &
              radiance2,      &
              calcemis,       &
              emissivity,     &
              calcrefl,       &
              reflectance,    &
              aer_opt_param,  &
              cld_opt_param,  &
              traj,           &
              traj_dyn,       &
              traj_sta,       &
              pccomp,         &
              channels_rec,   &
              lbl_check)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :    &
         rttov_coefs,        &
         rttov_pccomp,       &
         rttov_profile,      &
         rttov_transmission, &
         rttov_radiance,     &
         rttov_radiance2,    &
         rttov_options,      &
         rttov_chanprof,     &
         rttov_emissivity,   &
         rttov_reflectance,  &
         rttov_opt_param,    &
         rttov_traj,         &
         rttov_traj_dyn,     &
         rttov_traj_sta,     &
         rttov_lbl_check
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
         sensor_id_mw,     &
         sensor_id_hi,     &
         sensor_id_po,     &
         gas_unit_ppmvdry, &
         vis_scatt_dom,    &
         ir_scatt_dom,     &
         pi_r
  USE rttov_htfrtc_interface_mod, ONLY : htfrtc_interface
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),       INTENT(OUT)                     :: errorstatus     ! return flag
  TYPE(rttov_chanprof),     INTENT(IN)                      :: chanprof(:)     ! Profile/channel indices (nchanprof)
  TYPE(rttov_options),      INTENT(IN)                      :: opts
  TYPE(rttov_profile),      INTENT(IN)                      :: profiles(:)     ! Atmospheric profiles supplied on
                                                                               ! user levels (nprofiles)
  TYPE(rttov_coefs),        INTENT(IN)   , TARGET           :: coefs           ! TARGET attribute required
  TYPE(rttov_transmission), INTENT(INOUT)                   :: transmission    ! transmittances
  TYPE(rttov_radiance),     INTENT(INOUT)                   :: radiance        ! radiances (mw/cm-1/ster/sq.m),
                                                                               ! BTs (degK) and reflectances (BRF)
  TYPE(rttov_radiance2),    INTENT(INOUT), OPTIONAL         :: radiance2
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))    ! switches for emis calcs
  TYPE(rttov_emissivity),   INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))  ! surface emis
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))    ! switches for refl calcs
  TYPE(rttov_reflectance),  INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof)) ! surface refl
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_traj),         INTENT(INOUT), OPTIONAL, TARGET :: traj            ! TARGET attribute required (see rttov_check_traj)
  TYPE(rttov_traj_dyn),     INTENT(INOUT), OPTIONAL, TARGET :: traj_dyn
  TYPE(rttov_traj_sta),     INTENT(INOUT), OPTIONAL, TARGET :: traj_sta
  TYPE(rttov_pccomp),       INTENT(INOUT), OPTIONAL         :: pccomp
  INTEGER(KIND=jpim),       INTENT(IN)   , OPTIONAL         :: channels_rec(:)
  TYPE(rttov_lbl_check),    INTENT(IN)   , OPTIONAL         :: lbl_check
!INTF_END
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_apply_pc_aer_reg_lims.interface"
#include "rttov_apply_reg_limits.interface"
#include "rttov_calc_nearest_lev.interface"
#include "rttov_calcbt.interface"
#include "rttov_calcbt_pc.interface"
#include "rttov_calcemis_ir.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_calcsatrefl.interface"
#include "rttov_calcsurfrefl.interface"
#include "rttov_check_options.interface"
#include "rttov_check_profiles.interface"
#include "rttov_check_reg_limits.interface"
#include "rttov_check_traj.interface"
#include "rttov_checkpcchan.interface"
#include "rttov_cloud_overlap.interface"
#include "rttov_convert_profile_units.interface"
#include "rttov_copy_aux_prof.interface"
#include "rttov_copy_opdp_path.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_raytracing.interface"
#include "rttov_dom.interface"
#include "rttov_dom_setup_profile.interface"
#include "rttov_errorreport.interface"
#include "rttov_fresnel.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_intavg_chan.interface"
#include "rttov_intavg_prof.interface"
#include "rttov_integrate.interface"
#include "rttov_locpat.interface"
#include "rttov_mfasis.interface"
#include "rttov_mw_clw_absorption.interface"
#include "rttov_nlte_bias_correction.interface"
#include "rttov_opdep_13.interface"
#include "rttov_opdep_9.interface"
#include "rttov_opdep_78.interface"
#include "rttov_opdpscattir.interface"
#include "rttov_pcscores.interface"
#include "rttov_predictor_precalc_13.interface"
#include "rttov_predictor_precalc_789.interface"
#include "rttov_profaux_cldaer.interface"
#include "rttov_rayleigh_extinction.interface"
#include "rttov_reconstruct.interface"
#include "rttov_refsun.interface"
#include "rttov_setgeometry.interface"
#include "rttov_setpredictors_13.interface"
#include "rttov_setpredictors_789.interface"
#include "rttov_transmit.interface"
#include "rttov_transmit_solar.interface"

  INTEGER(KIND=jpim) :: err
  INTEGER(KIND=jpim) :: nlevels, nprofiles, nchanprof, npcscores
  LOGICAL(KIND=jplm) :: sensor_mw, sensor_ir

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj_dyn), TARGET  :: traj1_dyn
  TYPE(rttov_traj_dyn), POINTER :: traj0_dyn
  TYPE(rttov_traj_sta), TARGET  :: traj1_sta
  TYPE(rttov_traj_sta), POINTER :: traj0_sta
  LOGICAL(KIND=jplm) :: ltraj_dyn_dealloc
  LOGICAL(KIND=jplm) :: ltraj_sta_dealloc

  REAL(KIND=jprb) :: ZHOOK_HANDLE
!- End of header ------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 0_jpim, ZHOOK_HANDLE)

  errorstatus = errorstatus_success

  IF (opts%htfrtc_opts%htfrtc) THEN
    CALL htfrtc_interface(err, coefs, opts, profiles, pccomp, calcemis, emissivity)
    errorstatus = err
    THROWM(err.NE.0, 'Error in HTFRTC module')
    ! DAR: Since HTFRTC isn't RTTOV, skip to the end of rttov_direct and leave
    GOTO 998
  ENDIF

!-------------
! Initialize
!-------------
  nprofiles = SIZE(profiles)
  nchanprof = SIZE(chanprof)
  nlevels   = profiles(1)%nlevels
  sensor_mw = coefs%coef%id_sensor == sensor_id_mw .OR. &
              coefs%coef%id_sensor == sensor_id_po
  sensor_ir = .NOT. sensor_mw

  ltraj_dyn_dealloc = .FALSE.
  ltraj_sta_dealloc = .FALSE.
  NULLIFY (traj0, traj0_dyn, traj0_sta)

  CALL rttov_init_transmission(transmission)
  CALL rttov_init_rad(radiance, radiance2)

  CALL rttov_check_options(err, opts, coefs, chanprof, profiles, aer_opt_param, cld_opt_param)
  THROW(err.NE.0)

  IF (opts%rt_ir%pc%addpc) THEN
    IF (.NOT. PRESENT(pccomp)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "pccomp argument is mandatory for PC-RTTOV")
    ENDIF
    npcscores = opts%rt_ir%pc%npcscores * nprofiles
    IF (npcscores > SIZE(pccomp%total_pcscores)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "pccomp structure is too small for specified npcscores")
    ENDIF

    CALL rttov_checkpcchan( &
         nprofiles,         &
         nchanprof,         &
         opts,              &
         chanprof,          &
         coefs,             &
         err                )
    THROWM(err.NE.0, "rttov_checkpcchan fatal error")
  ENDIF

  IF (PRESENT(traj_sta)) THEN
    traj0_sta => traj_sta
  ELSE
    traj0_sta => traj1_sta
  ENDIF
  ltraj_sta_dealloc = .NOT. PRESENT(traj_sta)

  ! Allocate/initialise the static trajectory structure
  CALL rttov_alloc_traj_sta(err, traj0_sta, opts, coefs, chanprof, profiles, &
                            1_jpim, npcscores, channels_rec, calcemis, lbl_check)
  THROW(err.NE.0)

  ! Check surface emissivity and reflectance inputs
  IF (traj0_sta%dothermal) THEN
    IF (.NOT. (PRESENT(calcemis) .AND. PRESENT(emissivity))) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "calcemis and emissivity parameters required")
    ENDIF
    IF (ANY(traj0_sta%thermal .AND. .NOT. calcemis .AND. emissivity%emis_in < 0._jprb)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "invalid input surface emissivity (must be >=0)")
    ENDIF
    IF (ANY(traj0_sta%thermal .AND. .NOT. calcemis .AND. emissivity%emis_in > 1._jprb)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "invalid input surface emissivity (must be <=1)")
    ENDIF
    IF (opts%rt_all%do_lambertian) THEN
      IF (ANY(emissivity%specularity < 0._jprb)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "invalid input surface specularity (must be >=0)")
      ENDIF
      IF (ANY(emissivity%specularity > 1._jprb)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "invalid input surface specularity (must be <=1)")
      ENDIF
    ENDIF
  ENDIF
  IF (traj0_sta%dosolar) THEN
    IF (.NOT. (PRESENT(calcrefl) .AND. PRESENT(reflectance))) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "calcrefl and reflectance parameters required")
    ENDIF
    IF (ANY(traj0_sta%solar .AND. .NOT. calcrefl .AND. reflectance%refl_in < 0._jprb)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "invalid input surface BRDF (must be >=0)")
    ENDIF
  ENDIF

  ! Check/allocate trajectory structure
  CALL rttov_check_traj( &
          err,           &
          nprofiles,     &
          nchanprof,     &
          opts,          &
          nlevels,       &
          coefs,         &
          1_jpim,        &
          traj0 = traj0, &
          traj1 = traj1, &
          traj2 = traj)
  THROWM(err.NE.0, "rttov_check_traj fatal error")

!-----------------------------------------------------------------------
! Convert input profiles to units used internally
!-----------------------------------------------------------------------
  CALL rttov_convert_profile_units(opts, coefs, profiles, traj0%profiles_int)

!-----------------------------------------------------------------------
! Check input profiles for unphysical values
!-----------------------------------------------------------------------
  IF (opts%config%do_checkinput) THEN
    CALL rttov_check_profiles(err, opts, coefs, traj0_sta%dothermal, traj0_sta%do_mfasis, &
                              chanprof, profiles, traj0%profiles_int, radiance, aer_opt_param, &
                              cld_opt_param)
    THROW(err.NE.0)
  ENDIF

!-----------------------------------------------------------------------
! Calculate zenith and view angle geometry (direct model only)
!-----------------------------------------------------------------------
  CALL rttov_setgeometry(traj0_sta%plane_parallel, profiles, coefs%coef, &
                         traj0_sta%angles)

!-----------------------------------------------------------------------
! Calculate near surface and cloud top levels
!-----------------------------------------------------------------------
  CALL rttov_calc_nearest_lev(sensor_ir, profiles, traj0%aux_prof%s)

!------------------------------------------------------------------------
! Local path geometry calculations
!------------------------------------------------------------------------
  CALL rttov_locpat( &
      opts,                     &
      traj0_sta%dosolar,        &
      traj0_sta%plane_parallel, &
      profiles,                 &
      traj0%profiles_int,       &
      traj0%aux_prof%s,         &
      traj0_sta%angles,         &
      traj0%raytracing,         &
      chanprof,                 &
      radiance%geometric_height)

  ! If no channels are active (thermal and solar flags false for all channels)
  ! then skip to the end. In particular this can happen when using the parallel
  ! interface. Do this after we have populated the geometric_height output to
  ! ensure consistency in outputs via the parallel interface.
  IF (.NOT. (traj0_sta%dothermal .OR. traj0_sta%dosolar)) GOTO 998

!-----------------------------------------------------------------------
! RTTOV optical depth calculation: profiles are interpolated onto
! coef levels, ancillary data, predictor variables and finally optical
! depths are calculated on coef levels, and optical depths are
! interpolated back to user levels.
!-----------------------------------------------------------------------

  IF (traj0_sta%do_opdep_calc) THEN
!-----------------------------------------------------------------------
! Interpolator first call - input profiles from user levs to coef levs
!-----------------------------------------------------------------------
    CALL rttov_init_prof(traj0%profiles_coef, p = coefs%coef%ref_prfl_p)

    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_prof( &
              opts,                &
              nlevels,             &
              coefs%coef%nlevels,  &
              profiles,            &
              traj0%profiles_int,  &
              traj0%profiles_coef, &
              coefs%coef,          &
              coefs%coef_pccomp)
    ELSE
      ! If interpolator is off copy profile variables
      CALL rttov_copy_prof( &
              traj0%profiles_coef,    &
              profiles,               &
              larray = .TRUE._jplm,   &
              lscalar = .FALSE._jplm, &
              profiles_gas = traj0%profiles_int)
    ENDIF

    ! Copy non-interpolated variables (e.g. surface parameters)
    CALL rttov_copy_prof( &
            traj0%profiles_coef,   &
            profiles,              &
            larray = .FALSE._jplm, &
            lscalar = .TRUE._jplm, &
            profiles_gas = traj0%profiles_int)

    traj0%profiles_coef(:)%gas_units = gas_unit_ppmvdry

!------------------------------------------------------------------
! Check profile data against regression limits - coef levels
!------------------------------------------------------------------
    IF (opts%config%do_checkinput) THEN
      CALL rttov_check_reg_limits( &
               err,                 &
               opts,                &
               chanprof,            &
               traj0_sta%thermal,   &
               traj0_sta%solar,     &
               traj0%profiles_coef, &
               profiles,            &
               coefs%coef,          &
               coefs%coef_pccomp,   &
               radiance)
      THROW(err .NE. 0)
    ENDIF

    IF (opts%config%apply_reg_limits .OR. opts%interpolation%addinterp) THEN
      CALL rttov_copy_prof(traj0_sta%profiles_coef_ref, traj0%profiles_coef)

      CALL rttov_apply_reg_limits( &
               opts,                &
               profiles,            &
               traj0%profiles_coef, &
               coefs%coef,          &
               coefs%coef_pccomp)
    ENDIF

!--------------------------------------------------------------------------------
! Calculate near surface level and local path geometry calculations - coef levels
!--------------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_calc_nearest_lev(.FALSE._jplm, traj0%profiles_coef, traj0%aux_prof_coef%s)
      CALL rttov_locpat( &
          opts,                     &
          traj0_sta%dosolar,        &
          traj0_sta%plane_parallel, &
          traj0%profiles_coef,      &
          traj0%profiles_coef,      &
          traj0%aux_prof_coef%s,    &
          traj0_sta%angles,         &
          traj0%raytracing_coef)
    ELSE
      CALL rttov_copy_aux_prof(traj0%aux_prof_coef%s, traj0%aux_prof%s)
      CALL rttov_copy_raytracing(traj0_sta%dosolar, traj0%raytracing_coef, traj0%raytracing)
    ENDIF

!---------------------------------------------------
! Calculate predictors - coef levels
!---------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      CALL rttov_predictor_precalc_13( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0%raytracing_coef)

      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_13( &
          opts,                       &
          traj0%profiles_coef,        &
          coefs%coef,                 &
          traj0%aux_prof_coef,        &
          traj0%predictors%path1,     &
          traj0%raytracing_coef,      &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_13( &
          opts,                       &
          traj0%profiles_coef,        &
          coefs%coef,                 &
          traj0%aux_prof_coef,        &
          traj0%predictors%path2,     &
          traj0%raytracing_coef,      &
          raypath = 2_jpim) 
      ENDIF
    ELSE
      CALL rttov_predictor_precalc_789( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0_sta%angles,              &
          traj0%raytracing_coef)

      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_789( &
          opts,                       &
          traj0%profiles_coef,        &
          coefs%coef,                 &
          traj0%aux_prof_coef,        &
          traj0%predictors%path1,     &
          traj0%raytracing_coef,      &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_789( &
          opts,                       &
          traj0%profiles_coef,        &
          coefs%coef,                 &
          traj0%aux_prof_coef,        &
          traj0%predictors%path2,     &
          traj0%raytracing_coef,      &
          raypath = 2_jpim) 
      ENDIF
    ENDIF

!--------------------------------------------------------------------
! Predict atmospheric (emissive) and solar optical depths - coef levs
!--------------------------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for all thermal channels
        CALL rttov_opdep_13( &
                coefs%coef%nlayers,             &
                chanprof,                       &
                traj0_sta%thermal,              &
                traj0%predictors%path1,         &
                coefs%coef,                     &
                coefs%coef%thermal,             &
                coefs%coef%thermal_corr,        &
                traj0%opdp_path_coef%atm_level, &
                traj0_sta%thermal_path1%opdp_ref)

      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for all solar channels
        CALL rttov_opdep_13( &
                coefs%coef%nlayers,                   &
                chanprof,                             &
                traj0_sta%solar,                      &
                traj0%predictors%path2,               &
                coefs%coef,                           &
                coefs%coef%solar,                     &
                coefs%coef%solar_corr,                &
                traj0%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref)
      ENDIF
    ELSEIF (coefs%coef%fmv_model_ver == 9) THEN
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for all thermal channels
        CALL rttov_opdep_9( &
                coefs%coef%nlayers,                    &
                chanprof,                              &
                traj0_sta%thermal,                     &
                traj0%predictors%path1,                &
                coefs%coef,                            &
                coefs%coef%thermal,                    &
                traj0%opdp_path_coef%atm_level,        &
                traj0_sta%thermal_path1%opdp_ref_coef)
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for all solar channels
        CALL rttov_opdep_9( &
                coefs%coef%nlayers,                   &
                chanprof,                             &
                traj0_sta%solar,                      &
                traj0%predictors%path2,               &
                coefs%coef,                           &
                coefs%coef%solar,                     &
                traj0%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref_coef)
      ENDIF
    ELSE
      CALL rttov_opdep_78( &
              coefs%coef%nlayers,        &
              chanprof,                  &
              traj0%predictors%path1,    &
              coefs%coef,                &
              coefs%coef%thermal,        &
              traj0%opdp_path_coef,      &
              traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF

!--------------------------------------------------------------------------
! Interpolator second  call - optical depths from coef levs to user levs
!--------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_chan( &
              opts,                           &
              traj0_sta%thermal,              &
              traj0_sta%solar,                &
              coefs%coef%nlevels,             &
              nlevels,                        &
              chanprof,                       &
              traj0%profiles_coef,            &
              profiles,                       &
              traj0%opdp_path_coef,           &
              traj0%opdp_path)
    ELSE
      CALL rttov_copy_opdp_path(opts, traj0%opdp_path, traj0%opdp_path_coef)
    ENDIF

!--------------------------------------------------------------------------
! Set top-level optical depth to zero if spacetop is true
!--------------------------------------------------------------------------
    IF (opts%interpolation%spacetop) THEN
      traj0%opdp_path%atm_level(1,:) = 0._jprb
      IF (traj0_sta%do_solar_opdep_calc) THEN
        traj0%opdp_path%sun_level_path2(1,:) = 0._jprb
      ENDIF
    ENDIF

  ENDIF ! do_opdep_calc

  ! This is for coefficient testing; we replace the optical depths
  ! with those from the line-by-line integrated over the ISRF
  IF (PRESENT(lbl_check)) THEN
    IF (ASSOCIATED(lbl_check%atm_layer)) THEN
      traj0%opdp_path%atm_level(1,:) = 0._jprb
      traj0%opdp_path%atm_level(2:,:) = lbl_check%atm_layer(:,:)
      IF (traj0_sta%dosolar) THEN
        traj0%opdp_path%sun_level_path2(1,:) = 0._jprb
        traj0%opdp_path%sun_level_path2(2:,:) = lbl_check%atm_layer_path2(:,:)
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------
! Calculations from this point onwards are on user levels
!------------------------------------------------------------------------

!--------------------------------------------------------------------------
! MW CLW absorption optical depths
!--------------------------------------------------------------------------
  IF (sensor_mw .AND. opts%rt_mw%clw_data) THEN
    CALL rttov_mw_clw_absorption( &
            opts,             &
            coefs%coef,       &
            chanprof,         &
            traj0%raytracing, &
            profiles,         &
            traj0%aux_prof,   &
            traj0%opdp_path)
  ENDIF

!--------------------------------------------------------------------------
! Rayleigh extinction optical depths
!--------------------------------------------------------------------------
  IF (traj0_sta%do_rayleigh_param) THEN
    CALL rttov_rayleigh_extinction( &
            opts,                            &
            chanprof,                        &
            traj0_sta%do_rayleigh_dom,       &
            traj0_sta%solar,                 &
            profiles,                        &
            coefs%coef,                      &
            traj0%raytracing,                &
            traj0%opdp_path%sun_level_path2, &
            traj0%transmission_scatt_ir)
  ENDIF

!------------------------------------------------------
! Calculate the number of cloud columns and the cloud distribution in each column
!------------------------------------------------------
  IF (PRESENT(traj_dyn)) THEN
    traj0_dyn => traj_dyn
  ELSE
    traj0_dyn => traj1_dyn
  ENDIF
  ltraj_dyn_dealloc = .NOT. PRESENT(traj_dyn)

  traj0_dyn%ncolumns = 0_jpim

  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cloud_overlap( &
            opts%rt_ir,         &
            profiles,           &
            traj0%profiles_int, &
            traj0%ircld,        &
            traj0_dyn%ncolumns)
  ELSE
    traj0%ircld%xcolclr = 1._jprb
    traj0%ircld%ncolumn = 0_jpim
    traj0%ircld%icldarr = 0_jpim
  ENDIF

  CALL rttov_alloc_traj_dyn(err, traj0_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                            traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                            traj0_sta%thermal, traj0_sta%solar, &
                            traj0_sta%dothermal, traj0_sta%do_mfasis, 1_jpim)
  THROW(err.NE.0)

!-----------------------------------------------------------------------
! Apply PC-RTTOV aerosol regression limits
!-----------------------------------------------------------------------
  IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addaerosl) THEN
    CALL rttov_apply_pc_aer_reg_lims( &
          opts,                 &
          chanprof,             &
          coefs%coef_pccomp,    &
          profiles,             &
          traj0%profiles_int,   &
          radiance,             &
          traj0_sta%pc_aer_ref, &
          traj0_sta%pc_aer_min, &
          traj0_sta%pc_aer_max)
  ENDIF

!----------------------------------------------------------------------------
! Calculate optical depths of aerosols and/or clouds
!----------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_profaux_cldaer( &
            opts,                &
            profiles,            &
            traj0%profiles_int,  &
            traj0%aux_prof)

    CALL rttov_opdpscattir( &
            profiles(1)%nlayers,                  &
            chanprof,                             &
            opts,                                 &
            traj0_sta%dom_nstreams,               &
            traj0%aux_prof,                       &
            traj0%ircld,                          &
            profiles,                             &
            traj0%profiles_int,                   &
            aer_opt_param,                        &
            cld_opt_param,                        &
            traj0_sta%dothermal,                  &
            traj0_sta%thermal,                    &
            traj0_sta%dosolar,                    &
            traj0_sta%solar,                      &
            traj0_sta%do_rayleigh_dom,            &
            coefs%coef,                           &
            coefs%coef_scatt,                     &
            coefs%coef_mfasis_cld,                &
            traj0_sta%angles,                     &
            traj0%raytracing,                     &
            traj0%transmission_scatt_ir,          &
            traj0_dyn%transmission_scatt_ir_dyn,  &
            transmission)
  ENDIF

!-------------------------------------------------
! Calculate transmittances
!-------------------------------------------------
  IF (traj0_sta%dothermal) THEN
    ! Thermal path1 transmittances for all thermal channels
    CALL rttov_transmit( &
            opts,                                     &
            traj0_sta%do_lambertian,                  &
            profiles(1)%nlayers,                      &
            chanprof,                                 &
            traj0_sta%thermal,                        &
            traj0%aux_prof,                           &
            coefs%coef,                               &
            traj0%ircld,                              &
            traj0_sta%angles,                         &
            traj0%opdp_path%atm_level,                &
            traj0_sta%thermal_path1%od_level,         &
            transmission,                             &
            traj0_dyn%transmission_aux,               &
            traj0_dyn%transmission_aux%thermal_path1, &
            traj0%transmission_scatt_ir,              &
            traj0_dyn%transmission_scatt_ir_dyn,      &
            traj0_sta%thermal_path1%tau_surf,         &
            traj0_sta%thermal_path1%tau_level)
  ENDIF

  IF (traj0_sta%do_solar_opdep_calc) THEN
    ! Solar path2 transmittances for all solar channels
    CALL rttov_transmit_solar( &
            opts,                                   &
            profiles(1)%nlayers,                    &
            chanprof,                               &
            traj0_sta%solar,                        &
            traj0%aux_prof,                         &
            coefs%coef,                             &
            traj0%raytracing,                       &
            traj0%ircld,                            &
            traj0%opdp_path,                        &
            traj0_sta%solar_path2,                  &
            traj0_sta%solar_path1,                  &
            transmission,                           &
            traj0_dyn%transmission_aux,             &
            traj0%transmission_scatt_ir,            &
            traj0_dyn%transmission_scatt_ir_dyn)
  ENDIF

!-------------------------------------------------
! Calculate surface emissivity
!-------------------------------------------------
  IF (PRESENT(emissivity)) THEN
    ! Specify output emissivity values for all thermal channels and zeroes elsewhere.
    ! Emissivity models below will overwrite values where appropriate.
    ! This ensures consistent outputs when using the parallel interface
    ! which can result in some threads processing only solar or only
    ! thermal channels.
    WHERE (traj0_sta%thermal)
      emissivity%emis_out = emissivity%emis_in
    ELSEWHERE
      emissivity%emis_out = 0._jprb
    ENDWHERE
  ENDIF

  traj0%diffuse_refl(:) = 0._jprb

  IF (traj0_sta%dothermal) THEN

    WHERE (traj0_sta%thermal .AND. .NOT. calcemis)
      traj0%diffuse_refl = 1._jprb - emissivity%emis_out
    ENDWHERE

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%diffuse_refl for selected channels

      IF (sensor_mw) THEN
        ! Microwave
        CALL rttov_calcemis_mw( &
                opts,                          &
                profiles,                      &
                traj0_sta%angles,              &
                coefs%coef,                    &
                chanprof,                      &
                traj0_dyn%transmission_aux,    &
                calcemis,                      &
                emissivity%emis_out,           &
                traj0%diffuse_refl,            &
                err)
        THROWM(err.NE.0, "calcemis_mw")

      ELSE
        ! Infrared
        CALL rttov_calcemis_ir( &
                err,                  &
                opts,                 &
                chanprof,             &
                profiles,             &
                traj0_sta%angles,     &
                coefs,                &
                traj0_sta%thermal,    &
                calcemis,             &
                emissivity%emis_out)
        THROWM(err.NE.0, "calcemis_ir")

        WHERE (traj0_sta%thermal .AND. calcemis)
          traj0%diffuse_refl(:) = 1._jprb - emissivity(:)%emis_out
        ENDWHERE
      ENDIF

    ENDIF ! calcemis

  ENDIF ! dothermal

!-------------------------------------------------------
! Calculate surface reflectance
!-------------------------------------------------------
  IF (PRESENT(reflectance)) THEN
    ! Specify output BRDF values for all solar channels and zeroes elsewhere.
    ! BRDF models below will overwrite values where appropriate.
    ! This ensures consistent outputs when using the parallel interface
    ! which can result in some threads processing only solar or only
    ! thermal channels.
    WHERE (traj0_sta%solar)
      reflectance%refl_out = reflectance%refl_in
    ELSEWHERE
      reflectance%refl_out = 0._jprb
    ENDWHERE
  ENDIF

  IF (traj0_sta%dosolar) THEN

    IF (ANY(calcrefl)) THEN
      CALL rttov_refsun( &
              opts,             &
              profiles,         &
              coefs%coef,       &
              traj0%aux_prof,   &
              traj0%sunglint,   &
              traj0%raytracing)

      CALL rttov_fresnel( &
              chanprof,        &
              calcrefl,        &
              profiles,        &
              traj0_sta%solar, &
              coefs%coef,      &
              traj0%sunglint,  &
              traj0%fresnrefl)
    ENDIF

    ! rttov_calcsurfrefl populates refl_norm and diffuse_refl so is called regardless of calcrefl
    CALL rttov_calcsurfrefl(     &
            coefs%coef,          &
            profiles,            &
            traj0%sunglint,      &
            traj0%fresnrefl,     &
            traj0_sta%solar,     &
            chanprof,            &
            traj0_sta%refl_norm, &
            calcrefl,            &
            emissivity,          &
            reflectance,         &
            traj0%diffuse_refl)
    THROWM(err.NE.0, "calcsurfrefl")

  ENDIF

  IF (PRESENT(reflectance)) THEN
    reflectance(:)%diffuse_refl_out = traj0%diffuse_refl(:) * pi_r
  ENDIF

!---------------------------------------------------------
! Solve the radiative transfer equation
!---------------------------------------------------------
  IF (traj0_sta%do_mfasis) THEN
    CALL rttov_mfasis( &
            err,                           &
            chanprof,                      &
            traj0_sta%solar,               &
            opts,                          &
            profiles,                      &
            traj0%profiles_int,            &
            coefs,                         &
            traj0%ircld,                   &
            traj0%aux_prof,                &
            reflectance,                   &
            traj0_sta%solar_spec_esd,      &
            traj0%transmission_scatt_ir,   &
            radiance,                      &
            traj0_dyn%mfasis_refl)
    THROW(err.NE.0)
  ENDIF

  IF (traj0_sta%dothermal .OR. .NOT. traj0_sta%do_mfasis) THEN
    CALL rttov_integrate( &
            sensor_mw,                              &
            opts,                                   &
            traj0_dyn%ncolumns,                     &
            chanprof,                               &
            emissivity,                             &
            reflectance,                            &
            traj0_sta%refl_norm,                    &
            traj0%diffuse_refl,                     &
            traj0_sta%do_lambertian,                &
            traj0_sta%do_mfasis,                    &
            traj0_sta%thermal,                      &
            traj0_sta%dothermal,                    &
            traj0_sta%solar,                        &
            traj0_sta%dosolar,                      &
            traj0_sta%do_rayleigh_ss,               &
            traj0_sta%solar_spec_esd,               &
            traj0_dyn%transmission_aux,             &
            traj0%transmission_scatt_ir,            &
            profiles,                               &
            traj0%profiles_int,                     &
            traj0%aux_prof,                         &
            coefs%coef,                             &
            traj0%raytracing,                       &
            traj0%ircld,                            &
            radiance,                               &
            radiance2,                              &
            traj0%auxrad,                           &
            traj0_dyn%auxrad_column)

    IF (traj0_sta%dosolar .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
      CALL rttov_dom_setup_profile(                 &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .TRUE._jplm,                            &  ! TRUE => solar source term
            traj0_sta%solar,                        &
            traj0_sta%do_rayleigh_dom,              &
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0%opdp_path%sun_level_path2,        &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_dyn%profiles_dom_solar)

      CALL rttov_dom(                                     &
            err,                                          &
            opts,                                         &
            chanprof,                                     &
            .TRUE._jplm,                                  &  ! TRUE => solar source term
            traj0_sta%solar,                              &
            traj0_sta%dom_nstreams,                       &
            traj0_sta%dom_nstreams - 1_jpim,              &
            profiles,                                     &
            traj0_dyn%profiles_dom_solar,                 &
            MAXVAL(traj0_dyn%profiles_dom_solar%nlayers), &
            traj0%auxrad,                                 &
            traj0%transmission_scatt_ir,                  &
            traj0_dyn%transmission_scatt_ir_dyn,          &
            emissivity,                                   &
            reflectance,                                  &
            traj0%diffuse_refl,                           &
            traj0_sta%solar_spec_esd,                     &
            traj0%raytracing,                             &
            traj0%ircld,                                  &
            radiance,                                     &
            traj0_dyn%dom_state_solar)
      THROW(err.NE.0)
    ENDIF

    IF (traj0_sta%dothermal .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
      CALL rttov_dom_setup_profile(                 &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .FALSE._jplm,                           &  ! FALSE => thermal source term
            traj0_sta%thermal,                      &
            .FALSE._jplm,                           &  ! FALSE => do_rayleigh_dom
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0%opdp_path%atm_level,              &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_dyn%profiles_dom_thermal)

      CALL rttov_dom(                                       &
            err,                                            &
            opts,                                           &
            chanprof,                                       &
            .FALSE._jplm,                                   &  ! FALSE => thermal source term
            traj0_sta%thermal,                              &
            traj0_sta%dom_nstreams,                         &
            0_jpim,                                         & ! maxnaz
            profiles,                                       &
            traj0_dyn%profiles_dom_thermal,                 &
            MAXVAL(traj0_dyn%profiles_dom_thermal%nlayers), &
            traj0%auxrad,                                   &
            traj0%transmission_scatt_ir,                    &
            traj0_dyn%transmission_scatt_ir_dyn,            &
            emissivity,                                     &
            reflectance,                                    &
            traj0%diffuse_refl,                             &
            traj0_sta%solar_spec_esd,                       &
            traj0%raytracing,                               &
            traj0%ircld,                                    &
            radiance,                                       &
            traj0_dyn%dom_state_thermal)
      THROW(err.NE.0)
    ENDIF
  ENDIF ! dothermal .OR. .NOT. do_mfasis

!---------------------------------------------------------
! Do NLTE bias correction for hyperspectral instruments
!---------------------------------------------------------
  IF (coefs%coef%nltecoef) THEN
    IF (opts%rt_ir%do_nlte_correction .AND. coefs%coef%id_sensor == sensor_id_hi) THEN
      CALL rttov_nlte_bias_correction(opts, coefs%coef, profiles, chanprof, radiance)
    ENDIF
  ENDIF

!---------------------------------------------------------
! Calculate output BTs/reflectances/PCscores
!---------------------------------------------------------
  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_pcscores( &
            opts,                        &
            chanprof,                    &
            traj0_sta%chanprof_pc,       &
            pccomp,                      &
            coefs%coef_pccomp,           &
            radiance)

    IF (opts%rt_ir%pc%addradrec) THEN
      CALL rttov_reconstruct( &
              opts,                        &
              traj0_sta%chanprof_in,       &
              traj0_sta%chanprof_pc,       &
              pccomp,                      &
              coefs%coef_pccomp)

      CALL rttov_calcbt_pc(traj0_sta%chanprof_in, coefs%coef_pccomp, pccomp)
    ENDIF
  ELSE
    IF (traj0_sta%dothermal) THEN
      CALL rttov_calcbt(chanprof, coefs%coef, traj0_sta%thermal, radiance)
    ENDIF
    IF (traj0_sta%do_solar_opdep_calc) THEN
      CALL rttov_calcsatrefl(chanprof, profiles, traj0_sta%solar_spec_esd, traj0_sta%solar, radiance)
    ENDIF
  ENDIF

  radiance%plane_parallel = traj0_sta%plane_parallel

!---------------------
! Deallocate memory
!---------------------
998 CONTINUE ! from HTFRTC or if no channels are active
  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_DIRECT', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (opts%htfrtc_opts%htfrtc) RETURN

    IF (ltraj_dyn_dealloc) THEN
      IF (ASSOCIATED(traj0_dyn)) THEN
        CALL rttov_alloc_traj_dyn(error, traj0_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                                  traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                                  traj0_sta%thermal, traj0_sta%solar, &
                                  traj0_sta%dothermal, traj0_sta%do_mfasis, 0_jpim)
      ENDIF
    ENDIF

    IF (ltraj_sta_dealloc) THEN
      IF (ASSOCIATED(traj0_sta)) THEN
        CALL rttov_alloc_traj_sta(error, traj0_sta, opts, coefs, chanprof, profiles, &
                                  0_jpim, npcscores, channels_rec)
      ENDIF
    ENDIF

    IF (ASSOCIATED(traj0)) THEN
      CALL rttov_check_traj( &
              error,         &
              nprofiles,     &
              nchanprof,     &
              opts,          &
              nlevels,       &
              coefs,         &
              0_jpim,        &
              traj0 = traj0, &
              traj1 = traj1, &
              traj2 = traj)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_direct
