! Description:
!> @file
!!   Runs RTTOV tangent linear (TL) model
!
!> @brief
!!   Runs RTTOV tangent linear (TL) model
!!
!! @details
!!   Given an input profile and profile perturbation, computes
!!   the corresponding radiance perturbation from the tangent
!!   linear (TL) to the direct model evaluated at the input profile.
!!
!! @param[out]    errorstatus       status on exit
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     opts              options to configure the simulations
!! @param[in]     profiles          input atmospheric profiles and surface variables
!! @param[in]     profiles_tl       input atmospheric profile and surface variable perturbations
!! @param[in]     coefs             coefficients structure for instrument to simulate
!! @param[in,out] transmission      output transmittances
!! @param[in,out] transmission_tl   output transmittance perturbations
!! @param[in,out] radiance          output radiances and corresponding BTs and BRFs
!! @param[in,out] radiance_tl       output radiance, BT and BRF perturbations
!! @param[in,out] radiance2         secondary output radiances, optional
!! @param[in]     calcemis          flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity        input/output surface emissivities, optional
!! @param[in,out] emissivity_tl     input/output surface emissivity perturbations, optional
!! @param[in]     calcrefl          flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance       input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_tl    input/output surface BRDF perturbations, optional
!! @param[in]     aer_opt_param     input aerosol optical parameters, optional
!! @param[in]     aer_opt_param_tl  input aerosol optical parameter pertubations, optional
!! @param[in]     cld_opt_param     input cloud optical parameters, optional
!! @param[in]     cld_opt_param_tl  input cloud optical parameter pertubations, optional
!! @param[in,out] traj              RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_tl           RTTOV TL internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp            output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_tl         output PC score and radiance perturbations from PC-RTTOV, optional
!! @param[in]     channels_rec      list of channels for which to calculate reconstructed radiances, optional
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
SUBROUTINE rttov_tl( &
              errorstatus,       &
              chanprof,          &
              opts,              &
              profiles,          &
              profiles_tl,       &
              coefs,             &
              transmission,      &
              transmission_tl,   &
              radiance,          &
              radiance_tl,       &
              radiance2,         &
              calcemis,          &
              emissivity,        &
              emissivity_tl,     &
              calcrefl,          &
              reflectance,       &
              reflectance_tl,    &
              aer_opt_param,     &
              aer_opt_param_tl,  &
              cld_opt_param,     &
              cld_opt_param_tl,  &
              traj,              &
              traj_tl,           &
              pccomp,            &
              pccomp_tl,         &
              channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
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
         rttov_traj
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
         sensor_id_mw,     &
         sensor_id_hi,     &
         sensor_id_po,     &
         vis_scatt_dom,    &
         ir_scatt_dom,     &
         pi_r
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_types, ONLY : &
         rttov_traj_dyn,  &
         rttov_traj_sta
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),       INTENT(OUT)                     :: errorstatus
  TYPE(rttov_chanprof),     INTENT(IN)                      :: chanprof(:)
  TYPE(rttov_options),      INTENT(IN)                      :: opts
  TYPE(rttov_profile),      INTENT(IN)                      :: profiles(:)
  TYPE(rttov_profile),      INTENT(IN)                      :: profiles_tl(SIZE(profiles))
  TYPE(rttov_coefs),        INTENT(IN)   , TARGET           :: coefs
  TYPE(rttov_transmission), INTENT(INOUT)                   :: transmission
  TYPE(rttov_transmission), INTENT(INOUT)                   :: transmission_tl
  TYPE(rttov_radiance),     INTENT(INOUT)                   :: radiance
  TYPE(rttov_radiance),     INTENT(INOUT)                   :: radiance_tl
  TYPE(rttov_radiance2),    INTENT(INOUT), OPTIONAL         :: radiance2
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity),   INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity),   INTENT(INOUT), OPTIONAL         :: emissivity_tl(SIZE(chanprof))
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance),  INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),  INTENT(INOUT), OPTIONAL         :: reflectance_tl(SIZE(chanprof))
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: aer_opt_param_tl
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: cld_opt_param_tl
  TYPE(rttov_traj),         INTENT(INOUT), OPTIONAL, TARGET :: traj, traj_tl    ! target is needed here
  TYPE(rttov_pccomp),       INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_pccomp),       INTENT(INOUT), OPTIONAL         :: pccomp_tl
  INTEGER(KIND=jpim),       INTENT(IN)   , OPTIONAL         :: channels_rec(:)
!INTF_END
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_apply_pc_aer_reg_lims_tl.interface"
#include "rttov_apply_reg_limits_tl.interface"
#include "rttov_calc_nearest_lev_tl.interface"
#include "rttov_calcbt_pc_tl.interface"
#include "rttov_calcbt_tl.interface"
#include "rttov_calcemis_ir_tl.interface"
#include "rttov_calcemis_mw_tl.interface"
#include "rttov_calcsatrefl_tl.interface"
#include "rttov_calcsurfrefl_tl.interface"
#include "rttov_check_traj.interface"
#include "rttov_cloud_overlap_tl.interface"
#include "rttov_convert_profile_units_tl.interface"
#include "rttov_copy_aux_prof.interface"
#include "rttov_copy_opdp_path.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_raytracing.interface"
#include "rttov_direct.interface"
#include "rttov_dom_setup_profile_tl.interface"
#include "rttov_dom_tl.interface"
#include "rttov_errorreport.interface"
#include "rttov_fresnel_tl.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_intavg_chan_tl.interface"
#include "rttov_intavg_prof_tl.interface"
#include "rttov_integrate_tl.interface"
#include "rttov_locpat_tl.interface"
#include "rttov_mfasis_tl.interface"
#include "rttov_mw_clw_absorption_tl.interface"
#include "rttov_nlte_bias_correction_tl.interface"
#include "rttov_opdep_13_tl.interface"
#include "rttov_opdep_9_tl.interface"
#include "rttov_opdep_78_tl.interface"
#include "rttov_opdpscattir_tl.interface"
#include "rttov_pcscores_tl.interface"
#include "rttov_predictor_precalc_13_tl.interface"
#include "rttov_predictor_precalc_789_tl.interface"
#include "rttov_profaux_cldaer_tl.interface"
#include "rttov_rayleigh_extinction_tl.interface"
#include "rttov_reconstruct_tl.interface"
#include "rttov_refsun_tl.interface"
#include "rttov_setpredictors_13_tl.interface"
#include "rttov_setpredictors_789_tl.interface"
#include "rttov_transmit_solar_tl.interface"
#include "rttov_transmit_tl.interface"

  INTEGER(KIND=jpim) :: err
  INTEGER(KIND=jpim) :: nlevels, nprofiles, nchanprof, npcscores
  LOGICAL(KIND=jplm) :: sensor_mw, sensor_ir

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_tl
  TYPE(rttov_traj), POINTER :: traj0_tl
  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_tl_dyn
  LOGICAL(KIND=jplm)   :: ltraj_tl_dyn_dealloc

  REAL(KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 0_jpim, ZHOOK_HANDLE)

  errorstatus = errorstatus_success

  IF (opts%htfrtc_opts%htfrtc) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "HT-FRTC: TL model not available")
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

  ltraj_tl_dyn_dealloc = .FALSE.
  traj0_dyn%from_tladk = .TRUE.

  CALL rttov_init_transmission(transmission_tl)
  CALL rttov_init_rad(radiance_tl)

  IF (opts%rt_ir%pc%addpc) THEN
    IF (.NOT. PRESENT(pccomp_tl)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "pccomp_tl argument is mandatory for PC-RTTOV")
    ENDIF
    npcscores = opts%rt_ir%pc%npcscores * nprofiles
  ENDIF

  NULLIFY (traj0, traj0_tl)
  CALL rttov_check_traj( &
          err,                 &
          nprofiles,           &
          nchanprof,           &
          opts,                &
          nlevels,             &
          coefs,               &
          1_jpim,              &
          traj0 = traj0,       &
          traj0_tl = traj0_tl, &
          traj1 = traj1,       &
          traj1_tl = traj1_tl, &
          traj2 = traj,        &
          traj2_tl = traj_tl)
  THROWM(err.NE.0, "rttov_check_traj fatal error")

!-----------------------------------------------------------------------
! Call direct model
!-----------------------------------------------------------------------
  CALL rttov_direct( &
              err,           &
              chanprof,      &
              opts,          &
              profiles,      &
              coefs,         &
              transmission,  &
              radiance,      &
              radiance2,     &
              calcemis,      &
              emissivity,    &
              calcrefl,      &
              reflectance,   &
              aer_opt_param, &
              cld_opt_param, &
              traj0,         &
              traj0_dyn,     &
              traj0_sta,     &
              pccomp,        &
              channels_rec)
  THROW(err.NE.0)

  radiance_tl%plane_parallel = traj0_sta%plane_parallel
  radiance_tl%quality        = radiance%quality

  ! If no channels are active (thermal and solar flags false for all channels) then
  ! skip to the end. In particular this can happen when using the parallel interface.
  IF (.NOT. (traj0_sta%dothermal .OR. traj0_sta%dosolar)) GOTO 998

  IF (traj0_sta%dothermal .AND. .NOT. PRESENT(emissivity_tl)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "emissivity_tl parameter required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. PRESENT(reflectance_tl)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "reflectance_tl parameter required")
  END IF

  CALL rttov_alloc_traj_dyn(err, traj0_tl_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                            traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                            traj0_sta%thermal, traj0_sta%solar, &
                            traj0_sta%dothermal, traj0_sta%do_mfasis, 1_jpim, traj0_dyn)
  THROW(err.NE.0)
  ltraj_tl_dyn_dealloc = .TRUE.


!----------------
! Tangent Linear
!----------------

!-----------------------------------------------------------------------
! Convert input profiles to units used internally
!-----------------------------------------------------------------------
  CALL rttov_convert_profile_units_tl(opts, coefs, profiles, profiles_tl, &
                                      traj0%profiles_int, traj0_tl%profiles_int)

!-----------------------------------------------------------------------
! Calculate near surface and cloud top levels
!-----------------------------------------------------------------------
  CALL rttov_calc_nearest_lev_tl(sensor_ir, opts, profiles, profiles_tl, &
                                 traj0%aux_prof%s, traj0_tl%aux_prof%s)

!------------------------------------------------------------------------
! Local path geometry calculations
!------------------------------------------------------------------------
  CALL rttov_locpat_tl( &
      opts,                     &
      traj0_sta%dosolar,        &
      traj0_sta%plane_parallel, &
      profiles,                 &
      profiles_tl,              &
      traj0%profiles_int,       &
      traj0_tl%profiles_int,    &
      traj0%aux_prof%s,         &
      traj0_sta%angles,         &
      traj0%raytracing,         &
      traj0_tl%raytracing)

!-----------------------------------------------------------------------
! RTTOV optical depth calculation: on coef levels
!-----------------------------------------------------------------------

  IF (traj0_sta%do_opdep_calc) THEN
!-----------------------------------------------------------------------
! Interpolator first call - input profiles from user levs to coef levs
!-----------------------------------------------------------------------
    CALL rttov_init_prof(traj0_tl%profiles_coef)

    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_prof_tl( &
              opts,                           &
              nlevels,                        &
              coefs%coef%nlevels,             &
              profiles,                       &
              profiles_tl,                    &
              traj0%profiles_int,             &
              traj0_tl%profiles_int,          &
              traj0%profiles_coef,            &
              traj0_tl%profiles_coef,         &
              coefs%coef,                     &
              coefs%coef_pccomp)
    ELSE
      ! If interpolator is off copy profile variables
      CALL rttov_copy_prof( &
              traj0_tl%profiles_coef, &
              profiles_tl,            &
              larray = .TRUE._jplm,   &
              lscalar = .FALSE._jplm, &
              profiles_gas = traj0_tl%profiles_int)
    ENDIF

    ! Copy non-interpolated variables (e.g. surface parameters)
    CALL rttov_copy_prof( &
            traj0_tl%profiles_coef, &
            profiles_tl,            &
            larray = .FALSE._jplm,  &
            lscalar = .TRUE._jplm,  &
            profiles_gas = traj0_tl%profiles_int)

!------------------------------------------------------------------
! Check input data is within suitable physical limits - coef levels
!------------------------------------------------------------------
    IF (opts%config%apply_reg_limits .OR. opts%interpolation%addinterp) THEN
      CALL rttov_apply_reg_limits_tl( &
              opts,                        &
              profiles,                    &
              traj0_sta%profiles_coef_ref, &
              traj0_tl%profiles_coef,      &
              coefs%coef,                  &
              coefs%coef_pccomp)
    ENDIF

!--------------------------------------------------------------------------------
! Calculate near surface level and local path geometry calculations - coef levels
!--------------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_calc_nearest_lev_tl(.FALSE._jplm, opts, traj0%profiles_coef, &
            traj0_tl%profiles_coef, traj0%aux_prof_coef%s, traj0_tl%aux_prof_coef%s)
      CALL rttov_locpat_tl( &
          opts,                     &
          traj0_sta%dosolar,        &
          traj0_sta%plane_parallel, &
          traj0%profiles_coef,      &
          traj0_tl%profiles_coef,   &
          traj0%profiles_coef,      &
          traj0_tl%profiles_coef,   &
          traj0%aux_prof_coef%s,    &
          traj0_sta%angles,         &
          traj0%raytracing_coef,    &
          traj0_tl%raytracing_coef)
    ELSE
      CALL rttov_copy_aux_prof(traj0_tl%aux_prof_coef%s, traj0_tl%aux_prof%s)
      CALL rttov_copy_raytracing(traj0_sta%dosolar, traj0_tl%raytracing_coef, traj0_tl%raytracing)
    ENDIF

!---------------------------------------------------
! Calculate predictors - coef levels
!---------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      CALL rttov_predictor_precalc_13_tl( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0_sta%plane_parallel,      &
          traj0%profiles_coef,           &
          traj0_tl%profiles_coef,        &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef)

      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_13_tl( &
          opts,                          &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0%predictors%path1,        &
          traj0_tl%predictors%path1,     &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef,      &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_13_tl( &
          opts,                          &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0%predictors%path2,        &
          traj0_tl%predictors%path2,     &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef,      &
          raypath = 2_jpim) 
      ENDIF
    ELSE
      CALL rttov_predictor_precalc_789_tl( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0_sta%plane_parallel,      &
          traj0%profiles_coef,           &
          traj0_tl%profiles_coef,        &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0_sta%angles,              &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef)

      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_789_tl( &
          opts,                          &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0%predictors%path1,        &
          traj0_tl%predictors%path1,     &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef,      &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_789_tl( &
          opts,                          &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_tl%aux_prof_coef,        &
          traj0%predictors%path2,        &
          traj0_tl%predictors%path2,     &
          traj0%raytracing_coef,         &
          traj0_tl%raytracing_coef,      &
          raypath = 2_jpim) 
      ENDIF
    ENDIF

!--------------------------------------------------------------------
! Predict atmospheric (emissive) and solar optical depths - coef levs
!--------------------------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for all thermal channels
        CALL rttov_opdep_13_tl( &
                coefs%coef%nlayers,                &
                chanprof,                          &
                traj0_sta%thermal,                 &
                traj0_tl%predictors%path1,         &
                coefs%coef,                        &
                coefs%coef%thermal,                &
                coefs%coef%thermal_corr,           &
                traj0_tl%opdp_path_coef%atm_level, &
                traj0_sta%thermal_path1%opdp_ref)
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for all solar channels
        CALL rttov_opdep_13_tl( &
                coefs%coef%nlayers,                      &
                chanprof,                                &
                traj0_sta%solar,                         &
                traj0_tl%predictors%path2,               &
                coefs%coef,                              &
                coefs%coef%solar,                        &
                coefs%coef%solar_corr,                   &
                traj0_tl%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref)
      ENDIF
    ELSEIF (coefs%coef%fmv_model_ver == 9) THEN
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for all thermal channels
        CALL rttov_opdep_9_tl( &
                coefs%coef%nlayers,                    &
                chanprof,                              &
                traj0_sta%thermal,                     &
                traj0%predictors%path1,                &
                traj0_tl%predictors%path1,             &
                coefs%coef,                            &
                coefs%coef%thermal,                    &
                traj0_tl%opdp_path_coef%atm_level,     &
                traj0_sta%thermal_path1%opdp_ref_coef)
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for all solar channels
        CALL rttov_opdep_9_tl( &
                coefs%coef%nlayers,                      &
                chanprof,                                &
                traj0_sta%solar,                         &
                traj0%predictors%path2,                  &
                traj0_tl%predictors%path2,               &
                coefs%coef,                              &
                coefs%coef%solar,                        &
                traj0_tl%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref_coef)
      ENDIF
    ELSE
      CALL rttov_opdep_78_tl( &
              coefs%coef%nlayers,           &
              chanprof,                     &
              traj0_tl%predictors%path1,    &
              coefs%coef,                   &
              coefs%coef%thermal,           &
              traj0_tl%opdp_path_coef,      &
              traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF

!--------------------------------------------------------------------------
! Interpolator second  call - optical depths from coef levs to user levs
!--------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_chan_tl( &
              opts,                    &
              traj0_sta%thermal,       &
              traj0_sta%solar,         &
              coefs%coef%nlevels,      &
              nlevels,                 &
              chanprof,                &
              traj0%profiles_coef,     &
              profiles,                &
              profiles_tl,             &
              traj0%opdp_path_coef,    &
              traj0_tl%opdp_path_coef, &
              traj0_tl%opdp_path)
    ELSE
      CALL rttov_copy_opdp_path(opts, traj0_tl%opdp_path, traj0_tl%opdp_path_coef)
    ENDIF

!--------------------------------------------------------------------------
! Move top level to space boundary (opdep=0) if spacetop flag is set
!--------------------------------------------------------------------------
    IF (opts%interpolation%spacetop) THEN
      traj0_tl%opdp_path%atm_level(1,:) = 0._jprb
      IF (traj0_sta%do_solar_opdep_calc) THEN
        traj0_tl%opdp_path%sun_level_path2(1,:) = 0._jprb
      ENDIF
    ENDIF

  ENDIF ! do_opdep_calc

!--------------------------------------------------------------------------
! MW CLW absorption optical depths
!--------------------------------------------------------------------------
  IF (sensor_mw .AND. opts%rt_mw%clw_data) THEN
    CALL rttov_mw_clw_absorption_tl( &
            opts,                &
            coefs%coef,          &
            chanprof,            &
            traj0%raytracing,    &
            traj0_tl%raytracing, &
            profiles,            &
            profiles_tl,         &
            traj0%aux_prof,      &
            traj0_tl%aux_prof,   &
            traj0_tl%opdp_path)
  ENDIF

!--------------------------------------------------------------------------
! Rayleigh extinction optical depths
!--------------------------------------------------------------------------
  IF (traj0_sta%do_rayleigh_param) THEN
    CALL rttov_rayleigh_extinction_tl( &
            opts,                               &
            chanprof,                           &
            traj0_sta%do_rayleigh_dom,          &
            traj0_sta%solar,                    &
            profiles,                           &
            profiles_tl,                        &
            coefs%coef,                         &
            traj0%raytracing,                   &
            traj0_tl%raytracing,                &
            traj0_tl%opdp_path%sun_level_path2, &
            traj0_tl%transmission_scatt_ir)
  ENDIF

!------------------------------------------------------
! Calculate the number of cloud columns and
! the cloud distribution in each column
!------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cloud_overlap_tl( &
            opts%rt_ir,            &
            profiles,              &
            profiles_tl,           &
            traj0%profiles_int,    &
            traj0_tl%profiles_int, &
            traj0%ircld,           &
            traj0_tl%ircld)
  ELSE
    traj0_tl%ircld%xcolclr = 0._jprb
  ENDIF

!-----------------------------------------------------------------------
! Apply PC-RTTOV aerosol regression limits
!-----------------------------------------------------------------------
  IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addaerosl) THEN
    CALL rttov_apply_pc_aer_reg_lims_tl( &
          opts,                  &
          coefs%coef_pccomp,     &
          profiles,              &
          traj0_tl%profiles_int, &
          traj0_sta%pc_aer_ref,  &
          traj0_sta%pc_aer_min,  &
          traj0_sta%pc_aer_max)
  ENDIF

!----------------------------------------------------------------------------
! Calculate optical depths of aerosols and/or clouds
!----------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_profaux_cldaer_tl( &
            opts,                  &
            profiles,              &
            profiles_tl,           &
            traj0%profiles_int,    &
            traj0_tl%profiles_int, &
            traj0%aux_prof,        &
            traj0_tl%aux_prof)

    CALL rttov_opdpscattir_tl( &
            profiles(1)%nlayers,                     &
            chanprof,                                &
            opts,                                    &
            traj0_sta%dom_nstreams,                  &
            traj0%aux_prof,                          &
            traj0_tl%aux_prof,                       &
            traj0%ircld,                             &
            profiles,                                &
            profiles_tl,                             &
            traj0%profiles_int,                      &
            traj0_tl%profiles_int,                   &
            aer_opt_param,                           &
            aer_opt_param_tl,                        &
            cld_opt_param,                           &
            cld_opt_param_tl,                        &
            traj0_sta%dothermal,                     &
            traj0_sta%thermal,                       &
            traj0_sta%dosolar,                       &
            traj0_sta%solar,                         &
            traj0_sta%do_rayleigh_dom,               &
            coefs%coef,                              &
            coefs%coef_scatt,                        &
            coefs%coef_mfasis_cld,                   &
            traj0_sta%angles,                        &
            traj0%raytracing,                        &
            traj0_tl%raytracing,                     &
            traj0%transmission_scatt_ir,             &
            traj0_tl%transmission_scatt_ir,          &
            traj0_dyn%transmission_scatt_ir_dyn,     &
            traj0_tl_dyn%transmission_scatt_ir_dyn)
  ENDIF

!-------------------------------------------------
! Calculate transmittances
!-------------------------------------------------
  IF (traj0_sta%dothermal) THEN
    CALL rttov_transmit_tl( &
            opts,                                           &
            traj0_sta%do_lambertian,                        &
            profiles(1)%nlayers,                            &
            chanprof,                                       &
            traj0_sta%thermal,                              &
            traj0%aux_prof,                                 &
            traj0_tl%aux_prof,                              &
            coefs%coef,                                     &
            traj0%ircld,                                    &
            traj0_sta%angles,                               &
            traj0_tl%opdp_path%atm_level,                   &
            traj0_sta%thermal_path1%od_level,               &
            transmission_tl,                                &
            traj0_dyn%transmission_aux%thermal_path1,       &
            traj0_tl_dyn%transmission_aux%thermal_path1,    &
            traj0_tl%transmission_scatt_ir,                 &
            traj0_dyn%transmission_scatt_ir_dyn,            &
            traj0_tl_dyn%transmission_scatt_ir_dyn,         &
            traj0_sta%thermal_path1%tau_surf,               &
            traj0_sta%thermal_path1%tau_level)
  ENDIF

  IF (traj0_sta%do_solar_opdep_calc) THEN
    CALL rttov_transmit_solar_tl( &
            opts,                                          &
            profiles(1)%nlayers,                           &
            nprofiles,                                     &
            chanprof,                                      &
            traj0_sta%solar,                               &
            traj0%aux_prof,                                &
            traj0_tl%aux_prof,                             &
            coefs%coef,                                    &
            traj0%raytracing,                              &
            traj0_tl%raytracing,                           &
            traj0%ircld,                                   &
            traj0_tl%opdp_path,                            &
            traj0_sta%solar_path2,                         &
            traj0_sta%solar_path1,                         &
            transmission_tl,                               &
            traj0_dyn%transmission_aux,                    &
            traj0_tl_dyn%transmission_aux,                 &
            traj0%transmission_scatt_ir,                   &
            traj0_tl%transmission_scatt_ir,                &
            traj0_dyn%transmission_scatt_ir_dyn,           &
            traj0_tl_dyn%transmission_scatt_ir_dyn)

  ENDIF

!-------------------------------------------------
! Calculate surface emissivity
!-------------------------------------------------
  IF (PRESENT(emissivity)) THEN
    WHERE (traj0_sta%thermal)
      emissivity_tl%emis_out = emissivity_tl%emis_in
    ELSEWHERE
      emissivity_tl%emis_out = 0._jprb
    ENDWHERE
  ENDIF

  traj0_tl%diffuse_refl(:) = 0._jprb

  IF (traj0_sta%dothermal) THEN

    WHERE (traj0_sta%thermal .AND. .NOT. calcemis)
      traj0_tl%diffuse_refl =  - emissivity_tl%emis_out
    ENDWHERE

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%diffuse_refl for selected channels

      IF (sensor_mw) THEN
        CALL rttov_calcemis_mw_tl( &
                opts,                             &
                profiles,                         &
                profiles_tl,                      &
                traj0_sta%angles,                 &
                coefs%coef,                       &
                chanprof,                         &
                traj0_dyn%transmission_aux,       &
                traj0_tl_dyn%transmission_aux,    &
                calcemis,                         &
                emissivity_tl%emis_out,           &
                traj0_tl%diffuse_refl)
      ELSE
        CALL rttov_calcemis_ir_tl( &
                err,               &
                opts,              &
                chanprof,          &
                profiles,          &
                profiles_tl,       &
                coefs,             &
                traj0_sta%thermal, &
                calcemis,          &
                emissivity_tl%emis_out)
        THROWM(err.NE.0, "calcemis_ir_tl")

        WHERE (traj0_sta%thermal .AND. calcemis)
          traj0_tl%diffuse_refl =  - emissivity_tl%emis_out
        ENDWHERE
      ENDIF

    ENDIF

  ENDIF

!-------------------------------------------------------
! Calculate surface reflectance
!-------------------------------------------------------
  IF (PRESENT(reflectance)) THEN
    WHERE (traj0_sta%solar)
      reflectance_tl%refl_out = reflectance_tl%refl_in
    ELSEWHERE
      reflectance_tl%refl_out = 0._jprb
    ENDWHERE
  ENDIF

  IF (traj0_sta%dosolar) THEN

    IF (ANY(calcrefl)) THEN
      CALL rttov_refsun_tl( &
              opts,                &
              profiles,            &
              profiles_tl,         &
              coefs%coef,          &
              traj0%aux_prof,      &
              traj0%sunglint,      &
              traj0_tl%sunglint,   &
              traj0%raytracing,    &
              traj0_tl%raytracing)

      CALL rttov_fresnel_tl( &
              chanprof,           &
              calcrefl,           &
              profiles,           &
              traj0_sta%solar,    &
              coefs%coef,         &
              traj0%sunglint,     &
              traj0_tl%sunglint,  &
              traj0_tl%fresnrefl)
    ENDIF

    CALL rttov_calcsurfrefl_tl(  &
            coefs%coef,          &
            profiles,            &
            traj0%sunglint,      &
            traj0_tl%sunglint,   &
            traj0%fresnrefl,     &
            traj0_tl%fresnrefl,  &
            traj0_sta%solar,     &
            chanprof,            &
            traj0_sta%refl_norm, &
            calcrefl,            &
            emissivity_tl,       &
            reflectance,         &
            reflectance_tl,      &
            traj0_tl%diffuse_refl)

  ENDIF

  IF (PRESENT(reflectance_tl)) THEN
    reflectance_tl(:)%diffuse_refl_out = traj0_tl%diffuse_refl(:) * pi_r
  ENDIF

!---------------------------------------------------------
! Solve the radiative transfer equation
!---------------------------------------------------------
  IF (traj0_sta%do_mfasis) THEN
    CALL rttov_mfasis_tl( &
            err,                            &
            chanprof,                       &
            traj0_sta%solar,                &
            opts,                           &
            profiles,                       &
            profiles_tl,                    &
            traj0%profiles_int,             &
            traj0_tl%profiles_int,          &
            coefs,                          &
            traj0%ircld,                    &
            traj0_tl%ircld,                 &
            traj0%aux_prof,                 &
            traj0_tl%aux_prof,              &
            reflectance,                    &
            reflectance_tl,                 &
            traj0_sta%solar_spec_esd,       &
            traj0%transmission_scatt_ir,    &
            traj0_tl%transmission_scatt_ir, &
            traj0_dyn%mfasis_refl,          &
            radiance_tl)
    THROW(err.NE.0)
  ENDIF

  IF (traj0_sta%dothermal .OR. .NOT. traj0_sta%do_mfasis) THEN
    CALL rttov_integrate_tl( &
            opts,                                         &
            traj0_dyn%ncolumns,                           &
            chanprof,                                     &
            emissivity,                                   &
            emissivity_tl,                                &
            reflectance,                                  &
            reflectance_tl,                               &
            traj0_sta%refl_norm,                          &
            traj0%diffuse_refl,                           &
            traj0_tl%diffuse_refl,                        &
            traj0_sta%do_lambertian,                      &
            traj0_sta%do_mfasis,                          &
            traj0_sta%thermal,                            &
            traj0_sta%dothermal,                          &
            traj0_sta%solar,                              &
            traj0_sta%dosolar,                            &
            traj0_sta%do_rayleigh_ss,                     &
            traj0_sta%solar_spec_esd,                     &
            traj0_dyn%transmission_aux,                   &
            traj0_tl_dyn%transmission_aux,                &
            traj0%transmission_scatt_ir,                  &
            traj0_tl%transmission_scatt_ir,               &
            profiles,                                     &
            profiles_tl,                                  &
            traj0%profiles_int,                           &
            traj0_tl%profiles_int,                        &
            traj0%aux_prof,                               &
            traj0_tl%aux_prof,                            &
            coefs%coef,                                   &
            traj0%raytracing,                             &
            traj0_tl%raytracing,                          &
            traj0%ircld,                                  &
            traj0_tl%ircld,                               &
            radiance,                                     &
            traj0%auxrad,                                 &
            traj0_tl%auxrad,                              &
            traj0_dyn%auxrad_column,                      &
            traj0_tl_dyn%auxrad_column,                   &
            radiance_tl)

    IF (traj0_sta%dosolar .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
      CALL rttov_dom_setup_profile_tl(              &
            err,                                    &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .TRUE._jplm,                            &  ! TRUE => solar source term
            traj0_sta%solar,                        &
            traj0_sta%do_rayleigh_dom,              &
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0_tl%aux_prof,                      &
            traj0%opdp_path%sun_level_path2,        &
            traj0_tl%opdp_path%sun_level_path2,     &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_tl%transmission_scatt_ir,         &
            traj0_dyn%profiles_dom_solar,           &
            traj0_tl_dyn%profiles_dom_solar)
      THROW(err.NE.0)

      CALL rttov_dom_tl(                                  &
            err,                                          &
            chanprof,                                     &
            .TRUE._jplm,                                  &  ! TRUE => solar source term
            traj0_sta%solar,                              &
            traj0_sta%dom_nstreams,                       &
            profiles,                                     &
            traj0_dyn%profiles_dom_solar,                 &
            traj0_tl_dyn%profiles_dom_solar,              &
            MAXVAL(traj0_dyn%profiles_dom_solar%nlayers), &
            traj0%auxrad,                                 &
            traj0_tl%auxrad,                              &
            traj0%transmission_scatt_ir,                  &
            traj0_tl%transmission_scatt_ir,               &
            traj0_dyn%transmission_scatt_ir_dyn,          &
            traj0_tl_dyn%transmission_scatt_ir_dyn,       &
            emissivity,                                   &
            emissivity_tl,                                &
            reflectance,                                  &
            reflectance_tl,                               &
            traj0%diffuse_refl,                           &
            traj0_tl%diffuse_refl,                        &
            traj0_sta%solar_spec_esd,                     &
            traj0%raytracing,                             &
            traj0%ircld,                                  &
            traj0_tl%ircld,                               &
            radiance_tl,                                  &
            traj0_dyn%dom_state_solar)
      THROW(err.NE.0)
    ENDIF

    IF (traj0_sta%dothermal .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
      CALL rttov_dom_setup_profile_tl(              &
            err,                                    &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .FALSE._jplm,                           &  ! FALSE => thermal source term
            traj0_sta%thermal,                      &
            .FALSE._jplm,                           &  ! FALSE => do_rayleigh_dom
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0_tl%aux_prof,                      &
            traj0%opdp_path%atm_level,              &
            traj0_tl%opdp_path%atm_level,           &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_tl%transmission_scatt_ir,         &
            traj0_dyn%profiles_dom_thermal,         &
            traj0_tl_dyn%profiles_dom_thermal)
      THROW(err.NE.0)

      CALL rttov_dom_tl(                                    &
            err,                                            &
            chanprof,                                       &
            .FALSE._jplm,                                   &  ! FALSE => thermal source term
            traj0_sta%thermal,                              &
            traj0_sta%dom_nstreams,                         &
            profiles,                                       &
            traj0_dyn%profiles_dom_thermal,                 &
            traj0_tl_dyn%profiles_dom_thermal,              &
            MAXVAL(traj0_dyn%profiles_dom_thermal%nlayers), &
            traj0%auxrad,                                   &
            traj0_tl%auxrad,                                &
            traj0%transmission_scatt_ir,                    &
            traj0_tl%transmission_scatt_ir,                 &
            traj0_dyn%transmission_scatt_ir_dyn,            &
            traj0_tl_dyn%transmission_scatt_ir_dyn,         &
            emissivity,                                     &
            emissivity_tl,                                  &
            reflectance,                                    &
            reflectance_tl,                                 &
            traj0%diffuse_refl,                             &
            traj0_tl%diffuse_refl,                          &
            traj0_sta%solar_spec_esd,                       &
            traj0%raytracing,                               &
            traj0%ircld,                                    &
            traj0_tl%ircld,                                 &
            radiance_tl,                                    &
            traj0_dyn%dom_state_thermal)
      THROW(err.NE.0)
    ENDIF
  ENDIF ! dothermal .OR. .NOT. do_mfasis

!---------------------------------------------------------
! Do NLTE bias correction for hyperspectral instruments
!---------------------------------------------------------
  IF (coefs%coef%nltecoef) THEN
    IF (opts%rt_ir%do_nlte_correction .AND. coefs%coef%id_sensor == sensor_id_hi) THEN
      CALL rttov_nlte_bias_correction_tl(opts, coefs%coef, profiles, profiles_tl, chanprof, radiance_tl)
    ENDIF
  ENDIF

!---------------------------------------------------------
! Calculate output BTs/reflectances/PCscores
!---------------------------------------------------------
  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_pcscores_tl(          &
            opts,                    &
            chanprof,                &
            traj0_sta%chanprof_pc,   &
            pccomp_tl,               &
            coefs%coef_pccomp,       &
            radiance_tl)

    IF (opts%rt_ir%pc%addradrec) THEN
      CALL rttov_reconstruct_tl(     &
              opts,                  &
              traj0_sta%chanprof_in, &
              traj0_sta%chanprof_pc, &
              pccomp_tl,             &
              coefs%coef_pccomp)

      CALL rttov_calcbt_pc_tl(&
           traj0_sta%chanprof_in, &
           coefs%coef_pccomp,     &
           pccomp,                &
           pccomp_tl)
    ENDIF
  ELSE
    IF (traj0_sta%dothermal) THEN
      CALL rttov_calcbt_tl(chanprof, coefs%coef, traj0_sta%thermal, radiance, radiance_tl)
    ENDIF
    IF (traj0_sta%do_solar_opdep_calc) THEN
      CALL rttov_calcsatrefl_tl(chanprof, profiles, traj0_sta%solar_spec_esd, traj0_sta%solar, radiance_tl)
    ENDIF
  ENDIF

!---------------------
! Deallocate memory
!---------------------
998 CONTINUE ! no channels are active
  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (opts%htfrtc_opts%htfrtc) RETURN

    IF (ltraj_tl_dyn_dealloc) THEN
      CALL rttov_alloc_traj_dyn(error, traj0_tl_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                                traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                                traj0_sta%thermal, traj0_sta%solar, &
                                traj0_sta%dothermal, traj0_sta%do_mfasis, 0_jpim, traj0_dyn)
    ENDIF

    IF (traj0_dyn%ncolumns >= 0) THEN
      CALL rttov_alloc_traj_dyn(error, traj0_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                                traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                                traj0_sta%thermal, traj0_sta%solar, &
                                traj0_sta%dothermal, traj0_sta%do_mfasis, 0_jpim)
    ENDIF

    IF (ASSOCIATED(traj0_sta%thermal)) THEN
      CALL rttov_alloc_traj_sta(error, traj0_sta, opts, coefs, chanprof, profiles, &
                                0_jpim, npcscores, channels_rec)
    ENDIF

    IF (ASSOCIATED(traj0)) THEN
      CALL rttov_check_traj( &
              error,               &
              nprofiles,           &
              nchanprof,           &
              opts,                &
              nlevels,             &
              coefs,               &
              0_jpim,              &
              traj0 = traj0,       &
              traj0_tl = traj0_tl, &
              traj1 = traj1,       &
              traj1_tl = traj1_tl, &
              traj2 = traj,        &
              traj2_tl = traj_tl)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_tl
