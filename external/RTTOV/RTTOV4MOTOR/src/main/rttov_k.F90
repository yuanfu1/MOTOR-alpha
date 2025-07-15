! Description:
!> @file
!!   Runs RTTOV Jacobian (K) model
!
!> @brief
!!   Runs RTTOV Jacobian (K) model
!!
!! @details
!!   Given an input profile, calculates the Jacobian of the direct
!!   model with the outputs being the derivatives of the channel
!!   radiances or BTs with respect to each profile variable.
!!
!!   All K variables should be initialised to zero before calling.
!!   An input radiance or BT perturbation should be specified. This is
!!   usually 1.0: the value linearly scales the corresponding output
!!   Jacobians.
!!
!!   The input perturbations should then be written to the radiance_k
!!   structure. If opts\%rt_all\%switchrad is false the perturbations are
!!   in radiance in radiance_k%total, otherwise they are in BT in
!!   radiance_k\%bt. For visible/near-IR channels (at wavelengths less
!!   than 3 microns) the input perturbations are always in
!!   radiance_k\%total.
!!
!!   The output Jacobians are contained in profiles_k.
!!
!!   For PC-RTTOV the input perturbations are specified in pccomp_k instead.
!!   If opts\%rt_ir\%addpc\%addradrec is false the perturbations are in PC
!!   scores in pccomp_k\%total_pcscores. Otherwise they are in
!!   pccomp_k\%total_pccomp or pccomp_k\%bt_pccomp according to the
!!   setting of opts\%rt_all\%switchrad (false/true respectively).
!!
!!   For PC-RTTOV profiles_k contains the Jacobians of the PC predictor
!!   channels. The Jacobians of PC scores wrt profile variables are in
!!   profiles_k_pc while the Jacobians of reconstructed radiances or BTs
!!   are contained in profiles_k_rec.
!!
!! @param[out]    errorstatus       status on exit
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     opts              options to configure the simulations
!! @param[in]     profiles          input atmospheric profiles and surface variables
!! @param[in,out] profiles_k        output Jacobians wrt atmospheric profile and surface variables
!! @param[in]     coefs             coefficients structure for instrument to simulate
!! @param[in,out] transmission      output transmittances
!! @param[in,out] transmission_k    input perturbations wrt transmittances (usually zero)
!! @param[in,out] radiance          output radiances and corresponding BTs and BRFs
!! @param[in,out] radiance_k        input perturbations wrt radiances or BTs
!! @param[in,out] radiance2         secondary output radiances, optional
!! @param[in]     calcemis          flags for internal RTTOV surface emissivity calculation, optional
!! @param[in,out] emissivity        input/output surface emissivities, optional
!! @param[in,out] emissivity_k      output Jacobians wrt surface emissivities, optional
!! @param[in]     calcrefl          flags for internal RTTOV surface BRDF calculation, optional
!! @param[in,out] reflectance       input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param[in,out] reflectance_k     output Jacobians wrt surface BRDFs, optional
!! @param[in]     aer_opt_param     input aerosol optical parameters, optional
!! @param[in,out] aer_opt_param_k   output Jacobians wrt aerosol optical parameters, optional
!! @param[in]     cld_opt_param     input cloud optical parameters, optional
!! @param[in,out] cld_opt_param_k   output Jacobians wrt cloud optical parameters, optional
!! @param[in,out] traj              RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_k            RTTOV K internal state, can be initialised outside RTTOV, optional
!! @param[in,out] pccomp            output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_k          input gradient in PC scores, radiances or BTs for PC-RTTOV, optional
!! @param[in,out] profiles_k_pc     output PC score Jacobians wrt atmospheric profile and surface variables, optional
!! @param[in,out] profiles_k_rec    output reconstructed radiance/BT Jacobians wrt atmospheric profile and surface
!!                                    variables, optional
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
SUBROUTINE rttov_k( &
              errorstatus,      &
              chanprof,         &
              opts,             &
              profiles,         &
              profiles_k,       &
              coefs,            &
              transmission,     &
              transmission_k,   &
              radiance,         &
              radiance_k,       &
              radiance2,        &
              calcemis,         &
              emissivity,       &
              emissivity_k,     &
              calcrefl,         &
              reflectance,      &
              reflectance_k,    &
              aer_opt_param,    &
              aer_opt_param_k,  &
              cld_opt_param,    &
              cld_opt_param_k,  &
              traj,             &
              traj_k,           &
              pccomp,           &
              pccomp_k,         &
              profiles_k_pc,    &
              profiles_k_rec,   &
              channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
         rttov_options,      &
         rttov_coefs,        &
         rttov_pccomp,       &
         rttov_transmission, &
         rttov_profile,      &
         rttov_radiance,     &
         rttov_radiance2,    &
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
         rttov_traj_sta,  &
         rttov_traj_dyn
  USE rttov_htfrtc_interface_mod, ONLY : htfrtc_interface
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),       INTENT(OUT)                     :: errorstatus
  TYPE(rttov_chanprof),     INTENT(IN)                      :: chanprof(:)
  TYPE(rttov_options),      INTENT(IN)                      :: opts
  TYPE(rttov_profile),      INTENT(IN)                      :: profiles(:)
  TYPE(rttov_profile),      INTENT(INOUT)                   :: profiles_k(SIZE(chanprof))
  TYPE(rttov_coefs),        INTENT(IN)   , TARGET           :: coefs
  TYPE(rttov_transmission), INTENT(INOUT)                   :: transmission
  TYPE(rttov_transmission), INTENT(INOUT)                   :: transmission_k
  TYPE(rttov_radiance),     INTENT(INOUT)                   :: radiance
  TYPE(rttov_radiance),     INTENT(INOUT)                   :: radiance_k
  TYPE(rttov_radiance2),    INTENT(INOUT), OPTIONAL         :: radiance2
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity),   INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity),   INTENT(INOUT), OPTIONAL         :: emissivity_k(SIZE(chanprof))
  LOGICAL(KIND=jplm),       INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance),  INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),  INTENT(INOUT), OPTIONAL         :: reflectance_k(SIZE(chanprof))
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param),    INTENT(INOUT), OPTIONAL         :: aer_opt_param_k
  TYPE(rttov_opt_param),    INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_opt_param),    INTENT(INOUT), OPTIONAL         :: cld_opt_param_k
  TYPE(rttov_traj),         INTENT(INOUT), OPTIONAL, TARGET :: traj, traj_k        ! Target is needed here
  TYPE(rttov_pccomp),       INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_pccomp),       INTENT(INOUT), OPTIONAL         :: pccomp_k
  TYPE(rttov_profile),      INTENT(INOUT), OPTIONAL         :: profiles_k_pc(:)
  TYPE(rttov_profile),      INTENT(INOUT), OPTIONAL         :: profiles_k_rec(:)
  INTEGER(KIND=jpim),       INTENT(IN)   , OPTIONAL         :: channels_rec(:)
!INTF_END
#include "rttov_add_aux_prof.interface"
#include "rttov_add_opdp_path.interface"
#include "rttov_add_prof.interface"
#include "rttov_add_raytracing.interface"
#include "rttov_alloc_traj_dyn.interface"
#include "rttov_alloc_traj_sta.interface"
#include "rttov_apply_pc_aer_reg_lims_k.interface"
#include "rttov_apply_reg_limits_k.interface"
#include "rttov_calc_nearest_lev_k.interface"
#include "rttov_calcbt_ad.interface"
#include "rttov_calcbt_pc_ad.interface"
#include "rttov_calcemis_ir_k.interface"
#include "rttov_calcemis_mw_k.interface"
!#include "rttov_calcsatrefl_ad.interface"
#include "rttov_calcsurfrefl_k.interface"
#include "rttov_check_traj.interface"
#include "rttov_cloud_overlap_k.interface"
#include "rttov_convert_profile_units_k.interface"
#include "rttov_direct.interface"
#include "rttov_dom_k.interface"
#include "rttov_dom_setup_profile_k.interface"
#include "rttov_errorreport.interface"
#include "rttov_fresnel_k.interface"
#include "rttov_init_aux_prof.interface"
#include "rttov_init_auxrad_column.interface"
!#include "rttov_init_ircld.interface"
#include "rttov_init_opdp_path.interface"
#include "rttov_init_predictor.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_prof_internal.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_raytracing.interface"
#include "rttov_init_sunglint.interface"
#include "rttov_init_trans_scatt_ir.interface"
#include "rttov_init_transmission_aux.interface"
#include "rttov_intavg_chan_k.interface"
#include "rttov_intavg_prof_k.interface"
#include "rttov_integrate_k.interface"
#include "rttov_locpat_k.interface"
#include "rttov_mfasis_k.interface"
#include "rttov_mult_profiles_k.interface"
#include "rttov_mw_clw_absorption_k.interface"
#include "rttov_nlte_bias_correction_k.interface"
#include "rttov_opdep_13_ad.interface"
#include "rttov_opdep_9_k.interface"
#include "rttov_opdep_78_k.interface"
#include "rttov_opdpscattir_k.interface"
#include "rttov_pcscores_k.interface"
#include "rttov_pcscores_rec_k.interface"
#include "rttov_predictor_precalc_13_k.interface"
#include "rttov_predictor_precalc_789_k.interface"
#include "rttov_profaux_cldaer_k.interface"
#include "rttov_rayleigh_extinction_k.interface"
#include "rttov_reconstruct_k.interface"
#include "rttov_refsun_k.interface"
#include "rttov_setpredictors_13_k.interface"
#include "rttov_setpredictors_789_k.interface"
#include "rttov_transmit_solar_k.interface"
#include "rttov_transmit_k.interface"

  INTEGER(KIND=jpim) :: err
  INTEGER(KIND=jpim) :: nlevels, nprofiles, nchanprof, npcscores, nchannels_rec
  LOGICAL(KIND=jplm) :: sensor_mw, sensor_ir

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_k
  TYPE(rttov_traj), POINTER :: traj0_k
  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_k_dyn
  LOGICAL(KIND=jplm)   :: ltraj_k_dyn_dealloc

  REAL(KIND=jprb), ALLOCATABLE :: total_k_pc (:,:,:), pcscores_k(:,:,:)
  REAL(KIND=jprb) :: ZHOOK_HANDLE
!- End of header ------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 0_jpim, ZHOOK_HANDLE)

  errorstatus = errorstatus_success

  IF (opts%htfrtc_opts%htfrtc) THEN
    CALL htfrtc_interface(err, coefs, opts, profiles, pccomp, &
                          calcemis, emissivity,               &
                          emissivity_k = emissivity_k,        &
                          profiles_k_pc = profiles_k_pc,      &
                          profiles_k_rec = profiles_k_rec)
    errorstatus = err
    THROWM(err.NE.0, "Error in HTFRTC module")
    ! Since HTFRTC isn't RTTOV, skip to the end of rttov_k and leave
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

  ltraj_k_dyn_dealloc  = .FALSE.
  traj0_dyn%from_tladk = .TRUE.

  IF (opts%rt_ir%pc%addpc) THEN
    IF (.NOT. PRESENT(pccomp_k)) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "pccomp_k argument is mandatory for PC-RTTOV")
    ENDIF
    IF (opts%rt_ir%pc%addradrec) THEN
      IF (.NOT. PRESENT(profiles_k_rec)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "profiles_k_rec argument is mandatory for PC-RTTOV with addradrec true")
      ENDIF
    ELSE
      IF (.NOT. PRESENT(profiles_k_pc)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "profiles_k_pc argument is mandatory for PC-RTTOV with addradrec false")
      ENDIF
    ENDIF
    npcscores = opts%rt_ir%pc%npcscores * nprofiles
  ENDIF

  NULLIFY (traj0, traj0_k)
  CALL rttov_check_traj( &
          err,               &
          nprofiles,         &
          nchanprof,         &
          opts,              &
          nlevels,           &
          coefs,             &
          1_jpim,            &
          traj0 = traj0,     &
          traj0_k = traj0_k, &
          traj1 = traj1,     &
          traj1_k = traj1_k, &
          traj2 = traj,      &
          traj2_k = traj_k)
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

  radiance_k%plane_parallel = traj0_sta%plane_parallel
  radiance_k%quality        = radiance%quality

  ! If no channels are active (thermal and solar flags false for all channels) then
  ! skip to the end. In particular this can happen when using the parallel interface.
  IF (.NOT. (traj0_sta%dothermal .OR. traj0_sta%dosolar)) GOTO 998

  IF (traj0_sta%dothermal .AND. .NOT. PRESENT(emissivity_k)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "emissivity_k parameter required")
  END IF
  IF (traj0_sta%dosolar .AND. .NOT. PRESENT(reflectance_k)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "reflectance_k parameter required")
  END IF

  CALL rttov_alloc_traj_dyn(err, traj0_k_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
                            traj0_dyn%ncolumns, traj0_sta%dom_nstreams, &
                            traj0_sta%thermal, traj0_sta%solar, &
                            traj0_sta%dothermal, traj0_sta%do_mfasis, 1_jpim, traj0_dyn)
  THROW(err.NE.0)
  ltraj_k_dyn_dealloc = .TRUE.

!----------------
! K matrix
!----------------

!---------------------------------------------
! Initialize K variables
!---------------------------------------------
  ! xcolclr and xcol are the only ircld variables used outside rttov_cloud_overlap:
  ! if traj_k is not passed in by the user some time can be saved
  ! by avoiding unnecessary initialisations in traj0_k%ircld
!   CALL rttov_init_ircld(traj0_k%ircld)
  traj0_k%ircld%xcolclr = 0._jprb
  IF (opts%rt_ir%addclouds) THEN
    traj0_k%ircld%xcol = 0._jprb
  ENDIF
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_init_trans_scatt_ir(traj0_k%transmission_scatt_ir)
    CALL rttov_init_trans_scatt_ir(traj0_k_dyn%transmission_scatt_ir_dyn)
  ENDIF

  IF (traj0_sta%do_opdep_calc) THEN
    CALL rttov_init_opdp_path(opts, traj0_k%opdp_path)
    CALL rttov_init_opdp_path(opts, traj0_k%opdp_path_coef)
! DARFIX: TEST removing init for rttov9 predictors
    IF (coefs%coef%fmv_model_ver >= 9) &
      CALL rttov_init_predictor(traj0_sta%dosolar, traj0_k%predictors)

    CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_k%raytracing_coef)
    CALL rttov_init_prof(traj0_k%profiles_coef)
    CALL rttov_init_aux_prof_coef(traj0_k%aux_prof_coef)
  ENDIF

  CALL rttov_init_raytracing(traj0_sta%dosolar, traj0_k%raytracing)
  CALL rttov_init_prof_internal(traj0_k%profiles_int)
  CALL rttov_init_aux_prof(traj0_k%aux_prof)

  IF (traj0_sta%dothermal .OR. .NOT. traj0_sta%do_mfasis) THEN
    CALL rttov_init_transmission_aux(opts, traj0_k_dyn%transmission_aux)
    CALL rttov_init_auxrad_column(traj0_k_dyn%auxrad_column)
  ENDIF

  traj0_k%diffuse_refl(:) = 0._jprb
  IF (traj0_sta%dosolar) THEN
    traj0_k%fresnrefl(:) = 0._jprb
    CALL rttov_init_sunglint(traj0_k%sunglint)
  ENDIF

  traj0_k%auxrad%surfair(:) = 0._jprb
  traj0_k%auxrad%skin(:)    = 0._jprb
  traj0_k%auxrad%air(:,:)   = 0._jprb

  ! Copy units from direct model profile (not strictly essential for K code)
  profiles_k(:)%gas_units = profiles(1)%gas_units
  profiles_k(:)%mmr_cldaer = profiles(1)%mmr_cldaer

!--------------------------------------------------
! Set up radiance/BT inputs
!--------------------------------------------------
  IF (opts%rt_ir%pc%addpc) THEN

    CALL rttov_init_rad(radiance_k)
    radiance_k%total = 1._jprb

  ELSE

    IF (opts%rt_all%switchrad) THEN
      IF (traj0_sta%dothermal) &
        CALL rttov_calcbt_ad(chanprof, coefs%coef, traj0_sta%thermal, radiance, radiance_k)
!       Input K perturbation is always in radiance for pure-solar channels, never in reflectance
!       IF (traj0_sta%do_solar_opdep_calc) &
!         CALL rttov_calcsatrefl_ad(chanprof, profiles, traj0_sta%solar_spec_esd, &
!                                   traj0_sta%thermal, traj0_sta%solar, radiance_k)
    ELSE
! AJGDB this prevents RTTOV-SCATT using the clear radiance as adjoint input
!         radiance_k%clear(:) = 0._jprb
    ENDIF

  ENDIF

!---------------------------------------------------------
! Do NLTE bias correction for hyperspectral instruments
!---------------------------------------------------------
  IF (coefs%coef%nltecoef) THEN
    IF (opts%rt_ir%do_nlte_correction .AND. coefs%coef%id_sensor == sensor_id_hi) THEN
      CALL rttov_nlte_bias_correction_k(coefs%coef, profiles, profiles_k, chanprof, radiance_k)
    ENDIF
  ENDIF

!---------------------------------------------------------
! Solve the radiative transfer equation
!---------------------------------------------------------
  IF (traj0_sta%do_mfasis) THEN
    CALL rttov_mfasis_k( &
            err,                           &
            chanprof,                      &
            traj0_sta%solar,               &
            opts,                          &
            profiles,                      &
            profiles_k,                    &
            traj0%profiles_int,            &
            traj0_k%profiles_int,          &
            coefs,                         &
            traj0%ircld,                   &
            traj0_k%ircld,                 &
            traj0%aux_prof,                &
            traj0_k%aux_prof,              &
            reflectance,                   &
            reflectance_k,                 &
            traj0_sta%solar_spec_esd,      &
            traj0%transmission_scatt_ir,   &
            traj0_k%transmission_scatt_ir, &
            traj0_dyn%mfasis_refl,         &
            radiance_k)
    THROW(err.NE.0)
  ENDIF

  IF (traj0_sta%dothermal .OR. .NOT. traj0_sta%do_mfasis) THEN
    IF (traj0_sta%dothermal .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%ir_scatt_model == ir_scatt_dom) THEN
      CALL rttov_dom_k(                                     &
            err,                                            &
            chanprof,                                       &
            .FALSE._jplm,                                   &  ! FALSE => thermal source term
            traj0_sta%thermal,                              &
            traj0_sta%dom_nstreams,                         &
            profiles,                                       &
            traj0_dyn%profiles_dom_thermal,                 &
            traj0_k_dyn%profiles_dom_thermal,               &
            MAXVAL(traj0_dyn%profiles_dom_thermal%nlayers), &
            traj0%auxrad,                                   &
            traj0_k%auxrad,                                 &
            traj0%transmission_scatt_ir,                    &
            traj0_k%transmission_scatt_ir,                  &
            traj0_dyn%transmission_scatt_ir_dyn,            &
            traj0_k_dyn%transmission_scatt_ir_dyn,          &
            emissivity,                                     &
            emissivity_k,                                   &
            reflectance,                                    &
            reflectance_k,                                  &
            traj0%diffuse_refl,                             &
            traj0_k%diffuse_refl,                           &
            traj0_sta%solar_spec_esd,                       &
            traj0%raytracing,                               &
            traj0%ircld,                                    &
            traj0_k%ircld,                                  &
            radiance_k,                                     &
            traj0_dyn%dom_state_thermal)
      THROW(err.NE.0)

      CALL rttov_dom_setup_profile_k(               &
            err,                                    &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .FALSE._jplm,                           &  ! FALSE => thermal source term
            traj0_sta%thermal,                      &
            .FALSE._jplm,                           &  ! FALSE => do_rayleigh_dom
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0_k%aux_prof,                       &
            traj0%opdp_path%atm_level,              &
            traj0_k%opdp_path%atm_level,            &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_k%transmission_scatt_ir,          &
            traj0_dyn%profiles_dom_thermal,         &
            traj0_k_dyn%profiles_dom_thermal)
      THROW(err.NE.0)
    ENDIF

    IF (traj0_sta%dosolar .AND. (opts%rt_ir%addclouds .OR. opts%rt_ir%addaerosl) .AND. &
        opts%rt_ir%vis_scatt_model == vis_scatt_dom) THEN
      CALL rttov_dom_k(                                   &
            err,                                          &
            chanprof,                                     &
            .TRUE._jplm,                                  &  ! TRUE => solar source term
            traj0_sta%solar,                              &
            traj0_sta%dom_nstreams,                       &
            profiles,                                     &
            traj0_dyn%profiles_dom_solar,                 &
            traj0_k_dyn%profiles_dom_solar,               &
            MAXVAL(traj0_dyn%profiles_dom_solar%nlayers), &
            traj0%auxrad,                                 &
            traj0_k%auxrad,                               &
            traj0%transmission_scatt_ir,                  &
            traj0_k%transmission_scatt_ir,                &
            traj0_dyn%transmission_scatt_ir_dyn,          &
            traj0_k_dyn%transmission_scatt_ir_dyn,        &
            emissivity,                                   &
            emissivity_k,                                 &
            reflectance,                                  &
            reflectance_k,                                &
            traj0%diffuse_refl,                           &
            traj0_k%diffuse_refl,                         &
            traj0_sta%solar_spec_esd,                     &
            traj0%raytracing,                             &
            traj0%ircld,                                  &
            traj0_k%ircld,                                &
            radiance_k,                                   &
            traj0_dyn%dom_state_solar)
      THROW(err.NE.0)

      CALL rttov_dom_setup_profile_k(               &
            err,                                    &
            opts,                                   &
            coefs%coef,                             &
            chanprof,                               &
            .TRUE._jplm,                            &  ! TRUE => solar source term
            traj0_sta%solar,                        &
            traj0_sta%do_rayleigh_dom,              &
            profiles(1)%nlayers,                    &
            traj0%aux_prof,                         &
            traj0_k%aux_prof,                       &
            traj0%opdp_path%sun_level_path2,        &
            traj0_k%opdp_path%sun_level_path2,      &
            traj0_sta%angles,                       &
            traj0%ircld,                            &
            traj0%transmission_scatt_ir,            &
            traj0_k%transmission_scatt_ir,          &
            traj0_dyn%profiles_dom_solar,           &
            traj0_k_dyn%profiles_dom_solar)
      THROW(err.NE.0)
    ENDIF

    CALL rttov_integrate_k( &
            opts,                                     &
            traj0_dyn%ncolumns,                       &
            chanprof,                                 &
            emissivity,                               &
            emissivity_k,                             &
            reflectance,                              &
            reflectance_k,                            &
            traj0_sta%refl_norm,                      &
            traj0%diffuse_refl,                       &
            traj0_k%diffuse_refl,                     &
            traj0_sta%do_lambertian,                  &
            traj0_sta%do_mfasis,                      &
            traj0_sta%thermal,                        &
            traj0_sta%dothermal,                      &
            traj0_sta%solar,                          &
            traj0_sta%dosolar,                        &
            traj0_sta%do_rayleigh_ss,                 &
            traj0_sta%solar_spec_esd,                 &
            traj0_dyn%transmission_aux,               &
            traj0_k_dyn%transmission_aux,             &
            traj0%transmission_scatt_ir,              &
            traj0_k%transmission_scatt_ir,            &
            profiles,                                 &
            profiles_k,                               &
            traj0%profiles_int,                       &
            traj0_k%profiles_int,                     &
            traj0%aux_prof,                           &
            traj0_k%aux_prof,                         &
            coefs%coef,                               &
            traj0%raytracing,                         &
            traj0_k%raytracing,                       &
            traj0%ircld,                              &
            traj0_k%ircld,                            &
            radiance,                                 &
            traj0%auxrad,                             &
            traj0_k%auxrad,                           &
            traj0_dyn%auxrad_column,                  &
            traj0_k_dyn%auxrad_column,                &
            radiance_k)
  ENDIF ! dothermal .OR. .NOT. do_mfasis

!-------------------------------------------------------
! Calculate surface reflectance
!-------------------------------------------------------
  IF (PRESENT(reflectance_k)) THEN
    traj0_k%diffuse_refl(1:nchanprof) = traj0_k%diffuse_refl(1:nchanprof) + &
                                        reflectance_k(1:nchanprof)%diffuse_refl_out * pi_r
  ENDIF

  IF (traj0_sta%dosolar) THEN

    CALL rttov_calcsurfrefl_k( &
            coefs%coef,          &
            profiles,            &
            traj0%sunglint,      &
            traj0_k%sunglint,    &
            traj0%fresnrefl,     &
            traj0_k%fresnrefl,   &
            traj0_sta%solar,     &
            chanprof,            &
            traj0_sta%refl_norm, &
            calcrefl,            &
            emissivity_k,        &
            reflectance,         &
            reflectance_k,       &
            traj0_k%diffuse_refl)

    IF (ANY(calcrefl)) THEN
      CALL rttov_fresnel_k( &
              chanprof,          &
              calcrefl,          &
              profiles,          &
              traj0_sta%solar,   &
              coefs%coef,        &
              traj0%sunglint,    &
              traj0_k%sunglint,  &
              traj0_k%fresnrefl)

      CALL rttov_refsun_k( &
              opts,               &
              chanprof,           &
              traj0_sta%solar,    &
              calcrefl,           &
              profiles,           &
              profiles_k,         &
              coefs%coef,         &
              traj0%aux_prof,     &
              traj0%sunglint,     &
              traj0_k%sunglint,   &
              traj0%raytracing,   &
              traj0_k%raytracing)
    ENDIF

    reflectance_k%refl_in = reflectance_k%refl_in + reflectance_k%refl_out

  ENDIF

!-------------------------------------------------
! Calculate surface emissivity
!-------------------------------------------------
  IF (traj0_sta%dothermal) THEN

    IF (ANY(calcemis)) THEN
      ! Calculate surface emissivity and traj0%diffuse_refl for selected channels

      IF (sensor_mw) THEN
        CALL rttov_calcemis_mw_k( &
                opts,                           &
                profiles,                       &
                profiles_k,                     &
                traj0_sta%angles,               &
                coefs%coef,                     &
                chanprof,                       &
                traj0_dyn%transmission_aux,     &
                traj0_k_dyn%transmission_aux,   &
                calcemis,                       &
                emissivity_k%emis_out,          &
                traj0_k%diffuse_refl)
      ELSE
        WHERE (traj0_sta%thermal .AND. calcemis)
          emissivity_k(:)%emis_out = emissivity_k(:)%emis_out - traj0_k%diffuse_refl(:)
        ENDWHERE

        CALL rttov_calcemis_ir_k( &
                err,               &
                opts,              &
                chanprof,          &
                profiles,          &
                profiles_k,        &
                coefs,             &
                traj0_sta%thermal, &
                calcemis,          &
                emissivity_k%emis_out)
        THROWM(err.NE.0, "calcemis_ir_k")
      ENDIF

    ENDIF

    WHERE (traj0_sta%thermal .AND. .NOT. calcemis)
      emissivity_k%emis_out = emissivity_k%emis_out - traj0_k%diffuse_refl
    ENDWHERE

    emissivity_k%emis_in = emissivity_k%emis_in + emissivity_k%emis_out

  ENDIF

  ! traj0_k%diffuse_refl(:) = 0._jprb

!-------------------------------------------------
! Calculate transmittances
!-------------------------------------------------
  IF (traj0_sta%do_solar_opdep_calc) THEN
    CALL rttov_transmit_solar_k( &
            opts,                                     &
            profiles(1)%nlayers,                      &
            chanprof,                                 &
            traj0_sta%solar,                          &
            traj0%aux_prof,                           &
            traj0_k%aux_prof,                         &
            coefs%coef,                               &
            traj0%raytracing,                         &
            traj0_k%raytracing,                       &
            traj0%ircld,                              &
            traj0_k%opdp_path,                        &
            traj0_sta%solar_path2,                    &
            traj0_sta%solar_path1,                    &
            transmission_k,                           &
            traj0_dyn%transmission_aux,               &
            traj0_k_dyn%transmission_aux,             &
            traj0%transmission_scatt_ir,              &
            traj0_k%transmission_scatt_ir,            &
            traj0_dyn%transmission_scatt_ir_dyn,      &
            traj0_k_dyn%transmission_scatt_ir_dyn)
  ENDIF

  IF (traj0_sta%dothermal) THEN
    CALL rttov_transmit_k( &
            opts,                                       &
            traj0_sta%do_lambertian,                    &
            profiles(1)%nlayers,                        &
            chanprof,                                   &
            traj0_sta%thermal,                          &
            traj0%aux_prof,                             &
            traj0_k%aux_prof,                           &
            coefs%coef,                                 &
            traj0%ircld,                                &
            traj0_sta%angles,                           &
            traj0_k%opdp_path%atm_level,                &
            traj0_sta%thermal_path1%od_level,           &
            transmission_k,                             &
            traj0_dyn%transmission_aux%thermal_path1,   &
            traj0_k_dyn%transmission_aux%thermal_path1, &
            traj0_k%transmission_scatt_ir,              &
            traj0_dyn%transmission_scatt_ir_dyn,        &
            traj0_k_dyn%transmission_scatt_ir_dyn,      &
            traj0_sta%thermal_path1%tau_surf,           &
            traj0_sta%thermal_path1%tau_level)
  ENDIF

!----------------------------------------------------------------------------
! Calculate optical depths of aerosols and/or clouds
!----------------------------------------------------------------------------
  IF (opts%rt_ir%addaerosl .OR. opts%rt_ir%addclouds) THEN
    CALL rttov_opdpscattir_k( &
            profiles(1)%nlayers,                    &
            chanprof,                               &
            opts,                                   &
            traj0_sta%dom_nstreams,                 &
            traj0%aux_prof,                         &
            traj0_k%aux_prof,                       &
            traj0%ircld,                            &
            profiles,                               &
            profiles_k,                             &
            traj0%profiles_int,                     &
            traj0_k%profiles_int,                   &
            aer_opt_param,                          &
            aer_opt_param_k,                        &
            cld_opt_param,                          &
            cld_opt_param_k,                        &
            traj0_sta%dothermal,                    &
            traj0_sta%thermal,                      &
            traj0_sta%dosolar,                      &
            traj0_sta%solar,                        &
            traj0_sta%do_rayleigh_dom,              &
            coefs%coef,                             &
            coefs%coef_scatt,                       &
            coefs%coef_mfasis_cld,                  &
            traj0_sta%angles,                       &
            traj0%raytracing,                       &
            traj0_k%raytracing,                     &
            traj0%transmission_scatt_ir,            &
            traj0_k%transmission_scatt_ir,          &
            traj0_dyn%transmission_scatt_ir_dyn,    &
            traj0_k_dyn%transmission_scatt_ir_dyn)

    CALL rttov_profaux_cldaer_k( &
            opts,                  &
            chanprof,              &
            profiles,              &
            profiles_k,            &
            traj0%profiles_int,    &
            traj0_k%profiles_int,  &
            traj0%aux_prof,        &
            traj0_k%aux_prof)
  ENDIF

!-----------------------------------------------------------------------
! Apply PC-RTTOV aerosol regression limits
!-----------------------------------------------------------------------
  IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addaerosl) THEN
    CALL rttov_apply_pc_aer_reg_lims_k( &
          opts,                 &
          chanprof,             &
          coefs%coef_pccomp,    &
          profiles,             &
          traj0_k%profiles_int, &
          traj0_sta%pc_aer_ref, &
          traj0_sta%pc_aer_min, &
          traj0_sta%pc_aer_max)
  ENDIF

!------------------------------------------------------
! Calculate the number of cloud columns and
! the cloud distribution in each column
!------------------------------------------------------
  IF (opts%rt_ir%addclouds) THEN
    CALL rttov_cloud_overlap_k( &
            opts%rt_ir,           &
            chanprof,             &
            profiles,             &
            profiles_k,           &
            traj0%profiles_int,   &
            traj0_k%profiles_int, &
            traj0%ircld,          &
            traj0_k%ircld)
  ENDIF

!--------------------------------------------------------------------------
! Rayleigh extinction optical depths
!--------------------------------------------------------------------------
  IF (traj0_sta%do_rayleigh_param) THEN
    CALL rttov_rayleigh_extinction_k( &
            opts,                              &
            chanprof,                          &
            traj0_sta%do_rayleigh_dom,         &
            traj0_sta%solar,                   &
            profiles,                          &
            profiles_k,                        &
            coefs%coef,                        &
            traj0%raytracing,                  &
            traj0_k%raytracing,                &
            traj0_k%opdp_path%sun_level_path2, &
            traj0_k%transmission_scatt_ir)
  ENDIF

!--------------------------------------------------------------------------
! MW CLW absorption optical depths
!--------------------------------------------------------------------------
  IF (sensor_mw .AND. opts%rt_mw%clw_data) THEN
    CALL rttov_mw_clw_absorption_k( &
            opts,               &
            coefs%coef,         &
            chanprof,           &
            traj0%raytracing,   &
            traj0_k%raytracing, &
            profiles,           &
            profiles_k,         &
            traj0%aux_prof,     &
            traj0_k%aux_prof,   &
            traj0_k%opdp_path)
  ENDIF

!-----------------------------------------------------------------------
! RTTOV optical depth calculation: on coef levels
!-----------------------------------------------------------------------

  IF (traj0_sta%do_opdep_calc) THEN
!--------------------------------------------------------------------------
! Move top level to space boundary (opdep=0) if spacetop flag is set
!--------------------------------------------------------------------------
    IF (opts%interpolation%spacetop) THEN
      traj0_k%opdp_path%atm_level(1,:) = 0._jprb
      IF (traj0_sta%do_solar_opdep_calc) THEN
        traj0_k%opdp_path%sun_level_path2(1,:) = 0._jprb
      ENDIF
    ENDIF

!--------------------------------------------------------------------------
! Interpolator second  call - optical depths from coef levs to user levs
!--------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_chan_k( &
              opts,                           &
              traj0_sta%thermal,              &
              traj0_sta%solar,                &
              coefs%coef%nlevels,             &
              nlevels,                        &
              chanprof,                       &
              traj0%profiles_coef,            &
              profiles,                       &
              profiles_k,                     &
              traj0%opdp_path_coef,           &
              traj0_k%opdp_path_coef,         &
              traj0_k%opdp_path)
    ELSE
      CALL rttov_add_opdp_path(opts, traj0_k%opdp_path_coef, traj0_k%opdp_path_coef, traj0_k%opdp_path)
    ENDIF

!--------------------------------------------------------------------
! Predict atmospheric (emissive) and solar optical depths - coef levs
!--------------------------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for solar channels
        CALL rttov_opdep_13_ad( &
                coefs%coef%nlayers,                     &
                chanprof,                               &
                traj0_sta%solar,                        &
                traj0_k%predictors%path2,               &
                coefs%coef,                             &
                coefs%coef%solar,                       &
                coefs%coef%solar_corr,                  &
                traj0_k%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref)
      ENDIF
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for thermal channels
        CALL rttov_opdep_13_ad( &
                coefs%coef%nlayers,               &
                chanprof,                         &
                traj0_sta%thermal,                &
                traj0_k%predictors%path1,         &
                coefs%coef,                       &
                coefs%coef%thermal,               &
                coefs%coef%thermal_corr,          &
                traj0_k%opdp_path_coef%atm_level, &
                traj0_sta%thermal_path1%opdp_ref)
      ENDIF
    ELSEIF (coefs%coef%fmv_model_ver == 9) THEN
      IF (traj0_sta%do_solar_opdep_calc) THEN
        ! Calculate solar path2 optical depths for solar channels
        CALL rttov_opdep_9_k( &
                coefs%coef%nlayers,                     &
                chanprof,                               &
                traj0_sta%solar,                        &
                traj0%predictors%path2,                 &
                traj0_k%predictors%path2,               &
                coefs%coef,                             &
                coefs%coef%solar,                       &
                traj0_k%opdp_path_coef%sun_level_path2, &
                traj0_sta%solar_path2%opdp_ref_coef)
      ENDIF
      IF (traj0_sta%dothermal) THEN
        ! Calculate thermal path1 optical depths for thermal channels
        CALL rttov_opdep_9_k( &
                coefs%coef%nlayers,                    &
                chanprof,                              &
                traj0_sta%thermal,                     &
                traj0%predictors%path1,                &
                traj0_k%predictors%path1,              &
                coefs%coef,                            &
                coefs%coef%thermal,                    &
                traj0_k%opdp_path_coef%atm_level,      &
                traj0_sta%thermal_path1%opdp_ref_coef)
      ENDIF
    ELSE
      CALL rttov_opdep_78_k( &
              coefs%coef%nlayers,          &
              chanprof,                    &
              traj0_k%predictors%path1,    &
              coefs%coef,                  &
              coefs%coef%thermal,          &
              traj0_k%opdp_path_coef,      &
              traj0_sta%thermal_path1%opdp_ref_coef)
    ENDIF

!---------------------------------------------------
! Calculate predictors - coef levels
!---------------------------------------------------
    IF (coefs%coef%fmv_model_ver == 13) THEN
      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_13_k(  &
          opts,                          &
          chanprof,                      &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0%predictors%path1,        &
          traj0_k%predictors%path1,      &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef,       &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_13_k(  &
          opts,                          &
          chanprof,                      &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0%predictors%path2,        &
          traj0_k%predictors%path2,      &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef,       &
          raypath = 2_jpim) 
      ENDIF

      CALL rttov_predictor_precalc_13_k( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0_sta%plane_parallel,      &
          chanprof,                      &
          traj0%profiles_coef,           &
          traj0_k%profiles_coef,         &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef)
    ELSE
      IF (traj0_sta%dothermal) THEN
        CALL rttov_setpredictors_789_k(  &
          opts,                          &
          chanprof,                      &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0%predictors%path1,        &
          traj0_k%predictors%path1,      &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef,       &
          raypath = 1_jpim) 
      ENDIF
      IF (traj0_sta%do_solar_opdep_calc) THEN
        CALL rttov_setpredictors_789_k(  &
          opts,                          &
          chanprof,                      &
          traj0%profiles_coef,           &
          coefs%coef,                    &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0%predictors%path2,        &
          traj0_k%predictors%path2,      &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef,       &
          raypath = 2_jpim) 
      ENDIF

      CALL rttov_predictor_precalc_789_k( &
          opts,                          &
          traj0_sta%do_solar_opdep_calc, &
          traj0_sta%plane_parallel,      &
          chanprof,                      &
          traj0%profiles_coef,           &
          traj0_k%profiles_coef,         &
          coefs%coef,                    &
          coefs%coef_pccomp,             &
          traj0%aux_prof_coef,           &
          traj0_k%aux_prof_coef,         &
          traj0_sta%angles,              &
          traj0%raytracing_coef,         &
          traj0_k%raytracing_coef)
    ENDIF
!--------------------------------------------------------------------------------
! Calculate near surface level and local path geometry calculations - coef levels
!--------------------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_locpat_k( &
          opts,                     &
          traj0_sta%dosolar,        &
          traj0_sta%plane_parallel, &
          chanprof,                 &
          traj0%profiles_coef,      &
          traj0_k%profiles_coef,    &
          traj0%profiles_coef,      &
          traj0_k%profiles_coef,    &
          traj0%aux_prof_coef%s,    &
          traj0_sta%angles,         &
          traj0%raytracing_coef,    &
          traj0_k%raytracing_coef)
      CALL rttov_calc_nearest_lev_k(.FALSE._jplm, opts, chanprof, traj0%profiles_coef, &
            traj0_k%profiles_coef, traj0%aux_prof_coef%s, traj0_k%aux_prof_coef%s)
    ELSE
      CALL rttov_add_raytracing(traj0_sta%dosolar, traj0_k%raytracing, &
                                traj0_k%raytracing, traj0_k%raytracing_coef)
      CALL rttov_add_aux_prof(traj0_k%aux_prof%s, traj0_k%aux_prof%s, traj0_k%aux_prof_coef%s)
    ENDIF

!------------------------------------------------------------------
! Check input data is within suitable physical limits - coef levels
!------------------------------------------------------------------
    IF (opts%config%apply_reg_limits .OR. opts%interpolation%addinterp) THEN
      CALL rttov_apply_reg_limits_k( &
              opts,                         &
              chanprof,                     &
              profiles,                     &
              traj0_sta%profiles_coef_ref,  &
              traj0_k%profiles_coef,        &
              coefs%coef,                   &
              coefs%coef_pccomp)
    ENDIF

!-----------------------------------------------------------------------
! Interpolator first call - input profiles from user levs to coef levs
!-----------------------------------------------------------------------
    IF (opts%interpolation%addinterp) THEN
      CALL rttov_intavg_prof_k( &
              opts,                           &
              chanprof,                       &
              nlevels,                        &
              coefs%coef%nlevels,             &
              profiles,                       &
              profiles_k,                     &
              traj0%profiles_int,             &
              traj0_k%profiles_int,           &
              traj0%profiles_coef,            &
              traj0_k%profiles_coef,          &
              coefs%coef,                     &
              coefs%coef_pccomp)
    ELSE
      ! If interpolator is off copy profile variables
      CALL rttov_add_prof( &
              profiles_k,             &
              profiles_k,             &
              traj0_k%profiles_coef,  &
              larray = .TRUE._jplm,   &
              lscalar = .FALSE._jplm, &
              profiles_gas = traj0_k%profiles_int)
    ENDIF

    ! Copy non-interpolated variables (e.g. surface parameters)
    CALL rttov_add_prof( &
            profiles_k,            &
            profiles_k,            &
            traj0_k%profiles_coef, &
            larray = .FALSE._jplm, &
            lscalar = .TRUE._jplm, &
            profiles_gas = traj0_k%profiles_int)

  ENDIF ! do_opdep_calc

!------------------------------------------------------------------------
! Local path geometry calculations
!------------------------------------------------------------------------
  CALL rttov_locpat_k( &
      opts,                     &
      traj0_sta%dosolar,        &
      traj0_sta%plane_parallel, &
      chanprof,                 &
      profiles,                 &
      profiles_k,               &
      traj0%profiles_int,       &
      traj0_k%profiles_int,     &
      traj0%aux_prof%s,         &
      traj0_sta%angles,         &
      traj0%raytracing,         &
      traj0_k%raytracing)

!-----------------------------------------------------------------------
! Calculate near surface and cloud top levels
!-----------------------------------------------------------------------
  CALL rttov_calc_nearest_lev_k(sensor_ir, opts, chanprof, profiles, &
                                profiles_k, traj0%aux_prof%s, traj0_k%aux_prof%s)

!-----------------------------------------------------------------------
! Convert input profiles to units used internally
!-----------------------------------------------------------------------
  CALL rttov_convert_profile_units_k(opts, chanprof, coefs, profiles, profiles_k, &
                                     traj0%profiles_int, traj0_k%profiles_int)

!-----------------------------------------------------------------------
! For PC-RTTOV complete the K calculation
!-----------------------------------------------------------------------
  IF (opts%rt_ir%pc%addpc) THEN
    npcscores     = SIZE(traj0_sta%chanprof_pc)

    IF (opts%rt_ir%pc%addradrec) THEN
      nchannels_rec = SIZE(traj0_sta%chanprof_in)

      ALLOCATE ( &
        pcscores_k(npcscores / nprofiles, nchannels_rec / nprofiles, nprofiles),&
        total_k_pc(nchanprof / nprofiles, nchannels_rec / nprofiles, nprofiles),&
        STAT = err)
      THROWM(err.NE.0, "allocation pc k arrays")

      IF (opts%rt_all%switchrad) THEN
        CALL rttov_calcbt_pc_ad( &
            traj0_sta%chanprof_in,       &
            coefs%coef_pccomp,           &
            pccomp,                      &
            pccomp_k)
      ENDIF

      CALL rttov_reconstruct_k( &
          opts,                        &
          traj0_sta%chanprof_in,       &
          traj0_sta%chanprof_pc,       &
          pccomp_k,                    &
          pcscores_k,                  &
          coefs%coef_pccomp)

      CALL rttov_pcscores_rec_k(       &
          opts,                        &
          chanprof,                    &
          pcscores_k,                  &
          coefs%coef_pccomp,           &
          total_k_pc)

! DAR: Don't need to initialise profiles if they're being overwritten immediately
!     CALL rttov_init_prof(profiles_k_rec)
      CALL rttov_mult_profiles_k(profiles_k_rec, profiles_k, total_k_pc, opts, coefs%coef_pccomp)

      ! Copy units from direct model profile (not strictly essential for K code)
      profiles_k_rec(:)%gas_units = profiles(1)%gas_units
      profiles_k_rec(:)%mmr_cldaer = profiles(1)%mmr_cldaer

      DEALLOCATE (pcscores_k, total_k_pc, STAT = err)
      THROWM(err.NE.0, "deallocation pc k arrays")

    ELSE
      ALLOCATE ( &
        total_k_pc(nchanprof / nprofiles, npcscores / nprofiles, nprofiles), &
        STAT = err)
      THROWM(err.NE.0, "allocation pc k arrays")

      CALL rttov_pcscores_k( &
              opts,                        &
              chanprof,                    &
              traj0_sta%chanprof_pc,       &
              pccomp_k,                    &
              coefs%coef_pccomp,           &
              total_k_pc)

! DAR: Don't need to initialise profiles if they're being overwritten immediately
!     CALL rttov_init_prof(profiles_k_pc)
      CALL rttov_mult_profiles_k(profiles_k_pc, profiles_k, total_k_pc, opts, coefs%coef_pccomp)

      ! Copy units from direct model profile (not strictly essential for K code)
      profiles_k_pc(:)%gas_units = profiles(1)%gas_units
      profiles_k_pc(:)%mmr_cldaer = profiles(1)%mmr_cldaer

      DEALLOCATE (total_k_pc, STAT = err)
      THROWM(err.NE.0, "deallocation pc k arrays")

    ENDIF
  ENDIF

!---------------------
! Deallocate memory
!---------------------
998 CONTINUE ! from HTFRTC or if no channels are active
  CALL cleanup()

  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)
  CATCH_C
  errorstatus = err
  IF (LHOOK) CALL DR_HOOK('RTTOV_K', 1_jpim, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE cleanup()
    INTEGER(KIND=jpim) :: error

    IF (opts%htfrtc_opts%htfrtc) RETURN

    IF (ltraj_k_dyn_dealloc) THEN
      CALL rttov_alloc_traj_dyn(error, traj0_k_dyn, opts, coefs, nchanprof, profiles(1)%nlayers, &
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
              error,             &
              nprofiles,         &
              nchanprof,         &
              opts,              &
              nlevels,           &
              coefs,             &
              0_jpim,            &
              traj0 = traj0,     &
              traj0_k = traj0_k, &
              traj1 = traj1,     &
              traj1_k = traj1_k, &
              traj2 = traj,      &
              traj2_k = traj_k)
    ENDIF
  END SUBROUTINE cleanup

END SUBROUTINE rttov_k
