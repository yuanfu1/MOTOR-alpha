! Description:
!> @file
!!   Compute brute force Jacobian using direct model
!
!> @brief
!!   Compute brute force Jacobian using direct model
!!
!! @param[out]    err              status on exit
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     frequencies      RTTOV-SCATT frequency list
!! @param[in]     opts             options to configure the simulations
!! @param[in]     opts_scatt       RTTOV-SCATT options to configure the simulations
!! @param[in]     profiles         atmospheric profiles and surface variables
!! @param[in,out] profiles_k       computed profile Jacobian
!! @param[in]     cld_profiles     RTTOV-SCATT cloud profiles
!! @param[in,out] cld_profiles_k   computed RTTOV-SCATT cloud profile Jacobian
!! @param[in]     coefs            coefficients structure
!! @param[in]     coefs_scatt      RTTOV-SCTT coefficients structure
!! @param[in,out] transmission     computed direct model transmittances
!! @param[in,out] radiance         computed direct model radiances
!! @param[in,out] reflectivity     computed direct model radar reflectivities
!! @param[in]     calcemis         flags for internal RTTOV surface emissivity calculation
!! @param[in,out] emissivity       input/output surface emissivities
!! @param[in,out] emissivity_k     computed emissivity Jacobian
!! @param[in]     calcrefl         flags for internal RTTOV surface BRDF calculation
!! @param[in,out] reflectance      input/output surface BRDFs, input cloud top BRDF for simple cloud
!! @param[in,out] reflectance_k    computed BRDF Jacobian
!! @param[in]     aer_opt_param    input aerosol optical parameters
!! @param[in,out] aer_opt_param_k  computed aerosol optical parameter Jacobians
!! @param[in]     cld_opt_param    input cloud optical parameters
!! @param[in,out] cld_opt_param_k  computed cloud optical parameter Jacobians
!! @param[in,out] pccomp           computed direct model PC scores and radiances from PC-RTTOV
!! @param[in,out] profiles_k_pc    computed PC score profile Jacobian
!! @param[in,out] profiles_k_rec   computed reconstructed radiance/BT profile Jacobian
!! @param[in]     channels_rec     reconstructed radiance channels for PC-RTTOV simulations with addradrec true
!! @param[in]     nthreads         number of threads (if >0 the parallel interface is used which can be much faster)
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
SUBROUTINE rttov_k_bf( &
              err,              &
              chanprof,         &
              frequencies,      &
              opts,             &
              opts_scatt,       &
              profiles,         &
              profiles_k,       &
              cld_profiles,     &
              cld_profiles_k,   &
              coefs,            &
              coefs_scatt,      &
              transmission,     &
              radiance,         &
              reflectivity,     &
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
              pccomp,           &
              profiles_k_pc,    &
              profiles_k_rec,   &
              channels_rec,     &
              nthreads)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY :  &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_transmission,  &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance,   &
         rttov_opt_param,     &
         rttov_radiance,      &
         rttov_reflectivity,  &
         rttov_pccomp
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_traj
  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi
  USE rttov_chain, ONLY :  &
         chain,              &
         pchain,             &
         size_chain,         &
         get_pointer_chain,  &
         advance_chain
  USE rttov_test_k_mod, ONLY : &
         make_chain_profile,     &
         make_chain_cld_profile, &
         free_chain_profile,     &
         assign_chain_profile,   &
         radar_k_lev
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),        INTENT(OUT)           :: err
  TYPE(rttov_chanprof),      INTENT(IN)            :: chanprof(:)
  INTEGER(KIND=jpim),        INTENT(IN)            :: frequencies(:)
  TYPE(rttov_options),       INTENT(IN)            :: opts
  TYPE(rttov_options_scatt), INTENT(IN)            :: opts_scatt
  TYPE(rttov_profile),       INTENT(IN)   , TARGET :: profiles(:)
  TYPE(rttov_profile),       INTENT(INOUT), TARGET :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile_cloud), POINTER               :: cld_profiles(:)
  TYPE(rttov_profile_cloud), POINTER               :: cld_profiles_k(:)
  TYPE(rttov_coefs),         INTENT(IN)   , TARGET :: coefs
  TYPE(rttov_scatt_coef),    INTENT(IN)   , TARGET :: coefs_scatt
  TYPE(rttov_transmission),  INTENT(INOUT)         :: transmission
  TYPE(rttov_radiance),      INTENT(INOUT)         :: radiance
  TYPE(rttov_reflectivity),  INTENT(INOUT)         :: reflectivity
  LOGICAL(KIND=jplm),        INTENT(IN)   , OPTIONAL         :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity),    INTENT(INOUT), OPTIONAL         :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity),    INTENT(INOUT), OPTIONAL         :: emissivity_k(SIZE(chanprof))
  LOGICAL(KIND=jplm),        INTENT(IN)   , OPTIONAL         :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance),   INTENT(INOUT), OPTIONAL         :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),   INTENT(INOUT), OPTIONAL         :: reflectance_k(SIZE(chanprof))
  TYPE(rttov_opt_param),     INTENT(IN)   , OPTIONAL, TARGET :: aer_opt_param
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL, TARGET :: aer_opt_param_k
  TYPE(rttov_opt_param),     INTENT(IN)   , OPTIONAL, TARGET :: cld_opt_param
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL, TARGET :: cld_opt_param_k
  TYPE(rttov_pccomp),        INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_profile),       INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_pc(:)
  TYPE(rttov_profile),       INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_rec(:)
  INTEGER(KIND=jpim),        INTENT(IN)   , OPTIONAL         :: channels_rec(:)
  INTEGER(KIND=jpim),        INTENT(IN)                      :: nthreads
!INTF_END

#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_scatt.interface"
#include "rttov_parallel_scatt.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_alloc_pc_dimensions.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_scatt_prof.interface"
#include "rttov_copy_opt_param.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_make_profile_inc.interface"
#include "rttov_make_cld_profile_inc.interface"
#include "rttov_make_opt_param_inc.interface"
#include "rttov_errorreport.interface"

  REAL(KIND=jprb), PARAMETER :: demis = -0.01_jprb
  REAL(KIND=jprb), PARAMETER :: dspec = 0.01_jprb
  REAL(KIND=jprb), PARAMETER :: drefl = 0.01_jprb

  TYPE(rttov_traj)       :: traj
  TYPE(rttov_pccomp)     :: pccomp_bf
  TYPE(rttov_chanprof), POINTER :: chanprof_pc(:)
  TYPE(rttov_chanprof), POINTER :: chanprof_in(:)
  REAL(KIND=jprb),      POINTER :: w(:)
  REAL(KIND=jprb),      POINTER :: w_pc(:)
  REAL(KIND=jprb),      POINTER :: w_in(:)
  TYPE(rttov_profile),       TARGET :: profiles_bf(SIZE(profiles))
  TYPE(rttov_profile),       TARGET :: profiles_inc(SIZE(profiles))
  TYPE(rttov_profile_cloud), TARGET :: cld_profiles_bf (SIZE(profiles))
  TYPE(rttov_profile_cloud), TARGET :: cld_profiles_inc(SIZE(profiles))
  TYPE(rttov_opt_param), TARGET :: aer_opt_param_bf
  TYPE(rttov_opt_param), TARGET :: aer_opt_param_inc
  TYPE(rttov_opt_param), TARGET :: cld_opt_param_bf
  TYPE(rttov_opt_param), TARGET :: cld_opt_param_inc
  TYPE(rttov_opt_param), POINTER :: opt_param, opt_param_bf, opt_param_inc, opt_param_k
  LOGICAL(KIND=jplm)       :: calcemis_bf(SIZE(chanprof))
  TYPE(rttov_emissivity)   :: emissivity_bf(SIZE(chanprof))
  LOGICAL(KIND=jplm)       :: calcrefl_bf(SIZE(chanprof))
  TYPE(rttov_reflectance)  :: reflectance_bf(SIZE(chanprof))
  TYPE(rttov_transmission) :: transmission_bf
  TYPE(rttov_radiance)     :: radiance_bf
  TYPE(rttov_reflectivity), POINTER :: reflectivity_bf => NULL()
  REAL(KIND=jprb)          :: out_bf(SIZE(chanprof))
  TYPE(chain),  POINTER :: chain_profiles(:)
  TYPE(chain),  POINTER :: chain_profiles_bf(:)
  TYPE(chain),  POINTER :: chain_profiles_k(:)
  TYPE(chain),  POINTER :: chain_cld_profiles(:)
  TYPE(chain),  POINTER :: chain_cld_profiles_bf(:)
  TYPE(chain),  POINTER :: chain_cld_profiles_k(:)
  TYPE(chain),  POINTER :: chain_profiles_k_pc(:)
  TYPE(chain),  POINTER :: chain_profiles_k_rec(:)
  TYPE(chain),  POINTER :: chain_profiles_inc(:)
  TYPE(chain),  POINTER :: chain_cld_profiles_inc(:)
  TYPE(pchain), POINTER :: c(:), c_bf(:), c_k(:), c_inc(:)
  TYPE(pchain), POINTER :: c_cld(:), c_cld_bf(:), c_cld_k(:), c_cld_inc(:)
  TYPE(pchain), POINTER :: c_k_pc(:), c_k_in(:)
  INTEGER(KIND=jpim) :: iprof, ichan, i, aercld, var, lay
  REAL(KIND=jprb), POINTER :: x_profiles, x_profiles_inc, x_profiles_bf
  REAL(KIND=jprb), POINTER :: var_in(:), var_bf(:), var_inc(:), var_out(:)
  INTEGER(KIND=jpim) :: nsize
  REAL   (KIND=jprb) :: emisout(SIZE(chanprof))
  REAL   (KIND=jprb) :: reflout(SIZE(chanprof))
  REAL   (KIND=jprb) :: diffusereflout(SIZE(chanprof))
  REAL   (KIND=jprb) :: inc(SIZE(profiles))
  INTEGER(KIND=jpim) :: nprofiles, nchannels, npcscores, nchannels_rec, nlevels, nlayers, nphangle
  LOGICAL(KIND=jplm) :: do_pc, do_radrec, do_rttovscatt, do_radar

  TRY

  ! ---------------------------------------------------------------------------
  ! Set up, allocations, etc
  ! ---------------------------------------------------------------------------

  nprofiles = SIZE(profiles)
  nchannels = SIZE(chanprof)
  nlevels = profiles(1)%nlevels
  nlayers = profiles(1)%nlevels - 1

  do_rttovscatt = ASSOCIATED(cld_profiles)
  do_radar = do_rttovscatt .AND. ASSOCIATED(reflectivity%zef)

  do_pc = opts%rt_ir%pc%addpc .OR. opts%htfrtc_opts%htfrtc
  do_radrec = (opts%rt_ir%pc%addpc .AND. opts%rt_ir%pc%addradrec) .OR. &
              (opts%htfrtc_opts%htfrtc .AND. opts%htfrtc_opts%reconstruct)

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_bf,         &
          profiles(1)%nlevels, &
          opts,                &
          1_jpim,              &
          coefs = coefs,       &
          init = .TRUE._jplm)
  THROW(err.NE.0)

  CALL rttov_copy_prof(profiles_bf, profiles)

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_inc,        &
          profiles(1)%nlevels, &
          opts,                &
          1_jpim,              &
          coefs = coefs,       &
          init = .TRUE._jplm)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL rttov_alloc_scatt_prof( &
          err,                          &
          nprofiles,                    &
          cld_profiles_bf,              &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          1_jpim,                       &
          init=.TRUE._jplm)
    THROW(err.NE.0)

    CALL rttov_copy_scatt_prof(cld_profiles_bf, cld_profiles)

    CALL rttov_alloc_scatt_prof( &
          err,                          &
          nprofiles,                    &
          cld_profiles_inc,             &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          1_jpim,                       &
          init=.TRUE._jplm)
    THROW(err.NE.0)

    IF (do_radar) THEN
      ALLOCATE(reflectivity_bf)
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_bf,         &
            nlevels,                 &
            1_jpim,                  &
            init = .TRUE._jplm)
      THROW(err.NE.0)
    ENDIF
  ENDIF

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_bf,         &
          profiles(1)%nlevels, &
          1_jpim, init=.TRUE._jplm)
  THROW(err.NE.0)

  CALL rttov_alloc_transmission(err, transmission_bf, profiles(1)%nlevels, nchannels, 1_jpim)
  THROW(err.NE.0)

  IF (opts%rt_ir%user_aer_opt_param) THEN
    nphangle = SIZE(aer_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, aer_opt_param_bf, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)

    CALL rttov_copy_opt_param(aer_opt_param_bf, aer_opt_param)
    CALL rttov_init_opt_param(err, opts, aer_opt_param_bf)

    CALL rttov_alloc_opt_param(err, aer_opt_param_inc, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)

    CALL rttov_make_opt_param_inc(aer_opt_param_inc, aer_opt_param)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    nphangle = SIZE(cld_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, cld_opt_param_bf, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)

    CALL rttov_copy_opt_param(cld_opt_param_bf, cld_opt_param)
    CALL rttov_init_opt_param(err, opts, cld_opt_param_bf)

    CALL rttov_alloc_opt_param(err, cld_opt_param_inc, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)

    CALL rttov_make_opt_param_inc(cld_opt_param_inc, cld_opt_param)
  ENDIF

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                 &
          err,                  nprofiles, &
          nchannels,            opts,      &
          profiles(1)%nlevels,  coefs,     &
          1_jpim,               traj = traj)
    THROW(err.NE.0)
  ENDIF

  ! ---------------------------------------------------------------------------
  ! Generate the increments in all profile and cloud profile variables
  ! ---------------------------------------------------------------------------
  CALL rttov_make_profile_inc(profiles_inc, profiles, opts)
  IF (do_rttovscatt) CALL rttov_make_cld_profile_inc(cld_profiles_inc, cld_profiles)

  ! ---------------------------------------------------------------------------
  ! Initial run of direct model for given profile
  ! ---------------------------------------------------------------------------
  IF (nthreads <= 0) THEN
    IF (do_rttovscatt) THEN
      IF (do_radar) THEN
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance,       &
            reflectivity = reflectivity)
      ELSE
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance)
      ENDIF
    ELSE
      CALL rttov_direct(                     &
              err,                           &
              chanprof,                      &
              opts,                          &
              profiles,                      &
              coefs,                         &
              transmission,                  &
              radiance,                      &
              calcemis = calcemis,           &
              emissivity = emissivity,       &
              calcrefl = calcrefl,           &
              reflectance = reflectance,     &
              aer_opt_param = aer_opt_param, &
              cld_opt_param = cld_opt_param, &
              traj         = traj,           &
              pccomp       = pccomp,         &
              channels_rec = channels_rec)
    ENDIF
  ELSE
    IF (do_rttovscatt) THEN
      IF (do_radar) THEN
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance,       &
            reflectivity = reflectivity,      &
            nthreads = nthreads)
      ELSE
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance,       &
            nthreads = nthreads)
      ENDIF
    ELSE
      CALL rttov_parallel_direct(            &
              err,                           &
              chanprof,                      &
              opts,                          &
              profiles,                      &
              coefs,                         &
              transmission,                  &
              radiance,                      &
              calcemis = calcemis,           &
              emissivity = emissivity,       &
              calcrefl = calcrefl,           &
              reflectance = reflectance,     &
              aer_opt_param = aer_opt_param, &
              cld_opt_param = cld_opt_param, &
              pccomp       = pccomp,         &
              channels_rec = channels_rec,   &
              nthreads = nthreads)
    ENDIF
  ENDIF
  THROW(err.NE.0)

  emisout(:)        = emissivity(:)%emis_out
  reflout(:)        = reflectance(:)%refl_out
  diffusereflout(:) = reflectance(:)%diffuse_refl_out

  CALL make_chain_profile(err, chain_profiles, c, "PROFILES", profiles, .FALSE._jplm)
  THROW(err.NE.0)

  CALL make_chain_profile(err, chain_profiles_inc, c_inc, "PROFILES_INC", profiles_inc, .FALSE._jplm)
  THROW(err.NE.0)

  CALL make_chain_profile(err, chain_profiles_bf, c_bf, "PROFILES_BF", profiles_bf, .FALSE._jplm)
  THROW(err.NE.0)

  CALL make_chain_profile(err, chain_profiles_k, c_k, "PROFILES_K", profiles_k, .TRUE._jplm)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL make_chain_cld_profile(err, chain_cld_profiles, c_cld, "CLD_PROFILES", cld_profiles, .FALSE._jplm)
    THROW(err.NE.0)

    CALL make_chain_cld_profile(err, chain_cld_profiles_inc, c_cld_inc, "CLD_PROFILES_INC", cld_profiles_inc, .FALSE._jplm)
    THROW(err.NE.0)

    CALL make_chain_cld_profile(err, chain_cld_profiles_bf, c_cld_bf, "CLD_PROFILES_BF", cld_profiles_bf, .FALSE._jplm)
    THROW(err.NE.0)

    CALL make_chain_cld_profile(err, chain_cld_profiles_k, c_cld_k, "CLD_PROFILES_K", cld_profiles_k, .TRUE._jplm)
    THROW(err.NE.0)
  ENDIF

  NULLIFY (chain_profiles_k_rec, c_k_in)
  nchannels_rec = -1_jpim

  NULLIFY (chain_profiles_k_pc, c_k_pc)
  npcscores = -1_jpim

  NULLIFY (chanprof_in, chanprof_pc)
  NULLIFY (w, w_in, w_pc)

  IF (do_pc) THEN

    npcscores = SIZE(pccomp%total_pcscores)

    IF (do_radrec) THEN

      ! Reconstructed radiances Jacobians
      nchannels_rec = SIZE(pccomp%bt_pccomp)
      CALL make_chain_profile(err, chain_profiles_k_rec, c_k_in, "PROFILES_K_REC", profiles_k_rec, .true._jplm)
      THROW(err.NE.0)

      ALLOCATE (w_in(nchannels_rec), STAT = err)
      THROWM(err.NE.0,"Cannot allocate w_in")

    ELSE

      ! PC-scores Jacobians
      CALL make_chain_profile(err, chain_profiles_k_pc, c_k_pc, "PROFILES_K_PC", profiles_k_pc, .true._jplm)
      THROW(err.NE.0)

      ALLOCATE (w_pc(npcscores), STAT = err)
      THROWM(err.NE.0,"Cannot allocate w_pc")

    ENDIF

    CALL rttov_alloc_pc_dimensions(err, opts, npcscores, nprofiles, chanprof_in, chanprof_pc, &
                                   1_jpim, channels_rec = channels_rec)
    THROW(err.NE.0)

    CALL rttov_alloc_pccomp(err, pccomp_bf, npcscores, 1_jpim, nchannels_rec = nchannels_rec, &
                            opts = opts, nlevels = nlevels)
    THROW(err.NE.0)

  ENDIF

  ALLOCATE (w(nchannels), STAT = err)
  THROWM(err.NE.0,"Cannot allocate w")

  ! ---------------------------------------------------------------------------
  ! Perturb RTTOV profile variables
  ! ---------------------------------------------------------------------------
  CALL size_chain(nsize, chain_profiles(1))
  DO i = 1, nsize

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c(iprof)%p, x_profiles)
      CALL get_pointer_chain(c_bf(iprof)%p, x_profiles_bf)
      CALL get_pointer_chain(c_inc(iprof)%p, x_profiles_inc)
      x_profiles_bf = x_profiles + x_profiles_inc
      inc(iprof)    = x_profiles_inc
      IF (do_rttovscatt) cld_profiles_bf(iprof)%ph(nlevels+1) = profiles_bf(iprof)%s2m%p
    ENDDO

    IF (nthreads <= 0) THEN
      IF (do_rttovscatt) THEN
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles_bf,    &
            cld_profiles_bf,  coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance_bf,    &
            reflectivity = reflectivity_bf)
      ELSE
        CALL rttov_direct(                     &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles_bf,                   &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity,       &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                traj         = traj,           &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec)
      ENDIF
    ELSE
      IF (do_rttovscatt) THEN
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles_bf,    &
            cld_profiles_bf,  coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance_bf,    &
            reflectivity = reflectivity_bf,   &
            nthreads = nthreads)
      ELSE
        CALL rttov_parallel_direct(            &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles_bf,                   &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity,       &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec,   &
                nthreads     = nthreads)
      ENDIF
    ENDIF
    THROW(err.NE.0)

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c(iprof)%p, x_profiles)
      CALL get_pointer_chain(c_bf(iprof)%p, x_profiles_bf)
      x_profiles_bf = x_profiles
    ENDDO

    DO iprof = 1, nprofiles
      CALL advance_chain(c(iprof)%p)
      CALL advance_chain(c_bf(iprof)%p)
      CALL advance_chain(c_inc(iprof)%p)
    ENDDO

    IF (do_pc) THEN
      IF (do_radrec) THEN
        WHERE (ABS(inc(chanprof_in(:)%prof)) > 0._jprb)
          w_in = 1._jprb / inc(chanprof_in(:)%prof)
        ELSEWHERE
          w_in = 0._jprb
        ENDWHERE
      ELSE
        WHERE (ABS(inc(chanprof_pc(:)%prof)) > 0._jprb)
          w_pc = 1._jprb / inc(chanprof_pc(:)%prof)
        ELSEWHERE
          w_pc = 0._jprb
        ENDWHERE
      ENDIF
    END IF
    IF (.NOT. opts%htfrtc_opts%htfrtc) THEN
      WHERE (ABS(inc(chanprof(:)%prof)) > 0._jprb)
        w = 1._jprb / inc(chanprof(:)%prof)
      ELSEWHERE
        w = 0._jprb
      ENDWHERE
    ENDIF

    IF (do_rttovscatt .OR. opts%rt_all%switchrad) THEN
      IF (do_radar) THEN
        CALL assign_chain_profile (c_k, (reflectivity_bf%azef(radar_k_lev,:)-reflectivity%azef(radar_k_lev,:))*w)
      ELSEIF (do_pc) THEN
        IF (.NOT. opts%htfrtc_opts%htfrtc) CALL assign_chain_profile (c_k, (radiance_bf%total-radiance%total)*w)
        IF (do_radrec) THEN
          CALL assign_chain_profile (c_k_in, (pccomp_bf%bt_pccomp-pccomp%bt_pccomp)*w_in)
        ELSE
          CALL assign_chain_profile (c_k_pc, (pccomp_bf%total_pcscores-pccomp%total_pcscores)*w_pc)
        ENDIF
      ELSEIF (.NOT. opts%htfrtc_opts%htfrtc) THEN
        DO ichan = 1, nchannels
          IF (do_rttovscatt .OR. coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
            out_bf(ichan) = (radiance_bf%bt(ichan)-radiance%bt(ichan))*w(ichan)
          ELSE
            out_bf(ichan) = (radiance_bf%total(ichan)-radiance%total(ichan))*w(ichan)
          ENDIF
        ENDDO
        CALL assign_chain_profile (c_k, out_bf)
      ENDIF
    ELSE
      IF (.NOT. opts%htfrtc_opts%htfrtc) CALL assign_chain_profile (c_k, (radiance_bf%total-radiance%total)*w)
      IF (do_pc) THEN
        IF (do_radrec) THEN
          CALL assign_chain_profile (c_k_in, (pccomp_bf%total_pccomp-pccomp%total_pccomp)*w_in)
        ELSE
          CALL assign_chain_profile (c_k_pc, (pccomp_bf%total_pcscores-pccomp%total_pcscores)*w_pc)
        ENDIF
      ENDIF
    ENDIF

  ENDDO ! chain loop

  ! ---------------------------------------------------------------------------
  ! Perturb RTTOV-SCATT cloud profile variables
  ! ---------------------------------------------------------------------------
  IF (do_rttovscatt) THEN
    CALL rttov_copy_prof(profiles_bf, profiles)
    CALL rttov_copy_scatt_prof(cld_profiles_bf, cld_profiles)
    CALL size_chain(nsize, chain_cld_profiles(1))
    DO i = 1, nsize

      DO iprof = 1, nprofiles
        CALL get_pointer_chain(c_cld(iprof)%p, x_profiles)
        CALL get_pointer_chain(c_cld_bf(iprof)%p, x_profiles_bf)
        CALL get_pointer_chain(c_cld_inc(iprof)%p, x_profiles_inc)
        x_profiles_bf = x_profiles + x_profiles_inc
        inc(iprof)    = x_profiles_inc
        profiles_bf(iprof)%s2m%p = cld_profiles_bf(iprof)%ph(nlevels+1)
      ENDDO

      IF (nthreads <= 0) THEN
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles_bf,    &
            cld_profiles_bf,  coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance_bf,    &
            reflectivity = reflectivity_bf)
      ELSE
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles_bf,    &
            cld_profiles_bf,  coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity,       radiance_bf,    &
            reflectivity = reflectivity_bf,   &
            nthreads = nthreads)
      ENDIF
      THROW(err.NE.0)

      DO iprof = 1, nprofiles
        CALL get_pointer_chain(c_cld(iprof)%p, x_profiles)
        CALL get_pointer_chain(c_cld_bf(iprof)%p, x_profiles_bf)
        x_profiles_bf = x_profiles
      ENDDO

      DO iprof = 1, nprofiles
        CALL advance_chain(c_cld(iprof)%p)
        CALL advance_chain(c_cld_bf(iprof)%p)
        CALL advance_chain(c_cld_inc(iprof)%p)
      ENDDO


      WHERE (inc(chanprof(:)%prof) .NE. 0._jprb)
        w = 1._jprb / inc(chanprof(:)%prof)
      ELSEWHERE
        w = 0._jprb
      ENDWHERE

      IF (do_radar) THEN
        out_bf = (reflectivity_bf%azef(radar_k_lev,:) - reflectivity%azef(radar_k_lev,:)) * w
      ELSE
        out_bf = (radiance_bf%bt - radiance%bt) * w
      ENDIF
      CALL assign_chain_profile (c_cld_k, out_bf)

    ENDDO ! cld_profiles chain loop
  ENDIF

  CALL free_chain_profile (err, chain_profiles, c)
  THROW(err.NE.0)

  CALL free_chain_profile (err, chain_profiles_bf, c_bf)
  THROW(err.NE.0)

  CALL free_chain_profile (err, chain_profiles_inc, c_inc)
  THROW(err.NE.0)

  CALL free_chain_profile (err, chain_profiles_k, c_k)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL free_chain_profile (err, chain_cld_profiles, c_cld)
    THROW(err.NE.0)

    CALL free_chain_profile (err, chain_cld_profiles_bf, c_cld_bf)
    THROW(err.NE.0)

    CALL free_chain_profile (err, chain_cld_profiles_inc, c_cld_inc)
    THROW(err.NE.0)

    CALL free_chain_profile (err, chain_cld_profiles_k, c_cld_k)
    THROW(err.NE.0)
  ENDIF

  IF (do_pc) THEN
    IF (do_radrec) THEN
      CALL free_chain_profile (err, chain_profiles_k_rec, c_k_in)
      THROW(err.NE.0)
    ELSE
      CALL free_chain_profile (err, chain_profiles_k_pc, c_k_pc)
      THROW(err.NE.0)
    ENDIF
  ENDIF


  ! ---------------------------------------------------------------------------
  ! Perturb explicit aerosol/cloud optical properties
  ! ---------------------------------------------------------------------------
  ! The chain approach doesn't work for the optical properties so we perturb each variable in turn.
  ! All channels are perturbed together for each variable in each layer.
  ! For efficiency we avoid perturbing non-scattering layers since Jacobians are zero anyway in this case.

  IF ((opts%rt_ir%user_aer_opt_param .OR. opts%rt_ir%user_cld_opt_param) .AND. &
      .NOT. opts%dev%no_opt_param_tladk) THEN
    DO aercld = 1, 2

      IF (aercld == 1) THEN
        IF (.NOT. opts%rt_ir%user_aer_opt_param) CYCLE
        opt_param => aer_opt_param
        opt_param_bf => aer_opt_param_bf
        opt_param_inc => aer_opt_param_inc
        opt_param_k => aer_opt_param_k

        IF (opts%rt_ir%user_cld_opt_param) CALL rttov_copy_opt_param(cld_opt_param_bf, cld_opt_param)
      ELSE
        IF (.NOT. opts%rt_ir%user_cld_opt_param) CYCLE
        opt_param => cld_opt_param
        opt_param_bf => cld_opt_param_bf
        opt_param_inc => cld_opt_param_inc
        opt_param_k => cld_opt_param_k

        IF (opts%rt_ir%user_aer_opt_param) CALL rttov_copy_opt_param(aer_opt_param_bf, aer_opt_param)
      ENDIF

      CALL rttov_init_opt_param(err, opts, opt_param_k, .TRUE._jplm)

      DO var = 1, 5
        DO lay = 1, nlayers
          IF (SUM(opt_param%abs(lay,:)) + SUM(opt_param%sca(lay,:)) == 0._jprb) CYCLE

          IF (var <= 3) THEN
            nsize = 1
          ELSEIF (var == 4) THEN
            nsize = opt_param%nmom + 1
          ELSEIF (var == 5) THEN
            nsize = SIZE(opt_param%phangle)
          ENDIF

          DO i = 1, nsize
            IF (var == 1) THEN
              var_in => opt_param%abs(lay,:)
              var_bf => opt_param_bf%abs(lay,:)
              var_inc => opt_param_inc%abs(lay,:)
              var_out => opt_param_k%abs(lay,:)
            ELSEIF (var == 2) THEN
              var_in => opt_param%sca(lay,:)
              var_bf => opt_param_bf%sca(lay,:)
              var_inc => opt_param_inc%sca(lay,:)
              var_out => opt_param_k%sca(lay,:)
            ELSEIF (var == 3) THEN
              var_in => opt_param%bpr(lay,:)
              var_bf => opt_param_bf%bpr(lay,:)
              var_inc => opt_param_inc%bpr(lay,:)
              var_out => opt_param_k%bpr(lay,:)
            ELSEIF (var == 4) THEN
              var_in => opt_param%legcoef(i,lay,:)
              var_bf => opt_param_bf%legcoef(i,lay,:)
              var_inc => opt_param_inc%legcoef(i,lay,:)
              var_out => opt_param_k%legcoef(i,lay,:)
            ELSEIF (var == 5) THEN
              var_in => opt_param%pha(i,lay,:)
              var_bf => opt_param_bf%pha(i,lay,:)
              var_inc => opt_param_inc%pha(i,lay,:)
              var_out => opt_param_k%pha(i,lay,:)
            ENDIF
            var_bf(:) = var_in(:) + var_inc(:)

            IF (nthreads <= 0) THEN
              CALL rttov_direct( &
                      err,                               &
                      chanprof,                          &
                      opts,                              &
                      profiles,                          &
                      coefs,                             &
                      transmission_bf,                   &
                      radiance_bf,                       &
                      calcemis      = calcemis,          &
                      emissivity    = emissivity,        &
                      calcrefl      = calcrefl,          &
                      reflectance   = reflectance,       &
                      aer_opt_param = aer_opt_param_bf,  &
                      cld_opt_param = cld_opt_param_bf,  &
                      traj          = traj,              &
                      channels_rec  = channels_rec)
            ELSE
              CALL rttov_parallel_direct( &
                      err,                               &
                      chanprof,                          &
                      opts,                              &
                      profiles,                          &
                      coefs,                             &
                      transmission_bf,                   &
                      radiance_bf,                       &
                      calcemis      = calcemis,          &
                      emissivity    = emissivity,        &
                      calcrefl      = calcrefl,          &
                      reflectance   = reflectance,       &
                      aer_opt_param = aer_opt_param_bf,  &
                      cld_opt_param = cld_opt_param_bf,  &
                      channels_rec  = channels_rec,      &
                      nthreads      = nthreads)
            ENDIF
            THROW(err.NE.0)

            var_bf(:) = var_in(:)

            DO ichan = 1, nchannels
              IF (ABS(var_inc(ichan)) > 0._jprb) THEN
                IF (opts%rt_all%switchrad .AND. coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
                  var_out(ichan) = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / var_inc(ichan)
                ELSE
                  var_out(ichan) = (radiance_bf%total(ichan) - radiance%total(ichan)) / var_inc(ichan)
                ENDIF
              ELSE
                var_out(ichan) = 0._jprb
              ENDIF
            ENDDO
          ENDDO ! nsize
        ENDDO ! nlayers
      ENDDO ! nvar
    ENDDO ! aer/cld
  ENDIF


  ! ---------------------------------------------------------------------------
  ! Perturb surface emissivities
  ! ---------------------------------------------------------------------------
  ! All channels can be perturbed together with just one call to RTTOV
  ! which is much more efficient than doing each one individually.

  ! For PC models emissivities are computed on predictor/centroid frequencies.
  ! This means the Jacobians are expressed in terms of predictor radiances, not
  ! PC scores or reconstructed radiances. Since HTFRTC doesn't explicitly output
  ! the centroid radiances we can't compute these Jacobians.

  IF (.NOT. opts%htfrtc_opts%htfrtc) THEN

    emissivity_bf%emis_in = emisout + demis
    emissivity_bf%specularity = emissivity%specularity
    DO ichan = 1, nchannels
      ! Set calcemis_bf false for IR sensors with calcemis true over land, sea-ice or
      ! with ISEM since the Jacobian is non-zero so we need to manually perturb the
      ! original emissivity. For MW there is no issue.
      IF (coefs%coef%id_sensor == sensor_id_ir .OR. coefs%coef%id_sensor == sensor_id_hi) THEN
        calcemis_bf(ichan) = calcemis(ichan) .AND. &
                             profiles(chanprof(ichan)%prof)%skin%surftype == 1 .AND. &
                             opts%rt_ir%ir_sea_emis_model /= 1
      ELSE
        calcemis_bf(ichan) = calcemis(ichan)
      ENDIF
    ENDDO

    IF (nthreads <= 0) THEN
      IF (do_rttovscatt) THEN
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis_bf,    &
            emissivity_bf,    radiance_bf,    &
            reflectivity = reflectivity_bf)
      ELSE
        CALL rttov_direct(                     &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis_bf,        &
                emissivity = emissivity_bf,    &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                traj         = traj,           &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec)
      ENDIF
    ELSE
      IF (do_rttovscatt) THEN
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis_bf,    &
            emissivity_bf,    radiance_bf,    &
            reflectivity = reflectivity_bf,   &
            nthreads = nthreads)
      ELSE
        CALL rttov_parallel_direct(            &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis_bf,        &
                emissivity = emissivity_bf,    &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec,   &
                nthreads     = nthreads)
      ENDIF
    ENDIF
    THROW(err.NE.0)

    DO ichan = 1, nchannels
      IF (do_rttovscatt .OR. &
          (opts%rt_all%switchrad .AND. .NOT. do_pc .AND. &
           coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2)) THEN
        IF (do_radar) THEN
          emissivity_k(ichan)%emis_in = (reflectivity_bf%azef(radar_k_lev,ichan) - &
                                         reflectivity%azef(radar_k_lev,ichan)) / demis
        ELSE
          emissivity_k(ichan)%emis_in = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / demis
        ENDIF
      ELSE
        emissivity_k(ichan)%emis_in = (radiance_bf%total(ichan) - radiance%total(ichan)) / demis
      ENDIF
    ENDDO

  ENDIF ! not HTFRTC

  ! ---------------------------------------------------------------------------
  ! Perturb surface specularity parameters
  ! ---------------------------------------------------------------------------
  ! All channels can be perturbed together with just one call to RTTOV
  ! which is much more efficient than doing each one individually.

  ! See emissivity section above for notes on PC models.

  IF (opts%rt_all%do_lambertian .AND. .NOT. opts%htfrtc_opts%htfrtc) THEN

    emissivity_bf%emis_in = emisout
    ! Ensure the perturbed specularity remains in the interval [0,1]
    WHERE (emissivity%specularity < 0.5_jprb)
      emissivity_bf%specularity = emissivity%specularity + dspec
    ELSEWHERE
      emissivity_bf%specularity = emissivity%specularity - dspec
    ENDWHERE

    IF (nthreads <= 0) THEN
      IF (do_rttovscatt) THEN
        CALL rttov_scatt(                     &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity_bf,    radiance_bf,    &
            reflectivity = reflectivity_bf)
      ELSE
        CALL rttov_direct(                     &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity_bf,    &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                traj         = traj,           &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec)
      ENDIF
    ELSE
      IF (do_rttovscatt) THEN
        CALL rttov_parallel_scatt(            &
            err,              opts_scatt,     &
            nlevels,          chanprof,       &
            frequencies,      profiles,       &
            cld_profiles,     coefs,          &
            coefs_scatt,      calcemis,       &
            emissivity_bf,    radiance_bf,    &
            reflectivity = reflectivity_bf,   &
            nthreads = nthreads)
      ELSE
        CALL rttov_parallel_direct(            &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity_bf,    &
                calcrefl = calcrefl,           &
                reflectance = reflectance,     &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec,   &
                nthreads     = nthreads)
      ENDIF
    ENDIF
    THROW(err.NE.0)

    DO ichan = 1, nchannels
      IF (do_rttovscatt .OR. &
          (opts%rt_all%switchrad .AND. .NOT. do_pc .AND. &
           coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2)) THEN
        IF (do_radar) THEN
          IF (emissivity(ichan)%specularity < 0.5_jprb) THEN
            emissivity_k(ichan)%specularity = (reflectivity_bf%azef(radar_k_lev,ichan) - &
                                               reflectivity%azef(radar_k_lev,ichan)) / dspec
          ELSE
            emissivity_k(ichan)%specularity = (reflectivity_bf%azef(radar_k_lev,ichan) - &
                                               reflectivity%azef(radar_k_lev,ichan)) / (-dspec)
          ENDIF
        ELSE
          IF (emissivity(ichan)%specularity < 0.5_jprb) THEN
            emissivity_k(ichan)%specularity = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / dspec
          ELSE
            emissivity_k(ichan)%specularity = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / (-dspec)
          ENDIF
        ENDIF
      ELSE
        IF (emissivity(ichan)%specularity < 0.5_jprb) THEN
          emissivity_k(ichan)%specularity = (radiance_bf%total(ichan) - radiance%total(ichan)) / dspec
        ELSE
          emissivity_k(ichan)%specularity = (radiance_bf%total(ichan) - radiance%total(ichan)) / (-dspec)
        ENDIF
      ENDIF
    ENDDO

  ENDIF ! do_lambertian and not HTFRTC

  ! ---------------------------------------------------------------------------
  ! Perturb surface reflectances
  ! ---------------------------------------------------------------------------
  ! All channels can be perturbed together with just one call to RTTOV
  ! which is much more efficient than doing each one individually.
  IF (opts%rt_ir%addsolar) THEN

    ! Surface BRDF

    ! Reflectances may depend on profile/emis values when calcrefl is true. Therefore we
    ! use the saved reflectances from the original call to rttov_direct.

    reflectance_bf%refl_cloud_top = reflectance%refl_cloud_top
    reflectance_bf%refl_in = reflout + drefl
    reflectance_bf%diffuse_refl_in = diffusereflout
    DO ichan = 1, nchannels
      ! Set calcrefl_bf false with calcrefl true over land or sea-ice for pure
      ! solar channels since the Jacobian is non-zero so we need to manually
      ! perturb the original BRDF.
      calcrefl_bf(ichan) = calcrefl(ichan) .AND. &
                           (profiles(chanprof(ichan)%prof)%skin%surftype == 1 .OR. &
                            coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2)
    ENDDO

    IF (nthreads <= 0) THEN
      CALL rttov_direct(                     &
              err,                           &
              chanprof,                      &
              opts,                          &
              profiles,                      &
              coefs,                         &
              transmission_bf,               &
              radiance_bf,                   &
              calcemis = calcemis,           &
              emissivity = emissivity,       &
              calcrefl = calcrefl_bf,        &
              reflectance = reflectance_bf,  &
              aer_opt_param = aer_opt_param, &
              cld_opt_param = cld_opt_param, &
              traj         = traj,           &
              pccomp       = pccomp_bf,      &
              channels_rec = channels_rec)
    ELSE
      CALL rttov_parallel_direct(            &
              err,                           &
              chanprof,                      &
              opts,                          &
              profiles,                      &
              coefs,                         &
              transmission_bf,               &
              radiance_bf,                   &
              calcemis = calcemis,           &
              emissivity = emissivity,       &
              calcrefl = calcrefl_bf,        &
              reflectance = reflectance_bf,  &
              aer_opt_param = aer_opt_param, &
              cld_opt_param = cld_opt_param, &
              pccomp       = pccomp_bf,      &
              channels_rec = channels_rec,   &
              nthreads     = nthreads)
    ENDIF
    THROW(err.NE.0)

    DO ichan = 1, nchannels
      IF (opts%rt_all%switchrad .AND. &
        coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
        reflectance_k(ichan)%refl_in = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / drefl
      ELSE
        reflectance_k(ichan)%refl_in = (radiance_bf%total(ichan) - radiance%total(ichan)) / drefl
      ENDIF
    ENDDO


    ! Surface diffuse reflectance
    ! Only relevant if some diffuse_refl_in > 0

    IF (ANY(reflectance%diffuse_refl_in > 0._jprb)) THEN
      ! Reflectances may depend on profile/emis values when calcrefl is true. Therefore we
      ! use the saved reflectances from the original call to rttov_direct.

      ! We can perturb the diffuse reflectance for all channels: it is only used
      ! for pure solar channels where calcrefl is false, and is ignored in all
      ! other cases.

      reflectance_bf%refl_cloud_top = reflectance%refl_cloud_top
      reflectance_bf%refl_in = reflout
      reflectance_bf%diffuse_refl_in = diffusereflout + drefl
      calcrefl_bf = calcrefl

      IF (nthreads <= 0) THEN
        CALL rttov_direct(                     &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity,       &
                calcrefl = calcrefl_bf,        &
                reflectance = reflectance_bf,  &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                traj         = traj,           &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec)
      ELSE
        CALL rttov_parallel_direct(            &
                err,                           &
                chanprof,                      &
                opts,                          &
                profiles,                      &
                coefs,                         &
                transmission_bf,               &
                radiance_bf,                   &
                calcemis = calcemis,           &
                emissivity = emissivity,       &
                calcrefl = calcrefl_bf,        &
                reflectance = reflectance_bf,  &
                aer_opt_param = aer_opt_param, &
                cld_opt_param = cld_opt_param, &
                pccomp       = pccomp_bf,      &
                channels_rec = channels_rec,   &
                nthreads     = nthreads)
      ENDIF
      THROW(err.NE.0)

      DO ichan = 1, nchannels
        IF (opts%rt_all%switchrad .AND. &
          coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
          reflectance_k(ichan)%diffuse_refl_in = (radiance_bf%bt(ichan) - radiance%bt(ichan)) / drefl
        ELSE
          reflectance_k(ichan)%diffuse_refl_in = (radiance_bf%total(ichan) - radiance%total(ichan)) / drefl
        ENDIF
      ENDDO
    ENDIF

  ENDIF

  ! ---------------------------------------------------------------------------
  ! Deallocations
  ! ---------------------------------------------------------------------------
  IF (do_pc) THEN

    CALL rttov_alloc_pccomp (err, pccomp_bf, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.NE.0)

    CALL rttov_alloc_pc_dimensions(err, opts, npcscores, nprofiles, chanprof_in, chanprof_pc, &
                                   0_jpim)
    THROW(err.NE.0)

    IF (do_radrec) THEN
      DEALLOCATE (w_in, STAT = err)
      THROWM(err.NE.0,"Cannot deallocate w_in")
    ELSE
      DEALLOCATE (w_pc, STAT = err)
      THROWM(err.NE.0,"Cannot deallocate w_pc")
    ENDIF

  ENDIF

  DEALLOCATE (w, STAT = err)
  THROWM(err.NE.0,"Cannot deallocate w")

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                          &
          err,                  nprofiles,          &
          nchannels,            opts,               &
          profiles(1)%nlevels,  coefs,              &
          0_jpim,               traj = traj)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_aer_opt_param) THEN
    CALL rttov_alloc_opt_param(err, aer_opt_param_bf, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)

    CALL rttov_alloc_opt_param(err, aer_opt_param_inc, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    CALL rttov_alloc_opt_param(err, cld_opt_param_bf, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)

    CALL rttov_alloc_opt_param(err, cld_opt_param_inc, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_bf, profiles(1)%nlevels, nchannels, 0_jpim)
  THROW(err.NE.0)

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_bf,         &
          profiles(1)%nlevels, &
          0_jpim)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    IF (do_radar) THEN
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_bf,         &
            nlevels,                 &
            0_jpim)
      THROW(err.NE.0)
      ALLOCATE(reflectivity_bf)
    ENDIF

    CALL rttov_alloc_scatt_prof( &
          err,                          &
          nprofiles,                    &
          cld_profiles_inc,             &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          0_jpim)
    THROW(err.NE.0)

    CALL rttov_alloc_scatt_prof( &
          err,                          &
          nprofiles,                    &
          cld_profiles_bf,              &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_inc,        &
          profiles(1)%nlevels, &
          opts,                &
          0_jpim)
  THROW(err.NE.0)

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_bf,         &
          profiles(1)%nlevels, &
          opts,                &
          0_jpim)
  THROW(err.NE.0)

  CATCH

END SUBROUTINE rttov_k_bf
