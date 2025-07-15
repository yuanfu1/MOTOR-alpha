! Description:
!> @file
!!   Compute Jacobian using AD model
!
!> @brief
!!   Compute Jacobian using AD model
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
SUBROUTINE rttov_k_ad( &
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
  USE rttov_test_k_mod, ONLY : radar_k_lev
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
  TYPE(rttov_opt_param),     INTENT(IN)   , OPTIONAL         :: aer_opt_param
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL         :: aer_opt_param_k
  TYPE(rttov_opt_param),     INTENT(IN)   , OPTIONAL         :: cld_opt_param
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL         :: cld_opt_param_k
  TYPE(rttov_pccomp),        INTENT(INOUT), OPTIONAL         :: pccomp
  TYPE(rttov_profile),       INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_pc(:)
  TYPE(rttov_profile),       INTENT(INOUT), OPTIONAL, TARGET :: profiles_k_rec(:)
  INTEGER(KIND=jpim),        INTENT(IN)   , OPTIONAL         :: channels_rec(:)
  INTEGER(KIND=jpim),        INTENT(IN)                      :: nthreads
!INTF_END

#include "rttov_ad.interface"
#include "rttov_parallel_ad.interface"
#include "rttov_scatt_ad.interface"
#include "rttov_parallel_scatt_ad.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_scatt_prof.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_init_pccomp.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_init_reflectivity.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_scatt_prof.interface"
#include "rttov_copy_opt_param.interface"
#include "rttov_errorreport.interface"

  TYPE(rttov_traj)                  :: traj, traj_ad
  TYPE(rttov_pccomp)                :: pccomp_ad
  TYPE(rttov_profile)               :: profiles_ad(SIZE(profiles))
  TYPE(rttov_profile_cloud)         :: cld_profiles_ad(SIZE(profiles))
  TYPE(rttov_opt_param)             :: aer_opt_param_ad
  TYPE(rttov_opt_param)             :: cld_opt_param_ad
  TYPE(rttov_emissivity)            :: emissivity_ad(SIZE(chanprof))
  TYPE(rttov_reflectance)           :: reflectance_ad(SIZE(chanprof))
  TYPE(rttov_transmission)          :: transmission_ad
  TYPE(rttov_radiance)              :: radiance_ad
  TYPE(rttov_reflectivity), POINTER :: reflectivity_ad => NULL()
  INTEGER(KIND=jpim)                :: ichan, iprof, ipcscore, ichannels_rec
  INTEGER(KIND=jpim)                :: nprofiles, nchannels, nlevels, npcscores, nchannels_rec, nphangle
  LOGICAL(KIND=jplm)                :: do_rttovscatt, do_radar

  TRY

  nprofiles = SIZE(profiles)
  nchannels = SIZE(chanprof)
  nlevels = profiles(1)%nlevels

  do_rttovscatt = ASSOCIATED(cld_profiles)
  do_radar = do_rttovscatt .AND. ASSOCIATED(reflectivity%zef)

  npcscores = -1_jpim
  nchannels_rec = -1_jpim

  IF (opts%rt_ir%pc%addpc) THEN
    npcscores = SIZE(pccomp%total_pcscores)
    IF (opts%rt_ir%pc%addradrec) THEN
      nchannels_rec = SIZE(pccomp%bt_pccomp)
    ENDIF
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp(err, pccomp_ad, npcscores, 1_jpim, nchannels_rec = nchannels_rec)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_ad,         &
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
          cld_profiles_ad,              &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          1_jpim,                       &
          init=.TRUE._jplm)
    THROW(err.NE.0)
    IF (do_radar) THEN
      ALLOCATE(reflectivity_ad)
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_ad,         &
            nlevels,                 &
            1_jpim,                  &
            init = .TRUE._jplm)
      THROW(err.NE.0)
    ENDIF
  ENDIF

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_ad,         &
          profiles(1)%nlevels, &
          1_jpim)
  THROW(err.NE.0)

  CALL rttov_alloc_transmission(err, transmission_ad, profiles(1)%nlevels, nchannels, 1_jpim)
  THROW(err.NE.0)

  IF (opts%rt_ir%user_aer_opt_param) THEN
    nphangle = SIZE(aer_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, aer_opt_param_ad, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    nphangle = SIZE(cld_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, cld_opt_param_ad, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                   &
          err,                  nprofiles,   &
          nchannels,            opts,        &
          profiles(1)%nlevels,  coefs,       &
          1_jpim,               traj = traj, &
          traj_ad = traj_ad )
    THROW(err.NE.0)
  ENDIF

!
! Regular channels
!
  DO ichan = 1, nchannels

    CALL rttov_init_emis_refl(emissivity_ad, reflectance_ad)
    CALL rttov_init_transmission(transmission_ad)
    CALL rttov_init_rad(radiance_ad)
    IF (do_radar) CALL rttov_init_reflectivity(reflectivity_ad)
    IF (opts%rt_ir%pc%addpc) CALL rttov_init_pccomp(pccomp_ad)

    IF (do_radar) THEN
      reflectivity_ad%azef(radar_k_lev,ichan) = 1._jprb
    ELSE IF (opts%rt_ir%pc%addpc) THEN
      radiance_ad%total(ichan) = 1._jprb
    ELSE IF (opts%rt_all%switchrad .AND. &
             coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
      radiance_ad%bt(ichan) = 1._jprb
    ELSE
      radiance_ad%total(ichan) = 1._jprb
    ENDIF

    CALL rttov_init_prof(profiles_ad)
    IF (do_rttovscatt) CALL rttov_init_scatt_prof(cld_profiles_ad)
    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_init_opt_param(err, opts, aer_opt_param_ad, .TRUE._jplm)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_init_opt_param(err, opts, cld_opt_param_ad, .TRUE._jplm)

    IF (nthreads <= 0) THEN
      IF (do_rttovscatt) THEN
        IF (do_radar) THEN
          CALL rttov_scatt_ad(                  &
              err,              opts_scatt,     &
              nlevels,          chanprof,       &
              frequencies,      profiles,       &
              cld_profiles,     coefs,          &
              coefs_scatt,      calcemis,       &
              emissivity,       profiles_ad,    &
              cld_profiles_ad,  emissivity_ad,  &
              radiance,         radiance_ad,    &
              reflectivity,     reflectivity_ad)
        ELSE
          CALL rttov_scatt_ad(                  &
              err,              opts_scatt,     &
              nlevels,          chanprof,       &
              frequencies,      profiles,       &
              cld_profiles,     coefs,          &
              coefs_scatt,      calcemis,       &
              emissivity,       profiles_ad,    &
              cld_profiles_ad,  emissivity_ad,  &
              radiance,         radiance_ad)
        ENDIF
      ELSE
        CALL rttov_ad(                                &
                err,                                  &
                chanprof,                             &
                opts,                                 &
                profiles,                             &
                profiles_ad,                          &
                coefs,                                &
                transmission,                         &
                transmission_ad,                      &
                radiance,                             &
                radiance_ad,                          &
                calcemis         = calcemis,          &
                emissivity       = emissivity,        &
                emissivity_ad    = emissivity_ad,     &
                calcrefl         = calcrefl,          &
                reflectance      = reflectance,       &
                reflectance_ad   = reflectance_ad,    &
                aer_opt_param    = aer_opt_param,     &
                aer_opt_param_ad = aer_opt_param_ad,  &
                cld_opt_param    = cld_opt_param,     &
                cld_opt_param_ad = cld_opt_param_ad,  &
                traj             = traj,              &
                traj_ad          = traj_ad,           &
                pccomp           = pccomp,            &
                pccomp_ad        = pccomp_ad,         &
                channels_rec     = channels_rec)
      ENDIF
    ELSE
      IF (do_rttovscatt) THEN
        IF (do_radar) THEN
          CALL rttov_parallel_scatt_ad(          &
              err,              opts_scatt,      &
              nlevels,          chanprof,        &
              frequencies,      profiles,        &
              cld_profiles,     coefs,           &
              coefs_scatt,      calcemis,        &
              emissivity,       profiles_ad,     &
              cld_profiles_ad,  emissivity_ad,   &
              radiance,         radiance_ad,     &
              reflectivity,     reflectivity_ad, &
              nthreads = nthreads)
        ELSE
          CALL rttov_parallel_scatt_ad(          &
              err,              opts_scatt,      &
              nlevels,          chanprof,        &
              frequencies,      profiles,        &
              cld_profiles,     coefs,           &
              coefs_scatt,      calcemis,        &
              emissivity,       profiles_ad,     &
              cld_profiles_ad,  emissivity_ad,   &
              radiance,         radiance_ad,     &
              nthreads = nthreads)
        ENDIF
      ELSE
        CALL rttov_parallel_ad(                       &
                err,                                  &
                chanprof,                             &
                opts,                                 &
                profiles,                             &
                profiles_ad,                          &
                coefs,                                &
                transmission,                         &
                transmission_ad,                      &
                radiance,                             &
                radiance_ad,                          &
                calcemis         = calcemis,          &
                emissivity       = emissivity,        &
                emissivity_ad    = emissivity_ad,     &
                calcrefl         = calcrefl,          &
                reflectance      = reflectance,       &
                reflectance_ad   = reflectance_ad,    &
                aer_opt_param    = aer_opt_param,     &
                aer_opt_param_ad = aer_opt_param_ad,  &
                cld_opt_param    = cld_opt_param,     &
                cld_opt_param_ad = cld_opt_param_ad,  &
                pccomp           = pccomp,            &
                pccomp_ad        = pccomp_ad,         &
                channels_rec     = channels_rec,      &
                nthreads         = nthreads)
      ENDIF
    ENDIF
    THROW(err.NE.0)

    iprof = chanprof(ichan)%prof
    CALL rttov_copy_prof(profiles_k(ichan:ichan), profiles_ad(iprof:iprof))
    IF (do_rttovscatt) CALL rttov_copy_scatt_prof(cld_profiles_k(ichan:ichan), cld_profiles_ad(iprof:iprof))
    emissivity_k(ichan)%emis_in = emissivity_ad(ichan)%emis_in
    emissivity_k(ichan)%specularity = emissivity_ad(ichan)%specularity
    reflectance_k(ichan)%refl_in = reflectance_ad(ichan)%refl_in
    reflectance_k(ichan)%diffuse_refl_in = reflectance_ad(ichan)%diffuse_refl_in

  ENDDO

!
! Explicit optical properties: AD is the same as the K because the state vector contains
! per-channel profiles of optical properties, so only need to make one call to rttov_ad
!
  IF ((opts%rt_ir%user_aer_opt_param .OR. opts%rt_ir%user_cld_opt_param) .AND. &
      .NOT. opts%dev%no_opt_param_tladk) THEN

    CALL rttov_init_emis_refl(emissivity_ad, reflectance_ad)
    CALL rttov_init_transmission(transmission_ad)
    CALL rttov_init_rad(radiance_ad)

    radiance_ad%bt(:) = 1._jprb
    radiance_ad%total(:) = 1._jprb

    CALL rttov_init_prof(profiles_ad)
    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_init_opt_param(err, opts, aer_opt_param_ad, .TRUE._jplm)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_init_opt_param(err, opts, cld_opt_param_ad, .TRUE._jplm)

    IF (nthreads <= 0) THEN
      CALL rttov_ad(                                &
              err,                                  &
              chanprof,                             &
              opts,                                 &
              profiles,                             &
              profiles_ad,                          &
              coefs,                                &
              transmission,                         &
              transmission_ad,                      &
              radiance,                             &
              radiance_ad,                          &
              calcemis         = calcemis,          &
              emissivity       = emissivity,        &
              emissivity_ad    = emissivity_ad,     &
              calcrefl         = calcrefl,          &
              reflectance      = reflectance,       &
              reflectance_ad   = reflectance_ad,    &
              aer_opt_param    = aer_opt_param,     &
              aer_opt_param_ad = aer_opt_param_ad,  &
              cld_opt_param    = cld_opt_param,     &
              cld_opt_param_ad = cld_opt_param_ad,  &
              traj             = traj,              &
              traj_ad          = traj_ad,           &
              channels_rec     = channels_rec)
    ELSE
      CALL rttov_parallel_ad(                       &
              err,                                  &
              chanprof,                             &
              opts,                                 &
              profiles,                             &
              profiles_ad,                          &
              coefs,                                &
              transmission,                         &
              transmission_ad,                      &
              radiance,                             &
              radiance_ad,                          &
              calcemis         = calcemis,          &
              emissivity       = emissivity,        &
              emissivity_ad    = emissivity_ad,     &
              calcrefl         = calcrefl,          &
              reflectance      = reflectance,       &
              reflectance_ad   = reflectance_ad,    &
              aer_opt_param    = aer_opt_param,     &
              aer_opt_param_ad = aer_opt_param_ad,  &
              cld_opt_param    = cld_opt_param,     &
              cld_opt_param_ad = cld_opt_param_ad,  &
              channels_rec     = channels_rec,      &
              nthreads         = nthreads)
    ENDIF
    THROW(err.NE.0)

    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_copy_opt_param(aer_opt_param_k, aer_opt_param_ad)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_copy_opt_param(cld_opt_param_k, cld_opt_param_ad)

  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN

    IF (opts%rt_ir%pc%addradrec) THEN

!
! Reconstructed radiance
!

      DO ichannels_rec = 1, nchannels_rec

        CALL rttov_init_emis_refl(emissivity_ad, reflectance_ad)
        CALL rttov_init_transmission(transmission_ad)
        CALL rttov_init_rad(radiance_ad)
        CALL rttov_init_pccomp(pccomp_ad)

        IF (opts%rt_all%switchrad) THEN
          pccomp_ad%bt_pccomp(ichannels_rec) = 1._jprb
        ELSE
          pccomp_ad%total_pccomp(ichannels_rec) = 1._jprb
        ENDIF

        CALL rttov_init_prof(profiles_ad)

        IF (nthreads <= 0) THEN
          CALL rttov_ad(                             &
                  err,                               &
                  chanprof,                          &
                  opts,                              &
                  profiles,                          &
                  profiles_ad,                       &
                  coefs,                             &
                  transmission,                      &
                  transmission_ad,                   &
                  radiance,                          &
                  radiance_ad,                       &
                  calcemis         = calcemis,       &
                  emissivity       = emissivity,     &
                  emissivity_ad    = emissivity_ad,  &
                  calcrefl         = calcrefl,       &
                  reflectance      = reflectance,    &
                  reflectance_ad   = reflectance_ad, &
                  traj             = traj,           &
                  traj_ad          = traj_ad,        &
                  pccomp           = pccomp,         &
                  pccomp_ad        = pccomp_ad,      &
                  channels_rec     = channels_rec)
        ELSE
          CALL rttov_parallel_ad(                    &
                  err,                               &
                  chanprof,                          &
                  opts,                              &
                  profiles,                          &
                  profiles_ad,                       &
                  coefs,                             &
                  transmission,                      &
                  transmission_ad,                   &
                  radiance,                          &
                  radiance_ad,                       &
                  calcemis         = calcemis,       &
                  emissivity       = emissivity,     &
                  emissivity_ad    = emissivity_ad,  &
                  calcrefl         = calcrefl,       &
                  reflectance      = reflectance,    &
                  reflectance_ad   = reflectance_ad, &
                  pccomp           = pccomp,         &
                  pccomp_ad        = pccomp_ad,      &
                  channels_rec     = channels_rec,   &
                  nthreads         = nthreads)
        ENDIF
        THROW(err.NE.0)

        iprof = 1 + ((nprofiles * (ichannels_rec-1))/nchannels_rec)
        CALL rttov_copy_prof(profiles_k_rec(ichannels_rec:ichannels_rec), profiles_ad(iprof:iprof))

      ENDDO

    ELSE
!
! PC-scores
!
      DO ipcscore = 1, npcscores

        CALL rttov_init_emis_refl(emissivity_ad, reflectance_ad)
        CALL rttov_init_transmission(transmission_ad)
        CALL rttov_init_rad(radiance_ad)
        CALL rttov_init_pccomp(pccomp_ad)

        pccomp_ad%total_pcscores(ipcscore) = 1._jprb

        CALL rttov_init_prof(profiles_ad)

        IF (nthreads <= 0) THEN
          CALL rttov_ad(                             &
                  err,                               &
                  chanprof,                          &
                  opts,                              &
                  profiles,                          &
                  profiles_ad,                       &
                  coefs,                             &
                  transmission,                      &
                  transmission_ad,                   &
                  radiance,                          &
                  radiance_ad,                       &
                  calcemis         = calcemis,       &
                  emissivity       = emissivity,     &
                  emissivity_ad    = emissivity_ad,  &
                  calcrefl         = calcrefl,       &
                  reflectance      = reflectance,    &
                  reflectance_ad   = reflectance_ad, &
                  traj             = traj,           &
                  traj_ad          = traj_ad,        &
                  pccomp           = pccomp,         &
                  pccomp_ad        = pccomp_ad,      &
                  channels_rec     = channels_rec)
        ELSE
          CALL rttov_parallel_ad(                    &
                  err,                               &
                  chanprof,                          &
                  opts,                              &
                  profiles,                          &
                  profiles_ad,                       &
                  coefs,                             &
                  transmission,                      &
                  transmission_ad,                   &
                  radiance,                          &
                  radiance_ad,                       &
                  calcemis         = calcemis,       &
                  emissivity       = emissivity,     &
                  emissivity_ad    = emissivity_ad,  &
                  calcrefl         = calcrefl,       &
                  reflectance      = reflectance,    &
                  reflectance_ad   = reflectance_ad, &
                  pccomp           = pccomp,         &
                  pccomp_ad        = pccomp_ad,      &
                  channels_rec     = channels_rec,   &
                  nthreads         = nthreads)
        ENDIF
        THROW(err.NE.0)

        iprof = 1 + ((nprofiles * (ipcscore-1))/npcscores)
        CALL rttov_copy_prof(profiles_k_pc(ipcscore:ipcscore), profiles_ad(iprof:iprof))

      ENDDO

    ENDIF

  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp(err, pccomp_ad, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.NE.0)
  ENDIF

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                   &
          err,                  nprofiles,   &
          nchannels,            opts,        &
          profiles(1)%nlevels,  coefs,       &
          0_jpim,               traj = traj, &
          traj_ad = traj_ad )
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    CALL rttov_alloc_opt_param(err, cld_opt_param_ad, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_aer_opt_param) THEN
    CALL rttov_alloc_opt_param(err, aer_opt_param_ad, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_ad, profiles(1)%nlevels, nchannels, 0_jpim)
  THROW(err.NE.0)

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_ad,         &
          profiles(1)%nlevels, &
          0_jpim)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    IF (do_radar) THEN
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_ad,         &
            nlevels,                 &
            0_jpim)
      THROW(err.NE.0)
      DEALLOCATE(reflectivity_ad)
    ENDIF
    CALL rttov_alloc_scatt_prof( &
          err,                          &
          nprofiles,                    &
          cld_profiles_ad,              &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_ad,         &
          profiles(1)%nlevels, &
          opts,                &
          0_jpim)
  THROW(err.NE.0)

  CATCH

END SUBROUTINE rttov_k_ad
