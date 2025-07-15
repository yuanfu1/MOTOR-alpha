! Description:
!> @file
!!   Run Taylor test
!
!> @brief
!!   Run Taylor test
!!
!! @details
!!   This subroutine runs the direct and TL models for the given profiles, and
!!   then repeatedly runs the direct model on perturbed profiles. The perturbed
!!   profiles are obtained by adding TL perturbations to the original profiles.
!!   At each iteration the TL perturbations are reduced in size by a factor of
!!   10. The test records the ratio of the difference between the perturbed
!!   and original direct model radiances and the output from the TL for each
!!   iteration. The ratio should approach 1. if the direct and TL models are
!!   consistent and the direct model is "sufficiently" linear given the size of
!!   the TL perturbation. In practice at some point the perturbations usually
!!   grow so small that numerical noise causes the ratio to diverge from 1 in
!!   the final iterations.
!!
!!   If the taylor_by_chan flag is false the test is applied per profile by
!!   summing output for all channels for each profile, otherwise the test is
!!   applied to each channel individually.
!!
!!   If the taylor_on_btrefl is false the test is applied to radiances,
!!   otherwise it is applied to BTs where they are non-zero and reflectances
!!   otherwise. If the simulation involves both visible and IR channels then
!!   this option should be used with taylor_by_chan (the code does not enforce
!!   this).
!!
!!   For PC-RTTOV simulations the test is always applied per profile and is
!!   computed on BTs of reconstructed radiances if addradrec is true, otherwise
!!   it is computed on PC scores.
!!
!!   For RTTOV-SCATT the test is always calculated on BTs for passive
!!   calculations. For radar reflectivity simulations the test is computed on
!!   both zef and azef for every level (though the output does not yet
!!   distinguish clearly the levels/variables).
!!
!!   Output is written to taylor_test.log. The rttov_test.pl script checks this
!!   output for results which look anomalous. This subroutine simply runs the
!!   test and creates the output file: it only reports an error if there is a
!!   run-time failure, not if there is an apparent inconsistency between the
!!   direct and TL according to the test results.
!!
!!
!! @param[out]    err              status on exit
!! @param[in]     path             directory in which to write the test output file
!! @param[in]     opts             options to configure the simulations
!! @param[in]     opts_scatt       RTTOV-SCATT options to configure the simulations
!! @param[in]     coefs            coefficients structure
!! @param[in]     coefs_scatt      RTTOV-SCTT coefficients structure
!! @param[in]     chanprof         specifies channels and profiles to simulate
!! @param[in]     frequencies      RTTOV-SCATT frequency list
!! @param[in]     npcscores        number of PC scores for PC-RTTOV simulations
!! @param[in]     channels_rec     reconstructed radiance channels for PC-RTTOV simulations with addradrec true
!! @param[in]     profiles         atmospheric profiles and surface variables
!! @param[in]     cld_profiles     RTTOV-SCATT cloud profiles
!! @param[in]     aer_opt_param    input aerosol optical parameters
!! @param[in]     cld_opt_param    input cloud optical parameters
!! @param[in]     calcemis_in      flags indicating whether RTTOV or user should supply emissivities
!! @param[in]     calcrefl_in      flags indicating whether RTTOV or user should supply BRDFs
!! @param[in]     emissivity_in    input surface emissivities
!! @param[in]     reflectance_in   input surface BRDFs
!! @param[in]     radar            flag to indicate RTTOV-SCATT radar reflectivity calculations
!! @param[in]     taylor_by_chan   flag to select whether to apply test per channel or per profile
!! @param[in]     taylor_on_btrefl flag to select whether to apply test to BTs/refls or radiances
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
SUBROUTINE rttov_taylor_test(err,             &
                             path,            &
                             opts,            &
                             opts_scatt,      &
                             coefs,           &
                             coefs_scatt,     &
                             chanprof,        &
                             frequencies,     &
                             npcscores,       &
                             channels_rec,    &
                             profiles,        &
                             cld_profiles,    &
                             aer_opt_param,   &
                             cld_opt_param,   &
                             calcemis_in,     &
                             calcrefl_in,     &
                             emissivity_in,   &
                             reflectance_in,  &
                             radar,           &
                             taylor_by_chan,  &
                             taylor_on_btrefl)
#include "throw.h"

  USE parkind1, ONLY : jpim, jplm

  USE rttov_types, ONLY :  &
      rttov_options,       &
      rttov_options_scatt, &
      rttov_chanprof,      &
      rttov_coefs,         &
      rttov_scatt_coef,    &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_emissivity,    &
      rttov_reflectance,   &
      rttov_opt_param
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_ir, sensor_id_hi, surftype_sea

  USE parkind1, ONLY : jprb

  USE rttov_types, ONLY : &
      rttov_transmission, &
      rttov_radiance,     &
      rttov_reflectivity, &
      rttov_pccomp

  USE rttov_lun, ONLY : rttov_get_lun, rttov_put_lun
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),             INTENT(OUT)   :: err
  CHARACTER(LEN=*),          INTENT(IN)    :: path
  TYPE(rttov_options),       INTENT(IN)    :: opts
  TYPE(rttov_options_scatt), INTENT(IN)    :: opts_scatt
  TYPE(rttov_coefs),         INTENT(IN)    :: coefs
  TYPE(rttov_scatt_coef),    INTENT(IN)    :: coefs_scatt
  TYPE(rttov_chanprof),      INTENT(IN)    :: chanprof(:)
  INTEGER(jpim),             INTENT(IN)    :: frequencies(:)
  INTEGER(jpim),             INTENT(IN)    :: npcscores
  INTEGER(jpim),             INTENT(IN)    :: channels_rec(:)
  TYPE(rttov_profile),       INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile_cloud), POINTER       :: cld_profiles(:)
  TYPE(rttov_opt_param),     INTENT(IN)    :: aer_opt_param
  TYPE(rttov_opt_param),     INTENT(IN)    :: cld_opt_param
  LOGICAL(jplm),             INTENT(IN)    :: calcemis_in(SIZE(chanprof))
  LOGICAL(jplm),             INTENT(IN)    :: calcrefl_in(SIZE(chanprof))
  TYPE(rttov_emissivity),    INTENT(IN)    :: emissivity_in(SIZE(chanprof))
  TYPE(rttov_reflectance),   INTENT(IN)    :: reflectance_in(SIZE(chanprof))
  LOGICAL(jplm),             INTENT(IN)    :: radar
  LOGICAL(jplm),             INTENT(IN)    :: taylor_by_chan
  LOGICAL(jplm),             INTENT(IN)    :: taylor_on_btrefl
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_scale_profile_inc.interface"
#include "rttov_scale_cld_profile_inc.interface"
#include "rttov_scale_opt_param_inc.interface"
#include "rttov_make_profile_inc.interface"
#include "rttov_make_cld_profile_inc.interface"
#include "rttov_make_opt_param_inc.interface"
#include "rttov_make_emisrefl_inc.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_alloc_tl.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_copy_prof.interface"
#include "rttov_copy_scatt_prof.interface"
#include "rttov_copy_opt_param.interface"
#include "rttov_add_prof.interface"
#include "rttov_add_scatt_prof.interface"
#include "rttov_add_opt_param.interface"
#include "rttov_direct.interface"
#include "rttov_tl.interface"
#include "rttov_scatt.interface"
#include "rttov_scatt_tl.interface"

  CHARACTER(LEN=*), PARAMETER        :: taylor_filename = 'taylor_test.log'
  INTEGER(jpim),    PARAMETER        :: niter = 8

  LOGICAL(jplm)                      :: calcemis(SIZE(chanprof))
  TYPE(rttov_emissivity)             :: emissivity(SIZE(chanprof))
  TYPE(rttov_emissivity)             :: emissivity_tl(SIZE(chanprof))
  TYPE(rttov_emissivity)             :: emissivity_pert(SIZE(chanprof))

  LOGICAL(jplm)                      :: calcrefl(SIZE(chanprof))
  TYPE(rttov_reflectance)            :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance)            :: reflectance_tl(SIZE(chanprof))
  TYPE(rttov_reflectance)            :: reflectance_pert(SIZE(chanprof))

  TYPE(rttov_profile),       POINTER :: profiles_pert(:)
  TYPE(rttov_profile),       POINTER :: profiles_tl(:)
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles_pert(:)
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles_tl(:)
  TYPE(rttov_opt_param)              :: aer_opt_param_pert
  TYPE(rttov_opt_param)              :: aer_opt_param_tl
  TYPE(rttov_opt_param)              :: cld_opt_param_pert
  TYPE(rttov_opt_param)              :: cld_opt_param_tl
  TYPE(rttov_transmission)           :: transmission
  TYPE(rttov_transmission)           :: transmission_tl
  TYPE(rttov_radiance)               :: radiance_saved
  TYPE(rttov_radiance)               :: radiance
  TYPE(rttov_radiance)               :: radiance_tl
  TYPE(rttov_reflectivity), POINTER  :: reflectivity_saved => NULL()
  TYPE(rttov_reflectivity), POINTER  :: reflectivity       => NULL()
  TYPE(rttov_reflectivity), POINTER  :: reflectivity_tl    => NULL()
  TYPE(rttov_pccomp)                 :: pccomp_saved
  TYPE(rttov_pccomp)                 :: pccomp
  TYPE(rttov_pccomp)                 :: pccomp_tl

  CHARACTER(LEN=256)                 :: output_filename
  REAL(jprb), ALLOCATABLE            :: taylor_direct(:,:), taylor_tl(:,:)
  REAL(jprb)                         :: ratio
  INTEGER(jpim)                      :: nchanprof, nprofiles, nlevels, nchannels_rec, nphangle, ntaylor
  INTEGER(jpim)                      :: i, j, k, lev, chan, prof, iter, file_id
  LOGICAL(jplm)                      :: do_rttovscatt

  TRY

  nchanprof     = SIZE(chanprof)
  nprofiles     = SIZE(profiles)
  nlevels       = profiles(1)%nlevels
  nchannels_rec = SIZE(channels_rec)

  do_rttovscatt = ASSOCIATED(cld_profiles)
  IF (do_rttovscatt .AND. radar) ALLOCATE(reflectivity_saved, reflectivity, reflectivity_tl)
  IF (taylor_by_chan .AND. .NOT. opts%rt_ir%pc%addpc) THEN
    ntaylor = nchanprof
  ELSE
    ntaylor = nprofiles
  ENDIF

  CALL rttov_alloc_direct(                               &
              err,                                       &
              1_jpim,                                    &
              nprofiles,                                 &
              nchanprof,                                 &
              nlevels,                                   &
              opts = opts,                               &
              coefs = coefs,                             &
              radiance = radiance_saved,                 &
              npcscores = npcscores * nprofiles,         &
              nchannels_rec = nchannels_rec * nprofiles, &
              pccomp = pccomp_saved,                     &
              reflectivity = reflectivity_saved,         &
              init = .TRUE._jplm)
  THROWM(err.NE.0, "Allocation for Taylor test")

  CALL rttov_alloc_tl(                                   &
              err,                                       &
              1_jpim,                                    &
              nprofiles,                                 &
              nchanprof,                                 &
              nlevels,                                   &
              opts = opts,                               &
              profiles = profiles_pert,                  &
              profiles_tl = profiles_tl,                 &
              coefs = coefs,                             &
              transmission = transmission,               &
              transmission_tl = transmission_tl,         &
              radiance = radiance,                       &
              radiance_tl = radiance_tl,                 &
              npcscores = npcscores * nprofiles,         &
              nchannels_rec = nchannels_rec * nprofiles, &
              pccomp = pccomp,                           &
              pccomp_tl = pccomp_tl,                     &
              reflectivity = reflectivity,               &
              reflectivity_tl = reflectivity_tl,         &
              init = .TRUE._jplm)
  THROWM(err.NE.0, "Allocation for Taylor test")

  IF (opts%rt_ir%user_aer_opt_param) THEN
    nphangle = SIZE(aer_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, aer_opt_param_pert, nchanprof, nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROWM(err.NE.0, "Allocation for Taylor test")

    aer_opt_param_pert%phangle = aer_opt_param%phangle
    CALL rttov_init_opt_param(err, opts, aer_opt_param_pert)
    THROWM(err.NE.0, "Allocation for Taylor test")

    CALL rttov_alloc_opt_param(err, aer_opt_param_tl, nchanprof, nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROWM(err.NE.0, "Allocation for Taylor test")
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    nphangle = SIZE(cld_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, cld_opt_param_pert, nchanprof, nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROWM(err.NE.0, "Allocation for Taylor test")

    cld_opt_param_pert%phangle = cld_opt_param%phangle
    CALL rttov_init_opt_param(err, opts, cld_opt_param_pert)
    THROWM(err.NE.0, "Allocation for Taylor test")

    CALL rttov_alloc_opt_param(err, cld_opt_param_tl, nchanprof, nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROWM(err.NE.0, "Allocation for Taylor test")
  ENDIF

  IF (do_rttovscatt) THEN
    ALLOCATE(cld_profiles_pert(nprofiles), cld_profiles_tl(nprofiles), STAT=err)
    THROWM(err.NE.0, "Allocation for Taylor test")

    CALL rttov_alloc_scatt_prof(err, nprofiles, cld_profiles_pert, nlevels, &
        cld_profiles(1)%nhydro, cld_profiles(1)%nhydro_frac, 1_jpim, init=.TRUE._jplm)
    THROWM(err.NE.0, "Allocation for Taylor test")

    CALL rttov_alloc_scatt_prof(err, nprofiles, cld_profiles_tl, nlevels, &
        cld_profiles(1)%nhydro, cld_profiles(1)%nhydro_frac, 1_jpim, init=.TRUE._jplm)
    THROWM(err.NE.0, "Allocation for Taylor test")
  ENDIF

  calcemis = calcemis_in
  calcrefl = calcrefl_in
  emissivity = emissivity_in
  reflectance = reflectance_in

  ! Direct call on given profiles
  IF (do_rttovscatt) THEN
    CALL rttov_scatt( &
        err,               &
        opts_scatt,        &
        nlevels,           &
        chanprof,          &
        frequencies,       &
        profiles,          &
        cld_profiles,      &
        coefs,             &
        coefs_scatt,       &
        calcemis,          &
        emissivity,        &
        radiance_saved,    &
        reflectivity = reflectivity_saved)
  ELSE
    CALL rttov_direct(                  &
        err,                            &
        chanprof,                       &
        opts,                           &
        profiles,                       &
        coefs,                          &
        transmission,                   &
        radiance_saved,                 &
        calcemis = calcemis,            &
        emissivity = emissivity,        &
        calcrefl = calcrefl,            &
        reflectance = reflectance,      &
        aer_opt_param = aer_opt_param,  &
        cld_opt_param = cld_opt_param,  &
        pccomp       = pccomp_saved,    &
        channels_rec = channels_rec)
  ENDIF
  THROWM(err.NE.0, "Error calling RTTOV in Taylor test")

  ! Create TL perturbations
  CALL rttov_make_profile_inc(profiles_tl, profiles, opts)
  IF (.NOT. opts%dev%no_opt_param_tladk) THEN
    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_make_opt_param_inc(aer_opt_param_tl, aer_opt_param)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_make_opt_param_inc(cld_opt_param_tl, cld_opt_param)
  ELSE
    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_init_opt_param(err, opts, aer_opt_param_tl, zero_only=.TRUE._jplm)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_init_opt_param(err, opts, cld_opt_param_tl, zero_only=.TRUE._jplm)
  ENDIF
  IF (do_rttovscatt) CALL rttov_make_cld_profile_inc(cld_profiles_tl, cld_profiles)

  CALL rttov_make_emisrefl_inc(opts, coefs, profiles, &
        chanprof, calcemis, &
        emissivity_tl%emis_in, reflectance_tl%refl_in, &
        reflectance%diffuse_refl_in, reflectance_tl%diffuse_refl_in, &
        emissivity%specularity, emissivity_tl%specularity)

  ! TL call on given profiles
  IF (do_rttovscatt) THEN
    CALL rttov_scatt_tl( &
        err,               &
        opts_scatt,        &
        nlevels,           &
        chanprof,          &
        frequencies,       &
        profiles,          &
        cld_profiles,      &
        coefs,             &
        coefs_scatt,       &
        calcemis,          &
        emissivity,        &
        profiles_tl,       &
        cld_profiles_tl,   &
        emissivity_tl,     &
        radiance,          &
        radiance_tl,       &
        reflectivity,      &
        reflectivity_tl)
  ELSE
    CALL rttov_tl(                           &
        err,                                 &
        chanprof,                            &
        opts,                                &
        profiles,                            &
        profiles_tl,                         &
        coefs,                               &
        transmission,                        &
        transmission_tl,                     &
        radiance,                            &
        radiance_tl,                         &
        calcemis         = calcemis,         &
        emissivity       = emissivity,       &
        emissivity_tl    = emissivity_tl,    &
        calcrefl         = calcrefl,         &
        reflectance      = reflectance,      &
        reflectance_tl   = reflectance_tl,   &
        aer_opt_param    = aer_opt_param,    &
        aer_opt_param_tl = aer_opt_param_tl, &
        cld_opt_param    = cld_opt_param,    &
        cld_opt_param_tl = cld_opt_param_tl, &
        pccomp           = pccomp,           &
        pccomp_tl        = pccomp_tl,        &
        channels_rec     = channels_rec)
  ENDIF
  THROWM(err.NE.0, "Error calling RTTOV in Taylor test")


  ! Sort out emissivity/reflectance for Taylor iterations

  ! Where calcemis is true, but the emissivity does not depend on any active
  ! profile variables it is necessary to perturb the emissivities manually
  ! by recording the direct model emissivity value and setting calcemis to false.
  ! The same applies to surface reflectance.

  ! With the do_lambertian option we have an issue with ISEM over sea-surfaces:
  ! do_lambertian will be false for sea profiles with calcemis true, but if we
  ! set calcemis false for the direct model iterations then do_lambertian is
  ! activated. The simplest option is to zero the emissivity perturbation when
  ! ISEM and do_lambertian are selected together: this is done in
  ! rttov_make_emisrefl_inc. Then we avoid setting calcemis to false for this
  ! particular case.

  IF (coefs%coef%id_sensor == sensor_id_ir .OR. &
      coefs%coef%id_sensor == sensor_id_hi) THEN

    DO i = 1, nchanprof
      prof = chanprof(i)%prof
      chan = chanprof(i)%chan

      ! For land/sea-ice profiles or when ISEM is used
      IF (profiles(prof)%skin%surftype /= surftype_sea .OR. &
          (opts%rt_ir%ir_sea_emis_model == 1_jpim .AND. .NOT. &
           (opts%rt_all%do_lambertian .AND. &
            profiles(prof)%skin%surftype == surftype_sea))) THEN
        ! Save the calculated emissivity for this channel and turn calcemis
        ! off so we can supply the perturbed emissivity to rttov_direct
        emissivity(i)%emis_in = emissivity(i)%emis_out
        calcemis(i) = .FALSE.
      ENDIF

      ! For VIS/NIR channels over land/sea-ice
      IF (coefs%coef%ss_val_chn(chan) == 2_jpim .AND. &
          profiles(prof)%skin%surftype /= surftype_sea) THEN
        ! Save the calculated reflectance for this channel and turn calcrefl
        ! off so we can supply the perturbed reflectance to rttov_direct
        reflectance(i)%refl_in = reflectance(i)%refl_out
        calcrefl(i) = .FALSE.

        ! In this case, for the diffuse reflectance we must use the calculated
        ! value only if the input reflectance was non-zero (otherwise the
        ! diffuse reflectance is computed internally)
        IF (reflectance(i)%diffuse_refl_in > 0._jprb) &
          reflectance(i)%diffuse_refl_in = reflectance(i)%diffuse_refl_out
      ENDIF
    ENDDO

  ENDIF

  ! If TAYLOR_BY_CHAN is set Taylor ratios are calculated for every channel individually.
  ! Otherwise the ratios are calculated based on radiances/BTs/refls summed for each profile.
  ! For PC-RTTOV the Taylor test is always calculated per profile.
  ! For radar reflectivities the Taylor ratios are computed for each level for zef and azef.

  IF (do_rttovscatt .AND. radar) THEN
    ALLOCATE(taylor_direct(niter,2*nlevels*ntaylor), taylor_tl(niter,2*nlevels*ntaylor))
  ELSE
    ALLOCATE(taylor_direct(niter,ntaylor), taylor_tl(niter,ntaylor))
  ENDIF


  ! Run Taylor test: perturb input profile, emissivity, reflectance,
  ! call direct model and store results

  DO iter = 1, niter

    ! Scale the profile increments in place by a factor of 0.1 at each iteration
    IF (iter > 1) THEN
      CALL rttov_scale_profile_inc(profiles_tl, 0.1_jprb)
      IF (.NOT. opts%dev%no_opt_param_tladk) THEN
        IF (opts%rt_ir%user_aer_opt_param) CALL rttov_scale_opt_param_inc(aer_opt_param_tl, 0.1_jprb)
        IF (opts%rt_ir%user_cld_opt_param) CALL rttov_scale_opt_param_inc(cld_opt_param_tl, 0.1_jprb)
      ENDIF
      IF (do_rttovscatt) CALL rttov_scale_cld_profile_inc(cld_profiles_tl, 0.1_jprb)
      emissivity_tl(:)%emis_in = emissivity_tl(:)%emis_in * 0.1_jprb
      emissivity_tl(:)%specularity = emissivity_tl(:)%specularity * 0.1_jprb
      reflectance_tl(:)%refl_in = reflectance_tl(:)%refl_in * 0.1_jprb
      reflectance_tl(:)%diffuse_refl_in = reflectance_tl(:)%diffuse_refl_in * 0.1_jprb
    ENDIF

    CALL rttov_copy_prof(profiles_pert, profiles)
    CALL rttov_add_prof(profiles_pert, profiles, profiles_tl)
    IF (opts%rt_ir%user_aer_opt_param) THEN
      CALL rttov_copy_opt_param(aer_opt_param_pert, aer_opt_param)
      CALL rttov_add_opt_param(aer_opt_param_pert, aer_opt_param, aer_opt_param_tl)
    ENDIF
    IF (opts%rt_ir%user_cld_opt_param) THEN
      CALL rttov_copy_opt_param(cld_opt_param_pert, cld_opt_param)
      CALL rttov_add_opt_param(cld_opt_param_pert, cld_opt_param, cld_opt_param_tl)
    ENDIF
    IF (do_rttovscatt) THEN
      CALL rttov_copy_scatt_prof(cld_profiles_pert, cld_profiles)
      CALL rttov_add_scatt_prof(cld_profiles_pert, cld_profiles, cld_profiles_tl)
    ENDIF

    emissivity_pert(:)%emis_in = emissivity(:)%emis_in + emissivity_tl(:)%emis_in
    emissivity_pert(:)%specularity = emissivity(:)%specularity + emissivity_tl(:)%specularity
    reflectance_pert(:)%refl_in = reflectance(:)%refl_in + reflectance_tl(:)%refl_in
    reflectance_pert(:)%diffuse_refl_in = reflectance(:)%diffuse_refl_in + reflectance_tl(:)%diffuse_refl_in
    reflectance_pert(:)%refl_cloud_top = reflectance(:)%refl_cloud_top

    IF (do_rttovscatt) THEN
      CALL rttov_scatt( &
          err,               &
          opts_scatt,        &
          nlevels,           &
          chanprof,          &
          frequencies,       &
          profiles_pert,     &
          cld_profiles_pert, &
          coefs,             &
          coefs_scatt,       &
          calcemis,          &
          emissivity_pert,   &
          radiance,          &
          reflectivity = reflectivity)
    ELSE
      CALL rttov_direct(                       &
          err,                                 &
          chanprof,                            &
          opts,                                &
          profiles_pert,                       &
          coefs,                               &
          transmission,                        &
          radiance,                            &
          calcemis      = calcemis,            &
          emissivity    = emissivity_pert,     &
          calcrefl      = calcrefl,            &
          reflectance   = reflectance_pert,    &
          aer_opt_param = aer_opt_param_pert,  &
          cld_opt_param = cld_opt_param_pert,  &
          pccomp        = pccomp,              &
          channels_rec  = channels_rec)
    ENDIF
    THROWM(err.NE.0, "Error calling RTTOV in Taylor test")

    IF (opts%rt_ir%pc%addpc) THEN
      DO i = 1, nprofiles
        IF (opts%rt_ir%pc%addradrec) THEN
          j = (i - 1) * nchannels_rec + 1
          taylor_direct(iter,i) = SUM(pccomp%bt_pccomp(j:j+nchannels_rec-1) - &
                                      pccomp_saved%bt_pccomp(j:j+nchannels_rec-1))
          taylor_tl(iter,i) = SUM(pccomp_tl%bt_pccomp(j:j+nchannels_rec-1))
        ELSE
          j = (i - 1) * npcscores + 1
          taylor_direct(iter,i) = SUM(pccomp%total_pcscores(j:j+npcscores-1) - &
                                      pccomp_saved%total_pcscores(j:j+npcscores-1))
          taylor_tl(iter,i) = SUM(pccomp_tl%total_pcscores(j:j+npcscores-1))
        ENDIF
      ENDDO
    ELSE
      IF (taylor_by_chan) THEN
        DO i = 1, nchanprof
          IF (do_rttovscatt .AND. radar) THEN
            DO lev = 1, nlevels
              k = (i - 1) * nlevels + lev
              taylor_direct(iter,k) = reflectivity%zef(lev,i) - reflectivity_saved%zef(lev,i)
              taylor_tl(iter,k) = reflectivity_tl%zef(lev,i)
            ENDDO
            DO lev = 1, nlevels
              k = nlevels * nchanprof + (i - 1) * nlevels + lev
              taylor_direct(iter,k) = reflectivity%azef(lev,i) - reflectivity_saved%azef(lev,i)
              taylor_tl(iter,k) = reflectivity_tl%azef(lev,i)
            ENDDO
          ELSEIF (taylor_on_btrefl .OR. do_rttovscatt) THEN
            IF (radiance%bt(i) > 0._jprb) THEN
              taylor_direct(iter,i) = radiance%bt(i) - radiance_saved%bt(i)
              taylor_tl(iter,i) = radiance_tl%bt(i)
            ELSE !IF (radiance%refl(i) > 0._jprb) THEN
              taylor_direct(iter,i) = radiance%refl(i) - radiance_saved%refl(i)
              taylor_tl(iter,i) = radiance_tl%refl(i)
            ENDIF
          ELSE
            taylor_direct(iter,i) = radiance%total(i) - radiance_saved%total(i)
            taylor_tl(iter,i) = radiance_tl%total(i)
          ENDIF
        ENDDO
      ELSE
        j = 1
        taylor_direct(iter,:) = 0._jprb
        taylor_tl(iter,:) = 0._jprb
        DO i = 1, nprofiles
          DO WHILE (chanprof(j)%prof == i)
            IF (do_rttovscatt .AND. radar) THEN
              DO lev = 1, nlevels
                k = (i - 1) * nlevels + lev
                taylor_direct(iter,k) = taylor_direct(iter,k) + reflectivity%zef(lev,j) - &
                                                                reflectivity_saved%zef(lev,j)
                taylor_tl(iter,k) = taylor_tl(iter,k) + reflectivity_tl%zef(lev,j)
              ENDDO
              DO lev = 1, nlevels
                k = nlevels * nprofiles + (i - 1) * nlevels + lev
                taylor_direct(iter,k) = taylor_direct(iter,k) + reflectivity%azef(lev,j) - &
                                                                reflectivity_saved%azef(lev,j)
                taylor_tl(iter,k) = taylor_tl(iter,k) + reflectivity_tl%azef(lev,j)
              ENDDO
            ELSEIF (taylor_on_btrefl .OR. do_rttovscatt) THEN
              IF (radiance%bt(i) > 0._jprb) THEN
                taylor_direct(iter,i) = taylor_direct(iter,i) + radiance%bt(j) - &
                                                                radiance_saved%bt(j)
                taylor_tl(iter,i) = taylor_tl(iter,i) + radiance_tl%bt(j)
              ELSE !IF (radiance%refl(j) > 0._jprb) THEN
                taylor_direct(iter,i) = taylor_direct(iter,i) + radiance%refl(j) - &
                                                                radiance_saved%refl(j)
                taylor_tl(iter,i) = taylor_tl(iter,i) + radiance_tl%refl(j)
              ENDIF
            ELSE
              taylor_direct(iter,i) = taylor_direct(iter,i) + radiance%total(j) - &
                                                              radiance_saved%total(j)
              taylor_tl(iter,i) = taylor_tl(iter,i) + radiance_tl%total(j)
            ENDIF
            j = j + 1
            IF (j > SIZE(chanprof)) EXIT
          ENDDO
          IF (j > SIZE(chanprof)) EXIT
        ENDDO
      ENDIF
    ENDIF
  ENDDO


  ! Write output to file

  CALL rttov_get_lun(file_id)
  output_filename = TRIM(path)//'/'//TRIM(taylor_filename)
  OPEN(file_id, file = output_filename, form = 'formatted', iostat = err)
  THROWM(err.NE.0,"Cannot open "//TRIM(output_filename))

  WRITE(file_id, '(a,i6)') 'nprof: ', SIZE(taylor_direct(1,:))
  WRITE(file_id, '(a,i6)') 'niter: ', niter
  WRITE(file_id, '(a)') 'Taylor test: profile, channel, iteration, (H(x*)-H(x))/H(x*-x)'
  DO i = 1, SIZE(taylor_direct(1,:))

    IF (do_rttovscatt .AND. radar) THEN
      ! prof/chanprof index for radar case
      k = (i - 1) / nlevels + 1
      IF (i > nlevels * ntaylor) THEN
        k = (i - nlevels * ntaylor - 1) / nlevels + 1
      ENDIF
    ELSE
      k = i
    ENDIF

    DO iter = 1, niter
      IF (do_rttovscatt .AND. radar) THEN
        IF (ABS(taylor_tl(iter,i)) > 0._jprb) THEN
          ratio = taylor_direct(iter,i) * (10._jprb**(iter-1)) / taylor_tl(iter,i)
        ELSE
          ratio = 1._jprb
        ENDIF
      ELSE
        IF (ABS(taylor_tl(iter,i)) > 0._jprb) THEN
          ratio = taylor_direct(iter,i) * (10._jprb**(iter-1)) / taylor_tl(iter,i)
        ELSE
          ratio = 1._jprb
        ENDIF
      ENDIF
      IF (taylor_by_chan .AND. .NOT. opts%rt_ir%pc%addpc) THEN
        WRITE(file_id, '(14x,i5,4x,i5,4x,i6,5x,7f26.16)') &
            chanprof(k)%prof, chanprof(k)%chan, iter, ratio
      ELSE
        WRITE(file_id, '(14x,i5,4x,a9,i6,5x,7f26.16)') &
            k, '   all   ', iter, ratio
      ENDIF
    ENDDO
    WRITE(file_id, '(a)') ''
  ENDDO

  CLOSE(file_id, iostat = err)
  THROWM(err.NE.0,"Cannot close "//TRIM(output_filename))
  CALL rttov_put_lun(file_id)


  ! Tidy up

  DEALLOCATE(taylor_direct, taylor_tl)

  IF (ASSOCIATED(cld_profiles)) THEN
    CALL rttov_alloc_scatt_prof(err, nprofiles, cld_profiles_pert, nlevels, &
        cld_profiles(1)%nhydro, cld_profiles(1)%nhydro_frac, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")

    CALL rttov_alloc_scatt_prof(err, nprofiles, cld_profiles_tl, nlevels, &
        cld_profiles(1)%nhydro, cld_profiles(1)%nhydro_frac, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")

    DEALLOCATE(cld_profiles_pert, cld_profiles_tl, STAT=err)
    THROWM(err.NE.0, "Deallocation for Taylor test")
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    CALL rttov_alloc_opt_param(err, cld_opt_param_tl, nchanprof, nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")

    CALL rttov_alloc_opt_param(err, cld_opt_param_pert, nchanprof, nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")
  ENDIF

  IF (opts%rt_ir%user_aer_opt_param) THEN
    CALL rttov_alloc_opt_param(err, aer_opt_param_tl, nchanprof, nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")

    CALL rttov_alloc_opt_param(err, aer_opt_param_pert, nchanprof, nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROWM(err.NE.0, "Deallocation for Taylor test")
  ENDIF

  CALL rttov_alloc_tl(                                   &
              err,                                       &
              0_jpim,                                    &
              nprofiles,                                 &
              nchanprof,                                 &
              nlevels,                                   &
              opts = opts,                               &
              profiles = profiles_pert,                  &
              profiles_tl = profiles_tl,                 &
              coefs = coefs,                             &
              transmission = transmission,               &
              transmission_tl = transmission_tl,         &
              radiance = radiance,                       &
              radiance_tl = radiance_tl,                 &
              npcscores = npcscores * nprofiles,         &
              nchannels_rec = nchannels_rec * nprofiles, &
              pccomp = pccomp,                           &
              reflectivity = reflectivity,               &
              reflectivity_tl = reflectivity_tl,         &
              pccomp_tl = pccomp_tl)
  THROWM(err.NE.0, "Deallocation for Taylor test")

  CALL rttov_alloc_direct(                               &
              err,                                       &
              0_jpim,                                    &
              nprofiles,                                 &
              nchanprof,                                 &
              nlevels,                                   &
              opts = opts,                               &
              coefs = coefs,                             &
              radiance = radiance_saved,                 &
              npcscores = npcscores * nprofiles,         &
              nchannels_rec = nchannels_rec * nprofiles, &
              reflectivity = reflectivity_saved,         &
              pccomp = pccomp_saved)
  THROWM(err.NE.0, "Deallocation for Taylor test")

  IF (do_rttovscatt .AND. radar) DEALLOCATE(reflectivity_saved, reflectivity, reflectivity_tl)

  CATCH
END SUBROUTINE rttov_taylor_test
