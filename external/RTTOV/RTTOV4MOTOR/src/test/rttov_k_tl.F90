! Description:
!> @file
!!   Compute Jacobian using TL model
!
!> @brief
!!   Compute Jacobian using TL model
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
SUBROUTINE rttov_k_tl( &
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

#include "rttov_errorreport.interface"
#include "rttov_tl.interface"
#include "rttov_parallel_tl.interface"
#include "rttov_scatt_tl.interface"
#include "rttov_parallel_scatt_tl.interface"
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_alloc_opt_param.interface"
#include "rttov_alloc_pccomp.interface"
#include "rttov_alloc_traj.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_reflectivity.interface"
#include "rttov_init_emis_refl.interface"
#include "rttov_init_opt_param.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_reflectivity.interface"
#include "rttov_init_pccomp.interface"

  TYPE(rttov_traj)                  :: traj, traj_tl
  TYPE(rttov_pccomp)                :: pccomp_tl
  TYPE(rttov_profile),       TARGET :: profiles_tl(SIZE(profiles))
  TYPE(rttov_profile_cloud), TARGET :: cld_profiles_tl(SIZE(profiles))
  TYPE(rttov_opt_param),     TARGET :: aer_opt_param_tl
  TYPE(rttov_opt_param),     TARGET :: cld_opt_param_tl
  TYPE(rttov_opt_param),    POINTER :: opt_param, opt_param_tl, opt_param_k
  TYPE(rttov_emissivity)            :: emissivity_tl(SIZE(chanprof))
  TYPE(rttov_reflectance)           :: reflectance_tl(SIZE(chanprof))
  TYPE(rttov_transmission)          :: transmission_tl
  TYPE(rttov_radiance)              :: radiance_tl
  TYPE(rttov_reflectivity), POINTER :: reflectivity_tl => NULL()
  REAL(KIND=jprb)                   :: out_tl(SIZE(chanprof))
  TYPE(chain),  POINTER             :: chain_profiles_tl(:)
  TYPE(chain),  POINTER             :: chain_cld_profiles_tl(:)
  TYPE(chain),  POINTER             :: chain_profiles_k(:)
  TYPE(chain),  POINTER             :: chain_cld_profiles_k(:)
  TYPE(chain),  POINTER             :: chain_profiles_k_pc(:)
  TYPE(chain),  POINTER             :: chain_profiles_k_rec(:)
  TYPE(pchain), POINTER             :: c_tl(:)
  TYPE(pchain), POINTER             :: c_cld_tl(:)
  TYPE(pchain), POINTER             :: c_k(:)
  TYPE(pchain), POINTER             :: c_cld_k(:)
  TYPE(pchain), POINTER             :: c_k_pc(:)
  TYPE(pchain), POINTER             :: c_k_in(:)
  INTEGER(KIND=jpim)                :: i, nsize, ichan, iprof, lay, var, aercld
  REAL(KIND=jprb), POINTER          :: x
  REAL(KIND=jprb), POINTER          :: var_in(:), var_out(:)
  INTEGER(KIND=jpim)                :: nprofiles, nchannels, nchannels_rec, npcscores, nlevels, nlayers, nphangle
  LOGICAL(KIND=jplm)                :: do_rttovscatt, do_radar

  TRY

  nprofiles = SIZE(profiles)
  nchannels = SIZE(chanprof)
  nlevels = profiles(1)%nlevels
  nlayers = nlevels - 1

  do_rttovscatt = ASSOCIATED(cld_profiles)
  do_radar = do_rttovscatt .AND. ASSOCIATED(reflectivity%zef)

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_tl,         &
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
          cld_profiles_tl,              &
          nlevels,                      &
          cld_profiles(1)%nhydro,       &
          cld_profiles(1)%nhydro_frac,  &
          1_jpim,                       &
          init = .TRUE._jplm)
    THROW(err.NE.0)
    IF (do_radar) THEN
      ALLOCATE(reflectivity_tl)
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_tl,         &
            nlevels,                 &
            1_jpim,                  &
            init = .TRUE._jplm)
      THROW(err.NE.0)
    ENDIF
  ENDIF

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_tl,         &
          profiles(1)%nlevels, &
          1_jpim)
  THROW(err.NE.0)

  CALL rttov_alloc_transmission(err, transmission_tl, profiles(1)%nlevels, nchannels, 1_jpim)
  THROW(err.NE.0)

  IF (opts%rt_ir%user_aer_opt_param) THEN
    nphangle = SIZE(aer_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, aer_opt_param_tl, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    nphangle = SIZE(cld_opt_param%phangle)
    CALL rttov_alloc_opt_param(err, cld_opt_param_tl, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 1_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                   &
          err,                  nprofiles,   &
          nchannels,            opts,        &
          profiles(1)%nlevels,  coefs,       &
          1_jpim,               traj = traj, &
          traj_tl = traj_tl)
    THROW(err.NE.0)
  ENDIF

!
! Chain input profiles
!
  CALL make_chain_profile(err, chain_profiles_tl, c_tl, "PROFILES_TL", profiles_tl, .FALSE._jplm)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL make_chain_cld_profile(err, chain_cld_profiles_tl, c_cld_tl, "CLD_PROFILES_TL", cld_profiles_tl, .FALSE._jplm)
    THROW(err.NE.0)
  ENDIF

!
! Regular radiances Jacobians
!
  CALL make_chain_profile(err, chain_profiles_k, c_k, "PROFILES_K", profiles_k, .TRUE._jplm)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL make_chain_cld_profile(err, chain_cld_profiles_k, c_cld_k, "CLD_PROFILES_K", cld_profiles_k, .TRUE._jplm)
    THROW(err.NE.0)
  ENDIF

  NULLIFY (chain_profiles_k_rec, c_k_in)
  nchannels_rec = -1_jpim

  NULLIFY (chain_profiles_k_pc, c_k_pc)
  npcscores = -1_jpim

  IF (opts%rt_ir%pc%addpc) THEN

    npcscores = SIZE(pccomp%total_pcscores)

    IF (opts%rt_ir%pc%addradrec) THEN
!
! Reconstructed radiances Jacobians
!
      nchannels_rec = SIZE(pccomp%bt_pccomp)
      CALL make_chain_profile(err, chain_profiles_k_rec, c_k_in, "PROFILES_K_REC", profiles_k_rec, .TRUE._jplm)
      THROW(err.NE.0)
    ELSE
!
! PC-scores Jacobians
!
      CALL make_chain_profile(err, chain_profiles_k_pc, c_k_pc, "PROFILES_K_PC", profiles_k_pc, .TRUE._jplm)
      THROW(err.NE.0)

    ENDIF

!
! Allocate a pccomp_tl
!
    CALL rttov_alloc_pccomp(err, pccomp_tl, npcscores, 1_jpim, nchannels_rec = nchannels_rec)
    THROW(err.NE.0)
  ENDIF

  CALL size_chain(nsize, chain_profiles_tl(1))! assume all profiles have the same size

  DO i = 1, nsize
    CALL rttov_init_emis_refl(emissivity_tl, reflectance_tl)
    IF (opts%rt_ir%user_aer_opt_param) CALL rttov_init_opt_param(err, opts, aer_opt_param_tl, .TRUE._jplm)
    IF (opts%rt_ir%user_cld_opt_param) CALL rttov_init_opt_param(err, opts, cld_opt_param_tl, .TRUE._jplm)
    CALL rttov_init_rad(radiance_tl)
    IF (opts%rt_ir%pc%addpc) CALL rttov_init_pccomp(pccomp_tl)

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c_tl(iprof)%p, x)
      x = 1._jprb
    ENDDO

    CALL call_rttov_tl(err)
    THROW(err.NE.0)

    DO iprof = 1, nprofiles
      CALL get_pointer_chain(c_tl(iprof)%p, x)
      x = 0._jprb
      CALL advance_chain(c_tl(iprof)%p)
    ENDDO

    IF (do_radar) THEN
      CALL assign_chain_profile(c_k, reflectivity_tl%azef(radar_k_lev,:))
    ELSEIF (do_rttovscatt .OR. opts%rt_all%switchrad) THEN
      IF (opts%rt_ir%pc%addpc) THEN
        CALL assign_chain_profile(c_k, radiance_tl%total)
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile(c_k_in, pccomp_tl%bt_pccomp)
        ELSE
          CALL assign_chain_profile(c_k_pc, pccomp_tl%total_pcscores)
        ENDIF
      ELSE
        DO ichan = 1, nchannels
          IF (do_rttovscatt .OR. coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
            out_tl(ichan) = radiance_tl%bt(ichan)
          ELSE
            out_tl(ichan) = radiance_tl%total(ichan)
          ENDIF
        ENDDO
        CALL assign_chain_profile(c_k, out_tl)
      ENDIF
    ELSE
      CALL assign_chain_profile(c_k, radiance_tl%total)
      IF (opts%rt_ir%pc%addpc) THEN
        IF (opts%rt_ir%pc%addradrec) THEN
          CALL assign_chain_profile(c_k_in, pccomp_tl%total_pccomp)
        ELSE
          CALL assign_chain_profile(c_k_pc, pccomp_tl%total_pcscores)
        ENDIF
      ENDIF
    ENDIF

  ENDDO ! chain loop


  IF (do_rttovscatt) THEN
    CALL size_chain(nsize, chain_cld_profiles_tl(1))! assume all profiles have the same size

    DO i = 1, nsize
      CALL rttov_init_emis_refl(emissivity_tl, reflectance_tl)
      CALL rttov_init_rad(radiance_tl)
      IF (do_radar) CALL rttov_init_reflectivity(reflectivity_tl)

      DO iprof = 1, nprofiles
        CALL get_pointer_chain(c_cld_tl(iprof)%p, x)
        x = 1._jprb
      ENDDO

      CALL call_rttov_tl(err)
      THROW(err.NE.0)

      DO iprof = 1, nprofiles
        CALL get_pointer_chain(c_cld_tl(iprof)%p, x)
        x = 0._jprb
        CALL advance_chain(c_cld_tl(iprof)%p)
      ENDDO

      IF (do_radar) THEN
        out_tl = reflectivity_tl%azef(radar_k_lev,:)
      ELSE
        out_tl = radiance_tl%bt
      ENDIF
      CALL assign_chain_profile(c_cld_k, out_tl)

    ENDDO ! cld_profile chain loop
  ENDIF ! do_rttovscatt

  CALL free_chain_profile(err, chain_profiles_tl, c_tl)
  THROW(err.NE.0)

  CALL free_chain_profile(err, chain_profiles_k, c_k)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    CALL free_chain_profile(err, chain_cld_profiles_tl, c_cld_tl)
    THROW(err.NE.0)

    CALL free_chain_profile(err, chain_cld_profiles_k, c_cld_k)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%pc%addpc) THEN
    IF (opts%rt_ir%pc%addradrec) THEN
      CALL free_chain_profile(err, chain_profiles_k_rec, c_k_in)
      THROW(err.NE.0)
    ELSE
      CALL free_chain_profile(err, chain_profiles_k_pc, c_k_pc)
      THROW(err.NE.0)
    ENDIF
  ENDIF


  ! Explicit aerosol/cloud optical properties
  ! The chain approach doesn't work for the optical properties so we perturb each variable in turn.
  ! All channels are perturbed together for each variable in each layer.
  ! For efficiency we avoid perturbing non-scattering layers since Jacobians will be zero anyway in this case.

  IF ((opts%rt_ir%user_aer_opt_param .OR. opts%rt_ir%user_cld_opt_param) .AND. &
      .NOT. opts%dev%no_opt_param_tladk .AND. .NOT. do_rttovscatt) THEN
    DO aercld = 1, 2

      IF (aercld == 1) THEN
        IF (.NOT. opts%rt_ir%user_aer_opt_param) CYCLE
        opt_param => aer_opt_param
        opt_param_tl => aer_opt_param_tl
        opt_param_k => aer_opt_param_k
      ELSE
        IF (.NOT. opts%rt_ir%user_cld_opt_param) CYCLE
        opt_param => cld_opt_param
        opt_param_tl => cld_opt_param_tl
        opt_param_k => cld_opt_param_k
      ENDIF

      CALL rttov_init_opt_param(err, opts, opt_param_k, .TRUE._jplm)
      CALL rttov_init_emis_refl(emissivity_tl, reflectance_tl)

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
            CALL rttov_init_opt_param(err, opts, opt_param_tl, .TRUE._jplm)
            IF (var == 1) THEN
              var_in => opt_param_tl%abs(lay,:)
              var_out => opt_param_k%abs(lay,:)
            ELSEIF (var == 2) THEN
              var_in => opt_param_tl%sca(lay,:)
              var_out => opt_param_k%sca(lay,:)
            ELSEIF (var == 3) THEN
              var_in => opt_param_tl%bpr(lay,:)
              var_out => opt_param_k%bpr(lay,:)
            ELSEIF (var == 4) THEN
              var_in => opt_param_tl%legcoef(i,lay,:)
              var_out => opt_param_k%legcoef(i,lay,:)
            ELSEIF (var == 5) THEN
              var_in => opt_param_tl%pha(i,lay,:)
              var_out => opt_param_k%pha(i,lay,:)
            ENDIF
            var_in(:) = 1._jprb

            CALL call_rttov_tl(err)
            THROW(err.NE.0)

            IF (opts%rt_all%switchrad) THEN
              DO ichan = 1, nchannels
                IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
                  var_out(ichan) = radiance_tl%bt(ichan)
                ELSE
                  var_out(ichan) = radiance_tl%total(ichan)
                ENDIF
              ENDDO
            ELSE
              var_out(:) = radiance_tl%total(:)
            ENDIF
          ENDDO ! nsize
        ENDDO ! nlayers
      ENDDO ! nvar
      CALL rttov_init_opt_param(err, opts, opt_param_tl, .TRUE._jplm)
    ENDDO ! aer/cld
  ENDIF


  ! Surface emissivity
  ! All channels can be perturbed together with just one call to RTTOV
  ! which is much more efficient than doing each one individually.
  
  emissivity_tl%emis_in          = 1._jprb
  emissivity_tl%specularity      = 0._jprb
  reflectance_tl%refl_in         = 0._jprb
  reflectance_tl%diffuse_refl_in = 0._jprb
  
  CALL call_rttov_tl(err)
  THROW(err.NE.0)
  
  DO ichan = 1, nchannels
    IF (do_rttovscatt .OR. &
        (opts%rt_all%switchrad .AND. .NOT. opts%rt_ir%pc%addpc .AND. &
        coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2)) THEN
      IF (do_radar) THEN
        emissivity_k(ichan)%emis_in = reflectivity_tl%azef(radar_k_lev,ichan)
      ELSE
        emissivity_k(ichan)%emis_in = radiance_tl%bt(ichan)
      ENDIF
    ELSE
      emissivity_k(ichan)%emis_in = radiance_tl%total(ichan)
    ENDIF
  ENDDO

  IF (opts%rt_all%do_lambertian) THEN

    ! Surface specularity
    ! All channels can be perturbed together with just one call to RTTOV
    ! which is much more efficient than doing each one individually.

    emissivity_tl%emis_in          = 0._jprb
    emissivity_tl%specularity      = 1._jprb
    reflectance_tl%refl_in         = 0._jprb
    reflectance_tl%diffuse_refl_in = 0._jprb

    CALL call_rttov_tl(err)
    THROW(err.NE.0)

    DO ichan = 1, nchannels
      IF (do_rttovscatt .OR. &
          (opts%rt_all%switchrad .AND. .NOT. opts%rt_ir%pc%addpc .AND. &
          coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2)) THEN
          IF (do_radar) THEN
            emissivity_k(ichan)%specularity = reflectivity_tl%azef(radar_k_lev,ichan)
          ELSE
            emissivity_k(ichan)%specularity = radiance_tl%bt(ichan)
          ENDIF
      ELSE
        emissivity_k(ichan)%specularity = radiance_tl%total(ichan)
      ENDIF
    ENDDO

  ENDIF ! do_lambertian

  ! Surface reflectances
  ! All channels can be perturbed together with just one call to RTTOV
  ! which is much more efficient than doing each one individually.

  IF (opts%rt_ir%addsolar .AND. .NOT. do_rttovscatt) THEN

    ! Surface BRDF

    emissivity_tl%emis_in          = 0._jprb
    emissivity_tl%specularity      = 0._jprb
    reflectance_tl%refl_in         = 1._jprb
    reflectance_tl%diffuse_refl_in = 0._jprb

    CALL call_rttov_tl(err)
    THROW(err.NE.0)

    DO ichan = 1, nchannels
      IF (opts%rt_all%switchrad) THEN
        IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
          reflectance_k(ichan)%refl_in = radiance_tl%bt(ichan)
        ELSE
          reflectance_k(ichan)%refl_in = radiance_tl%total(ichan)
        ENDIF
      ELSE
        reflectance_k(ichan)%refl_in = radiance_tl%total(ichan)
      ENDIF
    ENDDO

    ! Surface diffuse reflectance
    ! Only relevant if some diffuse_refl_in > 0

    IF (ANY(reflectance%diffuse_refl_in > 0._jprb)) THEN
      emissivity_tl%emis_in          = 0._jprb
      emissivity_tl%specularity      = 0._jprb
      reflectance_tl%refl_in         = 0._jprb
      reflectance_tl%diffuse_refl_in = 1._jprb

      CALL call_rttov_tl(err)
      THROW(err.NE.0)

      DO ichan = 1, nchannels
        IF (opts%rt_all%switchrad) THEN
          IF (coefs%coef%ss_val_chn(chanprof(ichan)%chan) < 2) THEN
            reflectance_k(ichan)%diffuse_refl_in = radiance_tl%bt(ichan)
          ELSE
            reflectance_k(ichan)%diffuse_refl_in = radiance_tl%total(ichan)
          ENDIF
        ELSE
          reflectance_k(ichan)%diffuse_refl_in = radiance_tl%total(ichan)
        ENDIF
      ENDDO
    ENDIF
  ENDIF ! addsolar

  ! Deallocation

  IF (opts%rt_ir%pc%addpc) THEN
    CALL rttov_alloc_pccomp(err, pccomp_tl, npcscores, 0_jpim, nchannels_rec = nchannels_rec)
    THROW(err.NE.0)
  ENDIF

  IF (nthreads <= 0 .AND. .NOT. do_rttovscatt) THEN
    CALL rttov_alloc_traj(                   &
          err,                  nprofiles,   &
          nchannels,            opts,        &
          profiles(1)%nlevels,  coefs,       &
          0_jpim,               traj = traj, &
          traj_tl = traj_tl)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_cld_opt_param) THEN
    CALL rttov_alloc_opt_param(err, cld_opt_param_tl, nchannels, profiles(1)%nlevels, &
                               cld_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  IF (opts%rt_ir%user_aer_opt_param) THEN
    CALL rttov_alloc_opt_param(err, aer_opt_param_tl, nchannels, profiles(1)%nlevels, &
                               aer_opt_param%nmom, nphangle, 0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_transmission(err, transmission_tl, profiles(1)%nlevels, nchannels, 0_jpim)
  THROWM(err.NE.0,"Deallocation of transmission_tl failed")

  CALL rttov_alloc_rad( &
          err,                 &
          nchannels,           &
          radiance_tl,         &
          profiles(1)%nlevels, &
          0_jpim)
  THROW(err.NE.0)

  IF (do_rttovscatt) THEN
    IF (do_radar) THEN
      CALL rttov_alloc_reflectivity( &
            err,                     &
            nchannels,               &
            reflectivity_tl,         &
            nlevels,                 &
            0_jpim)
      THROW(err.NE.0)
      DEALLOCATE(reflectivity_tl)
    ENDIF
    CALL rttov_alloc_scatt_prof( &
          err,               &
          nprofiles,         &
          cld_profiles_tl,   &
          nlevels,           &
          cld_profiles(1)%nhydro, &
          cld_profiles(1)%nhydro_frac,  &
          0_jpim)
    THROW(err.NE.0)
  ENDIF

  CALL rttov_alloc_prof( &
          err,                 &
          nprofiles,           &
          profiles_tl,         &
          profiles(1)%nlevels, &
          opts,                &
          0_jpim)
  THROW(err.NE.0)

  CATCH

CONTAINS

  SUBROUTINE call_rttov_tl(err)
    INTEGER(jpim), INTENT(OUT) :: err

    IF (nthreads <= 0) THEN
      IF (do_rttovscatt) THEN
        IF (do_radar) THEN
          CALL rttov_scatt_tl(                  &
              err,              opts_scatt,     &
              nlevels,          chanprof,       &
              frequencies,      profiles,       &
              cld_profiles,     coefs,          &
              coefs_scatt,      calcemis,       &
              emissivity,       profiles_tl,    &
              cld_profiles_tl,  emissivity_tl,  &
              radiance,         radiance_tl,    &
              reflectivity,     reflectivity_tl)
        ELSE
          CALL rttov_scatt_tl(                  &
              err,              opts_scatt,     &
              nlevels,          chanprof,       &
              frequencies,      profiles,       &
              cld_profiles,     coefs,          &
              coefs_scatt,      calcemis,       &
              emissivity,       profiles_tl,    &
              cld_profiles_tl,  emissivity_tl,  &
              radiance,         radiance_tl)
        ENDIF
      ELSE
        CALL rttov_tl( &
                err,                                  &
                chanprof,                             &
                opts,                                 &
                profiles,                             &
                profiles_tl,                          &
                coefs,                                &
                transmission,                         &
                transmission_tl,                      &
                radiance,                             &
                radiance_tl,                          &
                calcemis         = calcemis,          &
                emissivity       = emissivity,        &
                emissivity_tl    = emissivity_tl,     &
                calcrefl         = calcrefl,          &
                reflectance      = reflectance,       &
                reflectance_tl   = reflectance_tl,    &
                aer_opt_param    = aer_opt_param,     &
                aer_opt_param_tl = aer_opt_param_tl,  &
                cld_opt_param    = cld_opt_param,     &
                cld_opt_param_tl = cld_opt_param_tl,  &
                traj             = traj,              &
                traj_tl          = traj_tl,           &
                pccomp           = pccomp,            &
                pccomp_tl        = pccomp_tl,         &
                channels_rec     = channels_rec)
      ENDIF
    ELSE
      IF (do_rttovscatt) THEN
        IF (do_radar) THEN
          CALL rttov_parallel_scatt_tl(          &
              err,              opts_scatt,      &
              nlevels,          chanprof,        &
              frequencies,      profiles,        &
              cld_profiles,     coefs,           &
              coefs_scatt,      calcemis,        &
              emissivity,       profiles_tl,     &
              cld_profiles_tl,  emissivity_tl,   &
              radiance,         radiance_tl,     &
              reflectivity,     reflectivity_tl, &
              nthreads = nthreads)
        ELSE
          CALL rttov_parallel_scatt_tl(          &
              err,              opts_scatt,      &
              nlevels,          chanprof,        &
              frequencies,      profiles,        &
              cld_profiles,     coefs,           &
              coefs_scatt,      calcemis,        &
              emissivity,       profiles_tl,     &
              cld_profiles_tl,  emissivity_tl,   &
              radiance,         radiance_tl,     &
              nthreads = nthreads)
        ENDIF
      ELSE
        CALL rttov_parallel_tl( &
                err,                                  &
                chanprof,                             &
                opts,                                 &
                profiles,                             &
                profiles_tl,                          &
                coefs,                                &
                transmission,                         &
                transmission_tl,                      &
                radiance,                             &
                radiance_tl,                          &
                calcemis         = calcemis,          &
                emissivity       = emissivity,        &
                emissivity_tl    = emissivity_tl,     &
                calcrefl         = calcrefl,          &
                reflectance      = reflectance,       &
                reflectance_tl   = reflectance_tl,    &
                aer_opt_param    = aer_opt_param,     &
                aer_opt_param_tl = aer_opt_param_tl,  &
                cld_opt_param    = cld_opt_param,     &
                cld_opt_param_tl = cld_opt_param_tl,  &
                pccomp           = pccomp,            &
                pccomp_tl        = pccomp_tl,         &
                channels_rec     = channels_rec,      &
                nthreads         = nthreads)
      ENDIF
    ENDIF

  END SUBROUTINE call_rttov_tl

END SUBROUTINE rttov_k_tl
