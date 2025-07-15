! Description:
!> @file
!!   Allocate/deallocate any/all arrays and structures for RTTOV or
!!   RTTOV-SCATT adjoint models.
!
!> @brief
!!   Allocate/deallocate any/all arrays and structures for RTTOV or
!!   RTTOV-SCATT adjoint models.
!!
!! @details
!!   Argument order follows that of rttov_ad where possible.
!!   This allows all necessary input/output arguments to rttov_ad to be
!!   allocated or deallocated in a single subroutine call: pass only those
!!   arguments you wish to allocate/deallocate.
!!
!!   All arguments to RTTOV-SCATT can also be allocated using this routine.
!!
!!   The coefficient file(s) should have been read in before calling this
!!   subroutine.
!!
!!   NB Array arguments must be pointers.
!!
!! @param[out]    err                 status on exit
!! @param[in]     asw                 1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     nprofiles           number of profiles being simulated
!! @param[in]     nchanprof           total number of channels being simulated over all profiles
!! @param[in]     nlevels             number of levels in input profiles
!! @param         chanprof            specifies channels and profiles to simulate, optional
!! @param[in]     opts                options to configure the simulations
!! @param         profiles            input atmospheric profiles and surface variables, optional
!! @param         profiles_ad         output gradient wrt atmospheric profile and surface variables, optional
!! @param[in]     coefs               coefficients structure for instrument to simulate
!! @param[in,out] transmission        output transmittances, optional
!! @param[in,out] transmission_ad     input gradient wrt transmittances (usually zero), optional
!! @param[in,out] radiance            output radiances and corresponding BTs and BRFs, optional
!! @param[in,out] radiance_ad         input gradient wrt radiances or BTs, optional
!! @param[in,out] radiance2           secondary output radiances, optional
!! @param         calcemis            flags for internal RTTOV surface emissivity calculation, optional
!! @param         emissivity          input/output surface emissivities, optional
!! @param         emissivity_ad       output gradient wrt surface emissivities, optional
!! @param         calcrefl            flags for internal RTTOV surface BRDF calculation, optional
!! @param         reflectance         input/output surface BRDFs, input cloud top BRDF for simple cloud, optional
!! @param         reflectance_ad      output gradient wrt surface BRDFs, optional
!! @param[in]     aer_maxnmom         max number of Legendre coefficients specified for aerosol phase functions, optional
!! @param[in]     aer_nphangle        size of angular grid for aerosol phase functions, optional
!! @param[in]     aer_opt_param       input aerosol optical parameters, optional
!! @param[in]     aer_opt_param_ad    output gradient wrt aerosol optical parameters, optional
!! @param[in]     cld_maxnmom         max number of Legendre coefficients specified for cloud phase functions, optional
!! @param[in]     cld_nphangle        size of angular grid for cloud phase functions, optional
!! @param[in]     cld_opt_param       input cloud optical parameters, optional
!! @param[in]     cld_opt_param_ad    output gradient wrt cloud optical parameters, optional
!! @param[in,out] traj                RTTOV direct internal state, can be initialised outside RTTOV, optional
!! @param[in,out] traj_ad             RTTOV AD internal state, can be initialised outside RTTOV, optional
!! @param[in]     npcscores           number of PC scores to calculate for PC-RTTOV, optional
!! @param[in]     nchannels_rec       number of channels for which to calculate reconstructed radiances, optional
!! @param[in,out] pccomp              output PC scores and radiances from PC-RTTOV, optional
!! @param[in,out] pccomp_ad           input gradient wrt PC scores, radiances or BTs for PC-RTTOV, optional
!! @param         channels_rec        array for channels for which to calculate reconstructed radiances, optional
!! @param         frequencies         RTTOV-SCATT frequency indices for each channel
!! @param[in]     coef_scatt          RTTOV-SCATT coefficients structure for instrument to simulate (already populated), optional
!! @param[in]     nhydro_frac         number of RTTOV-SCATT hydrometeor cloud fractions, optional
!! @param         cld_profiles        input hydrometeor profiles for RTTOV-SCATT, optional
!! @param         cld_profiles_ad     output gradient wrt hydrometeors for RTTOV-SCATT, optional
!! @param[in,out] reflectivity        output reflectivities from RTTOV-SCATT radar simulator, optional
!! @param[in,out] reflectivity_ad     input gradient wrt reflectivities for RTTOV-SCATT radar simulator, optional
!! @param[in]     init                set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_ad( &
              err,              &
              asw,              &
              nprofiles,        &
              nchanprof,        &
              nlevels,          &
              chanprof,         &
              opts,             &
              profiles,         &
              profiles_ad,      &
              coefs,            &
              transmission,     &
              transmission_ad,  &
              radiance,         &
              radiance_ad,      &
              radiance2,        &
              calcemis,         &
              emissivity,       &
              emissivity_ad,    &
              calcrefl,         &
              reflectance,      &
              reflectance_ad,   &
              aer_maxnmom,      &
              aer_nphangle,     &
              aer_opt_param,    &
              aer_opt_param_ad, &
              cld_maxnmom,      &
              cld_nphangle,     &
              cld_opt_param,    &
              cld_opt_param_ad, &
              traj,             &
              traj_ad,          &
              npcscores,        &
              nchannels_rec,    &
              pccomp,           &
              pccomp_ad,        &
              channels_rec,     &
              frequencies,      &
              coef_scatt,       &
              nhydro_frac,      &
              cld_profiles,     &
              cld_profiles_ad,  &
              reflectivity,     &
              reflectivity_ad,  &
              init)

!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm

  USE rttov_types, ONLY : &
      rttov_options,       &
      rttov_coefs,         &
      rttov_scatt_coef,    &
      rttov_chanprof,      &
      rttov_emissivity,    &
      rttov_reflectance,   &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_transmission,  &
      rttov_radiance,      &
      rttov_radiance2,     &
      rttov_reflectivity,  &
      rttov_traj,          &
      rttov_pccomp,        &
      rttov_opt_param

  IMPLICIT NONE

  ! Input arguments
  INTEGER(jpim),             INTENT(OUT)             :: err                ! return code
  INTEGER(jpim),             INTENT(IN)              :: asw                ! 1=allocate, 0=deallocate
  INTEGER(jpim),             INTENT(IN)              :: nprofiles          ! number of profiles
  INTEGER(jpim),             INTENT(IN)              :: nchanprof          ! size of chanprof array
  INTEGER(jpim),             INTENT(IN)              :: nlevels            ! number of levels
  TYPE(rttov_options),       INTENT(IN)              :: opts               ! RTTOV options
  TYPE(rttov_coefs),         INTENT(IN),    TARGET   :: coefs              ! coefficients

  ! Output arguments
  TYPE(rttov_profile),       POINTER,       OPTIONAL :: profiles(:)        ! input profiles
  TYPE(rttov_profile),       POINTER,       OPTIONAL :: profiles_ad(:)     ! output profile AD
  TYPE(rttov_transmission),  INTENT(INOUT), OPTIONAL :: transmission       ! output transmission
  TYPE(rttov_transmission),  INTENT(INOUT), OPTIONAL :: transmission_ad    ! transmission AD
  TYPE(rttov_radiance),      INTENT(INOUT), OPTIONAL :: radiance           ! output radiances
  TYPE(rttov_radiance),      INTENT(INOUT), OPTIONAL :: radiance_ad        ! input radiance perturbation
  TYPE(rttov_radiance2),     INTENT(INOUT), OPTIONAL :: radiance2          ! secondary output radiances

  TYPE(rttov_chanprof),      POINTER,       OPTIONAL :: chanprof(:)        ! channel and profile indices
  LOGICAL(jplm),             POINTER,       OPTIONAL :: calcemis(:)        ! switch for emissivity calculations
  TYPE(rttov_emissivity),    POINTER,       OPTIONAL :: emissivity(:)      ! input/output emissivity values
  TYPE(rttov_emissivity),    POINTER,       OPTIONAL :: emissivity_ad(:)   ! input/output emissivity values AD
  LOGICAL(jplm),             POINTER,       OPTIONAL :: calcrefl(:)        ! switch for surface BRDF calculations
  TYPE(rttov_reflectance),   POINTER,       OPTIONAL :: reflectance(:)     ! input/output surface BRDF values
  TYPE(rttov_reflectance),   POINTER,       OPTIONAL :: reflectance_ad(:)  ! input/output surface BRDF values AD

  TYPE(rttov_traj),          INTENT(INOUT), OPTIONAL :: traj               ! RTTOV internal data
  TYPE(rttov_traj),          INTENT(INOUT), OPTIONAL :: traj_ad            ! RTTOV internal data AD

  ! Aerosol/cloud optical parameter arguments
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: aer_maxnmom        ! max number of aerosol phase fn Legendre coefs
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: aer_nphangle       ! size of angular grid for aerosol phase fns
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL :: aer_opt_param      ! input aerosol optical parameters
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL :: aer_opt_param_ad
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: cld_maxnmom        ! max number of cloud phase fn Legendre coefs
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: cld_nphangle       ! size of angular grid for cloud phase fns
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL :: cld_opt_param      ! input cloud optical parameters
  TYPE(rttov_opt_param),     INTENT(INOUT), OPTIONAL :: cld_opt_param_ad

  ! PC-RTTOV output arguments
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: npcscores          ! number of PC scores
  INTEGER(jpim),             INTENT(IN),    OPTIONAL :: nchannels_rec      ! number of reconstructed radiances
  TYPE(rttov_pccomp),        INTENT(INOUT), OPTIONAL :: pccomp             ! output PC scores
  TYPE(rttov_pccomp),        INTENT(INOUT), OPTIONAL :: pccomp_ad          ! input PC scores perturbation
  INTEGER(jpim),             POINTER,       OPTIONAL :: channels_rec(:)    ! reconstructed channel list

  ! RTTOV-SCATT arguments
  INTEGER(jpim),             POINTER,       OPTIONAL :: frequencies(:)     ! frequency indices
  TYPE(rttov_scatt_coef),    INTENT(IN),    OPTIONAL :: coef_scatt         ! RTTOV-SCATT coefficients
  INTEGER(jpim),                            OPTIONAL :: nhydro_frac        ! number of hydrometeor cloud fractions for RTTOV-SCATT
  TYPE(rttov_profile_cloud), POINTER,       OPTIONAL :: cld_profiles(:)    ! input hydrometeor profiles
  TYPE(rttov_profile_cloud), POINTER,       OPTIONAL :: cld_profiles_ad(:) ! output hydrometeor profiles AD
  TYPE(rttov_reflectivity),  INTENT(INOUT), OPTIONAL :: reflectivity       ! output radar reflectivities
  TYPE(rttov_reflectivity),  INTENT(INOUT), OPTIONAL :: reflectivity_ad    ! input radar reflectivity perturbation

  ! Initialisation flag
  LOGICAL(jplm),             INTENT(IN),    OPTIONAL :: init               ! flag to select initialisation of structures
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_alloc_traj.interface"

!- End of header --------------------------------------------------------

  TRY

  ! Check AD arguments before allocating anything
  IF (PRESENT(pccomp_ad) .AND. .NOT. PRESENT(npcscores)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, 'npscores required when (de)allocating pccomp_ad')
  ENDIF

  CALL rttov_alloc_direct( &
          err,                        &
          asw,                        &
          nprofiles,                  &
          nchanprof,                  &
          nlevels,                    &
          chanprof,                   &
          opts,                       &
          profiles,                   &
          coefs,                      &
          transmission,               &
          radiance,                   &
          radiance2,                  &
          calcemis,                   &
          emissivity,                 &
          calcrefl,                   &
          reflectance,                &
          aer_maxnmom,                &
          aer_nphangle,               &
          aer_opt_param,              &
          cld_maxnmom,                &
          cld_nphangle,               &
          cld_opt_param,              &
          traj,                       &
          npcscores,                  &
          nchannels_rec,              &
          pccomp,                     &
          channels_rec,               &
          frequencies,                &
          coef_scatt,                 &
          nhydro_frac,                &
          cld_profiles,               &
          reflectivity=reflectivity,  &
          init=init)
  THROWM(err.NE.0, 'Error allocating direct structures')

  ! Use rttov_alloc_direct to allocate most of the AD variables
  CALL rttov_alloc_direct( &
          err,                            &
          asw,                            &
          nprofiles,                      &
          nchanprof,                      &
          nlevels,                        &
          opts=opts,                      &
          profiles=profiles_ad,           &
          coefs=coefs,                    &
          transmission=transmission_ad,   &
          radiance=radiance_ad,           &
          emissivity=emissivity_ad,       &
          reflectance=reflectance_ad,     &
          aer_maxnmom=aer_maxnmom,        &
          aer_nphangle=aer_nphangle,      &
          aer_opt_param=aer_opt_param_ad, &
          cld_maxnmom=cld_maxnmom,        &
          cld_nphangle=cld_nphangle,      &
          cld_opt_param=cld_opt_param_ad, &
          npcscores=npcscores,            &
          nchannels_rec=nchannels_rec,    &
          pccomp=pccomp_ad,               &
          coef_scatt=coef_scatt,          &
          nhydro_frac=nhydro_frac,        &
          cld_profiles=cld_profiles_ad,   &
          reflectivity=reflectivity_ad,   &
          init=init)
  THROWM(err.NE.0, 'Error allocating AD structures')

  ! traj_ad is (de)allocated here
  IF (PRESENT(traj_ad)) THEN
    CALL rttov_alloc_traj( &
            err,           &
            nprofiles,     &
            nchanprof,     &
            opts,          &
            nlevels,       &
            coefs,         &
            asw,           &
            traj_ad=traj_ad)
    THROWM(err.NE.0, 'Error (de)allocating traj_ad structure')
  ENDIF

  CATCH
END SUBROUTINE rttov_alloc_ad
