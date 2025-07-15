! Description:
!> @file
!!   Clip aerosol profiles to PC-RTTOV regression limits.
!
!> @brief
!!   Clip aerosol profiles to PC-RTTOV regression limits.
!!
!! @details
!!   PC-RTTOV aerosol regression limits are provided in units of number
!!   density (cm^-3) on coefficient layers. This subroutine interpolates
!!   the layer limits onto the input pressure levels and compares against
!!   the input profiles converted to RTTOV internal units.
!!
!!   The climatological OPAC profiles used for training PC-RTTOV specify
!!   no aerosols in some layers (where the min and max are both zero) and
!!   have non-zero values in other layers which do not vary (where the min
!!   and max are equal and non-zero). This subroutine always enforces the
!!   limits for these layers where the training data have no variability.
!!
!!   In the remaining layers the aerosol profiles are clipped to the limits
!!   if the apply_reg_limits option is true (as done for gases). Otherwise
!!   it will print out warnings if the verbose option is true. In any case
!!   it sets the qflag_pc_aer_reg_limits bit in the radiance%quality output
!!   for PC predictor channels associated with the affected profile (this
!!   flag is not set for those layers where there is no aerosol variability).
!!
!! @param[in]     opts           options to configure the simulations
!! @param[in]     chanprof       specifies channels and profiles to simulate
!! @param[in]     coef_pccomp    PC-RTTOV coefficient structure
!! @param[in]     profiles       user profiles
!! @param[in,out] profiles_int   profiles in RTTOV internal units
!! @param[in,out] radiance       radiance structure
!! @param[in,out] aer_ref        original aerosol concentrations
!! @param[in,out] aer_min        aerosol min limits on user layers
!! @param[in,out] aer_max        aerosol max limits on user layers
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_apply_pc_aer_reg_lims( &
        opts,         &
        chanprof,     &
        coef_pccomp,  &
        profiles,     &
        profiles_int, &
        radiance,     &
        aer_ref,      &
        aer_min,      &
        aer_max)

  USE rttov_types, ONLY : &
    rttov_options,     &
    rttov_chanprof,    &
    rttov_coef_pccomp, &
    rttov_profile,     &
    rttov_radiance

  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim, jplm

  USE rttov_const, ONLY : &
    errorstatus_success, &
    qflag_pc_aer_reg_limits

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int(:)
  TYPE(rttov_radiance),    INTENT(INOUT) :: radiance
  REAL(jprb),              INTENT(INOUT) :: aer_ref(:,:,:)
  REAL(jprb),              INTENT(INOUT) :: aer_min(:,:,:)
  REAL(jprb),              INTENT(INOUT) :: aer_max(:,:,:)
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim) :: t, lay, prof
  INTEGER(jpim) :: naer_coef, nlayers, nprofiles
  INTEGER(jpim) :: lay_user, lev_user_upper, lev_user_lower
  INTEGER(jpim) :: lev_coef_upper, lev_coef_lower
  LOGICAL(jplm) :: reg_lim_flag_min, reg_lim_flag_max
  CHARACTER(32) :: sprof
  REAL(jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS',0_jpim,ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers
  naer_coef = SIZE(profiles(1)%aerosols, DIM=1)

  aer_min = 0.
  aer_max = 0.

  DO prof = 1, nprofiles

    ! See also rttov_user_profile_checkinput.F90
    IF (opts%interpolation%addinterp) THEN
      ! Interpolate aerosol regression limits

      DO lay_user = 1, nlayers
        lev_user_upper = lay_user
        lev_user_lower = lay_user + 1

        ! Ignore user layers which are entirely above the top coef level
        IF (profiles(prof)%p(lev_user_lower) < coef_pccomp%ref_pc_prfl_p(1)) CYCLE

        ! Find coef layers which overlap user layer (bounded by lev_coef_upper and lev_coef_lower)
        DO lev_coef_upper = coef_pccomp%fmv_pc_nlev-1, 2, -1
          ! Starting at the bottom of the profile, the upper coef level is the first one
          ! above (or equal to) the upper user pressure level (or failing that it's just the
          ! top level)
          IF (coef_pccomp%ref_pc_prfl_p(lev_coef_upper) <= profiles(prof)%p(lev_user_upper)) EXIT
        ENDDO

        DO lev_coef_lower = lev_coef_upper+1, coef_pccomp%fmv_pc_nlev-1
          ! Starting at the top of the profile, the lower coef level is the first one
          ! below (or equal to) the lower user pressure level (or failing that it's just the
          ! bottom level)
          IF (coef_pccomp%ref_pc_prfl_p(lev_coef_lower) >= profiles(prof)%p(lev_user_lower)) EXIT
        ENDDO

        ! Set user-layer limits based on limits from overlapping coef layers
        DO t = 1, coef_pccomp%fmv_pc_naer_types
          ! Most generous limit from overlapping layers (smallest min, largest max)
          aer_min(t,lay_user,prof) = MINVAL(coef_pccomp%lim_pc_prfl_aermin(t,lev_coef_upper:lev_coef_lower-1))
          aer_max(t,lay_user,prof) = MAXVAL(coef_pccomp%lim_pc_prfl_aermax(t,lev_coef_upper:lev_coef_lower-1))
        ENDDO

      ENDDO ! user layers

    ELSE
      aer_min(:,:coef_pccomp%fmv_pc_naer_types,prof) = coef_pccomp%lim_pc_prfl_aermin
      aer_max(:,:coef_pccomp%fmv_pc_naer_types,prof) = coef_pccomp%lim_pc_prfl_aermax
    ENDIF

    ! Store original profiles for TL/AD/K
    aer_ref(:,:,prof) = profiles_int(prof)%aerosols

    reg_lim_flag_min = .FALSE.
    reg_lim_flag_max = .FALSE.

    DO lay = 1, nlayers
      ! If the level bounding the top of the layer is not above the surface then
      ! this layer does not contribute to TOA radiance so we ignore it
      IF (profiles(prof)%p(lay) >= profiles(prof)%s2m%p) CYCLE

      DO t = 1, naer_coef

        IF (t > coef_pccomp%fmv_pc_naer_types) THEN

          ! Aerosol type not included in training at all so ensure it is zero
          profiles_int(prof)%aerosols(t,lay) = 0._jprb

        ELSE

          IF (aer_min(t,lay,prof) == aer_max(t,lay,prof)) THEN

            ! Fix aerosol concentrations in layers with no variability in training data
            ! (this includes layers with no aerosols in the training at all i.e min=max=0.)
            profiles_int(prof)%aerosols(t,lay) = aer_min(t,lay,prof)

          ELSE IF (profiles_int(prof)%aerosols(t,lay) < aer_min(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int(prof)%aerosols(t,lay) = aer_min(t,lay,prof)
            reg_lim_flag_min = .TRUE.

          ELSE IF (profiles_int(prof)%aerosols(t,lay) > aer_max(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int(prof)%aerosols(t,lay) = aer_max(t,lay,prof)
            reg_lim_flag_max = .TRUE.

          ENDIF
        ENDIF

      ENDDO ! aerosol types
    ENDDO ! layers

    IF (reg_lim_flag_min .OR. reg_lim_flag_max) THEN
      WHERE (chanprof%prof == prof)
        radiance%quality(1:SIZE(chanprof)) = IBSET(radiance%quality(1:SIZE(chanprof)), qflag_pc_aer_reg_limits)
      ENDWHERE
    ENDIF

    IF (opts%config%verbose) THEN
! For OpenMP only: Fortran I/O is not thread-safe
!$OMP CRITICAL
      WRITE(sprof,'(" (profile number = ",I8,")")') prof
!$OMP END CRITICAL
      IF (reg_lim_flag_min) THEN
        CALL rttov_errorreport(errorstatus_success, TRIM("PC-RTTOV: some aerosol exceeded lower coef limit "//sprof))
      ENDIF
      IF (reg_lim_flag_max) THEN
        CALL rttov_errorreport(errorstatus_success, TRIM("PC-RTTOV: some aerosol exceeded upper coef limit "//sprof))
      ENDIF
    ENDIF

  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_apply_pc_aer_reg_lims
