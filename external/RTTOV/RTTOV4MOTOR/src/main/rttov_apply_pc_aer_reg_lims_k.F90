! Description:
!> @file
!!   Clip aerosol profiles to PC-RTTOV regression limits K.
!
!> @brief
!!   Clip aerosol profiles to PC-RTTOV regression limits K.
!!
!! @param[in]     opts              options to configure the simulations
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     coef_pccomp       PC-RTTOV coefficient structure
!! @param[in]     profiles          user profiles
!! @param[in,out] profiles_int_k    profiles in RTTOV internal units
!! @param[in]     aer_ref           original aerosol concentrations
!! @param[in]     aer_min           aerosol min limits on user layers
!! @param[in]     aer_max           aerosol max limits on user layers
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
SUBROUTINE rttov_apply_pc_aer_reg_lims_k( &
        opts,           &
        chanprof,       &
        coef_pccomp,    &
        profiles,       &
        profiles_int_k, &
        aer_ref,        &
        aer_min,        &
        aer_max)

  USE rttov_types, ONLY : &
    rttov_options,     &
    rttov_chanprof,    &
    rttov_coef_pccomp, &
    rttov_profile

  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_chanprof),    INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_k(:)
  REAL(jprb),              INTENT(IN)    :: aer_ref(:,:,:)
  REAL(jprb),              INTENT(IN)    :: aer_min(:,:,:)
  REAL(jprb),              INTENT(IN)    :: aer_max(:,:,:)
!INTF_END

  INTEGER(jpim) :: t, lay, prof, j
  INTEGER(jpim) :: naer_coef, nlayers, nchanprof
  REAL(jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS_K',0_jpim,ZHOOK_HANDLE)

  nchanprof = SIZE(chanprof)
  nlayers = profiles(1)%nlayers
  naer_coef = SIZE(profiles(1)%aerosols, DIM=1)

  DO j = 1, nchanprof
    prof = chanprof(j)%prof

    DO lay = 1, nlayers
      IF (profiles(prof)%p(lay) >= profiles(prof)%s2m%p) CYCLE

      DO t = 1, naer_coef

        IF (t > coef_pccomp%fmv_pc_naer_types) THEN

          profiles_int_k(j)%aerosols(t,lay) = 0._jprb

        ELSE

          IF (aer_min(t,lay,prof) == aer_max(t,lay,prof)) THEN

            profiles_int_k(j)%aerosols(t,lay) = 0._jprb

          ELSE IF (aer_ref(t,lay,prof) < aer_min(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int_k(j)%aerosols(t,lay) = 0._jprb

          ELSE IF (aer_ref(t,lay,prof) > aer_max(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int_k(j)%aerosols(t,lay) = 0._jprb

          ENDIF

        ENDIF

      ENDDO ! aerosol types
    ENDDO ! layers
  ENDDO ! chanprof

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS_K',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_apply_pc_aer_reg_lims_k
