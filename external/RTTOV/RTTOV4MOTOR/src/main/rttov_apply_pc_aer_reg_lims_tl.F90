! Description:
!> @file
!!   Clip aerosol profiles to PC-RTTOV regression limits TL.
!
!> @brief
!!   Clip aerosol profiles to PC-RTTOV regression limits TL.
!!
!! @param[in]     opts              options to configure the simulations
!! @param[in]     coef_pccomp       PC-RTTOV coefficient structure
!! @param[in]     profiles          user profiles
!! @param[in,out] profiles_int_tl   profiles in RTTOV internal units
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
SUBROUTINE rttov_apply_pc_aer_reg_lims_tl( &
        opts,            &
        coef_pccomp,     &
        profiles,        &
        profiles_int_tl, &
        aer_ref,         &
        aer_min,         &
        aer_max)

  USE rttov_types, ONLY : &
    rttov_options,     &
    rttov_coef_pccomp, &
    rttov_profile

  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim

  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(rttov_profile),     INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_int_tl(:)
  REAL(jprb),              INTENT(IN)    :: aer_ref(:,:,:)
  REAL(jprb),              INTENT(IN)    :: aer_min(:,:,:)
  REAL(jprb),              INTENT(IN)    :: aer_max(:,:,:)
!INTF_END

  INTEGER(jpim) :: t, lay, prof
  INTEGER(jpim) :: naer_coef, nlayers, nprofiles
  REAL(jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS_TL',0_jpim,ZHOOK_HANDLE)

  nprofiles = SIZE(profiles)
  nlayers = profiles(1)%nlayers
  naer_coef = SIZE(profiles(1)%aerosols, DIM=1)

  DO prof = 1, nprofiles

    DO lay = 1, nlayers
      IF (profiles(prof)%p(lay) >= profiles(prof)%s2m%p) CYCLE

      DO t = 1, naer_coef

        IF (t > coef_pccomp%fmv_pc_naer_types) THEN

          profiles_int_tl(prof)%aerosols(t,lay) = 0._jprb

        ELSE

          IF (aer_min(t,lay,prof) == aer_max(t,lay,prof)) THEN

            profiles_int_tl(prof)%aerosols(t,lay) = 0._jprb

          ELSE IF (aer_ref(t,lay,prof) < aer_min(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int_tl(prof)%aerosols(t,lay) = 0._jprb

          ELSE IF (aer_ref(t,lay,prof) > aer_max(t,lay,prof)) THEN

            IF (opts%config%apply_reg_limits) profiles_int_tl(prof)%aerosols(t,lay) = 0._jprb

          ENDIF

        ENDIF

      ENDDO ! aerosol types
    ENDDO ! layers
  ENDDO ! profiles

  IF (LHOOK) CALL DR_HOOK('RTTOV_APPLY_PC_AER_REG_LIMS_TL',1_jpim,ZHOOK_HANDLE)

END SUBROUTINE rttov_apply_pc_aer_reg_lims_tl
