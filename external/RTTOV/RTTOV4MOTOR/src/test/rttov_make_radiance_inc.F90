! Description:
!> @file
!!   Compute radiance/BT increments suitable for use in AD
!
!> @brief
!!   Compute radiance/BT increments suitable for use in AD
!!
!! @details
!!   For channels at wavelengths above 3µm the radiance increment is calculated
!!   as the delta-Planck radiance for each channel corresponding to a change in
!!   temperature of 0.1K at 280K.
!!
!!   For channels below 3µm (i.e. visible/near-IR channels) the radiance
!!   increment is set to 0.1 mW/m2/sr/cm-1.
!!
!!   If switchrad is true a fixed BT increment of 0.1K is also set for all
!!   channels.
!!
!! @param[in]     coef             optical depth coefficients structure
!! @param[in,out] radiance_inc     computed radiance/BT increments
!! @param[in]     channels         all channels being simulated (i.e. chanprof(:)%chan)
!! @param[in]     nchannels        number of channels being simulated (i.e. SIZE(chanprof))
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
SUBROUTINE rttov_make_radiance_inc(coef, radiance_inc, channels, nchannels)

  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef, rttov_radiance
!INTF_OFF
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : planck_c1, planck_c2
  USE rttov_math_mod, ONLY : planck, planck_tl
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_coef),     INTENT(IN)    :: coef
  TYPE(rttov_radiance), INTENT(INOUT) :: radiance_inc
  INTEGER(KIND=jpim),   INTENT(IN)    :: nchannels
  INTEGER(KIND=jpim),   INTENT(IN)    :: channels(nchannels)
!INTF_END
  INTEGER(KIND=jpim) :: ichan
  REAL(KIND=jprb)    :: b

  radiance_inc%bt = 0.1_jprb

  DO ichan = 1, nchannels
    IF (coef%ss_val_chn(channels(ichan)) < 2) THEN
      CALL planck(coef%ff_cwn(channels(ichan)), planck_c1, planck_c2, 280._jprb, b)
      CALL planck_tl(coef%ff_cwn(channels(ichan)), planck_c1, planck_c2, &
                     280._jprb, 0.1_jprb, b, radiance_inc%total(ichan))
    ELSE
      radiance_inc%total(ichan) = 0.1_jprb
    ENDIF
  ENDDO

END SUBROUTINE rttov_make_radiance_inc
