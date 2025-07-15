! Description:
!> @file
!!   Compute cloud/aerosol optical parameters increments suitable for use in TL
!
!> @brief
!!   Compute cloud/aerosol optical parameters increments suitable for use in TL
!!
!! @details
!!   The various members of the increment are computed as small fractions of
!!   the corresponding members of the input optical property structure.
!!
!! @param[in,out] opt_param_inc     computed optical parameter increments
!! @param[in]     opt_param         cloud/aerosol optical parameters
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
SUBROUTINE rttov_make_opt_param_inc(opt_param_inc, opt_param)

  USE rttov_types, ONLY : rttov_opt_param
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_opt_param), INTENT(INOUT) :: opt_param_inc
  TYPE(rttov_opt_param), INTENT(IN)    :: opt_param
!INTF_END

  IF (ASSOCIATED(opt_param_inc%abs)) opt_param_inc%abs = -1.E-3_jprb * opt_param%abs
  IF (ASSOCIATED(opt_param_inc%sca)) opt_param_inc%sca = -1.E-3_jprb * opt_param%sca
  IF (ASSOCIATED(opt_param_inc%bpr)) opt_param_inc%bpr = -1.E-3_jprb * opt_param%bpr
  IF (ASSOCIATED(opt_param_inc%pha)) opt_param_inc%pha = -1.E-5_jprb * opt_param%pha
  IF (ASSOCIATED(opt_param_inc%legcoef)) opt_param_inc%legcoef = -1.E-5_jprb * opt_param%legcoef

END SUBROUTINE rttov_make_opt_param_inc
