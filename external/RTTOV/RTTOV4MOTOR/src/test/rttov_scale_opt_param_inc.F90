! Description:
!> @file
!!   Scales cloud/aerosol optical parameter increments
!
!> @brief
!!   Scales cloud/aerosol optical parameter increments
!!
!! @param[in,out] opt_param_inc     cloud/aerosol optical parameter increments to scale
!! @param[in]     factor            scale factor
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
SUBROUTINE rttov_scale_opt_param_inc(opt_param_inc, factor)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_opt_param
  IMPLICIT NONE
  TYPE(rttov_opt_param), INTENT(INOUT) :: opt_param_inc
  REAL(jprb),            INTENT(IN)    :: factor
!INTF_END

  IF (ASSOCIATED(opt_param_inc%abs)) opt_param_inc%abs = opt_param_inc%abs * factor
  IF (ASSOCIATED(opt_param_inc%sca)) opt_param_inc%sca = opt_param_inc%sca * factor
  IF (ASSOCIATED(opt_param_inc%bpr)) opt_param_inc%bpr = opt_param_inc%bpr * factor
  IF (ASSOCIATED(opt_param_inc%pha)) opt_param_inc%pha = opt_param_inc%pha * factor
  IF (ASSOCIATED(opt_param_inc%legcoef)) opt_param_inc%legcoef = opt_param_inc%legcoef * factor

END SUBROUTINE rttov_scale_opt_param_inc
