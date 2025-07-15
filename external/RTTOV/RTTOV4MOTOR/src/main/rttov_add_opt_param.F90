! Description:
!> @file
!!   Add two cloud/aerosol optical parameter structures.
!
!> @brief
!!   Add two cloud/aerosol optical parameter structures.
!!
!!
!! @param[in,out] opt_param      output summed cloud/aerosol optical parameter structure
!! @param[in]     opt_param1     first input cloud/aerosol optical parameter structure
!! @param[in]     opt_param2     second input cloud/aerosol optical parameter structure
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
SUBROUTINE rttov_add_opt_param(opt_param, opt_param1, opt_param2)

  USE rttov_types, ONLY : rttov_opt_param
  IMPLICIT NONE
  TYPE(rttov_opt_param), INTENT(INOUT) :: opt_param
  TYPE(rttov_opt_param), INTENT(IN)    :: opt_param1
  TYPE(rttov_opt_param), INTENT(IN)    :: opt_param2
!INTF_END

  IF (ASSOCIATED(opt_param%abs) .AND. ASSOCIATED(opt_param1%abs) .AND. ASSOCIATED(opt_param2%abs)) &
      opt_param%abs = opt_param1%abs + opt_param2%abs
  IF (ASSOCIATED(opt_param%sca) .AND. ASSOCIATED(opt_param1%sca) .AND. ASSOCIATED(opt_param2%sca)) &
      opt_param%sca = opt_param1%sca + opt_param2%sca
  IF (ASSOCIATED(opt_param%bpr) .AND. ASSOCIATED(opt_param1%bpr) .AND. ASSOCIATED(opt_param2%bpr)) &
      opt_param%bpr = opt_param1%bpr + opt_param2%bpr
  IF (ASSOCIATED(opt_param%pha) .AND. ASSOCIATED(opt_param1%pha) .AND. ASSOCIATED(opt_param2%pha)) &
      opt_param%pha = opt_param1%pha + opt_param2%pha
  IF (ASSOCIATED(opt_param%legcoef) .AND. ASSOCIATED(opt_param1%legcoef) .AND. ASSOCIATED(opt_param2%legcoef)) &
      opt_param%legcoef = opt_param1%legcoef + opt_param2%legcoef

END SUBROUTINE rttov_add_opt_param
