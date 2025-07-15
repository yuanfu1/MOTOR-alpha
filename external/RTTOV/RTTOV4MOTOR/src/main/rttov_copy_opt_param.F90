! Description:
!> @file
!!   Copy a cloud/aerosol optical parameter structure.
!
!> @brief
!!   Copy a cloud/aerosol optical parameter structure.
!!
!!
!! @param[in,out] opt_param1     copy of cloud/aerosol optical parameter structure
!! @param[in]     opt_param2     input cloud/aerosol optical parameter structure
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
SUBROUTINE rttov_copy_opt_param(opt_param1, opt_param2)

  USE rttov_types, ONLY : rttov_opt_param
  IMPLICIT NONE
  TYPE(rttov_opt_param), INTENT(INOUT) :: opt_param1
  TYPE(rttov_opt_param), INTENT(IN)    :: opt_param2
!INTF_END

  opt_param1%nmom = opt_param2%nmom
  IF (ASSOCIATED(opt_param1%abs) .AND. ASSOCIATED(opt_param2%abs)) &
      opt_param1%abs = opt_param2%abs
  IF (ASSOCIATED(opt_param1%sca) .AND. ASSOCIATED(opt_param2%sca)) &
      opt_param1%sca = opt_param2%sca
  IF (ASSOCIATED(opt_param1%bpr) .AND. ASSOCIATED(opt_param2%bpr)) &
      opt_param1%bpr = opt_param2%bpr
  IF (ASSOCIATED(opt_param1%phangle) .AND. ASSOCIATED(opt_param2%phangle)) &
      opt_param1%phangle = opt_param2%phangle
  IF (ASSOCIATED(opt_param1%pha) .AND. ASSOCIATED(opt_param2%pha)) &
      opt_param1%pha = opt_param2%pha
  IF (ASSOCIATED(opt_param1%legcoef) .AND. ASSOCIATED(opt_param2%legcoef)) &
      opt_param1%legcoef = opt_param2%legcoef
  IF (ASSOCIATED(opt_param1%phasefn_int%cosphangle) .AND. ASSOCIATED(opt_param2%phasefn_int%cosphangle)) &
      opt_param1%phasefn_int%cosphangle = opt_param2%phasefn_int%cosphangle
  IF (ASSOCIATED(opt_param1%phasefn_int%iphangle) .AND. ASSOCIATED(opt_param2%phasefn_int%iphangle)) &
      opt_param1%phasefn_int%iphangle = opt_param2%phasefn_int%iphangle

END SUBROUTINE rttov_copy_opt_param
