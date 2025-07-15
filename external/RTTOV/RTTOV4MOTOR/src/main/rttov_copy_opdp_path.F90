! Description:
!> @file
!!   Copy an internal optical depth path structure.
!
!> @brief
!!   Copy an internal optical depth path structure.
!!
!! @param[in]     opts           RTTOV options structure
!! @param[in,out] opdp_path1     copy of opdep_path structure
!! @param[in]     opdp_path2     input opdep_path structure
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_copy_opdp_path(opts, opdp_path1, opdp_path2)

  USE rttov_types, ONLY : rttov_options, rttov_opdp_path
  IMPLICIT NONE

  TYPE(rttov_options),   INTENT(IN)    :: opts
  TYPE(rttov_opdp_path), INTENT(INOUT) :: opdp_path1
  TYPE(rttov_opdp_path), INTENT(IN)    :: opdp_path2
!INTF_END

  opdp_path1%atm_level = opdp_path2%atm_level
  IF (opts%rt_ir%addsolar) THEN
    opdp_path1%sun_level_path2 = opdp_path2%sun_level_path2
  ENDIF
END SUBROUTINE 
