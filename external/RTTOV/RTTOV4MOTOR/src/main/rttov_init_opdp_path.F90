! Description:
!> @file
!!   Initialise an internal optical depth path structure.
!
!> @brief
!!   Initialise an internal optical depth path structure.
!!
!! @param[in]     opts           RTTOV options structure
!! @param[in,out] opdp_path      opdep_path structure to initialise
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
SUBROUTINE rttov_init_opdp_path(opts, opdp_path)

  USE rttov_types, ONLY : rttov_options, rttov_opdp_path
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),   INTENT(IN)    :: opts
  TYPE(rttov_opdp_path), INTENT(INOUT) :: opdp_path
!INTF_END

  opdp_path%atm_level = 0._jprb
  IF (opts%rt_ir%addsolar) THEN
    opdp_path%sun_level_path2 = 0._jprb
  ENDIF
END SUBROUTINE 
