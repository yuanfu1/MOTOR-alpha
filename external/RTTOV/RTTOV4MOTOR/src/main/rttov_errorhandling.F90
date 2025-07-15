! Description:
!> @file
!!   Set logical unit for error messages.
!
!> @brief
!!   Set logical unit for error messages.
!!
!! @details
!!   This subroutine may optionally be called to set the
!!   output unit for error, warning and informational messages.
!!   The default unit is defined in rttov_const.F90.
!!
!! @param[in] err_unit  logical unit number, if negative the default is used
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_errorhandling(err_unit)

  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : default_err_unit
  USE rttov_global, ONLY : error_unit
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim), INTENT(IN) :: err_unit
!INTF_END
  !- End of header --------------------------------------------------------

  ! Default error message logical unit is defined in rttov_const
  IF (err_unit >= 0) THEN
     error_unit = err_unit
  ELSE
     error_unit = default_err_unit
  ENDIF

END SUBROUTINE rttov_errorhandling
