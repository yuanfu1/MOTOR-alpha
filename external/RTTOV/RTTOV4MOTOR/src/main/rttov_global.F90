! Description:
!> @file
!!   Definition of global variables for RTTOV
!
!> @brief
!!   Definition of global variables for RTTOV
!!
!! @details
!!   Currently this is used only to hold the logical unit used for writing
!!   error, warning, and informational messages. The unit can be set by the
!!   user by calling rttov_errorhandling. The unit is used by rttov_errorreport
!!   and a number of other subroutines.
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
MODULE rttov_global

  USE parkind1, ONLY : jpim
  USE rttov_const, ONLY : default_err_unit
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: error_unit

  !> Logical unit for output messages
  INTEGER(KIND=jpim) :: error_unit

  DATA error_unit /default_err_unit/

END MODULE rttov_global
