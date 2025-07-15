! Description:
!> @file
!!   Initialise a reflectivity structure.
!
!> @brief
!!   Initialise a reflectivity structure.
!!
!! @param[in,out]  reflectivity   Reflectivity structure to initialise
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
!    Copyright 2020, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_init_reflectivity(reflectivity)

  USE rttov_types, ONLY : rttov_reflectivity
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_reflectivity), INTENT(INOUT) :: reflectivity
!INTF_END

  reflectivity%zef  = 0._jprb
  reflectivity%azef = 0._jprb

END SUBROUTINE
