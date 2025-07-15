! Description:
!> @file
!!   Scales reflectivity increment
!
!> @brief
!!   Scales reflectivity increment
!!
!! @param[in,out] reflectivity_inc     reflectivity structure to scale
!! @param[in]     factor               scale factor
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
SUBROUTINE rttov_scale_reflectivity_inc(reflectivity_inc, factor)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_reflectivity
  IMPLICIT NONE
  TYPE(rttov_reflectivity), INTENT(INOUT) :: reflectivity_inc
  REAL(jprb),               INTENT(IN)    :: factor
!INTF_END

  reflectivity_inc%zef = reflectivity_inc%zef * factor
  reflectivity_inc%azef = reflectivity_inc%azef * factor

END SUBROUTINE rttov_scale_reflectivity_inc
