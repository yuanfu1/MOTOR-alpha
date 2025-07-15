! Description:
!> @file
!!   Compute reflectivity increments suitable for use in AD
!
!> @brief
!!   Compute reflectivity increments suitable for use in AD
!!
!! @param[in,out] reflectivity_inc     computed reflectivity increments
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
SUBROUTINE rttov_make_reflectivity_inc(reflectivity_inc)

  USE rttov_types, ONLY : rttov_reflectivity
!INTF_OFF
  USE rttov_test_k_mod, ONLY : radar_k_lev
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_reflectivity), INTENT(INOUT) :: reflectivity_inc
!INTF_END

!  reflectivity_inc%zef = 1._jprb
!  reflectivity_inc%azef = 1._jprb

  reflectivity_inc%zef = 0._jprb
  reflectivity_inc%azef = 0._jprb
  reflectivity_inc%azef(radar_k_lev,:) = 1._jprb

END SUBROUTINE rttov_make_reflectivity_inc
