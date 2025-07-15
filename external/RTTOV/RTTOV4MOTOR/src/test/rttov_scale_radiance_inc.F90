! Description:
!> @file
!!   Scales radiance or brightness temperature increment
!
!> @brief
!!   Scales radiance or brightness temperature increment
!!
!! @details
!!   The scale factor is applied to radiances and additionally to BTs if
!!   PC scores if switchrad is true.
!!
!! @param[in,out] radiance_inc     radiance structure to scale
!! @param[in]     factor           scale factor
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
SUBROUTINE rttov_scale_radiance_inc(radiance_inc, factor)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_radiance
  IMPLICIT NONE
  TYPE(rttov_radiance), INTENT(INOUT) :: radiance_inc
  REAL(jprb),           INTENT(IN)    :: factor
!INTF_END

  radiance_inc%bt = radiance_inc%bt * factor
  radiance_inc%total = radiance_inc%total * factor

END SUBROUTINE rttov_scale_radiance_inc
